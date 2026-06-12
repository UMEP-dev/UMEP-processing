try:
    from sklearn.cluster import KMeans
except:
    pass
import numpy as np
import torch


def _to_tensor(x, device, dtype=torch.float32):
    if isinstance(x, torch.Tensor):
        return x.to(device)
    if x is None:
        return None
    return torch.tensor(x, dtype=dtype, device=device)


from . import svf_functions_torch as svf
from . import wallalgorithms_torch as wa


def wallscheme_prepare(
    dsm, scale, pixel_resolution, feedback, device=torch.device("cpu")
):
    dsm = _to_tensor(dsm, device)

    # Existing UMEP wall and aspect calculations
    walls = wa.findwalls_sp(
        dsm,
        2,
        device,
        torch.tensor([[1, 1, 1], [1, 0, 1], [1, 1, 1]], device=dsm.device),
    )
    walls_copy = torch.clone(walls)
    aspect = wa.filter1Goodwin_as_aspect_v3(
        walls_copy, scale, dsm, feedback, 100, device
    )

    walls_exact = walls.clone()
    walls_round = torch.ceil(walls).int()

    # 1. Vectorized identification of wall pixels
    wall_rows, wall_cols = torch.where(walls_round > 0)
    num_walls = wall_rows.shape[0]

    # 2. Assign Wall IDs (1-indexed) vectorially
    uniqueWallIDs = torch.zeros((dsm.shape[0], dsm.shape[1]), device=device)
    wall_indices = torch.arange(1, num_walls + 1, device=device)
    uniqueWallIDs[wall_rows, wall_cols] = wall_indices.float()

    # 3. Calculate number of voxels for every wall pixel at once
    # Extract wall attributes only for the valid wall locations
    walls_round_extracted = walls_round[wall_rows, wall_cols]
    walls_exact_extracted = walls_exact[wall_rows, wall_cols]

    number_of_voxels = (walls_round_extracted / pixel_resolution).int()

    # 4. Repeat wall properties based on their respective voxel counts
    # This replaces both the 'i' and 'j' loops entirely
    wall2d_id_tensor = torch.repeat_interleave(wall_indices, number_of_voxels)
    wall_height_tensor = torch.repeat_interleave(
        walls_round_extracted, number_of_voxels
    )
    wall_height_exact_tensor = torch.repeat_interleave(
        walls_exact_extracted, number_of_voxels
    )
    y_position_tensor = torch.repeat_interleave(wall_rows, number_of_voxels)
    x_position_tensor = torch.repeat_interleave(wall_cols, number_of_voxels)

    # 5. Generate voxel heights using a cumulative sequence sequence
    # This handles the equivalent of `j * pixel_resolution` vectorially
    # We create a 1-based local sequence for each repeated segment
    # e.g., if number_of_voxels is [2, 3], we generate [1, 2, 1, 2, 3]
    v_ids = torch.arange(wall2d_id_tensor.shape[0], device=device)
    # Find the starting index of each repeated wall sequence
    cum_voxels = torch.cumsum(number_of_voxels, dim=0)
    start_indices = torch.cat(
        [torch.tensor([0], device=device), cum_voxels[:-1]]
    )
    shift = torch.repeat_interleave(start_indices, number_of_voxels)

    local_voxel_index = v_ids - shift + 1
    voxel_height_tensor = local_voxel_index * pixel_resolution

    # 6. Append the final trailing zero-row (mimicking original code logic)
    zero_val = torch.tensor([0], device=device)

    wall2d_id_tensor = torch.cat([wall2d_id_tensor, zero_val])
    voxel_height_tensor = torch.cat([voxel_height_tensor, zero_val.float()])
    wall_height_tensor = torch.cat([wall_height_tensor, zero_val])
    wall_height_exact_tensor = torch.cat(
        [wall_height_exact_tensor, zero_val.float()]
    )
    y_position_tensor = torch.cat([y_position_tensor, zero_val])
    x_position_tensor = torch.cat([x_position_tensor, zero_val])

    # 7. Construct outputs
    voxelId_list = torch.arange(
        1, wall2d_id_tensor.shape[0] + 1, device=device
    )

    voxelTable = torch.column_stack(
        [
            voxelId_list.float(),
            voxel_height_tensor,
            wall_height_tensor.float(),
            wall_height_exact_tensor,
            wall2d_id_tensor.float(),
            y_position_tensor.float(),
            x_position_tensor.float(),
        ]
    )

    # 8. Construct dictionary mapping (Kept native Python as requested by return type)
    # We do this from the extracted wall tensors to avoid looping over the voxel table
    wall_dict = dict(
        zip(wall_indices.tolist(), walls_exact_extracted.tolist())
    )
    wall_dict[0] = 0.0

    # Convert specific lists to match original return signature formats if downstream requires it

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return (
        voxelTable,
        voxelId_list,
        wall_dict,
        walls,
        aspect,
        uniqueWallIDs,
        wall2d_id_tensor.tolist(),
        voxel_height_tensor.tolist(),
    )


def svf_for_voxels(
    dsm,
    dem,
    vegdsm,
    vegdsm2,
    transVeg,
    scale,
    usevegdem,
    pixel_resolution,
    voxelTable,
    svf_height,
    svf_array,
    svfbu_array,
    svfveg_array,
    svfaveg_array,
    svf_height_array,
    feedback,
    device=torch.device("cpu"),
):
    """This function calculates sky view factor at all voxel levels"""
    with torch.no_grad():
        dsm = _to_tensor(dsm, device)
        dem = _to_tensor(dem, device)
        vegdsm = _to_tensor(vegdsm, device)
        vegdsm2 = _to_tensor(vegdsm2, device)
        voxelTable = _to_tensor(voxelTable, device)
        svf_array = _to_tensor(svf_array, device)
        svfbu_array = _to_tensor(svfbu_array, device)
        svfveg_array = _to_tensor(svfveg_array, device)
        svfaveg_array = _to_tensor(svfaveg_array, device)
        svf_height_array = _to_tensor(svf_height_array, device)

        # Calculate where there are buildings and not. Used to elevate dem.
        ground = dsm - dem
        # Ground == 1 = ground
        ground[ground < 2] = 1.0
        # Ground == 0 = buildings
        ground[ground >= 2] = 0.0

        # Find maximum wall height, used to estimate how many iterations of svf_calc that are required
        maxWallHeight = torch.max(voxelTable[:, 2]) - svf_height

        # Counter to feedback current iteration
        counter = 1
        # How many iterations are required to calculate svf for all voxels
        loop_range = torch.arange(
            svf_height, maxWallHeight + svf_height, svf_height
        )

        # Loop for svf calculations of all voxel heights
        for i in loop_range:

            feedback.setProgressText(
                "SVF calculation number "
                + str(int(counter))
                + " of "
                + str(int(loop_range.shape[0]))
            )

            feedback.setProgressText(
                "Increasing ground level with " + str(i) + " meters."
            )

            # Elevate ground in dsm
            temp_dsm = ((dsm + i) * ground) + (dsm * (1 - ground))

            if usevegdem == 1:
                # Subtract from cdsm
                temp_cdsm = vegdsm - i
                temp_cdsm[temp_cdsm < 0] = 0
                # Subtract from tdsm
                temp_cdsm2 = vegdsm2 - i
                temp_cdsm2[temp_cdsm2 < 0] = 0
            else:
                temp_cdsm = dsm * 0.0
                temp_cdsm2 = dsm * 0.0
            # Calculate svf. wallScheme set to 0 as only svf is estimated and nothing on the location of voxels, etc.
            wallScheme = 0
            ret_ = svf.svfForProcessing153(
                temp_dsm,
                temp_cdsm,
                temp_cdsm2,
                scale,
                usevegdem,
                pixel_resolution,
                wallScheme,
                dem,
                feedback,
                device=device,
            )
            svfbu = ret_["svf"]

            if usevegdem == 0:
                svftotal = svfbu
                svfveg = ret_["svfveg"]
                svfaveg = ret_["svfaveg"]
            else:
                svfveg = ret_["svfveg"]
                svfaveg = ret_["svfaveg"]
                trans = transVeg / 100.0
                svftotal = svfbu - (1 - svfveg) * (1 - trans)

            # Get svf for each voxel
            voxel_y = torch.where(
                voxelTable[:, 1] == i + svf_height
            )  # +svf_height)
            for temp_y in voxel_y[0]:
                svf_array[temp_y] = svftotal[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svfbu_array[temp_y] = svfbu[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ].double()
                svfveg_array[temp_y] = svfveg[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svfaveg_array[temp_y] = svfaveg[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svf_height_array[temp_y] = i + svf_height  # +svf_height

            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break

            counter += 1

        svf_05 = svf_array.clone()
        svf_05[svf_05 > 0.5] = 0.5

        # Add svf arrays as volumns to voxelTable
        voxelTable = torch.column_stack(
            [
                voxelTable,
                svf_height_array,
                svf_array,
                svf_05,
                svfbu_array,
                svfveg_array,
                svfaveg_array,
            ]
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        return voxelTable


def svf_kmeans(
    dsm,
    dem,
    vegdsm,
    vegdsm2,
    wallHeights,
    transVeg,
    scale,
    usevegdem,
    pixel_resolution,
    voxelTable,
    clusters,
    svf_height,
    svf_array,
    svfbu_array,
    svfveg_array,
    svfaveg_array,
    svf_height_array,
    feedback,
    device=torch.device("cpu"),
):
    try:
        from sklearn.cluster import KMeans
    except:
        raise ImportError(
            "[UMEP-processing Error] pleas install sklearn via pip install scikit-learn or via osgeo4w"
        )

    with torch.no_grad():

        dsm = _to_tensor(dsm, device)
        dem = _to_tensor(dem, device)
        vegdsm = _to_tensor(vegdsm, device)
        vegdsm2 = _to_tensor(vegdsm2, device)
        wallHeights = _to_tensor(wallHeights, device)
        voxelTable = _to_tensor(voxelTable, device)
        svf_array = _to_tensor(svf_array, device)
        svfbu_array = _to_tensor(svfbu_array, device)
        svfveg_array = _to_tensor(svfveg_array, device)
        svfaveg_array = _to_tensor(svfaveg_array, device)
        svf_height_array = _to_tensor(svf_height_array, device)

        # Calculate where there are buildings and not. Used to elevate dem.
        ground = dsm - dem
        # Ground == 1 = ground
        ground[ground < 2] = 1.0
        # Ground == 0 = buildings
        ground[ground >= 2] = 0.0

        # building_heights = dsm - dem

        # Reshape data for clustering
        # data_reshaped = building_heights.reshape(-1, 1)
        data_reshaped = wallHeights.detach().cpu().numpy().reshape(-1, 1)

        # Apply K-means clustering
        kmeans = KMeans(n_clusters=clusters, random_state=0)
        labels = kmeans.fit_predict(data_reshaped)

        # Reshape the labels back into a torch tensor on the selected device
        kmeans_clusters = torch.from_numpy(
            labels.reshape(dsm.shape[0], dsm.shape[1])
        ).to(device)

        # Remove cluster representing ground areas, i.e. where dsm - dem = 0
        cluster_range = torch.arange(clusters, device=device)
        # cluster_range = cluster_range[cluster_range != torch.unique(kmeans_clusters[ground == 1])]

        # Array to store mean heights of clusters
        cluster_heights = torch.zeros((cluster_range.shape[0]), device=device)

        counter = 0
        for i in cluster_range:
            # cluster_heights[counter] = torch.round(building_heights[kmeans_clusters == i].mean())
            cluster_heights[counter] = (
                torch.round(wallHeights[kmeans_clusters == i].mean())
                - svf_height
            )  # Remove svf_height which is the voxel size to be below the top of the wall
            counter += 1

        # Unique heights based on mean height of clusters, sorted from min to max
        cluster_heights = torch.unique(cluster_heights)
        cluster_heights = cluster_heights[cluster_heights > 0]

        # Counter to feedback current iteration
        counter = 0

        for i in cluster_heights:
            if cluster_heights.shape[0] > 1:
                feedback.setProgressText(
                    "SVF calculation based on K-means. Calculation "
                    + str(int(counter + 1))
                    + " of "
                    + str(int(cluster_heights.shape[0]))
                    + " clusters."
                )
                feedback.setProgressText(
                    "Mean wall height of cluster is "
                    + str(int(i + svf_height))
                    + " meters. Increasing ground level with "
                    + str(int(i))
                    + " meters."
                )

            # Elevate ground in dsm
            temp_dsm = ((dsm + i) * ground) + (dsm * (1 - ground))
            # temp_dsm = dsm[ground == 1] + temp_mean

            if usevegdem == 1:
                # Subtract from cdsm
                temp_cdsm = vegdsm - i
                temp_cdsm[temp_cdsm < 0] = 0
                # Subtract from tdsm
                temp_cdsm2 = vegdsm2 - i
                temp_cdsm2[temp_cdsm2 < 0] = 0
            else:
                temp_cdsm = dsm * 0.0
                temp_cdsm2 = dsm * 0.0

            # Calculate svf. wallScheme set to 0 as only svf is estimated and nothing on the location of voxels, etc.
            wallScheme = 0
            ret_ = svf.svfForProcessing153(
                temp_dsm,
                temp_cdsm,
                temp_cdsm2,
                scale,
                usevegdem,
                pixel_resolution,
                wallScheme,
                dem,
                feedback,
                device=device,
            )

            svfbu = ret_["svf"]

            if usevegdem == 0:
                svftotal = svfbu
                svfveg = ret_["svfveg"]
                svfaveg = ret_["svfaveg"]
            else:
                svfveg = ret_["svfveg"]
                svfaveg = ret_["svfaveg"]
                trans = transVeg / 100.0
                svftotal = svfbu - (1 - svfveg) * (1 - trans)

            # Get svf for each voxel
            voxel_y = torch.where(
                voxelTable[:, 1] == i + svf_height
            )  # +svf_height)
            for temp_y in voxel_y[0]:
                svf_array[temp_y] = svftotal[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svfbu_array[temp_y] = svfbu[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ].double()
                svfveg_array[temp_y] = svfveg[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svfaveg_array[temp_y] = svfaveg[
                    int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                ]
                svf_height_array[temp_y] = i + svf_height

            if i == cluster_heights[-1]:
                temp_data = voxelTable[
                    voxelTable[:, 2] > i, :
                ]  # Get all walls that are taller than the mean of the lowest cluster
                unique_walls = torch.unique(
                    temp_data[:, 4]
                )  # Get their unique wall ids for slicing
                for (
                    unique_wall
                ) in (
                    unique_walls
                ):  # Loop over all unique walls lower than lowest cluster
                    temp_wall = temp_data[temp_data[:, 4] == unique_wall, :][
                        :, 1
                    ].max()  # Max height of highest voxel in unique_wall
                    temp_y = torch.where(
                        (voxelTable[:, 4] == unique_wall)
                        & (voxelTable[:, 1] == temp_wall)
                    )[
                        0
                    ]  # Get row of unique_wall and highest voxel in voxelTable

                    svf_array[temp_y] = (
                        0.5  # Set svf to 0.5 as these are the highest voxels and nothing or little should obstruct it, i.e. svf = 0.5
                    )
                    svfbu_array[temp_y] = svfbu[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ].double()  # save svfbu for current wall pixel
                    svfveg_array[temp_y] = svfveg[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ]  # save svfveg for current wall pixel
                    svfaveg_array[temp_y] = svfaveg[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ]  # save svfaveg for current wall pixel
                    svf_height_array[temp_y] = (
                        temp_wall  # set svf_height to highest voxel for current wall
                    )
            # Add 0.5 to highest voxel on walls that are lower than cluster with lowest mean height
            else:
                if counter == 0:
                    temp_data = voxelTable[
                        voxelTable[:, 2] < i, :
                    ]  # Get all walls that are lower than the mean of the lowest cluster
                else:
                    temp_data = voxelTable[
                        (voxelTable[:, 2] > cluster_heights[counter - 1])
                        & (voxelTable[:, 2] < i),
                        :,
                    ]
                unique_walls = torch.unique(
                    temp_data[:, 4]
                )  # Get their unique wall ids for slicing
                for (
                    unique_wall
                ) in (
                    unique_walls
                ):  # Loop over all unique walls lower than lowest cluster
                    temp_wall = temp_data[temp_data[:, 4] == unique_wall, :][
                        :, 1
                    ].max()  # Max height of highest voxel in unique_wall
                    temp_y = torch.where(
                        (voxelTable[:, 4] == unique_wall)
                        & (voxelTable[:, 1] == temp_wall)
                    )[
                        0
                    ]  # Get row of unique_wall and highest voxel in voxelTable
                    temp_svf = svftotal[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ]  # get the calculated svf for the wall pixel to check if it is higher or lower than 0.5
                    if (
                        temp_svf < 0.5
                    ):  # if current wall pixel is lower than 0.5, although it is estimated above the wall, save it at the highest voxel and use for interpolation
                        svf_array[temp_y] = temp_svf
                    else:  # else, give highest voxel a value of 0.5
                        svf_array[temp_y] = 0.5
                    svfbu_array[temp_y] = svfbu[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ].double()  # save svfbu for current wall pixel
                    svfveg_array[temp_y] = svfveg[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ]  # save svfveg for current wall pixel
                    svfaveg_array[temp_y] = svfaveg[
                        int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])
                    ]  # save svfaveg for current wall pixel
                    svf_height_array[temp_y] = (
                        temp_wall  # set svf_height to highest voxel for current wall
                    )

            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break

            counter += 1

        svf_05 = svf_array.clone()
        svf_05[svf_05 > 0.5] = 0.5

        # Add svf arrays as volumns to voxelTable
        voxelTable = torch.column_stack(
            [
                voxelTable,
                svf_height_array,
                svf_array,
                svf_05,
                svfbu_array,
                svfveg_array,
                svfaveg_array,
            ]
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        return voxelTable, cluster_heights


def interpolate_svf(voxelTable):
    with torch.no_grad():
        device = voxelTable.device
        voxelTable_np = voxelTable.cpu().numpy()

        unique_wall_pixels = np.unique(voxelTable_np[:, 4])
        unique_wall_pixels = unique_wall_pixels[unique_wall_pixels != 0]

        for unique_wall in unique_wall_pixels:
            # Mask for the current wall pixel
            wall_mask = voxelTable_np[:, 4] == unique_wall
            temp_data = voxelTable_np[
                wall_mask, :
            ]  # All data for current wall pixel

            # Mask for valid SVF calculation heights
            valid_svf_mask = temp_data[:, -1] != 0
            temp_heights = temp_data[valid_svf_mask, 1]  # Voxel heights
            temp_svf = temp_data[valid_svf_mask, -4]  # SVF values

            if len(temp_heights) == 1:
                calculated_svf = temp_svf[0] if len(temp_svf) > 0 else 0.0
                new_svf = np.full(temp_data.shape[0], calculated_svf)

            elif len(temp_heights) > 1:  # Interpolate
                new_svf = np.interp(temp_data[:, 1], temp_heights, temp_svf)
            else:
                continue

            voxelTable_np[wall_mask, -4] = new_svf
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        return torch.from_numpy(voxelTable_np).to(device)
