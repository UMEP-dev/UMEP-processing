import numpy as np

try:
    from sklearn.cluster import KMeans
except:
    pass

import time
from ..functions import svf_functions as svf
from ..functions import wallalgorithms as wa

def wallscheme_prepare(dsm, scale, pixel_resolution, feedback):
    total = 100. / (int(dsm.shape[0] * dsm.shape[1]))
    walls = wa.findwalls_sp(dsm, 2, np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]]))
    walls_copy = np.copy(walls)
    aspect = wa.filter1Goodwin_as_aspect_v3(walls_copy, scale, dsm, feedback, 100)

    # Copy to keep exact height values
    walls_exact = walls.copy()

    # Rounding wall heights (ceil or round?)
    # walls_round = np.round(walls).astype(int)
    walls_round = np.ceil(walls).astype(int)

    # create wall IDs
    wall_rows, wall_cols = np.where(walls_round > 0)
    voxel_height = list()
    wall2d_id = list()
    wall_height = list()
    wall_height_exact = list()
    y_position = list()
    x_position = list()
    index = 1
    uniqueWallIDs = np.zeros((dsm.shape[0], dsm.shape[1]))
    for i in np.arange(wall_rows.shape[0]):
        uniqueWallIDs[wall_rows[i], wall_cols[i]] = index
        number_of_voxels = int(walls_round[wall_rows[i], wall_cols[i]] / pixel_resolution)
        #temp_aspect = wallAspect[wall_rows[i], wall_cols[i]]
        voxel_index = 1
        for j in range(1, number_of_voxels+1):
            # wall_id.append(voxel_index)
            wall2d_id.append(index)
            voxel_height.append(j * pixel_resolution)
            wall_height.append(walls_round[wall_rows[i], wall_cols[i]])
            wall_height_exact.append(walls_exact[wall_rows[i], wall_cols[i]])
            y_position.append(wall_rows[i])
            x_position.append(wall_cols[i])
            #wall_aspect.append(temp_aspect)

            voxel_index += 1

        index += 1    

    wall2d_id.append(0)
    voxel_height.append(0)
    wall_height.append(0)
    wall_height_exact.append(0)
    y_position.append(0)
    x_position.append(0)

    wall_dict = {}
    for A, B in zip(wall2d_id, wall_height_exact):
        wall_dict[A] = B

    # saveraster(dataSet, output_uniquewallid, uniqueWallIDs)

    # Unique IDs for each voxel
    voxelId_list = np.arange(1, wall2d_id.__len__()+1)

    # Table with unique voxel ID, height of voxel, total height of wall, unique ID of wall (based on 2D-location in raster) and y and x coordinates
    voxelTable = np.column_stack([voxelId_list, voxel_height, wall_height, wall_height_exact, wall2d_id, y_position, x_position])

    return voxelTable, voxelId_list, wall_dict, walls, aspect, uniqueWallIDs, wall2d_id, voxel_height

def svf_for_voxels(dsm, dem, vegdsm, vegdsm2, transVeg, scale, usevegdem, pixel_resolution, voxelTable, 
                   svf_height, svf_array, svfbu_array, svfveg_array, svfaveg_array, svf_height_array, feedback):

    '''This function calculates sky view factor at all voxel levels'''

    # Calculate where there are buildings and not. Used to elevate dem.
    ground = dsm - dem
    # Ground == 1 = ground
    ground[ground < 2] = 1.
    # Ground == 0 = buildings
    ground[ground >= 2] = 0.

    # Find maximum wall height, used to estimate how many iterations of svf_calc that are required
    maxWallHeight = np.max(voxelTable[:,2]) - svf_height

    # Counter to feedback current iteration
    counter = 1
    # How many iterations are required to calculate svf for all voxels
    loop_range = np.arange(svf_height, maxWallHeight + svf_height, svf_height)
    
    # Loop for svf calculations of all voxel heights
    for i in loop_range:
        
        feedback.setProgressText('SVF calculation number ' + str(int(counter)) + ' of ' + str(int(loop_range.shape[0])))
        
        feedback.setProgressText('Increasing ground level with ' + str(i) + ' meters.')

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
            temp_cdsm = dsm * 0.
            temp_cdsm2 = dsm * 0.
        # Calculate svf. wallScheme set to 0 as only svf is estimated and nothing on the location of voxels, etc.
        wallScheme = 0
        ret_ = svf.svfForProcessing153(temp_dsm, temp_cdsm, temp_cdsm2, scale, usevegdem, pixel_resolution, wallScheme, dem, feedback)
        svfbu = ret_["svf"]

        if usevegdem == 0:
            svftotal = svfbu
            svfveg = ret_['svfveg']
            svfaveg = ret_['svfaveg']
        else:
            svfveg = ret_["svfveg"]
            svfaveg = ret_["svfaveg"]
            trans = transVeg / 100.0
            svftotal = (svfbu - (1 - svfveg) * (1 - trans))                    

        # Get svf for each voxel
        voxel_y = np.where(voxelTable[:, 1] == i + svf_height)# +svf_height)
        for temp_y in voxel_y[0]:
            svf_array[temp_y] = svftotal[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfbu_array[temp_y] = svfbu[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfveg_array[temp_y] = svfveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfaveg_array[temp_y] = svfaveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svf_height_array[temp_y] = i + svf_height# +svf_height

        if feedback.isCanceled():
            feedback.setProgressText("Calculation cancelled")
            break

        counter += 1
    
    svf_05 = svf_array.copy()
    svf_05[svf_05 > 0.5] = 0.5

    # Add svf arrays as volumns to voxelTable
    voxelTable = np.column_stack([voxelTable, svf_height_array, svf_array, svf_05, svfbu_array, svfveg_array, svfaveg_array])

    return voxelTable

def svf_kmeans(dsm, dem, vegdsm, vegdsm2, wallHeights, transVeg, scale, usevegdem, pixel_resolution, voxelTable, clusters,
                svf_height, svf_array, svfbu_array, svfveg_array, svfaveg_array, svf_height_array, feedback):
    
    # Calculate where there are buildings and not. Used to elevate dem.
    ground = dsm - dem
    # Ground == 1 = ground
    ground[ground < 2] = 1.
    # Ground == 0 = buildings
    ground[ground >= 2] = 0.

    # building_heights = dsm - dem

    # Reshape data for clustering
    # data_reshaped = building_heights.reshape(-1, 1)
    data_reshaped = wallHeights.reshape(-1, 1)

    # Apply K-means clustering
    # clusters = 3 # Number of clusters
    kmeans = KMeans(n_clusters=clusters, random_state=0)
    labels = kmeans.fit_predict(data_reshaped)

    # Reshape the labels back to the original data shape
    kmeans_clusters = labels.reshape(dsm.shape[0], dsm.shape[1])

    # Remove cluster representing ground areas, i.e. where dsm - dem = 0
    cluster_range = np.arange(clusters)
    # cluster_range = cluster_range[cluster_range != np.unique(kmeans_clusters[ground == 1])]

    # Array to store mean heights of clusters
    cluster_heights = np.zeros((cluster_range.shape[0]))

    counter = 0
    for i in cluster_range:
        # cluster_heights[counter] = np.round(building_heights[kmeans_clusters == i].mean())
        cluster_heights[counter] = np.round(wallHeights[kmeans_clusters == i].mean()) - svf_height # Remove svf_height which is the voxel size to be below the top of the wall
        counter += 1
    
    # Unique heights based on mean height of clusters, sorted from min to max
    cluster_heights = np.unique(cluster_heights)
    cluster_heights = cluster_heights[cluster_heights > 0]

    # Counter to feedback current iteration
    counter = 0

    for i in cluster_heights:
        if cluster_heights.shape[0] > 1:
            feedback.setProgressText('SVF calculation based on K-means. Calculation ' + str(int(counter + 1)) + ' of ' + str(int(cluster_heights.shape[0])) + ' clusters.')
            feedback.setProgressText('Mean wall height of cluster is ' + str(int(i + svf_height)) + ' meters. Increasing ground level with ' + str(int(i)) + ' meters.')
        
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
            temp_cdsm = dsm * 0.
            temp_cdsm2 = dsm * 0.
        
        # Calculate svf. wallScheme set to 0 as only svf is estimated and nothing on the location of voxels, etc.
        wallScheme = 0
        ret_ = svf.svfForProcessing153(temp_dsm, temp_cdsm, temp_cdsm2, scale, usevegdem, pixel_resolution, wallScheme, dem, feedback)
        svfbu = ret_["svf"]

        if usevegdem == 0:
            svftotal = svfbu
            svfveg = ret_['svfveg']
            svfaveg = ret_['svfaveg']
        else:
            svfveg = ret_["svfveg"]
            svfaveg = ret_["svfaveg"]
            trans = transVeg / 100.0
            svftotal = (svfbu - (1 - svfveg) * (1 - trans))                    

        # Get svf for each voxel
        voxel_y = np.where(voxelTable[:, 1] == i + svf_height)# +svf_height)
        for temp_y in voxel_y[0]:
            svf_array[temp_y] = svftotal[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfbu_array[temp_y] = svfbu[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfveg_array[temp_y] = svfveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svfaveg_array[temp_y] = svfaveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])]
            svf_height_array[temp_y] = i + svf_height

        if i == cluster_heights[-1]:
            temp_data = voxelTable[voxelTable[:, 2] > i, :] # Get all walls that are taller than the mean of the lowest cluster
            unique_walls = np.unique(temp_data[:, 4]) # Get their unique wall ids for slicing            
            for unique_wall in unique_walls: # Loop over all unique walls lower than lowest cluster
                temp_wall = temp_data[temp_data[:, 4] == unique_wall, :][:, 1].max() # Max height of highest voxel in unique_wall
                temp_y = np.where((voxelTable[:, 4] == unique_wall) & (voxelTable[:, 1] == temp_wall))[0] # Get row of unique_wall and highest voxel in voxelTable

                svf_array[temp_y] = 0.5 # Set svf to 0.5 as these are the highest voxels and nothing or little should obstruct it, i.e. svf = 0.5
                svfbu_array[temp_y] = svfbu[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfbu for current wall pixel
                svfveg_array[temp_y] = svfveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfveg for current wall pixel
                svfaveg_array[temp_y] = svfaveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfaveg for current wall pixel
                svf_height_array[temp_y] = temp_wall # set svf_height to highest voxel for current wall
        # Add 0.5 to highest voxel on walls that are lower than cluster with lowest mean height
        else:
            if counter == 0:
                temp_data = voxelTable[voxelTable[:, 2] < i, :] # Get all walls that are lower than the mean of the lowest cluster
            else:
                temp_data = voxelTable[(voxelTable[:, 2] > cluster_heights[counter - 1]) & (voxelTable[:, 2] < i), :]    
            unique_walls = np.unique(temp_data[:, 4]) # Get their unique wall ids for slicing
            for unique_wall in unique_walls: # Loop over all unique walls lower than lowest cluster
                temp_wall = temp_data[temp_data[:, 4] == unique_wall, :][:, 1].max() # Max height of highest voxel in unique_wall
                temp_y = np.where((voxelTable[:, 4] == unique_wall) & (voxelTable[:, 1] == temp_wall))[0] # Get row of unique_wall and highest voxel in voxelTable
                temp_svf = svftotal[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # get the calculated svf for the wall pixel to check if it is higher or lower than 0.5
                if temp_svf < 0.5: # if current wall pixel is lower than 0.5, although it is estimated above the wall, save it at the highest voxel and use for interpolation
                    svf_array[temp_y] = temp_svf
                else: # else, give highest voxel a value of 0.5
                    svf_array[temp_y] = 0.5
                svfbu_array[temp_y] = svfbu[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfbu for current wall pixel
                svfveg_array[temp_y] = svfveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfveg for current wall pixel
                svfaveg_array[temp_y] = svfaveg[int(voxelTable[temp_y, 5]), int(voxelTable[temp_y, 6])] # save svfaveg for current wall pixel
                svf_height_array[temp_y] = temp_wall # set svf_height to highest voxel for current wall

        if feedback.isCanceled():
            feedback.setProgressText("Calculation cancelled")
            break

        counter += 1

    svf_05 = svf_array.copy()
    svf_05[svf_05 > 0.5] = 0.5

    # Add svf arrays as volumns to voxelTable
    voxelTable = np.column_stack([voxelTable, svf_height_array, svf_array, svf_05, svfbu_array, svfveg_array, svfaveg_array])

    return voxelTable, cluster_heights

def interpolate_svf(voxelTable, cluster_heights, kmeans):

    unique_wall_pixels = np.unique(voxelTable[:, 4])
    unique_wall_pixels = unique_wall_pixels[unique_wall_pixels != 0]
    for unique_wall in unique_wall_pixels:
        temp_data = voxelTable[voxelTable[:, 4] == unique_wall, :] # All data for current wall pixel
        temp_heights = temp_data[temp_data[:, -1] != 0, 1] # Voxel heights for current wall pixel where svf has been calculated
        temp_svf = temp_data[temp_data[:, -1] != 0, -4] # SVF at voxel heights where svf has been calculated
        if temp_heights.size == 1:
            new_svf = temp_data[temp_data[:, -4] != 0, -4]
            new_svf[new_svf == 0] = new_svf[new_svf != 0]
        elif temp_heights.size > 1: # Interpolate
            new_svf = np.interp(temp_data[:, 1], temp_heights, temp_svf) # SVF for all voxels from interpolated values of calculated SVF at different heights (depend on svf_height)
        
        voxelTable[voxelTable[:, 4] == unique_wall, -4] = new_svf # Add the new SVFs to table

    return voxelTable

