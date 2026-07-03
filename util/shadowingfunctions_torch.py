# -*- coding: utf-8 -*-
# Ready for python action!
from math import radians

try:
    import torch
    import torch.nn.functional as F
except:
    pass


def _to_tensor(x, device, dtype=torch.float32):
    if isinstance(x, torch.Tensor):
        return x.to(device)
    if x is None:
        return None
    return torch.tensor(x, dtype=dtype, device=device)


def shadowingfunctionglobalradiation(
    a, azimuth, altitude, scale, feedback, forsvf, device=torch.device("cpu")
):
    """
    Computes global radiation shadows on a DSM using PyTorch hardware acceleration.

    This function leverages `F.grid_sample` to perform global map translations,
    completely removing multi-index array slicing and trigonometric quadrant switches.
    It is heavily optimized for GPU matrix operations.

    Args:
        a (torch.Tensor or numpy.ndarray): Digital Surface Model (DSM) matrix.
        azimuth (float): Sun azimuth angle in degrees.
        altitude (float): Sun altitude angle in degrees.
        scale (float): Spatial resolution modifier (1 pixel = 1 meter -> 1.0).
        feedback (QgsProcessingFeedback or similar): Object used to report algorithm progress.
        forsvf (int): Flag to toggle progress reporting (0 = report progress, 1 = skip).
        device (torch.device, optional): Device context for execution. Defaults to CPU.

    Returns:
        torch.Tensor: Binary shadow mask (1.0 = in sun, 0.0 = in shadow) matching `a.dtype`.
    """
    # Ensure inputs are correctly formatted PyTorch tensors on the target device
    if not isinstance(a, torch.Tensor):
        a = torch.tensor(a, device=device)
    else:
        a = a.to(device=device)

    # Convert angular sun positions to radians
    degrees = torch.pi / 180.0
    azimuth = azimuth * degrees
    altitude = altitude * degrees

    sizex, sizey = a.shape

    # Track progress bounds if requested
    if forsvf == 0:
        barstep = max(sizex, sizey)
        total = 100.0 / barstep

    # Clone initial DSM layout to accumulate the upper shadow envelope
    f = a.clone()

    # --- COORDINATE GRID GENERATION FOR GLOBAL TRANSFORMATION ---
    y, x = torch.meshgrid(
        torch.linspace(-1, 1, sizex, device=device),
        torch.linspace(-1, 1, sizey, device=device),
        indexing="ij",
    )
    # Target PyTorch shape: (1, H, W, 2) with (X, Y) layout at dim=-1
    grille_base = torch.stack((x, y), dim=-1).unsqueeze(0)

    # Convert the 2D DSM to a 4D tensor chunk (1, 1, H, W) for grid sampling
    a_4d = a.unsqueeze(0).unsqueeze(0)

    # Pre-calculate sun path factors
    amaxvalue = torch.max(a).item()
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanaltitudebyscale = torch.tan(altitude) / scale

    # Dynamically estimate the safety limit for our vector operations
    max_steps = int(amaxvalue / tanaltitudebyscale.item()) + 2

    # Vectorized loop structure (Eliminates slice variables and block switches)
    for step_idx in range(1, max_steps):
        index = float(step_idx)

        if forsvf == 0:
            feedback.setProgress(int(index * total))

        # 1. Project ground physical shadow offset (in pixels)
        shift_x = index * sinazimuth / scale
        shift_y = index * cosazimuth / scale
        dz = index * tanaltitudebyscale

        # Safe breakout loop check if rays completely clear the canvas boundaries
        if abs(shift_x) >= sizey and abs(shift_y) >= sizex:
            break

        # 2. Shift the underlying layout grid space (-1 to 1 space boundaries)
        grille_deplacee = grille_base.clone()
        grille_deplacee[..., 0] -= (2.0 * shift_x) / (sizey - 1)
        grille_deplacee[..., 1] -= (2.0 * shift_y) / (sizex - 1)

        # 3. Retrieve translated elevation layout using bilinear interpolation
        temp = (
            F.grid_sample(
                a_4d, grille_deplacee, mode="bilinear", padding_mode="border"
            ).squeeze()
            - dz
        )

        # 4. Amalgamate upper terrain shadow bounds
        f = torch.maximum(f, temp)

    # Calculate local shadow gaps
    f = f - a

    # Binary conversion: Where f == 0, the ground matches the shadow envelope (it is in sun)
    sh = (f == 0).to(dtype=a.dtype)

    return sh


# @jit(nopython=True)
def shadowingfunction_20(
    a,
    vegdem,
    vegdem2,
    azimuth,
    altitude,
    scale,
    amaxvalue,
    bush,
    feedback,
    forsvf,
    device=torch.device("cpu"),
):

    a = _to_tensor(a, device)
    vegdem = _to_tensor(vegdem, device)
    vegdem2 = _to_tensor(vegdem2, device)
    bush = _to_tensor(bush, device)

    degrees = torch.pi / 180.0
    azimuth = azimuth * degrees
    altitude = altitude * degrees

    sizex = a.shape[0]
    sizey = a.shape[1]

    if forsvf == 0:
        barstep = max(sizex, sizey)
        total = 100.0 / barstep
        feedback.setProgress(0)

    dz = 0.0
    temp = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    tempvegdem = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    tempvegdem2 = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    templastfabovea = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    templastgabovea = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    bushplant = bush > 1.0
    sh = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    vbshvegsh = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    vegsh = bushplant.to(dtype=a.dtype)
    f = a

    shvoveg = vegdem.clone()
    y, x = torch.meshgrid(
        torch.linspace(-1, 1, sizex, device=device),
        torch.linspace(-1, 1, sizey, device=device),
        indexing="ij",
    )
    grille_base = torch.stack((x, y), dim=-1).unsqueeze(0)

    # Preparation of the layers packed for grid_sample
    a_4d = a.unsqueeze(0).unsqueeze(0).float()
    veg_4d = vegdem.unsqueeze(0).unsqueeze(0).float()
    veg2_4d = vegdem2.unsqueeze(0).unsqueeze(0).float()

    # Trigonometric params for the sun's step
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanaltitudebyscale = torch.tan(altitude) / scale

    # Calculation of the maximum of steps
    # We stop when the shadow exceed the max possible height
    max_steps = int(amaxvalue / tanaltitudebyscale) + 2

    # Facteur d'avancement du rayon (équivalent géométrique de ton ancien 'ds')
    # Plus besoin de gros blocs If/Else selon les quadrants du soleil !
    delta_z_par_pas = tanaltitudebyscale

    for index in range(1, max_steps):
        # 1. Calculations of the slide
        shift_x = index * sinazimuth / scale
        shift_y = index * cosazimuth / scale
        dz = index * delta_z_par_pas

        # Stop if the shadow is completely outside of the map
        if abs(shift_x) >= sizey and abs(shift_y) >= sizex:
            break

        # 2. Applying the offset to our normalized coordinate sheet (-1 to 1)
        grille_deplacee = grille_base.clone()
        grille_deplacee[..., 0] -= (2.0 * shift_x) / (sizey - 1)
        grille_deplacee[..., 1] -= (2.0 * shift_y) / (sizex - 1)

        # 3. Global Material Sampling (Instant Layer Swipe)
        temp = (
            F.grid_sample(
                a_4d, grille_deplacee, mode="bilinear", padding_mode="border"
            ).squeeze()
            - dz
        )
        tempvegdem = (
            F.grid_sample(
                veg_4d, grille_deplacee, mode="bilinear", padding_mode="border"
            ).squeeze()
            - dz
        )
        tempvegdem2 = (
            F.grid_sample(
                veg2_4d,
                grille_deplacee,
                mode="bilinear",
                padding_mode="border",
            ).squeeze()
            - dz
        )

        # 4. Update of cast shadow volumes
        f = torch.fmax(f, temp)
        shvoveg = torch.fmax(shvoveg, tempvegdem)

        sh = torch.where(
            f > a,
            torch.tensor(1.0, device=device),
            torch.tensor(0.0, device=device),
        )
        fabovea = tempvegdem > a
        gabovea = tempvegdem2 > a

        # 5. Pergola logic (Overlaying layers with the & operator)
        templastfabovea = tempvegdem + delta_z_par_pas
        templastgabovea = tempvegdem2 + delta_z_par_pas

        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a

        # If all 4 layers are True at the same time, it's a pergola (light passes through).
        is_pergola = fabovea & gabovea & lastfabovea & lastgabovea

        # The shadow exists if one of the layers is True, UNLESS it's a pergola.
        vegsh2 = (fabovea | gabovea | lastfabovea | lastgabovea) & (
            ~is_pergola
        )
        vegsh2 = vegsh2.float()

        # Accumulation of vegetation shadows
        vegsh = torch.maximum(vegsh, vegsh2)
        vegsh = torch.where(
            vegsh * sh > 0, torch.tensor(0.0, device=device), vegsh
        )
        vbshvegsh.add_(vegsh)

    sh = 1.0 - sh
    vbshvegsh = torch.where(
        vbshvegsh > 0.0, torch.tensor(1.0, device=device), vbshvegsh
    )
    vbshvegsh = vbshvegsh - vegsh
    vegsh = 1.0 - vegsh
    vbshvegsh = 1.0 - vbshvegsh

    shadowresult = {"sh": sh, "vegsh": vegsh, "vbshvegsh": vbshvegsh}
    return shadowresult


def shadowingfunction_findwallID(
    dsm,
    azimuth,
    altitude,
    scale,
    walls,
    uniqueWallIDs,
    dem,
    wall2d_id,
    voxel_height,
    voxelId_list,
    facesh,
    wall_dict,
    sh,
    device,
):
    """
    Identifies which wall ID and voxel height are visible along solar ray vectors from ground pixels.

    This version removes explicit pixel-slice shifting loops by utilizing PyTorch's hardware-accelerated
    `F.grid_sample` engine. It implements nearest-neighbor interpolation to prevent spatial distortion
    of discrete categorical Wall IDs during matrix transformation layers.

    Args:
        dsm (torch.Tensor): Digital Surface Model tensor.
        azimuth (float): Sun azimuth angle in degrees.
        altitude (float): Sun altitude angle in degrees.
        scale (float): Spatial resolution modifier (1 pixel = 1 meter -> 1.0).
        walls (torch.Tensor): Height profile of the pixels representing building walls.
        uniqueWallIDs (torch.Tensor): Matrix containing distinct structural identifiers for each wall.
        dem (torch.Tensor): Digital Elevation Model representing bare earth topography.
        wall2d_id (torch.Tensor or list): Reference map containing baseline wall components.
        voxel_height (torch.Tensor or list): Reference metric indicating structural slice heights.
        voxelId_list (torch.Tensor or list): Linear index tracker for specific 3D voxel spaces.
        facesh (torch.Tensor): Binary mask identifying building walls shaded by their own geometry.
        wall_dict (dict): Dictionary mapping categorical string/int Wall IDs to their true heights.
        sh (torch.Tensor): Baseline ground shadow mask layer (1 = sun, 0 = shadow).
        device (torch.device): Execution hardware targeting context (CPU/GPU).

    Returns:
        tuple: A tuple containing:
            - buildIDSeen (torch.Tensor): Categorical ID of the wall casting a shadow onto the pixel.
            - voxelHeight (torch.Tensor): Accumulated vertical shadow volume height profile.
            - voxelId (torch.Tensor): Specific structural voxel block identifier seen by the layout.
    """
    # 1. PRE-PROCESSING & FORMATTING
    dsm = dsm - dem
    dsm = torch.where(dsm < 0.5, torch.tensor(0.0, device=device), dsm)

    # Conversion degrees to radians
    azimuth_rad = radians(azimuth)
    altitude_rad = radians(altitude)

    rows, cols = dsm.shape

    # Initialise tracked arrays
    buildIDSeen = torch.zeros((rows, cols), device=device)
    voxelHeight = torch.zeros((rows, cols), device=device)
    temp3 = torch.ones((rows, cols), device=device, dtype=torch.bool)

    # Mask wall tracking based on facing attributes
    uniqueWallIDs_masked = uniqueWallIDs * facesh

    # Build a high-performance native tensor lookup table for wall heights
    max_wall_id = int(max(wall_dict.keys())) if wall_dict else 0
    wall_height_lookup = torch.zeros(max_wall_id + 1, device=device)
    for k, v in wall_dict.items():
        wall_height_lookup[int(k)] = v

    # --- COORDINATE GRID GENERATION FOR THE SAMPLING TRANSFORMS ---
    y, x = torch.meshgrid(
        torch.linspace(-1, 1, rows, device=device),
        torch.linspace(-1, 1, cols, device=device),
        indexing="ij",
    )
    grille_base = torch.stack((x, y), dim=-1).unsqueeze(0)

    # Emballer uniqueWallIDs in 4D (1, 1, H, W) for grid_sample
    # We keep it as float for grid_sample, and convert back to long for indexing
    wall_id_4d = uniqueWallIDs_masked.unsqueeze(0).unsqueeze(0).float()

    # Trig variables
    sinazimuth = torch.sin(torch.tensor(azimuth_rad, device=device))
    cosazimuth = torch.cos(torch.tensor(azimuth_rad, device=device))
    tanaltitudebyscale = (
        torch.tan(torch.tensor(altitude_rad, device=device)) / scale
    )

    amaxvalue = torch.max(dsm)
    max_steps = int(amaxvalue / tanaltitudebyscale) + 2

    # 2. CORE VECTORIZED GRID-SHIFTING LOOP
    for index in range(1, max_steps):
        # Project spatial offsets on pixel dimensions
        shift_x = index * sinazimuth / scale
        shift_y = index * cosazimuth / scale
        dz = index * tanaltitudebyscale

        # Break early if ray projections completely escape spatial boundaries
        if abs(shift_x) >= cols and abs(shift_y) >= rows:
            break

        # Displace the structural tracking sheet
        grille_deplacee = grille_base.clone()
        grille_deplacee[..., 0] -= (2.0 * shift_x) / (cols - 1)
        grille_deplacee[..., 1] -= (2.0 * shift_y) / (rows - 1)

        # CRITICAL OPTIMIZATION: mode="nearest" keeps Wall IDs integer-pure (no interpolation fuzziness)
        tempwallID = F.grid_sample(
            wall_id_4d, grille_deplacee, mode="nearest", padding_mode="zeros"
        ).squeeze()

        # Batch-map wall ID categories directly to structural heights via tensor indexing
        temp_wallHeight = wall_height_lookup[tempwallID.long()]

        # Track wall degradation down ray segments
        temp2 = temp_wallHeight - dz

        # Process logical flags globally on the GPU cores
        valid_mask = temp2 > 0
        active_pixels_mask = valid_mask & temp3

        # Amalgamate ray intersection metrics
        buildIDSeen = torch.where(active_pixels_mask, tempwallID, buildIDSeen)
        voxelHeight = torch.where(
            active_pixels_mask, temp_wallHeight - temp2, voxelHeight
        )

        # Update remaining target pixel profiles using fast tensor boolean operations
        temp3 = (temp2 <= 0) & (buildIDSeen == 0.0)

    # 3. POST-PROCESSING & VOXEL IDENTIFICATION
    voxelHeight_ceil = torch.ceil(voxelHeight)
    voxelId = torch.zeros((rows, cols), device=device)

    # Secure mapping arrays are local tensors
    if not isinstance(wall2d_id, torch.Tensor):
        wall2d_id = torch.tensor(wall2d_id, device=device)
    else:
        wall2d_id = wall2d_id.to(device)

    if not isinstance(voxel_height, torch.Tensor):
        voxel_height = torch.tensor(voxel_height, device=device)
    else:
        voxel_height = voxel_height.to(device)

    if not isinstance(voxelId_list, torch.Tensor):
        voxelId_list = torch.tensor(
            voxelId_list, dtype=torch.long, device=device
        )
    else:
        voxelId_list = voxelId_list.to(device=device, dtype=torch.long)

    # Find unique Wall ID & Voxel Height pairings found across the canvas
    stacked_profiles = torch.column_stack(
        [buildIDSeen.flatten(), voxelHeight_ceil.flatten()]
    )
    unique_combinations = torch.unique(stacked_profiles, dim=0)
    unique_combinations = unique_combinations[
        ~torch.all(unique_combinations == 0, dim=1)
    ]

    # Map the linear 3D voxel index references
    for temp_id, temp_height in unique_combinations:
        mask = (wall2d_id == temp_id) & (voxel_height == temp_height)
        temp_fill_id = voxelId_list[mask]

        pixel_mask = (buildIDSeen == temp_id) & (
            voxelHeight_ceil == temp_height
        )

        if temp_fill_id.numel() > 0:
            voxelId[pixel_mask] = temp_fill_id[0].float()
        else:
            buildIDSeen[pixel_mask] = 0.0
            voxelHeight_ceil[pixel_mask] = 0.0

    # Invert original shadow notation layout logic (1 - sh)
    # This aligns the mapping accurately to raw cast shadows
    shadow_correction = 1.0 - sh
    buildIDSeen *= shadow_correction
    voxelHeight *= shadow_correction
    voxelId *= shadow_correction

    return buildIDSeen, voxelHeight, voxelId
