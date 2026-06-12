# -*- coding: utf-8 -*-
from __future__ import division
from math import radians

try:
    import torch
    import torch.nn.functional as F
except:
    pass


def shade_on_walls(azimuth, aspect, walls, dsm, f, device):
    """
    Calculates shadow heights and lighting states on building facades using vector geometry.

    This function determines wall segments that are self-shadowed (facing away from the sun)
    using the cosine of the angle difference, tracks the down-ray building shadow volume,
    and isolates sections of the walls that remain illuminated.

    Args:
        azimuth (float): Sun azimuth angle in radians.
        aspect (torch.Tensor): Aspect orientation of the walls in radians.
        walls (torch.Tensor): Height profile of the pixels representing building walls.
        dsm (torch.Tensor): Digital Surface Model (ground + building heights).
        f (torch.Tensor): Maximum building shadow volume tracked during the ray-tracing loop.
        device (torch.device): The PyTorch device (CPU/GPU) where computations occur.

    Returns:
        tuple: A tuple containing:
            - sh (torch.Tensor): Binary ground and roof shadow mask (1 = sun, 0 = shadow).
            - wallsh (torch.Tensor): Height of the wall that is in shadow (meters).
            - wallsun (torch.Tensor): Height of the wall that is in direct sun (meters).
            - facesh (torch.Tensor): Binary mask of walls in self-shadow (1 = shadowed by self).
            - facesun (torch.Tensor): Binary mask of walls facing the sun and not occluded.
    """
    # If the cosine of the angle difference is <= 0, the facade faces away from the sun (self-shadow)
    cos_incidence = torch.cos(aspect - azimuth)
    facesh = (cos_incidence <= 0).float() * (walls > 0).float()

    sh = f - dsm  # Shadow volume (height of cast shadow)

    # Facades facing the sun (facesh == 0) that are valid walls
    facesun = ((facesh == 0) & (walls > 0)).float()

    # Calculate the height of the wall segment receiving direct sunlight
    wallsun = walls - sh
    wallsun = torch.clamp(wallsun, min=0.0)
    wallsun = torch.where(
        facesh == 1, torch.tensor(0.0, device=device), wallsun
    )

    # The shadowed wall height is the total wall height minus the sunlit segment
    wallsh = walls - wallsun

    # Transform the cast shadow volume into a binary mask (0 = shadow, 1 = sun)
    sh = (sh > 0).float()
    sh = 1.0 - sh

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return sh, wallsh, wallsun, facesh, facesun


def shadowingfunction_wallheight_13(
    a,
    azimuth,
    altitude,
    scale,
    walls,
    aspect,
    walls_scheme=False,
    aspect_scheme=False,
):
    """
    Computes ground/roof shadows and shadow heights on building walls using a DSM.

    This function leverages hardware-accelerated grid sampling (`F.grid_sample`) in PyTorch
    to simulate solar ray tracing. It avoids explicit multi-index pixel slicing or nested
    loops, making it optimal for high-performance GPU execution.

    Args:
        a (torch.Tensor): Digital Surface Model (DSM) matrix.
        azimuth (float): Sun azimuth angle in degrees.
        altitude (float): Sun altitude angle in degrees.
        scale (float): Map scale modifier (e.g., 1 pixel = 1 meter -> 1.0; 2 meters -> 0.5).
        walls (torch.Tensor): Extruded standalone wall pixel heights.
        aspect (torch.Tensor): Surface normal orientation of the walls in radians.
        walls_scheme (torch.Tensor, optional): Alternative building scheme layout. Defaults to False.
        aspect_scheme (torch.Tensor, optional): Alternative building scheme wall orientation. Defaults to False.

    Returns:
        tuple: Depending on whether `walls_scheme` is provided, returns:
            - Standard: (sh, wallsh, wallsun, facesh, facesun)
            - Scheme: (sh, wallsh, wallsun, facesh, facesun, shade_on_wall)
    """
    # Automatically detect the execution device (CPU or GPU) from the input tensor
    device = a.device if isinstance(a, torch.Tensor) else torch.device("cpu")

    # Convert angular sun coordinates to radians
    azimuth = radians(azimuth)
    altitude = radians(altitude)

    # Extract map grid spatial dimensions
    sizex, sizey = a.shape

    # Clone the DSM matrix to track shifting horizon shadows
    f = a.clone()

    # --- COORDINATE GRID GENERATION FOR GLOBAL TRANSFORMATION ---
    y, x = torch.meshgrid(
        torch.linspace(-1, 1, sizex, device=device),
        torch.linspace(-1, 1, sizey, device=device),
        indexing="ij",
    )
    # Reshape grid layout to PyTorch standard: (1, H, W, 2) with (X, Y) ordered at dim=-1
    grille_base = torch.stack((x, y), dim=-1).unsqueeze(0)

    # Prepare 4D-batched representation of the DSM for the sampling engine (1, 1, H, W)
    a_4d = a.unsqueeze(0).unsqueeze(0)

    # Extract geometric limits and spatial sun vector components
    amaxvalue = torch.max(a)
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanaltitudebyscale = torch.tan(altitude) / scale

    # Dynamically compute the maximum loop range needed based on terrain peak height
    max_steps = int(amaxvalue / tanaltitudebyscale) + 2

    # Core Ray-Tracing Loop (Accelerated using structural translation layers)
    for index in range(1, max_steps):
        # 1. Project physical shadow displacement on the horizontal plane (in pixels)
        shift_x = index * sinazimuth / scale
        shift_y = index * cosazimuth / scale
        dz = index * tanaltitudebyscale

        # Early termination checkpoint if the casting rays clear the map canvas bounds
        if abs(shift_x) >= sizey and abs(shift_y) >= sizex:
            break

        # 2. Shift the normalized coordinate sheet (-1 to 1 space)
        grille_deplacee = grille_base.clone()
        grille_deplacee[..., 0] -= (2.0 * shift_x) / (sizey - 1)
        grille_deplacee[..., 1] -= (2.0 * shift_y) / (sizex - 1)

        # 3. Retrieve the shifted layer map natively using bilinear interpolation
        temp = (
            F.grid_sample(
                a_4d, grille_deplacee, mode="bilinear", padding_mode="border"
            ).squeeze()
            - dz
        )

        # 4. Amalgamate shadows via upper envelope tracking
        f = torch.maximum(f, temp)

    # --- FACADE EVALUATION & SHADOW METRICS POST-PROCESSING ---
    sh, wallsh, wallsun, facesh, facesun = shade_on_walls(
        azimuth, aspect, walls, a, f, device
    )

    if walls_scheme is not False:
        _, wallsh_, _, _, _ = shade_on_walls(
            azimuth, aspect_scheme, walls_scheme, a, f, device
        )
        shade_on_wall = wallsh_.clone()
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        return (sh, wallsh, wallsun, facesh, facesun, shade_on_wall)
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return (sh, wallsh, wallsun, facesh, facesun)
