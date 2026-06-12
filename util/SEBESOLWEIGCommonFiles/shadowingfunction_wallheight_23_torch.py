from __future__ import division

try:
    import torch
    import torch.nn.functional as F
except:
    pass


def shade_on_walls(azimuth, aspect, walls, dsm, f, shvoveg, device):
    cos_incidence = torch.cos(aspect - azimuth)
    facesh = (cos_incidence <= 0).float() * (walls > 0).float()

    shvo = f - dsm

    facesun = ((facesh == 0) & (walls > 0)).float()

    wallsun = walls - shvo
    wallsun = torch.clamp(wallsun, min=0.0)
    wallsun = torch.where(
        facesh == 1, torch.tensor(0.0, device=device), wallsun
    )

    wallsh = walls - wallsun

    wallshve = shvoveg * (walls > 0).float()
    wallshve = torch.clamp(wallshve - wallsh, min=0.0)

    mask_clip = wallshve > walls
    wallshve[mask_clip] = walls[mask_clip].to(dtype=wallshve.dtype)

    wallsun = wallsun - wallshve
    mask_neg = wallsun < 0
    wallshve[mask_neg] = 0
    wallsun[mask_neg] = 0

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()

    return wallsh, wallsun, wallshve, facesh, facesun


def shadowingfunction_wallheight_23(
    a,
    vegdem,
    vegdem2,
    azimuth,
    altitude,
    scale,
    amaxvalue,
    bush,
    walls,
    aspect,
    walls_scheme=False,
    aspect_scheme=False,
):
    """
    This function calculates shadows on a DSM and shadow height on building
    walls including both buildings and vegetion units.
    New functionallity to deal with pergolas, August 2021

    INPUTS:
    a = DSM
    vegdem = Vegetation canopy DSM (magl)
    vegdem2 = Trunkzone DSM (magl)
    azimuth and altitude = sun position
    scale= scale of DSM (1 meter pixels=1, 2 meter pixels=0.5)
    walls= pixel row 'outside' buildings. will be calculated if empty
    aspect=normal aspect of walls

    OUTPUT:
    sh=ground and roof shadow

    wallsh = height of wall that is in shadow
    wallsun = hieght of wall that is in sun

    original Matlab code:
    Fredrik Lindberg 2013-08-14
    fredrikl@gvc.gu.se

    :param a:
    :param vegdem:
    :param vegdem2:
    :param azimuth:
    :param altitude:
    :param scale:
    :param amaxvalue:
    :param bush:
    :param walls:
    :param aspect:
    :return:
    """

    # automatically get the device to use from the input
    device = a.device if isinstance(a, torch.Tensor) else torch.device("cpu")
    degrees = torch.pi / 180.0
    azimuth *= degrees
    altitude *= degrees

    sizex, sizey = a.shape
    bushplant = (bush > 1).float()

    # Tensor initialisation
    sh = torch.zeros((sizex, sizey), device=device)
    vbshvegsh = torch.zeros((sizex, sizey), device=device)
    vegsh = torch.zeros((sizex, sizey), device=device) + bushplant
    f = a.clone()
    shvoveg = vegdem.clone()

    # ---coordinates drag the layer  ---
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

    sh = 1 - sh
    vbshvegsh = torch.where(
        vbshvegsh > 0, torch.tensor(1.0, device=device), vbshvegsh
    )
    vbshvegsh = vbshvegsh - vegsh

    vegsh = torch.where(vegsh > 0, torch.tensor(1.0, device=device), vegsh)
    shvoveg = (shvoveg - a) * vegsh  # Vegetation shadow volume
    vegsh = 1 - vegsh
    vbshvegsh = 1 - vbshvegsh

    wallsh, wallsun, wallshve, facesh, facesun = shade_on_walls(
        azimuth, aspect, walls, a, f, shvoveg, device
    )

    if walls_scheme is not False:
        wallsh_, wallsun_, wallshve_, facesh_, facesun_ = shade_on_walls(
            azimuth, aspect_scheme, walls_scheme, a, f, shvoveg, device
        )

        shade_on_wall = wallsh_.clone()
        shade_on_wall[shade_on_wall < wallshve_] = wallshve_[
            shade_on_wall < wallshve_
        ]

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()

    return (
        (
            vegsh,
            sh,
            vbshvegsh,
            wallsh,
            wallsun,
            wallshve,
            facesh,
            facesun,
            shade_on_wall,
        )
        if walls_scheme is not False
        else (vegsh, sh, vbshvegsh, wallsh, wallsun, wallshve, facesh, facesun)
    )
