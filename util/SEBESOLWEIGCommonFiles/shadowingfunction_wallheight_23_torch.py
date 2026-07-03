from __future__ import division

try:
    import torch
    import torch.nn.functional as F
except:
    pass
author = "xlinfr and Lemap01"


def shade_on_walls(azimuth, aspect, walls, dsm, f, shvoveg):
    # wall shadows wall parameterization
    wallbol = (walls > 0).float()

    # Removing walls in shadow due to selfshadowing
    azilow = azimuth - torch.pi / 2
    azihigh = azimuth + torch.pi / 2
    if azilow >= 0 and azihigh < 2 * torch.pi:  # 90 to 270  (SHADOW)
        facesh = (
            torch.logical_or(aspect < azilow, aspect >= azihigh).float()
            - wallbol
            + 1
        )  # TODO check
    elif azilow < 0 and azihigh <= 2 * torch.pi:  # 0 to 90
        azilow = azilow + 2 * torch.pi
        facesh = (
            torch.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1
        )  # (SHADOW)
    elif azilow > 0 and azihigh >= 2 * torch.pi:  # 270 to 360
        azihigh -= 2 * torch.pi
        facesh = (
            torch.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1
        )  # (SHADOW)

    shvo = f - dsm  # building shadow volume
    facesun = torch.logical_and(
        facesh + (walls > 0).float() == 1, walls > 0
    ).float()
    wallsun = torch.clone(walls - shvo)
    wallsun[wallsun < 0] = 0
    wallsun[facesh == 1] = 0  # Removing walls in "self"-shadow
    wallsh = torch.clone(walls - wallsun)

    wallshve = shvoveg * wallbol
    wallshve = wallshve - wallsh
    wallshve[wallshve < 0] = 0
    id = torch.where(wallshve > walls)
    wallshve[id] = walls[id].to(dtype=wallshve.dtype)
    wallsun = wallsun - wallshve  # problem with wallshve only
    id = torch.where(wallsun < 0)
    wallshve[id] = 0
    wallsun[id] = 0

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
    device,
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

    # conversion
    degrees = torch.pi / 180.0
    azimuth *= degrees
    altitude *= degrees

    # measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]

    # initialise parameters
    dx = torch.tensor(0.0, device=device)
    dy = torch.tensor(0.0, device=device)
    dz = torch.tensor(0.0, device=device)
    temp = torch.zeros((sizex, sizey), device=device)
    tempvegdem = torch.zeros((sizex, sizey), device=device)
    tempvegdem2 = torch.zeros((sizex, sizey), device=device)
    templastfabovea = torch.zeros((sizex, sizey), device=device)
    templastgabovea = torch.zeros((sizex, sizey), device=device)
    bushplant = bush > 1
    sh = torch.zeros((sizex, sizey), device=device)  # shadows from buildings
    vbshvegsh = torch.clone(sh)  # vegetation blocking buildings
    vegsh = torch.add(
        torch.zeros((sizex, sizey), device=device), bushplant
    ).to(
        torch.float64
    )  # vegetation shadow
    f = torch.clone(a)
    shvoveg = torch.clone(vegdem)  # for vegetation shadowvolume

    # other loop parameters
    pibyfour = torch.pi / 4
    threetimespibyfour = 3 * pibyfour
    fivetimespibyfour = 5 * pibyfour
    seventimespibyfour = 7 * pibyfour
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanazimuth = torch.tan(azimuth)
    signsinazimuth = torch.sign(sinazimuth)
    signcosazimuth = torch.sign(cosazimuth)
    dssin = torch.abs(1 / sinazimuth)
    dscos = torch.abs(1 / cosazimuth)
    tanaltitudebyscale = torch.tan(altitude) / scale

    index = 0

    # new case with pergola (thin vertical layer of vegetation), August 2021
    dzprev = 0

    # main loop
    while (
        (amaxvalue >= dz)
        and (torch.abs(dx) < sizex)
        and (torch.abs(dy) < sizey)
    ):
        if ((pibyfour <= azimuth) and (azimuth < threetimespibyfour)) or (
            (fivetimespibyfour <= azimuth) and (azimuth < seventimespibyfour)
        ):
            dy = signsinazimuth * index
            dx = (
                -1
                * signcosazimuth
                * torch.abs(torch.round(index / tanazimuth))
            )
            ds = dssin
        else:
            dy = signsinazimuth * torch.abs(torch.round(index * tanazimuth))
            dx = -1 * signcosazimuth * index
            ds = dscos

        # note: dx and dy represent absolute values while ds is an incremental
        # value
        dz = (ds * index) * tanaltitudebyscale
        tempvegdem[0:sizex, 0:sizey] = 0
        tempvegdem2[0:sizex, 0:sizey] = 0
        temp[0:sizex, 0:sizey] = 0
        templastfabovea[0:sizex, 0:sizey] = 0.0
        templastgabovea[0:sizex, 0:sizey] = 0.0
        absdx = torch.abs(dx)
        absdy = torch.abs(dy)
        xc1 = int((dx + absdx) / 2)
        xc2 = int(sizex + (dx - absdx) / 2)
        yc1 = int((dy + absdy) / 2)
        yc2 = int(sizey + (dy - absdy) / 2)
        xp1 = -int((dx - absdx) / 2)
        xp2 = int(sizex - (dx + absdx) / 2)
        yp1 = -int((dy - absdy) / 2)
        yp2 = int(sizey - (dy + absdy) / 2)

        tempvegdem[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2] - dz

        f = torch.fmax(f, temp)  # Moving building shadow
        # moving vegetation shadow volume
        shvoveg = torch.fmax(shvoveg, tempvegdem)

        sh[f > a] = 1
        sh[f <= a] = 0
        fabovea = (tempvegdem > a).int()  # vegdem above DEM
        gabovea = (tempvegdem2 > a).int()  # vegdem2 above DEM

        # new pergola condition
        templastfabovea[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dzprev
        templastgabovea[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dzprev
        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a
        dzprev = dz
        vegsh2 = torch.add(
            torch.add(
                torch.add(fabovea, gabovea).to(torch.float64), lastfabovea
            ).to(torch.float64),
            lastgabovea,
        ).to(torch.float64)
        vegsh2[vegsh2 == 4] = 0.0
        # vegsh2[vegsh2 == 1] = 0. # This one is the ultimate question...
        vegsh2[vegsh2 > 0] = 1.0

        vegsh = torch.fmax(vegsh, vegsh2)
        vegsh[vegsh * sh > 0] = 0
        # removing shadows 'behind' buildings
        vbshvegsh = torch.clone(vegsh) + vbshvegsh

        index += 1

    sh = 1 - sh
    vbshvegsh[vbshvegsh > 0] = 1
    vbshvegsh = vbshvegsh - vegsh

    vegsh[vegsh > 0] = 1
    shvoveg = (shvoveg - a) * vegsh  # Vegetation shadow volume
    vegsh = 1 - vegsh
    vbshvegsh = 1 - vbshvegsh

    wallsh, wallsun, wallshve, facesh, facesun = shade_on_walls(
        azimuth, aspect, walls, a, f, shvoveg
    )

    if walls_scheme is not False:
        wallsh_, wallsun_, wallshve_, facesh_, facesun_ = shade_on_walls(
            azimuth, aspect_scheme, walls_scheme, a, f, shvoveg
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
