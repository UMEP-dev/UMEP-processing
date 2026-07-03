# -*- coding: utf-8 -*-
from __future__ import division
from math import radians
import torch


def shade_on_walls(azimuth, aspect, walls, dsm, f, device):
    # wall shadows wall parameterization
    wallbol = (walls > 0).astype(float)

    # Removing walls in shadow due to selfshadowing
    azilow = azimuth - torch.pi / 2
    azihigh = azimuth + torch.pi / 2

    if azilow >= 0 and azihigh < 2 * torch.pi:  # 90 to 270  (SHADOW)
        facesh = (
            torch.logical_or(aspect < azilow, aspect >= azihigh).astype(float)
            - wallbol
            + 1
        )
    elif azilow < 0 and azihigh <= 2 * torch.pi:  # 0 to 90
        azilow = azilow + 2 * torch.pi
        # (SHADOW)    # check for the -1
        facesh = torch.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1
    elif azilow > 0 and azihigh >= 2 * torch.pi:  # 270 to 360
        azihigh = azihigh - 2 * torch.pi
        facesh = (
            torch.logical_or(aspect > azilow, aspect <= azihigh) * -1 + 1
        )  # (SHADOW)

    sh = torch.clone(f - dsm)  # shadow volume
    facesun = torch.logical_and(
        facesh + (walls > 0).astype(float) == 1, walls > 0
    ).astype(float)
    wallsun = torch.clone(walls - sh)
    wallsun[wallsun < 0] = 0
    wallsun[facesh == 1] = 0
    wallsh = torch.clone(walls - wallsun)

    sh = torch.logical_not(torch.logical_not(sh)).astype(float)
    sh = sh * -1 + 1

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
    device,
    walls_scheme=False,
    aspect_scheme=False,
):
    """
    This m.file calculates shadows on a DSM and shadow height on building
    walls.

    INPUTS:
    a = DSM
    azimuth and altitude = sun position
    scale= scale of DSM (1 meter pixels=1, 2 meter pixels=0.5)
    walls= pixel row 'outside' buildings. will be calculated if empty
    aspect = normal aspect of buildings walls

    OUTPUT:
    sh=ground and roof shadow
    wallsh = height of wall that is in shadow
    wallsun = hieght of wall that is in sun

    Fredrik Lindberg 2012-03-19
    fredrikl@gvc.gu.se

     Utdate 2013-03-13 - bugfix for walls alinged with sun azimuths

    :param a:
    :param azimuth:
    :param altitude:
    :param scale:
    :param walls:
    :param aspect:
    :return:
    """

    # conversion
    # degrees = torch.pi/180
    azimuth = radians(azimuth)
    altitude = radians(altitude)

    # measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]

    # initialise parameters
    f = torch.clone(a)
    dx = 0
    dy = 0
    dz = 0
    temp = torch.zeros((sizex, sizey), device=device)

    # other loop parameters
    amaxvalue = torch.max(a)
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

    index = 1

    # main loop
    while (
        (amaxvalue >= dz)
        and (torch.abs(dx) < sizex)
        and (torch.abs(dy) < sizey)
    ):

        if (pibyfour <= azimuth and azimuth < threetimespibyfour) or (
            fivetimespibyfour <= azimuth and azimuth < seventimespibyfour
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
        dz = ds * index * tanaltitudebyscale
        temp[0:sizex, 0:sizey] = 0

        absdx = torch.abs(dx)
        absdy = torch.abs(dy)

        xc1 = int((dx + absdx) / 2)
        xc2 = int(sizex + (dx - absdx) / 2)
        yc1 = int((dy + absdy) / 2)
        yc2 = int(sizey + (dy - absdy) / 2)

        xp1 = int(-((dx - absdx) / 2))
        xp2 = int(sizex - (dx + absdx) / 2)
        yp1 = int(-((dy - absdy) / 2))
        yp2 = int(sizey - (dy + absdy) / 2)

        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2] - dz
        f = torch.fmax(f, temp)  # Moving building shadow

        index = index + 1

    sh, wallsh, wallsun, facesh, facesun = shade_on_walls(
        azimuth, aspect, walls, a, f
    )
    if walls_scheme is not False:
        sh_, wallsh_, wallsun_, facesh_, facesun_ = shade_on_walls(
            azimuth, aspect_scheme, walls_scheme, a, f
        )
        shade_on_wall = wallsh_.copy()

    return (
        (sh, wallsh, wallsun, facesh, facesun, shade_on_wall)
        if walls_scheme is not False
        else (sh, wallsh, wallsun, facesh, facesun)
    )
