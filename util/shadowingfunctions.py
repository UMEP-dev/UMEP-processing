# -*- coding: utf-8 -*-
# Ready for python action!
import torch
from math import radians
import matplotlib.pylab as plt
import numpy as np


def _to_tensor(x, device, dtype=torch.float32):
    if isinstance(x, torch.Tensor):
        return x.to(device)
    if x is None:
        return None
    return torch.tensor(x, dtype=dtype, device=device)


def shadowingfunctionglobalradiation(
    a, azimuth, altitude, scale, feedback, forsvf, device=torch.device("cpu")
):
    a = _to_tensor(a, device)

    degrees = torch.pi / 180.0
    azimuth = azimuth * degrees
    altitude = altitude * degrees

    sizex = a.shape[0]
    sizey = a.shape[1]
    if forsvf == 0:
        barstep = max(sizex, sizey)
        total = 100.0 / barstep

    f = a
    dx = 0.0
    dy = 0.0
    dz = 0.0
    temp = torch.zeros((sizex, sizey), device=device, dtype=a.dtype)
    index = 1.0

    amaxvalue = torch.max(a).item()
    pibyfour = torch.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanazimuth = torch.tan(azimuth)
    signsinazimuth = torch.sign(sinazimuth).item()
    signcosazimuth = torch.sign(cosazimuth).item()
    dssin = torch.abs(1.0 / sinazimuth).item()
    dscos = torch.abs(1.0 / cosazimuth).item()
    tanaltitudebyscale = torch.tan(altitude).item() / scale

    while amaxvalue >= dz and abs(dx) < sizex and abs(dy) < sizey:
        if forsvf == 0:
            feedback.setProgress(int(index * total))

        if (
            (pibyfour <= azimuth and azimuth < threetimespibyfour)
            or (fivetimespibyfour <= azimuth and azimuth < seventimespibyfour)
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * abs(torch.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * abs(torch.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos

        dz = ds * index * tanaltitudebyscale
        temp[:, :] = 0.0

        absdx = abs(dx)
        absdy = abs(dy)
        xc1 = int((dx + absdx) / 2.0 + 1.0)
        xc2 = int(sizex + (dx - absdx) / 2.0)
        yc1 = int((dy + absdy) / 2.0 + 1.0)
        yc2 = int(sizey + (dy - absdy) / 2.0)
        xp1 = int(-((dx - absdx) / 2.0) + 1.0)
        xp2 = int(sizex - (dx + absdx) / 2.0)
        yp1 = int(-((dy - absdy) / 2.0) + 1.0)
        yp2 = int(sizey - (dy + absdy) / 2.0)

        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2] - dz
        f = torch.maximum(f, temp)
        index += 1.0

    f = f - a
    f = f == 0
    sh = f.to(dtype=a.dtype)

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

    dx = 0.0
    dy = 0.0
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

    pibyfour = torch.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanazimuth = torch.tan(azimuth)
    signsinazimuth = torch.sign(sinazimuth).item()
    signcosazimuth = torch.sign(cosazimuth).item()
    dssin = torch.abs(1.0 / sinazimuth).item()
    dscos = torch.abs(1.0 / cosazimuth).item()
    tanaltitudebyscale = torch.tan(altitude).item() / scale
    index = 0.0

    dzprev = 0.0

    while (amaxvalue >= dz) and (abs(dx) < sizex) and (abs(dy) < sizey):
        if forsvf == 0:
            feedback.setProgress(int(index * total))

        if (
            (pibyfour <= azimuth)
            and (azimuth < threetimespibyfour)
            or (fivetimespibyfour <= azimuth)
            and (azimuth < seventimespibyfour)
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * abs(torch.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * abs(torch.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos

        dz = (ds * index) * tanaltitudebyscale
        tempvegdem[:, :] = 0.0
        tempvegdem2[:, :] = 0.0
        temp[:, :] = 0.0
        templastfabovea[:, :] = 0.0
        templastgabovea[:, :] = 0.0

        absdx = abs(dx)
        absdy = abs(dy)
        xc1 = int((dx + absdx) / 2.0)
        xc2 = int(sizex + (dx - absdx) / 2.0)
        yc1 = int((dy + absdy) / 2.0)
        yc2 = int(sizey + (dy - absdy) / 2.0)
        xp1 = int(-((dx - absdx) / 2.0))
        xp2 = int(sizex - (dx + absdx) / 2.0)
        yp1 = int(-((dy - absdy) / 2.0))
        yp2 = int(sizey - (dy + absdy) / 2.0)

        tempvegdem[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2] - dz

        f = torch.maximum(f, temp)
        sh = (f > a).to(dtype=a.dtype)
        vbshvegsh = vbshvegsh

        fabovea = tempvegdem > a
        gabovea = tempvegdem2 > a

        templastfabovea[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dzprev
        templastgabovea[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dzprev
        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a
        dzprev = dz

        vegsh2 = (
            fabovea.to(dtype=a.dtype)
            + gabovea.to(dtype=a.dtype)
            + lastfabovea.to(dtype=a.dtype)
            + lastgabovea.to(dtype=a.dtype)
        )
        vegsh2[vegsh2 == 4.0] = 0.0
        vegsh2[vegsh2 > 0.0] = 1.0

        vegsh = torch.maximum(vegsh, vegsh2)
        vegsh[(vegsh * sh) > 0.0] = 0.0
        vbshvegsh = vegsh + vbshvegsh

        index += 1.0

    sh = 1.0 - sh
    vbshvegsh[vbshvegsh > 0.0] = 1.0
    vbshvegsh = vbshvegsh - vegsh
    vegsh = 1.0 - vegsh
    vbshvegsh = 1.0 - vbshvegsh

    shadowresult = {"sh": sh, "vegsh": vegsh, "vbshvegsh": vbshvegsh}
    return shadowresult


def shadowingfunction_20_old(
    a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, dlg, forsvf
):
    a = _to_tensor(a, torch.device("cpu"))
    vegdem = _to_tensor(vegdem, torch.device("cpu"))
    vegdem2 = _to_tensor(vegdem2, torch.device("cpu"))
    bush = _to_tensor(bush, torch.device("cpu"))

    degrees = torch.pi / 180.0
    if azimuth == 0.0:
        azimuth = 1e-12
    azimuth = azimuth * degrees
    altitude = altitude * degrees

    sizex = a.shape[0]
    sizey = a.shape[1]
    if forsvf == 0:
        barstep = max(sizex, sizey)
        dlg.progressBar.setRange(0, barstep)
        dlg.progressBar.setValue(0)

    dx = 0.0
    dy = 0.0
    dz = 0.0
    temp = torch.zeros((sizex, sizey), dtype=a.dtype)
    tempvegdem = torch.zeros((sizex, sizey), dtype=a.dtype)
    tempvegdem2 = torch.zeros((sizex, sizey), dtype=a.dtype)
    sh = torch.zeros((sizex, sizey), dtype=a.dtype)
    vbshvegsh = torch.zeros((sizex, sizey), dtype=a.dtype)
    vegsh = torch.zeros((sizex, sizey), dtype=a.dtype)
    tempbush = torch.zeros((sizex, sizey), dtype=a.dtype)
    f = a
    g = torch.zeros((sizex, sizey), dtype=a.dtype)
    bushplant = bush > 1.0

    pibyfour = torch.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanazimuth = torch.tan(azimuth)
    signsinazimuth = torch.sign(sinazimuth).item()
    signcosazimuth = torch.sign(cosazimuth).item()
    dssin = torch.abs(1.0 / sinazimuth).item()
    dscos = torch.abs(1.0 / cosazimuth).item()
    tanaltitudebyscale = torch.tan(altitude).item() / scale
    index = 1.0

    while amaxvalue >= dz and abs(dx) < sizex and abs(dy) < sizey:
        if forsvf == 0:
            dlg.progressBar.setValue(int(index))
        if (
            pibyfour <= azimuth
            and azimuth < threetimespibyfour
            or fivetimespibyfour <= azimuth
            and azimuth < seventimespibyfour
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * abs(torch.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * abs(torch.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos

        dz = ds * index * tanaltitudebyscale
        tempvegdem[:, :] = 0.0
        tempvegdem2[:, :] = 0.0
        temp[:, :] = 0.0

        absdx = abs(dx)
        absdy = abs(dy)
        xc1 = int((dx + absdx) / 2.0 + 1.0)
        xc2 = int(sizex + (dx - absdx) / 2.0)
        yc1 = int((dy + absdy) / 2.0 + 1.0)
        yc2 = int(sizey + (dy - absdy) / 2.0)
        xp1 = int(-((dx - absdx) / 2.0) + 1.0)
        xp2 = int(sizex - (dx + absdx) / 2.0)
        yp1 = int(-((dy - absdy) / 2.0) + 1.0)
        yp2 = int(sizey - (dy + absdy) / 2.0)

        tempvegdem[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dz
        tempvegdem2[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dz
        temp[xp1:xp2, yp1:yp2] = a[xc1:xc2, yc1:yc2] - dz

        f = torch.maximum(f, temp)
        sh[(f > a)] = 1.0
        sh[(f <= a)] = 0.0
        fabovea = tempvegdem > a
        gabovea = tempvegdem2 > a
        vegsh2 = fabovea.to(dtype=a.dtype) - gabovea.to(dtype=a.dtype)
        vegsh = torch.maximum(vegsh, vegsh2)
        vegsh[(vegsh * sh) > 0.0] = 0.0
        vbshvegsh = vegsh + vbshvegsh

        if index == 1.0:
            firstvegdem = tempvegdem - temp
            firstvegdem[firstvegdem <= 0.0] = 1000.0
            vegsh[firstvegdem < dz] = 1.0
            vegsh = vegsh * (vegdem2 > a).to(dtype=a.dtype)
            vbshvegsh = torch.zeros((sizex, sizey), dtype=a.dtype)

        if torch.max(bush) > 0.0 and torch.max((fabovea.to(dtype=a.dtype) * bush)) > 0.0:
            tempbush[:, :] = 0.0
            tempbush[xp1:xp2, yp1:yp2] = bush[xc1:xc2, yc1:yc2] - dz
            g = torch.maximum(g, tempbush)
            g = g * bushplant

        index += 1.0

    sh = 1.0 - sh
    vbshvegsh[vbshvegsh > 0.0] = 1.0
    vbshvegsh = vbshvegsh - vegsh

    if torch.max(bush) > 0.0:
        g = g - bush
        g[g > 0.0] = 1.0
        g[g < 0.0] = 0.0
        vegsh = vegsh - bushplant.to(dtype=a.dtype) + g
        vegsh[vegsh < 0.0] = 0.0

    vegsh[vegsh > 0.0] = 1.0
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
    device
):
    """
    This function identifies what wall id and voxel height that is seen from a ground pixel

    INPUTS:
    dsm = Digital surface model
    azimuth and altitude = sun position in degrees
    scale= scale of DSM (1 meter pixels=1, 2 meter pixels=0.5)
    uniqueWallIDs = pixel row 'outside' buildings. will be calculated if empty
    walls = height of walls
    dem = Digital elevation model. (Should be excluded in future to incorporate ground elevation)

    OUTPUT:
    buildIDSeen = ID seen from ground pixel
    voxelHeight = Wall height shadow volume

    Fredrik Lindberg 2023-02-16
    fredrikl@gvc.gu.se

    """

    # Remove ground heights
    dsm = dsm - dem
    dsm[dsm < 0.5] = 0

    # conversion, degrees to radians
    azimuth = radians(azimuth)
    altitude = radians(altitude)

    # measure the size of the image
    rows = dsm.shape[0]
    cols = dsm.shape[1]

    # initialise parameters
    f = torch.clone(dsm)
    buildIDSeen = torch.zeros((rows, cols), device=device)

    dx = torch.tensor(0, device=device)
    dy = torch.tensor(0, device=device)
    dz = torch.tensor(0, device=device)
    temp = torch.zeros((rows, cols), device=device)
    temp2 = torch.zeros((rows, cols), device=device)  # walls
    tempwallID = torch.zeros((rows, cols), device=device)
    uniqueWallIDsOrig = torch.clone(uniqueWallIDs)

    voxelHeight = torch.zeros((rows, cols), device=device)
    temp3 = torch.ones((rows, cols), device=device)

    # create a fast PyTorch tensor lookup table for wall_dict
    max_wall_id = int(max(wall_dict.keys())) if wall_dict else 0
    wall_height_lookup = torch.zeros(max_wall_id + 1, device=device)
    for k, v in wall_dict.items():
        wall_height_lookup[int(k)] = v

    # other loop parameters
    amaxvalue = torch.max(dsm)
    pibyfour = torch.pi / 4
    threetimespibyfour = 3 * pibyfour
    fivetimespibyfour = 5 * pibyfour
    seventimespibyfour = 7 * pibyfour
    azimuth = torch.tensor(azimuth, device=device)
    sinazimuth = torch.sin(azimuth)
    cosazimuth = torch.cos(azimuth)
    tanazimuth = torch.tan(azimuth)
    signsinazimuth = torch.sign(sinazimuth)
    signcosazimuth = torch.sign(cosazimuth)
    dssin = torch.abs(1 / sinazimuth)
    dscos = torch.abs(1 / cosazimuth)
    altitude = torch.tensor(altitude, device=device)
    scale = torch.tensor(scale, device=device)

    tanaltitudebyscale = torch.tan(altitude) / scale

    index = 1

    # main loop
    while (amaxvalue >= dz) and (torch.abs(dx) < rows) and (torch.abs(dy) < cols):

        if (pibyfour <= azimuth and azimuth < threetimespibyfour) or (
            fivetimespibyfour <= azimuth and azimuth < seventimespibyfour
        ):
            dy = signsinazimuth * index
            dx = -1 * signcosazimuth * torch.abs(torch.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * torch.abs(torch.round(index * tanazimuth))
            dx = -1 * signcosazimuth * index
            ds = dscos

        dz = ds * index * tanaltitudebyscale
        temp[0:rows, 0:cols] = 0
        temp2[0:rows, 0:cols] = 0

        absdx = torch.abs(dx)
        absdy = torch.abs(dy)

        xc1 = int((dx + absdx) / 2)
        xc2 = int(rows + (dx - absdx) / 2)
        yc1 = int((dy + absdy) / 2)
        yc2 = int(cols + (dy - absdy) / 2)

        xp1 = int(-((dx - absdx) / 2))
        xp2 = int(rows - (dx + absdx) / 2)
        yp1 = int(-((dy - absdy) / 2))
        yp2 = int(cols - (dy + absdy) / 2)

        wallSeen = facesh
        uniqueWallIDs = uniqueWallIDs * wallSeen

        # Moving wall id
        tempwallID[xp1:xp2, yp1:yp2] = uniqueWallIDs[xc1:xc2, yc1:yc2]

        # FIX: Native PyTorch tensor lookup instead of np.vectorize
        temp_wallHeight = wall_height_lookup[tempwallID.long()]
        
        # Descending wall, how much of the wall that is still above ground level
        temp2 = temp_wallHeight - dz

        # buildIDSeen calculations
        buildIDSeen = (temp2 > 0) * temp3 * tempwallID + buildIDSeen

        # voxelHeight calculations
        voxelHeight = (temp2 > 0) * temp3 * (temp_wallHeight - temp2) + voxelHeight

        # Remember pixels previous iteration that walls have not progressed into yet.
        temp3 = torch.clone(temp2 <= 0) * (buildIDSeen == 0)

        index += 1

    # Ceil voxel height values to integers
    voxelHeight_ceil = torch.ceil(voxelHeight)

    # Empty raster to fill with voxel IDs
    voxelId = torch.zeros((rows, cols), device=device)
    
    # Ensure mapping references are native PyTorch tensors on the correct device
    wall2d_id = torch.tensor(wall2d_id, device=device)
    voxel_height = torch.tensor(voxel_height, device=device)
    voxelId_list = torch.tensor(voxelId_list, dtype=torch.long, device=device)

    # Flatten maps to find unique combinations
    a = buildIDSeen.flatten()
    b = voxelHeight_ceil.flatten()
    c = torch.column_stack([a, b])
    d = torch.unique(c, dim=0)
    d = d[~torch.all(d == 0, dim=1)]

    # Fill voxelId matrix with unique voxel IDs
    for temp_id, temp_height in d:
        # FIX: Safer element checking using .numel() and extraction via index [0]
        mask = (wall2d_id == temp_id) & (voxel_height == temp_height)
        temp_fill_id = voxelId_list[mask]
        
        pixel_mask = (buildIDSeen == temp_id) & (voxelHeight_ceil == temp_height)
        
        if temp_fill_id.numel() > 0:
            voxelId[pixel_mask] = temp_fill_id[0].float()
        else:
            buildIDSeen[pixel_mask] = 0
            voxelHeight_ceil[pixel_mask] = 0

    # Correct for shadows
    buildIDSeen = buildIDSeen * (1 - sh)
    voxelHeight = voxelHeight * (1 - sh)
    voxelId = voxelId * (1 - sh)

    return buildIDSeen, voxelHeight, voxelId