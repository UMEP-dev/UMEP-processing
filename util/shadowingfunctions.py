# -*- coding: utf-8 -*-
# Ready for python action!
import numpy as np
from math import radians


def shadowingfunctionglobalradiation(
    a, azimuth, altitude, scale, feedback, forsvf
):

    # %This m.file calculates shadows on a DEM
    # % conversion
    degrees = np.pi / 180.0
    # if azimuth == 0.0:
    # azimuth = 0.000000000001
    azimuth = np.dot(azimuth, degrees)
    altitude = np.dot(altitude, degrees)
    # % measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]
    if forsvf == 0:
        barstep = np.max([sizex, sizey])
        total = 100.0 / barstep  # dlg.progressBar.setRange(0, barstep)
    # % initialise parameters
    f = a
    dx = 0.0
    dy = 0.0
    dz = 0.0
    temp = np.zeros((sizex, sizey))
    index = 1.0
    # % other loop parameters
    amaxvalue = a.max()
    pibyfour = np.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1.0 / sinazimuth))
    dscos = np.abs((1.0 / cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    # % main loop
    while amaxvalue >= dz and np.abs(dx) < sizex and np.abs(dy) < sizey:
        if forsvf == 0:
            feedback.setProgress(int(index * total))

        if (
            pibyfour <= azimuth
            and azimuth < threetimespibyfour
            or fivetimespibyfour <= azimuth
            and azimuth < seventimespibyfour
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos

        # % note: dx and dy represent absolute values while ds is an incremental value
        dz = ds * index * tanaltitudebyscale
        temp[0:sizex, 0:sizey] = 0.0
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = (dx + absdx) / 2.0 + 1.0
        xc2 = sizex + (dx - absdx) / 2.0
        yc1 = (dy + absdy) / 2.0 + 1.0
        yc2 = sizey + (dy - absdy) / 2.0
        xp1 = -((dx - absdx) / 2.0) + 1.0
        xp2 = sizex - (dx + absdx) / 2.0
        yp1 = -((dy - absdy) / 2.0) + 1.0
        yp2 = sizey - (dy + absdy) / 2.0
        temp[int(xp1) - 1 : int(xp2), int(yp1) - 1 : int(yp2)] = (
            a[int(xc1) - 1 : int(xc2), int(yc1) - 1 : int(yc2)] - dz
        )
        f = np.fmax(f, temp)
        index += 1.0

    f = f - a
    f = np.logical_not(f)
    sh = np.double(f)

    return sh


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
):

    # This function casts shadows on buildings and vegetation units.
    # New capability to deal with pergolas 20210827

    # conversion
    degrees = np.pi / 180.0
    azimuth = azimuth * degrees
    altitude = altitude * degrees

    # measure the size of grid
    sizex = a.shape[0]
    sizey = a.shape[1]

    # progressbar for svf plugin
    if forsvf == 0:
        barstep = np.max([sizex, sizey])
        total = 100.0 / barstep
        feedback.setProgress(0)

    # initialise parameters
    dx = 0.0
    dy = 0.0
    dz = 0.0
    temp = np.zeros((sizex, sizey))
    tempvegdem = np.zeros((sizex, sizey))
    tempvegdem2 = np.zeros((sizex, sizey))
    templastfabovea = np.zeros((sizex, sizey))
    templastgabovea = np.zeros((sizex, sizey))
    bushplant = bush > 1.0
    sh = np.zeros((sizex, sizey))  # shadows from buildings
    vbshvegsh = np.zeros((sizex, sizey))  # vegetation blocking buildings
    vegsh = np.add(
        np.zeros((sizex, sizey)), bushplant, dtype=float
    )  # vegetation shadow
    f = a

    pibyfour = np.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1.0 / sinazimuth))
    dscos = np.abs((1.0 / cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    # index = 1
    index = 0

    # new case with pergola (thin vertical layer of vegetation), August 2021
    dzprev = 0

    # main loop
    while (amaxvalue >= dz) and (np.abs(dx) < sizex) and (np.abs(dy) < sizey):
        if forsvf == 0:
            # dlg.progressBar.setValue(index)
            feedback.setProgress(int(index * total))
        if (
            (pibyfour <= azimuth)
            and (azimuth < threetimespibyfour)
            or (fivetimespibyfour <= azimuth)
            and (azimuth < seventimespibyfour)
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos
        # note: dx and dy represent absolute values while ds is an incremental
        # value
        dz = (ds * index) * tanaltitudebyscale
        tempvegdem[0:sizex, 0:sizey] = 0.0
        tempvegdem2[0:sizex, 0:sizey] = 0.0
        temp[0:sizex, 0:sizey] = 0.0
        templastfabovea[0:sizex, 0:sizey] = 0.0
        templastgabovea[0:sizex, 0:sizey] = 0.0
        absdx = np.abs(dx)
        absdy = np.abs(dy)
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

        f = np.fmax(f, temp)  # Moving building shadow
        sh[(f > a)] = 1.0
        sh[(f <= a)] = 0.0
        fabovea = tempvegdem > a  # vegdem above DEM
        gabovea = tempvegdem2 > a  # vegdem2 above DEM

        # new pergola condition
        templastfabovea[xp1:xp2, yp1:yp2] = vegdem[xc1:xc2, yc1:yc2] - dzprev
        templastgabovea[xp1:xp2, yp1:yp2] = vegdem2[xc1:xc2, yc1:yc2] - dzprev
        lastfabovea = templastfabovea > a
        lastgabovea = templastgabovea > a
        dzprev = dz
        vegsh2 = np.add(
            np.add(
                np.add(fabovea, gabovea, dtype=float), lastfabovea, dtype=float
            ),
            lastgabovea,
            dtype=float,
        )
        vegsh2[vegsh2 == 4] = 0.0
        # vegsh2[vegsh2 == 1] = 0. # This one is the ultimate question...
        vegsh2[vegsh2 > 0] = 1.0

        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[(vegsh * sh > 0.0)] = 0.0
        vbshvegsh = vegsh + vbshvegsh  # removing shadows 'behind' buildings

        index += 1.0

    sh = 1.0 - sh
    vbshvegsh[(vbshvegsh > 0.0)] = 1.0
    vbshvegsh = vbshvegsh - vegsh
    vegsh = 1.0 - vegsh
    vbshvegsh = 1.0 - vbshvegsh

    shadowresult = {"sh": sh, "vegsh": vegsh, "vbshvegsh": vbshvegsh}

    return shadowresult


def shadowingfunction_20_old(
    a, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, dlg, forsvf
):

    # % This function casts shadows on buildings and vegetation units
    # % conversion
    degrees = np.pi / 180.0
    if azimuth == 0.0:
        azimuth = 0.000000000001
    azimuth = np.dot(azimuth, degrees)
    altitude = np.dot(altitude, degrees)
    # % measure the size of the image
    sizex = a.shape[0]
    sizey = a.shape[1]
    # % initialise parameters
    if forsvf == 0:
        barstep = np.max([sizex, sizey])
        dlg.progressBar.setRange(0, barstep)
        dlg.progressBar.setValue(0)

    dx = 0.0
    dy = 0.0
    dz = 0.0
    temp = np.zeros((sizex, sizey))
    tempvegdem = np.zeros((sizex, sizey))
    tempvegdem2 = np.zeros((sizex, sizey))
    sh = np.zeros((sizex, sizey))
    vbshvegsh = np.zeros((sizex, sizey))
    vegsh = np.zeros((sizex, sizey))
    tempbush = np.zeros((sizex, sizey))
    f = a
    g = np.zeros((sizex, sizey))
    bushplant = bush > 1.0
    pibyfour = np.pi / 4.0
    threetimespibyfour = 3.0 * pibyfour
    fivetimespibyfour = 5.0 * pibyfour
    seventimespibyfour = 7.0 * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs((1.0 / sinazimuth))
    dscos = np.abs((1.0 / cosazimuth))
    tanaltitudebyscale = np.tan(altitude) / scale
    index = 1

    # % main loop
    while amaxvalue >= dz and np.abs(dx) < sizex and np.abs(dy) < sizey:
        if forsvf == 0:
            dlg.progressBar.setValue(index)
        if (
            pibyfour <= azimuth
            and azimuth < threetimespibyfour
            or fivetimespibyfour <= azimuth
            and azimuth < seventimespibyfour
        ):
            dy = signsinazimuth * index
            dx = -1.0 * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1.0 * signcosazimuth * index
            ds = dscos
        # % note: dx and dy represent absolute values while ds is an incremental value
        dz = np.dot(np.dot(ds, index), tanaltitudebyscale)
        tempvegdem[0:sizex, 0:sizey] = 0.0
        tempvegdem2[0:sizex, 0:sizey] = 0.0
        temp[0:sizex, 0:sizey] = 0.0
        absdx = np.abs(dx)
        absdy = np.abs(dy)
        xc1 = (dx + absdx) / 2.0 + 1.0
        xc2 = sizex + (dx - absdx) / 2.0
        yc1 = (dy + absdy) / 2.0 + 1.0
        yc2 = sizey + (dy - absdy) / 2.0
        xp1 = -((dx - absdx) / 2.0) + 1.0
        xp2 = sizex - (dx + absdx) / 2.0
        yp1 = -((dy - absdy) / 2.0) + 1.0
        yp2 = sizey - (dy + absdy) / 2.0
        tempvegdem[int(xp1) - 1 : int(xp2), int(yp1) - 1 : int(yp2)] = (
            vegdem[int(xc1) - 1 : int(xc2), int(yc1) - 1 : int(yc2)] - dz
        )
        tempvegdem2[int(xp1) - 1 : int(xp2), int(yp1) - 1 : int(yp2)] = (
            vegdem2[int(xc1) - 1 : int(xc2), int(yc1) - 1 : int(yc2)] - dz
        )
        temp[int(xp1) - 1 : int(xp2), int(yp1) - 1 : int(yp2)] = (
            a[int(xc1) - 1 : int(xc2), int(yc1) - 1 : int(yc2)] - dz
        )
        # f = np.maximum(f, temp) # bad performance in python3. Replaced with
        # fmax
        f = np.fmax(f, temp)
        sh[(f > a)] = 1.0
        sh[(f <= a)] = 0.0
        # %Moving building shadow
        fabovea = tempvegdem > a
        # %vegdem above DEM
        gabovea = tempvegdem2 > a
        # %vegdem2 above DEM
        # vegsh2 = np.float(fabovea)-np.float(gabovea)
        vegsh2 = np.subtract(fabovea, gabovea, dtype=float)
        # vegsh = np.maximum(vegsh, vegsh2) # bad performance in python3.
        # Replaced with fmax
        vegsh = np.fmax(vegsh, vegsh2)
        vegsh[(vegsh * sh > 0.0)] = 0.0
        # % removing shadows 'behind' buildings
        vbshvegsh = vegsh + vbshvegsh
        # % vegsh at high sun altitudes
        if index == 1.0:
            firstvegdem = tempvegdem - temp
            firstvegdem[(firstvegdem <= 0.0)] = 1000.0
            vegsh[(firstvegdem < dz)] = 1.0
            vegsh = vegsh * (vegdem2 > a)
            vbshvegsh = np.zeros((sizex, sizey))

        # % Bush shadow on bush plant
        if np.logical_and(bush.max() > 0.0, np.max((fabovea * bush)) > 0.0):
            tempbush[0:sizex, 0:sizey] = 0.0
            tempbush[int(xp1) - 1 : int(xp2), int(yp1) - 1 : int(yp2)] = (
                bush[int(xc1) - 1 : int(xc2), int(yc1) - 1 : int(yc2)] - dz
            )
            # g = np.maximum(g, tempbush) # bad performance in python3.
            # Replaced with fmax
            g = np.fmax(g, tempbush)
            g *= bushplant
        index += 1.0

    sh = 1.0 - sh
    vbshvegsh[(vbshvegsh > 0.0)] = 1.0
    vbshvegsh = vbshvegsh - vegsh

    if bush.max() > 0.0:
        g = g - bush
        g[(g > 0.0)] = 1.0
        g[(g < 0.0)] = 0.0
        vegsh = vegsh - bushplant + g
        vegsh[(vegsh < 0.0)] = 0.0

    vegsh[(vegsh > 0.0)] = 1.0
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
    # buildings = 1 - ((dsm) > 0)
    dsm[dsm < 0.5] = 0

    # conversion, degrees to radians
    azimuth = radians(azimuth)
    altitude = radians(altitude)

    # measure the size of the image
    rows = np.shape(dsm)[0]
    cols = np.shape(dsm)[1]

    # initialise parameters
    f = np.copy(dsm)
    buildIDSeen = np.zeros((rows, cols))
    # h = np.zeros((rows, cols))

    dx = 0
    dy = 0
    dz = 0
    temp = np.zeros((rows, cols))
    temp2 = np.zeros((rows, cols))  # walls
    tempwallID = np.zeros((rows, cols))
    uniqueWallIDsOrig = np.copy(uniqueWallIDs)

    voxelHeight = np.zeros((rows, cols))
    temp3 = np.ones((rows, cols))

    # other loop parameters
    amaxvalue = np.max(dsm)
    pibyfour = np.pi / 4
    threetimespibyfour = 3 * pibyfour
    fivetimespibyfour = 5 * pibyfour
    seventimespibyfour = 7 * pibyfour
    sinazimuth = np.sin(azimuth)
    cosazimuth = np.cos(azimuth)
    tanazimuth = np.tan(azimuth)
    signsinazimuth = np.sign(sinazimuth)
    signcosazimuth = np.sign(cosazimuth)
    dssin = np.abs(1 / sinazimuth)
    dscos = np.abs(1 / cosazimuth)
    tanaltitudebyscale = np.tan(altitude) / scale

    index = 1

    # main loop
    while (amaxvalue >= dz) and (np.abs(dx) < rows) and (np.abs(dy) < cols):

        if (pibyfour <= azimuth and azimuth < threetimespibyfour) or (
            fivetimespibyfour <= azimuth and azimuth < seventimespibyfour
        ):
            dy = signsinazimuth * index
            dx = -1 * signcosazimuth * np.abs(np.round(index / tanazimuth))
            ds = dssin
        else:
            dy = signsinazimuth * np.abs(np.round(index * tanazimuth))
            dx = -1 * signcosazimuth * index
            ds = dscos

        # note: dx and dy represent absolute values while ds is an incremental
        # value
        dz = ds * index * tanaltitudebyscale
        temp[0:rows, 0:cols] = 0
        temp2[0:rows, 0:cols] = 0

        absdx = np.abs(dx)
        absdy = np.abs(dy)

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
        # uniqueWallIDs = ((uniqueWallIDs - firstMove) < 0) * uniqueWallIDsOrig + uniqueWallIDs # adding missing corner
        # wallSeenHeight = walls * wallSeen

        # temp2[xp1:xp2, yp1:yp2] = wallSeenHeight[xc1:xc2, yc1:yc2] - dz # Moving wall shadow
        # Moving wall id
        tempwallID[xp1:xp2, yp1:yp2] = uniqueWallIDs[xc1:xc2, yc1:yc2]

        # Get wall height from wall id
        temp_wallHeight = np.vectorize(wall_dict.__getitem__)(tempwallID)

        # Descending wall, how much of the wall that is still above ground
        # level
        temp2 = temp_wallHeight - dz

        # buildIDSeen = Wall pixels/voxels seen, i.e. only voxels that are positive (above ground level) (temp2 > 0).
        # temp3 indicates those pixels that the walls have not progressed into
        # yet (saved in previous iteration).
        buildIDSeen = (temp2 > 0) * temp3 * tempwallID + buildIDSeen

        # voxelHeight = the elevation on a wall that is seen from a pixel with the given altitude and azimuth (only above ground leve, i.e. (temp2 > 0)).
        # voxelHeight = wall height - descending wall, i.e. temp_wallHeight -
        # temp2. Only applicable to pixels where there is no value from
        # previous iterations (temp3).
        voxelHeight = (temp2 > 0) * temp3 * (
            temp_wallHeight - temp2
        ) + voxelHeight

        # Remember pixels previous iteration that walls have not progressed
        # into yet.
        temp3 = np.copy(temp2 <= 0) * (buildIDSeen == 0)

        index += 1

    # Ceil voxel height values to integers
    voxelHeight_ceil = np.ceil(voxelHeight)
    # voxelHeight_ceil = np.round(voxelHeight)

    # Empty raster to fill with voxel IDs
    voxelId = np.zeros((rows, cols))
    # Convert wall2d_id from list to numpy array
    wall2d_id = np.array(wall2d_id)
    # Convert voxel_height from list to numpy array
    voxel_height = np.array(voxel_height)
    # Convert voxelId_list from list to numpy array
    voxelId_list = np.array(voxelId_list, dtype=int)

    # Flatten buildIDseen from matrix to array
    a = buildIDSeen.flatten()
    # Flatten voxelHeight_ceil from matrix to array
    b = voxelHeight_ceil.flatten()
    # Combine the two above arrays into an n by 2 array
    c = np.column_stack([a, b])
    # Find unique values in c
    d = np.unique(c, axis=0)
    # Remove rows where both columns are zero
    d = d[~np.all(d == 0, axis=1)]

    not_in_list = 0
    in_list = 0
    # Fill voxelId matrix with unique voxel IDs
    for temp_id, temp_height in d:
        # print(str(temp_id) + ' ' + str(temp_height))
        temp_fill_id = voxelId_list[
            ((wall2d_id == temp_id) & (voxel_height == temp_height))
        ]
        if temp_fill_id.__len__() > 0:
            # print('temp_fill_id = ' + str(temp_fill_id))
            voxelId[
                (buildIDSeen == temp_id) & (voxelHeight_ceil == temp_height)
            ] = temp_fill_id
            in_list += 1
        else:
            not_in_list += 1
            buildIDSeen[
                (buildIDSeen == temp_id) & (voxelHeight_ceil == temp_height)
            ] = 0
            voxelHeight_ceil[
                (buildIDSeen == temp_id) & (voxelHeight_ceil == temp_height)
            ] = 0

    # Correct for shadows, i.e. remove weird pixels on top of buildings etc
    buildIDSeen = buildIDSeen * (1 - sh)
    voxelHeight = voxelHeight * (1 - sh)
    voxelId = voxelId * (1 - sh)

    return buildIDSeen, voxelHeight, voxelId
