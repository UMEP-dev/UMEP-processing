# -*- coding: utf-8 -*-
'''
Calculates morphometric parameters for an image based on prevailing wind
direction. Specify a dem on a square grid to load and averaging dimension

Date: 26 February 2004
Author:
   Offerle, B.
   Geovetarcentrum
   Goteborg University, Sweden
   Modified by Fredrik Lindberg 2010-01-09, fredrik.lindberg@kcl.ac.uk
   Translated to Python 20150108
   Extended to be albe to calculate on irregular grids, Fredrik 20230208
--------------------------------------------------------------------------
dsm = Digital Surface Model (Ground and building heights)
dem = Digital Elevation Model (Ground heights)
mid = Start from center of domain (1) or calculate thruogh whole grid (0)
scale = 1/pixel resolution (m)
dtheta = 5.  # degree interval
feedback = USed to communicate with QGIS
imp_point = used to communicate with QGIS


'''
import numpy as np
import scipy.ndimage.interpolation as sc
# import matplotlib as plt


def imagemorphparam_v2(dsm, dem, scale, mid, dtheta, feedback, imp_point):

    numPixels = len(dsm[np.where(dsm != -9999)]) # too deal with irregular grids

    build = dsm - dem
    build[(build < 2.)] = 0.  # building should be higher than 2 meter

    # new part (when did i write this?)
    buildvec = build[np.where(build > 0)]
    if buildvec.size > 0:
        zH_all = buildvec.mean()
        zHmax_all = buildvec.max()
        zH_sd_all = buildvec.std()
        pai_all = (buildvec.size * 1.0) / numPixels #(build.size * 1.0)
    else:
        zH_all = 0.0
        zHmax_all = 0.0
        zH_sd_all = 0.0
        pai_all = 0.0

    fai = np.zeros((int(360./dtheta), 1))
    zH = np.zeros((int(360./dtheta), 1))
    zHmax = np.zeros((int(360./dtheta), 1))
    zH_sd = np.zeros((int(360./dtheta), 1))
    pai = np.zeros((int(360./dtheta), 1))
    deg = np.zeros((int(360./dtheta), 1))
    test = np.zeros((int(360./dtheta), 1))
    #%subset and center (moved inside loop 20230208)
    # n = dsm.shape[0]
    # imid = np.floor((n/2.))
    # if mid == 1:
    #     dY = np.int16(np.arange(np.dot(1, imid)))  # half the length of the grid (y)
    # else:
    #     dY = np.int16(np.arange(np.dot(1, n)))  # the whole length of the grid (y)
    # # if imp_point == 1:
    #     # total = 100 / (360. / dtheta)
    # dX = np.int16(np.arange(imid, imid+1))
    # lx = dX.shape[0]
    # ly = dY.shape[0]
    # filt1 = np.ones((n, 1)) * -1.
    # filt2 = np.ones((n, 1))
    # filt = np.array(np.hstack((filt1, filt2))).conj().T
    j = int(0)
    for angle in np.arange(0, 360, dtheta):
        if imp_point == 1:
            feedback.setProgress(int(angle/3.6))

        # Rotating buildings
        # d = sc.rotate(build, angle, order=0, reshape=False, mode='nearest') #old
        a = sc.rotate(build, angle, order=0, reshape=True, mode='constant', cval=-99)
       
        #% convolve leading edge filter with domain
        c = a * 0.0

        n = c.shape[1]
        imid = np.floor((n/2.))

        # filter for fai. Moved inside loop since size change if grid is irregular
        filt1 = np.ones((n, 1)) * -1.
        filt2 = np.ones((n, 1))
        filt = np.array(np.hstack((filt1, filt2))).conj().T
        buildZero = np.copy(a)
        buildZero[buildZero == -99] = 0 # remove -99 to avoid one 99 meter tall building wall

        for i in np.arange(1, c.shape[0]):
            c[int(i)-1, :] = np.sum((filt*buildZero[int(i)-1:i+1, :]), 0)

        if mid == 1: # from center point
            ny = a.shape[0]
            imidy = np.floor((ny/2.)) # the mid (NtoS) line of the grid
            lineMid = a[0:imidy,imid]
            walltemp = c[0:imidy,imid]
        else: #whole grid
            lineMid = a[:,int(imid)] # whole center line
            walltemp = c[:,int(imid)]

        bld = lineMid[np.where(lineMid > -99)]
        wall = walltemp[np.where(lineMid > -99)]
        ly = bld.shape[0] #number of pixels to consider in NtoS
        lx = 1 #!TODO should this consider full length (EtoW) of grid and if so, how?
        
        wall = wall[np.where(wall > 2)]  # wall vector
        fai[j] = np.sum(wall)/((lx*ly)/scale)
        bld = bld[np.where(bld > 2)]  # building vector: change from 0 to 2  : 20150906
        pai[j] = np.float32(bld.shape[0]) / (lx*ly)
        deg[j] = angle
        if np.float32(bld.shape[0]) == 0:
            zH[j] = 0
            zHmax[j] = 0
            zH_sd[j] = 0
        else:
            zH[j] = bld.sum() / np.float32(bld.shape[0])
            zHmax[j] = bld.max()
            zH_sd[j] = bld.std()

        # if angle == 0:
            # test = wall
        test[j] = ly

        j += 1

    fai_all = np.mean(fai)

    immorphresult = {'fai': fai, 'pai': pai, 'zH': zH, 'deg': deg, 'zHmax': zHmax,'zH_sd': zH_sd, 'pai_all': pai_all,
                        'zH_all': zH_all, 'zHmax_all': zHmax_all, 'zH_sd_all': zH_sd_all, 'fai_all': fai_all,'test': test}

    return immorphresult


# def imagemorphparam_v1(dsm, dem, scale, mid, dtheta, dlg, imp_point):

#     build = dsm - dem
#     test = build.max()
#     build[(build < 2.)] = 0.  # building should be higher than 2 meter
#     test = build.max()
#     # new part
#     buildvec = build[np.where(build > 0)]
#     if buildvec.size > 0:
#         zH_all = buildvec.mean()
#         zHmax_all = buildvec.max()
#         zH_sd_all = buildvec.std()
#         pai_all = (buildvec.size * 1.0) / (build.size * 1.0)
#     else:
#         zH_all = 0
#         zHmax_all = 0
#         zH_sd_all = 0
#         pai_all = 0 

#     fai = np.zeros((int(360./dtheta), 1))
#     zH = np.zeros((int(360./dtheta), 1))
#     zHmax = np.zeros((int(360./dtheta), 1))
#     zH_sd = np.zeros((int(360./dtheta), 1))
#     pai = np.zeros((int(360./dtheta), 1))
#     deg = np.zeros((int(360./dtheta), 1))

#     #%subset and center
#     n = dsm.shape[0]
#     imid = np.floor((n/2.))
#     if mid == 1:
#         dY = np.int16(np.arange(np.dot(1, imid)))  # half the length of the grid (y)
#     else:
#         dY = np.int16(np.arange(np.dot(1, n)))  # the whole length of the grid (y)
#     if imp_point == 1:
#         dlg.progressBar.setRange(0., 360. / dtheta)
#     dX = np.int16(np.arange(imid, imid+1))
#     lx = dX.shape[0]
#     ly = dY.shape[0]
#     filt1 = np.ones((n, 1)) * -1.
#     filt2 = np.ones((n, 1))
#     filt = np.array(np.hstack((filt1, filt2))).conj().T
#     j = int(0)
#     for angle in np.arange(0, (360.-dtheta+0) + dtheta, dtheta):
#         if imp_point == 1:
#             dlg.progressBar.setValue(angle)

#         c = np.zeros((n, n))

#         # Rotating building
#         # d = sc.imrotate(build, angle, 'nearest')
#         d = sc.rotate(build, angle, reshape=False, mode='nearest')
#         b = ((build.max()-build.min())/d.max())*d+build.min()
#         a = b
#         if b.sum() != 0:  # ground heights
#             # d = sc.imrotate(dsm, angle, 'nearest')
#             d = sc.rotate(dsm, angle, reshape=False, mode='nearest')
#             a = ((dsm.max()-dsm.min())/d.max())*d+dsm.min()

#         #% convolve leading edge filter with domain
#         for i in np.arange(1, (n-1)+1):
#             c[int(i)-1, :] = np.sum((filt*a[int(i)-1:i+1, :]), 0)

#         wall = c[dY, dX]  # wall array
#         wall = wall[np.where(wall > 2)]  # wall vector
#         fai[j] = np.sum(wall)/((lx*ly)/scale)
#         bld = b[dY, dX]  # building array
#         bld = bld[np.where(bld > 0)]  # building vector
#         pai[j] = np.float32(bld.shape[0]) / (lx*ly)
#         deg[j] = angle
#         if np.float32(bld.shape[0]) == 0:
#             zH[j] = 0
#             zHmax[j] = 0
#             zH_sd[j] = 0
#         else:
#             zH[j] = bld.sum() / np.float32(bld.shape[0])
#             zHmax[j] = bld.max()
#             zH_sd[j] = bld.std()

#         if angle == 0:
#             test = wall

#         j += 1

#     fai_all = np.mean(fai)

#     immorphresult = {'fai': fai, 'pai': pai, 'zH': zH, 'deg': deg, 'zHmax': zHmax,'zH_sd': zH_sd, 'pai_all': pai_all,
#                      'zH_all': zH_all, 'zHmax_all': zHmax_all, 'zH_sd_all': zH_sd_all, 'fai_all': fai_all,'test': test}

#     return immorphresult
