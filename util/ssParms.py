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
from ..functions import wallalgorithms as wa

# import matplotlib as plt


def ss_calc(dsm, dem, cdsm, walls, numPixels, feedback):

    build = dsm - dem
    build[(build < 1.)] = 0.  # building should be higher than 1? meter

    #noOfPixels = int(dsm.shape[0] * dsm.shape[1])
    walllimit = 0.3 # 30 centimeters height variation identifies a vegetation edge pixel
    total = 100. / (int(dsm.shape[0] * dsm.shape[1]))

    if cdsm.max() > 0:
        vegEdges = wa.findwalls(cdsm, walllimit, feedback, total)
 
    buildvec = build[np.where(build > 0)]
    if buildvec.size > 0:
        zH_all = buildvec.mean()
        zHmax_all = buildvec.max()
        zH_sd_all = buildvec.std()
        pai_ground = (buildvec.size * 1.0) / numPixels
        iterHeights = int(np.floor(zHmax_all))
    else:
        zH_all = 0
        zHmax_all = 0
        zH_sd_all = 0
        pai_all = 0
        iterHeights = int(0)

    z = np.zeros((iterHeights, 1))
    paiZ_b = np.zeros((iterHeights, 1))
    waiZ_b = np.zeros((iterHeights, 1))
    paiZ_v = np.zeros((iterHeights, 1))
    waiZ_v = np.zeros((iterHeights, 1))
    
    for i in np.arange(0, iterHeights): 
        z[i] = i
        buildZ = build - i
        wallsZ = walls - i
        paiZ_b[i] = np.where(buildZ > 0)[0].shape[0] / numPixels
        waiZ_b[i] = np.where(wallsZ > 0)[0].shape[0] / numPixels

        if cdsm.max() > 0:
            vegZ = cdsm - i
            vegedgeZ = vegEdges - i
            paiZ_v[i] = np.where(vegZ > 0)[0].shape[0] / numPixels
            waiZ_v[i] = np.where(vegedgeZ > 0)[0].shape[0] / numPixels
        else:
            paiZ_v[i] = 0
            waiZ_v[i] = 0

    ssResult = {'z': z, 'paiZ_b': paiZ_b, 'waiZ_b': waiZ_b,'paiZ_v': paiZ_v, 'waiZ_v': waiZ_v}

    return ssResult


