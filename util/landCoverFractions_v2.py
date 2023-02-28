# -*- coding: utf-8 -*-
#%calculates morphometric parameters for an image based on prevailing wind
#%direction. Specify a dem on a square grid to load and averaging dimension
#%
#%Date: 26 February 2004
#%Author:
#%   Offerle, B.
#%   Geovetarcentrum
#%   Goteborg University, Sweden
#%   Modified by Fredrik Lindberg 2010-01-09, fredrik.lindberg@kcl.ac.uk
#%   Translated to Python 20150108
#    Extended to be albe to calculate on irregular grids, Fredrik 20230207
#%--------------------------------------------------------------------------

import numpy as np
import scipy.ndimage.interpolation as sc
# import matplotlib.pylab as plt


def landcover_v2(lc_grid, mid, dtheta, feedback, imp_point):

    # Isotropic (this is the same as before. Works on irregular grids)
    lc_frac_all = np.zeros((1, 7))
    for i in range(0, 7):
        lc_gridvec = lc_grid[np.where(lc_grid == i + 1)]
        if lc_gridvec.size > 0:
            # lc_frac_all[0, i] = round((lc_gridvec.size * 1.0) / (lc_grid.size * 1.0), 3)
            lc_frac_all[0, i] = round((lc_gridvec.size * 1.0) / (lc_grid.size - (lc_grid == 0).sum()),3) # ignoring NoData (0) pixels

    # Anisotropic (Adjusted for irregular grids)
    lc_frac = np.zeros((int(360./dtheta), 7))
    deg = np.zeros((int(360./dtheta), 1))

    #n = lc_grid.shape[0]
    #imid = np.floor((n/2.))
    # if mid == 1:
        # dY = np.int16(np.arange(np.dot(1, imid)))  # the half length of the grid (y)
    #else: #moved inside loop as it varies on an irregular grid
        #dY = np.int16(np.arange(np.dot(1, n)))  # the whole length of the grid (y)

        # dX = np.int16(np.arange(imid, imid+1))
        #lx = dX.shape[0]
        #ly = dY.shape[0]

    j = int(0)
    for angle in np.arange(0, 360, dtheta):
        if imp_point == 1:
            feedback.setProgress(int(angle/3.6))

        #d = sc.rotate(lc_grid, angle, order=0, reshape=False, mode='nearest') #old
        d = sc.rotate(lc_grid, angle, order=0, reshape=True, mode='constant', cval=-99)

        n = d.shape[1]
        imid = np.floor((n/2.)) # the mid (NtoS) line of the grid
        if mid == 1: # from center point
            ny = d.shape[0]
            imidy = np.floor((ny/2.)) # the mid (NtoS) line of the grid
            lineMid = d[0:imidy,imid]
        else: #whole grid
            lineMid = d[:,int(imid)] # whole center line
        bld = lineMid[np.where(lineMid > 0)] # line within grid only  

        #b = np.round(((lc_grid.max()-lc_grid.min())/d.max())*d+lc_grid.min(), 0) #not needed anymore
        #bld = b[dY, dX]  # lc array
        ly = bld.shape[0] #number of pixels to consider in NtoS
        lx = 1 #!TODO should this consider full length (EtoW) of grid and if so, how?

        for i in range(0, 7):
            bldtemp = bld[np.where(bld == i + 1)]  # lc vector
            lc_frac[j, i] = np.float32(bldtemp.shape[0]) / (lx*ly)

        deg[j] = angle
        j += 1

    landcoverresult = {'lc_frac_all': lc_frac_all, 'lc_frac': lc_frac, 'deg': deg}

    return landcoverresult

