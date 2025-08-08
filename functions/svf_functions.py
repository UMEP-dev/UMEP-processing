import numpy as np
from ..util import shadowingfunctions as shadow
from ..util.SEBESOLWEIGCommonFiles.create_patches import create_patches
# from ..functions.wallalgorithms import findwalls
from ..functions import wallalgorithms as wa
from ..functions import svf_for_voxels as svfv
from ..util.SEBESOLWEIGCommonFiles import shadowingfunction_wallheight_13 as shb
from ..util.SEBESOLWEIGCommonFiles import shadowingfunction_wallheight_23 as shbv
# remove
from ..util.misc import saveraster
from osgeo import gdal, osr
from osgeo.gdalconst import *

def annulus_weight(altitude, aziinterval):
    n = 90.
    steprad = (360./aziinterval) * (np.pi/180.)
    annulus = 91.-altitude
    w = (1./(2.*np.pi)) * np.sin(np.pi / (2.*n)) * np.sin((np.pi * (2. * annulus - 1.)) / (2. * n))
    weight = steprad * w

    return weight

def svf_angles_100121():
    azi1 = np.arange(1., 360., 360./16.)  #%22.5
    azi2 = np.arange(12., 360., 360./16.)  #%22.5
    azi3 = np.arange(5., 360., 360./32.)  #%11.25
    azi4 = np.arange(2., 360., 360./32.)  #%11.25
    azi5 = np.arange(4., 360., 360./40.)  #%9
    azi6 = np.arange(7., 360., 360./48.)  #%7.50
    azi7 = np.arange(6., 360., 360./48.)  #%7.50
    azi8 = np.arange(1., 360., 360./48.)  #%7.50
    azi9 = np.arange(4., 359., 360./52.)  #%6.9231
    azi10 = np.arange(5., 360., 360./52.)  #%6.9231
    azi11 = np.arange(1., 360., 360./48.)  #%7.50
    azi12 = np.arange(0., 359., 360./44.)  #%8.1818
    azi13 = np.arange(3., 360., 360./44.)  #%8.1818
    azi14 = np.arange(2., 360., 360./40.)  #%9
    azi15 = np.arange(7., 360., 360./32.)  #%10
    azi16 = np.arange(3., 360., 360./24.)  #%11.25
    azi17 = np.arange(10., 360., 360./16.)  #%15
    azi18 = np.arange(19., 360., 360./12.)  #%22.5
    azi19 = np.arange(17., 360., 360./8.)  #%45
    azi20 = 0.  #%360
    iazimuth = np.array(np.hstack((azi1, azi2, azi3, azi4, azi5, azi6, azi7, azi8, azi9, azi10, azi11, azi12, azi13,
                                    azi14, azi15, azi16, azi17, azi18, azi19, azi20)))
    aziinterval = np.array(np.hstack((16., 16., 32., 32., 40., 48., 48., 48., 52., 52., 48., 44., 44., 40., 32., 24.,
                                        16., 12., 8., 1.)))
    angleresult = {'iazimuth': iazimuth, 'aziinterval': aziinterval}

    return angleresult


def svfForProcessing153(dsm, vegdem, vegdem2, scale, usevegdem, pixel_resolution, wallScheme, demlayer, feedback):
    rows = dsm.shape[0]
    cols = dsm.shape[1]
    svf = np.zeros([rows, cols])
    svfE = np.zeros([rows, cols])
    svfS = np.zeros([rows, cols])
    svfW = np.zeros([rows, cols])
    svfN = np.zeros([rows, cols])
    svfveg = np.zeros((rows, cols))
    svfEveg = np.zeros((rows, cols))
    svfSveg = np.zeros((rows, cols))
    svfWveg = np.zeros((rows, cols))
    svfNveg = np.zeros((rows, cols))
    svfaveg = np.zeros((rows, cols))
    svfEaveg = np.zeros((rows, cols))
    svfSaveg = np.zeros((rows, cols))
    svfWaveg = np.zeros((rows, cols))
    svfNaveg = np.zeros((rows, cols))

    # % amaxvalue
    vegmax = vegdem.max()
    amaxvalue = dsm.max()
    amaxvalue = np.maximum(amaxvalue, vegmax)

    # % Elevation vegdems if buildingDSM inclused ground heights
    vegdem = vegdem + dsm
    vegdem[vegdem == dsm] = 0
    vegdem2 = vegdem2 + dsm
    vegdem2[vegdem2 == dsm] = 0
    # % Bush separation
    bush = np.logical_not((vegdem2 * vegdem)) * vegdem

    #index = int(0)

    # patch_option = 1 # 145 patches
    patch_option = 2 # 153 patches
    # patch_option = 3 # 306 patches
    # patch_option = 4 # 612 patches
    
    # Create patches based on patch_option
    skyvaultalt, skyvaultazi, annulino, skyvaultaltint, aziinterval, skyvaultaziint, azistart = create_patches(patch_option)

    skyvaultaziint = np.array([360/patches for patches in aziinterval])
    iazimuth = np.hstack(np.zeros((1, np.sum(aziinterval)))) # Nils

    shmat = np.zeros((rows, cols, np.sum(aziinterval)))
    vegshmat = np.zeros((rows, cols, np.sum(aziinterval)))
    vbshvegshmat = np.zeros((rows, cols, np.sum(aziinterval)))

    # Preparations for wall temperature scheme
    if wallScheme:
        feedback.setProgressText("Estimating view factors for wall voxels")
        voxelTable, voxelId_list, wall_dict, walls, aspect, uniqueWallIDs, wall2d_id, voxel_height = svfv.wallscheme_prepare(dsm, scale, pixel_resolution, feedback)

        # Rasters to fill with values in loop
        all_buildIDSeen = np.zeros((rows, cols, skyvaultalt.shape[0]))
        all_voxelHeight = np.zeros((rows, cols, skyvaultalt.shape[0]))
        all_voxelId = np.zeros((rows, cols, skyvaultalt.shape[0]))
    else:
        voxelTable = 0
        allbuildIDSeen = 0
        allvoxelHeight = 0
        all_voxelId = 0
        walls = 0

    index = 0
    for j in range(0, skyvaultaltint.shape[0]):
        for k in range(0, int(360 / skyvaultaziint[j])):
            iazimuth[index] = k * skyvaultaziint[j] + azistart[j]
            if iazimuth[index] > 360.:
                iazimuth[index] = iazimuth[index] - 360.
            index = index + 1
    aziintervalaniso = np.ceil(aziinterval / 2.0)
    index = int(0)
    for i in range(0, skyvaultaltint.shape[0]):
        for j in np.arange(0, (aziinterval[int(i)])):
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break
            altitude = skyvaultaltint[int(i)]
            azimuth = iazimuth[int(index)]

            # Casting shadow
            if wallScheme:
                if usevegdem == 1:
                    vegsh, sh, vbshvegsh, wallsh, wallsun, wallshve, facesh, facesun = shbv.shadowingfunction_wallheight_23(dsm, vegdem, vegdem2, azimuth, altitude, scale, amaxvalue, bush, walls, aspect * np.pi / 180)
                    vegshmat[:, :, index] = vegsh
                    vbshvegshmat[:, :, index] = vbshvegsh
                else:
                    sh, wallsh, wallsun, facesh, facesun = shb.shadowingfunction_wallheight_13(dsm, azimuth, altitude, scale, walls, aspect * np.pi / 180.)
                    vegsh = np.ones((sh.shape[0], sh.shape[1])).astype(float)
                    vbshvegsh = np.ones((sh.shape[0], sh.shape[1])).astype(float)
                    vegshmat[:, :, index] = vegsh
                    vbshvegshmat[:, :, index] = vbshvegsh
            else:
                if usevegdem == 1:
                    shadowresult = shadow.shadowingfunction_20(dsm, vegdem, vegdem2, azimuth, altitude,
                                                            scale, amaxvalue, bush, feedback, 1)
                    vegsh = shadowresult["vegsh"]
                    vbshvegsh = shadowresult["vbshvegsh"]
                    vegshmat[:, :, index] = vegsh
                    vbshvegshmat[:, :, index] = vbshvegsh
                    sh = shadowresult["sh"]
                else:
                    sh = shadow.shadowingfunctionglobalradiation(dsm, azimuth, altitude, scale, feedback, 1)

            shmat[:, :, index] = sh
            
            # Wall temperature scheme, i.e. finding out which voxel is seen from each pixel, where direction is patch azimuth and altitude
            if wallScheme:
                all_buildIDSeen[:,:, index], all_voxelHeight[:,:, index], all_voxelId[:,:, index] = shadow.shadowingfunction_findwallID(dsm, azimuth, altitude, scale, walls, uniqueWallIDs, demlayer, wall2d_id, voxel_height, voxelId_list, facesh, wall_dict, sh)

            # Calculate svfs
            for k in np.arange(annulino[int(i)]+1, (annulino[int(i+1.)])+1):
                weight = annulus_weight(k, aziinterval[i])*sh
                svf = svf + weight
                weight = annulus_weight(k, aziintervalaniso[i]) * sh
                if (azimuth >= 0) and (azimuth < 180):
                    svfE = svfE + weight
                if (azimuth >= 90) and (azimuth < 270):
                    svfS = svfS + weight
                if (azimuth >= 180) and (azimuth < 360):
                    svfW = svfW + weight
                if (azimuth >= 270) or (azimuth < 90):
                    svfN = svfN + weight

            if usevegdem == 1:
                for k in np.arange(annulino[int(i)] + 1, (annulino[int(i + 1.)]) + 1):
                    # % changed to include 90
                    weight = annulus_weight(k, aziinterval[i])
                    svfveg = svfveg + weight * vegsh
                    svfaveg = svfaveg + weight * vbshvegsh
                    weight = annulus_weight(k, aziintervalaniso[i])
                    if (azimuth >= 0) and (azimuth < 180):
                        svfEveg = svfEveg + weight * vegsh
                        svfEaveg = svfEaveg + weight * vbshvegsh
                    if (azimuth >= 90) and (azimuth < 270):
                        svfSveg = svfSveg + weight * vegsh
                        svfSaveg = svfSaveg + weight * vbshvegsh
                    if (azimuth >= 180) and (azimuth < 360):
                        svfWveg = svfWveg + weight * vegsh
                        svfWaveg = svfWaveg + weight * vbshvegsh
                    if (azimuth >= 270) or (azimuth < 90):
                        svfNveg = svfNveg + weight * vegsh
                        svfNaveg = svfNaveg + weight * vbshvegsh

            index += 1
            feedback.setProgress(int(index * (100. / np.sum(aziinterval))))

    svfS = svfS + 3.0459e-004
    svfW = svfW + 3.0459e-004
    # % Last azimuth is 90. Hence, manual add of last annuli for svfS and SVFW
    # %Forcing svf not be greater than 1 (some MATLAB crazyness)
    svf[(svf > 1.)] = 1.
    svfE[(svfE > 1.)] = 1.
    svfS[(svfS > 1.)] = 1.
    svfW[(svfW > 1.)] = 1.
    svfN[(svfN > 1.)] = 1.

    if usevegdem == 1:
        last = np.zeros((rows, cols))
        last[(vegdem2 == 0.)] = 3.0459e-004
        svfSveg = svfSveg + last
        svfWveg = svfWveg + last
        svfSaveg = svfSaveg + last
        svfWaveg = svfWaveg + last
        # %Forcing svf not be greater than 1 (some MATLAB crazyness)
        svfveg[(svfveg > 1.)] = 1.
        svfEveg[(svfEveg > 1.)] = 1.
        svfSveg[(svfSveg > 1.)] = 1.
        svfWveg[(svfWveg > 1.)] = 1.
        svfNveg[(svfNveg > 1.)] = 1.
        svfaveg[(svfaveg > 1.)] = 1.
        svfEaveg[(svfEaveg > 1.)] = 1.
        svfSaveg[(svfSaveg > 1.)] = 1.
        svfWaveg[(svfWaveg > 1.)] = 1.
        svfNaveg[(svfNaveg > 1.)] = 1.

    svfresult = {'svf': svf, 'svfE': svfE, 'svfS': svfS, 'svfW': svfW, 'svfN': svfN,
                    'svfveg': svfveg, 'svfEveg': svfEveg, 'svfSveg': svfSveg, 'svfWveg': svfWveg,
                    'svfNveg': svfNveg, 'svfaveg': svfaveg, 'svfEaveg': svfEaveg, 'svfSaveg': svfSaveg,
                    'svfWaveg': svfWaveg, 'svfNaveg': svfNaveg, 'shmat': shmat, 'vegshmat': vegshmat, 'vbshvegshmat': vbshvegshmat,
                    'voxelIds': all_voxelId, 'voxelTable': voxelTable, 'walls': walls}
                    # ,
                    # 'vbshvegshmat': vbshvegshmat, 'wallshmat': wallshmat, 'wallsunmat': wallsunmat,
                    # 'wallshvemat': wallshvemat, 'facesunmat': facesunmat}
    return svfresult


def svfForProcessing655(dsm, vegdem, vegdem2, scale, usevegdem, feedback):
    rows = dsm.shape[0]
    cols = dsm.shape[1]
    svf = np.zeros([rows, cols])
    svfE = np.zeros([rows, cols])
    svfS = np.zeros([rows, cols])
    svfW = np.zeros([rows, cols])
    svfN = np.zeros([rows, cols])
    svfveg = np.zeros((rows, cols))
    svfEveg = np.zeros((rows, cols))
    svfSveg = np.zeros((rows, cols))
    svfWveg = np.zeros((rows, cols))
    svfNveg = np.zeros((rows, cols))
    svfaveg = np.zeros((rows, cols))
    svfEaveg = np.zeros((rows, cols))
    svfSaveg = np.zeros((rows, cols))
    svfWaveg = np.zeros((rows, cols))
    svfNaveg = np.zeros((rows, cols))

    # % amaxvalue
    vegmax = vegdem.max()
    amaxvalue = dsm.max()
    amaxvalue = np.maximum(amaxvalue, vegmax)

    # % Elevation vegdems if buildingDSM inclused ground heights
    vegdem = vegdem + dsm
    vegdem[vegdem == dsm] = 0
    vegdem2 = vegdem2 + dsm
    vegdem2[vegdem2 == dsm] = 0
    # % Bush separation
    bush = np.logical_not((vegdem2 * vegdem)) * vegdem

    # shmat = np.zeros((rows, cols, 145))
    # vegshmat = np.zeros((rows, cols, 145))

    noa = 19.
    #% No. of anglesteps minus 1
    step = 89./noa
    iangle = np.array(np.hstack((np.arange(step/2., 89., step), 90.)))
    annulino = np.array(np.hstack((np.round(np.arange(0., 89., step)), 90.)))
    angleresult = svf_angles_100121()
    aziinterval = angleresult["aziinterval"]
    iazimuth = angleresult["iazimuth"]
    aziintervalaniso = np.ceil((aziinterval/2.))
    index = 1.

    for i in np.arange(0, iangle.shape[0]-1):
        for j in np.arange(0, (aziinterval[int(i)])):
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break
            altitude = iangle[int(i)]
            azimuth = iazimuth[int(index)-1]

            # Casting shadow
            if usevegdem == 1:
                shadowresult = shadow.shadowingfunction_20(dsm, vegdem, vegdem2, azimuth, altitude,
                                                            scale, amaxvalue, bush, feedback, 1)
                vegsh = shadowresult["vegsh"]
                vbshvegsh = shadowresult["vbshvegsh"]
                sh = shadowresult["sh"]
            else:
                sh = shadow.shadowingfunctionglobalradiation(dsm, azimuth, altitude, scale, feedback, 1)

            # Calculate svfs
            for k in np.arange(annulino[int(i)]+1, (annulino[int(i+1.)])+1):
                weight = annulus_weight(k, aziinterval[i])*sh
                svf = svf + weight
                weight = annulus_weight(k, aziintervalaniso[i]) * sh
                if (azimuth >= 0) and (azimuth < 180):
                    svfE = svfE + weight
                if (azimuth >= 90) and (azimuth < 270):
                    svfS = svfS + weight
                if (azimuth >= 180) and (azimuth < 360):
                    svfW = svfW + weight
                if (azimuth >= 270) or (azimuth < 90):
                    svfN = svfN + weight

            if usevegdem == 1:
                for k in np.arange(annulino[int(i)] + 1, (annulino[int(i + 1.)]) + 1):
                    # % changed to include 90
                    weight = annulus_weight(k, aziinterval[i])
                    svfveg = svfveg + weight * vegsh
                    svfaveg = svfaveg + weight * vbshvegsh
                    weight = annulus_weight(k, aziintervalaniso[i])
                    if (azimuth >= 0) and (azimuth < 180):
                        svfEveg = svfEveg + weight * vegsh
                        svfEaveg = svfEaveg + weight * vbshvegsh
                    if (azimuth >= 90) and (azimuth < 270):
                        svfSveg = svfSveg + weight * vegsh
                        svfSaveg = svfSaveg + weight * vbshvegsh
                    if (azimuth >= 180) and (azimuth < 360):
                        svfWveg = svfWveg + weight * vegsh
                        svfWaveg = svfWaveg + weight * vbshvegsh
                    if (azimuth >= 270) or (azimuth < 90):
                        svfNveg = svfNveg + weight * vegsh
                        svfNaveg = svfNaveg + weight * vbshvegsh

            index += 1
            feedback.setProgress(int(index * (100. / 655.)))

    svfS = svfS + 3.0459e-004
    svfW = svfW + 3.0459e-004
    # % Last azimuth is 90. Hence, manual add of last annuli for svfS and SVFW
    # %Forcing svf not be greater than 1 (some MATLAB crazyness)
    svf[(svf > 1.)] = 1.
    svfE[(svfE > 1.)] = 1.
    svfS[(svfS > 1.)] = 1.
    svfW[(svfW > 1.)] = 1.
    svfN[(svfN > 1.)] = 1.

    if usevegdem == 1:
        last = np.zeros((rows, cols))
        last[(vegdem2 == 0.)] = 3.0459e-004
        svfSveg = svfSveg + last
        svfWveg = svfWveg + last
        svfSaveg = svfSaveg + last
        svfWaveg = svfWaveg + last
        # %Forcing svf not be greater than 1 (some MATLAB crazyness)
        svfveg[(svfveg > 1.)] = 1.
        svfEveg[(svfEveg > 1.)] = 1.
        svfSveg[(svfSveg > 1.)] = 1.
        svfWveg[(svfWveg > 1.)] = 1.
        svfNveg[(svfNveg > 1.)] = 1.
        svfaveg[(svfaveg > 1.)] = 1.
        svfEaveg[(svfEaveg > 1.)] = 1.
        svfSaveg[(svfSaveg > 1.)] = 1.
        svfWaveg[(svfWaveg > 1.)] = 1.
        svfNaveg[(svfNaveg > 1.)] = 1.

    svfresult = {'svf': svf, 'svfE': svfE, 'svfS': svfS, 'svfW': svfW, 'svfN': svfN,
                    'svfveg': svfveg, 'svfEveg': svfEveg, 'svfSveg': svfSveg, 'svfWveg': svfWveg,
                    'svfNveg': svfNveg, 'svfaveg': svfaveg, 'svfEaveg': svfEaveg, 'svfSaveg': svfSaveg,
                    'svfWaveg': svfWaveg, 'svfNaveg': svfNaveg}

    return svfresult
