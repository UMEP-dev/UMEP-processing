from __future__ import absolute_import

import numpy as np
from .daylen import daylen
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import clearnessindex_2013b
from ...util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_13 import shadowingfunction_wallheight_13
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_23 import shadowingfunction_wallheight_23
from .gvf_2018a import gvf_2018a
from .cylindric_wedge import cylindric_wedge
from .TsWaveDelay_2015a import TsWaveDelay_2015a
from .Kup_veg_2015a import Kup_veg_2015a
# from .Lside_veg_v2015a import Lside_veg_v2015a
# from .Kside_veg_v2019a import Kside_veg_v2019a
from .Kside_veg_v2022a import Kside_veg_v2022a
from ...util.SEBESOLWEIGCommonFiles.Perez_v3 import Perez_v3
from ...util.SEBESOLWEIGCommonFiles.create_patches import create_patches

# Anisotropic longwave
from .Lcyl_v2022a import Lcyl_v2022a
from .Lside_veg_v2022a import Lside_veg_v2022a
from .anisotropic_sky import anisotropic_sky as ani_sky
from .patch_radiation import patch_steradians
from copy import deepcopy
import time

# Wall surface temperature scheme
from .wall_surface_temperature import wall_surface_temperature

def Solweig_2025a_calc(i, dsm, scale, rows, cols, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg, svfEveg, svfSveg,
                       svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, vegdem, vegdem2, albedo_b, absK, absL,
                       ewall, Fside, Fup, Fcyl, altitude, azimuth, zen, jday, usevegdem, onlyglobal, buildings, location, psi,
                       landcover, lc_grid, dectime, altmax, dirwalls, walls, cyl, elvis, Ta, RH, radG, radD, radI, P,
                       amaxvalue, bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
                       TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timestepdec, Tgmap1, 
                       Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI, TgOut1, diffsh, shmat, vegshmat, vbshvegshmat, anisotropic_sky, asvf, patch_option,
                       voxelMaps, voxelTable, ws, wallScheme, timeStep, steradians, walls_scheme, dirwalls_scheme):

#def Solweig_2021a_calc(i, dsm, scale, rows, cols, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg, svfEveg, svfSveg,
#                       svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, vegdem, vegdem2, albedo_b, absK, absL,
#                       ewall, Fside, Fup, Fcyl, altitude, azimuth, zen, jday, usevegdem, onlyglobal, buildings, location, psi,
#                       landcover, lc_grid, dectime, altmax, dirwalls, walls, cyl, elvis, Ta, RH, radG, radD, radI, P,
#                       amaxvalue, bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
#                       TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timestepdec, Tgmap1, 
#                       Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI, TgOut1, diffsh, ani):

    # This is the core function of the SOLWEIG model
    # 2016-Aug-28
    # Fredrik Lindberg, fredrikl@gvc.gu.se
    # Goteborg Urban Climate Group
    # Gothenburg University
    #
    # Input variables:
    # dsm = digital surface model
    # scale = height to pixel size (2m pixel gives scale = 0.5)
    # svf,svfN,svfW,svfE,svfS = SVFs for building and ground
    # svfveg,svfNveg,svfEveg,svfSveg,svfWveg = Veg SVFs blocking sky
    # svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg = Veg SVFs blocking buildings
    # vegdem = Vegetation canopy DSM
    # vegdem2 = Vegetation trunk zone DSM
    # albedo_b = building wall albedo
    # absK = human absorption coefficient for shortwave radiation
    # absL = human absorption coefficient for longwave radiation
    # ewall = Emissivity of building walls
    # Fside = The angular factors between a person and the surrounding surfaces
    # Fup = The angular factors between a person and the surrounding surfaces
    # Fcyl = The angular factors between a culidric person and the surrounding surfaces
    # altitude = Sun altitude (degree)
    # azimuth = Sun azimuth (degree)
    # zen = Sun zenith angle (radians)
    # jday = day of year
    # usevegdem = use vegetation scheme
    # onlyglobal = calculate dir and diff from global shortwave (Reindl et al. 1990)
    # buildings = Boolena grid to identify building pixels
    # location = geographic location
    # height = height of measurements point (center of gravity of human)
    # psi = 1 - Transmissivity of shortwave through vegetation
    # landcover = use landcover scheme !!!NEW IN 2015a!!!
    # lc_grid = grid with landcoverclasses
    # lc_class = table with landcover properties
    # dectime = decimal time
    # altmax = maximum sun altitude
    # dirwalls = aspect of walls
    # walls = one pixel row outside building footprint. height of building walls
    # cyl = consider man as cylinder instead of cude
    # elvis = dummy
    # Ta = air temp
    # RH
    # radG = global radiation
    # radD = diffuse
    # radI = direct
    # P = pressure
    # amaxvalue = max height of buildings
    # bush = grid representing bushes
    # Twater = temperature of water (daily)
    # TgK, Tstart, TgK_wall, Tstart_wall, TmaxLST,TmaxLST_wall, 
    # alb_grid, emis_grid = albedo and emmissivity on ground
    # first, second = conneted to old Ts model (source area based on Smidt et al.)
    # svfalfa = SVF recalculated to angle
    # svfbuveg = complete SVF 
    # firstdaytime, timeadd, timestepdec, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, 
    # CI = Clearness index
    # TgOut1 = old Ts model
    # diffsh, ani = Used in anisotrpic models (Wallenberg et al. 2019, 2022)

    # # # Core program start # # #
    # Instrument offset in degrees
    t = 0.

    # Stefan Bolzmans Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi/180

    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location['latitude'])

    # Vapor pressure
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH / 100.)

    # Determination of clear - sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (1 - (1 + msteg) * np.exp(-((1.2 + 3.0 * msteg) ** 0.5))) + elvis  # -0.04 old error from Jonsson et al.2006

    if altitude > 0: # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        #  factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt, I0et, CIuncorr = clearnessindex_2013b(zen, jday, Ta, RH / 100., radG, location, P)
        if (CI > 1) or (CI == np.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if onlyglobal == 1:
            I0, CI, Kt, I0et, CIuncorr = clearnessindex_2013b(zen, jday, Ta, RH / 100., radG, location, P)
            if (CI > 1) or (CI == np.inf):
                CI = 1

            radI, radD = diffusefraction(radG, altitude, Kt, Ta, RH)

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        if anisotropic_sky == 1:
            patchchoice = 1
            zenDeg = zen*(180/np.pi)
            # Relative luminance
            lv, pc_, pb_ = Perez_v3(zenDeg, azimuth, radD, radI, jday, patchchoice, patch_option)   
            # Total relative luminance from sky, i.e. from each patch, into each cell
            aniLum = np.zeros((rows, cols))
            for idx in range(lv.shape[0]):
                aniLum += diffsh[:,:,idx] * lv[idx,2]     

            dRad = aniLum * radD   # Total diffuse radiation from sky into each cell
        else:
            dRad = radD * svfbuveg
            patchchoice = 1
            lv = None

        # Shadow  images
        if usevegdem == 1:
            vegsh, sh, _, wallsh, wallsun, wallshve, _, facesun, wallsh_ = shadowingfunction_wallheight_23(dsm, vegdem, vegdem2,
                                        azimuth, altitude, scale, amaxvalue, bush, walls, dirwalls * np.pi / 180., walls_scheme, dirwalls_scheme * np.pi/180.)
            shadow = sh - (1 - vegsh) * (1 - psi)
        else:
            sh, wallsh, wallsun, facesh, facesun, wallsh_ = shadowingfunction_wallheight_13(dsm, azimuth, altitude, scale,
                                                                                   walls, dirwalls * np.pi / 180., walls_scheme, dirwalls_scheme * np.pi/180.)
            shadow = sh

        # # # Surface temperature parameterisation during daytime # # # #
        # new using max sun alt.instead of  dfm
        # Tgamp = (TgK * altmax - Tstart) + Tstart # Old
        Tgamp = TgK * altmax + Tstart # Fixed 2021
        # Tgampwall = (TgK_wall * altmax - (Tstart_wall)) + (Tstart_wall) # Old
        Tgampwall = TgK_wall * altmax + Tstart_wall
        Tg = Tgamp * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST / 24 - SNUP / 24)) * np.pi / 2) # 2015 a, based on max sun altitude
        Tgwall = Tgampwall * np.sin((((dectime - np.floor(dectime)) - SNUP / 24) / (TmaxLST_wall / 24 - SNUP / 24)) * np.pi / 2) # 2015a, based on max sun altitude

        if Tgwall < 0:  # temporary for removing low Tg during morning 20130205
            # Tg = 0
            Tgwall = 0

        # New estimation of Tg reduction for non - clear situation based on Reindl et al.1990
        radI0, _ = diffusefraction(I0, altitude, 1., Ta, RH)
        corr = 0.1473 * np.log(90 - (zen / np.pi * 180)) + 0.3454  # 20070329 correction of lat, Lindberg et al. 2008
        CI_Tg = (radG / radI0) + (1 - corr)
        if (CI_Tg > 1) or (CI_Tg == np.inf):
            CI_Tg = 1

        radG0 = radI0 * (np.sin(altitude * deg2rad)) + _
        CI_TgG = (radG / radG0) + (1 - corr)
        if (CI_TgG > 1) or (CI_TgG == np.inf):
            CI_TgG = 1
        
        # Tg = Tg * CI_Tg  # new estimation
        # Tgwall = Tgwall * CI_Tg
        Tg = Tg * CI_TgG  # new estimation
        Tgwall = Tgwall * CI_TgG
        if landcover == 1:
            Tg[Tg < 0] = 0  # temporary for removing low Tg during morning 20130205

        # # # # Ground View Factors # # # #
        gvfLup, gvfalb, gvfalbnosh, gvfLupE, gvfalbE, gvfalbnoshE, gvfLupS, gvfalbS, gvfalbnoshS, gvfLupW, gvfalbW,\
        gvfalbnoshW, gvfLupN, gvfalbN, gvfalbnoshN, gvfSum, gvfNorm = gvf_2018a(wallsun, walls, buildings, scale, shadow, first,
                second, dirwalls, Tg, Tgwall, Ta, emis_grid, ewall, alb_grid, SBC, albedo_b, rows, cols,
                                                                 Twater, lc_grid, landcover)

        # # # # Lup, daytime # # # #
        # Surface temperature wave delay - new as from 2014a
        Lup, timeaddnotused, Tgmap1 = TsWaveDelay_2015a(gvfLup, firstdaytime, timeadd, timestepdec, Tgmap1)
        LupE, timeaddnotused, Tgmap1E = TsWaveDelay_2015a(gvfLupE, firstdaytime, timeadd, timestepdec, Tgmap1E)
        LupS, timeaddnotused, Tgmap1S = TsWaveDelay_2015a(gvfLupS, firstdaytime, timeadd, timestepdec, Tgmap1S)
        LupW, timeaddnotused, Tgmap1W = TsWaveDelay_2015a(gvfLupW, firstdaytime, timeadd, timestepdec, Tgmap1W)
        LupN, timeaddnotused, Tgmap1N = TsWaveDelay_2015a(gvfLupN, firstdaytime, timeadd, timestepdec, Tgmap1N)
        
        # # For Tg output in POIs
        TgTemp = Tg * shadow + Ta
        TgOut, timeadd, TgOut1 = TsWaveDelay_2015a(TgTemp, firstdaytime, timeadd, timestepdec, TgOut1) #timeadd only here v2021a

        # Building height angle from svf
        F_sh = cylindric_wedge(zen, svfalfa, rows, cols)  # Fraction shadow on building walls based on sun alt and svf
        F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * shadow * np.sin(altitude * (np.pi / 180)) + dRad + albedo_b * (1 - svfbuveg) * \
                            (radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))

        Kup, KupE, KupS, KupW, KupN = Kup_veg_2015a(radI, radD, radG, altitude, svfbuveg, albedo_b, F_sh, gvfalb,
                    gvfalbE, gvfalbS, gvfalbW, gvfalbN, gvfalbnosh, gvfalbnoshE, gvfalbnoshS, gvfalbnoshW, gvfalbnoshN)

        Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside = Kside_veg_v2022a(radI, radD, radG, shadow, svfS, svfW, svfN, svfE,
                    svfEveg, svfSveg, svfWveg, svfNveg, azimuth, altitude, psi, t, albedo_b, F_sh, KupE, KupS, KupW,
                    KupN, cyl, lv, anisotropic_sky, diffsh, rows, cols, asvf, shmat, vegshmat, vbshvegshmat)
        
        firstdaytime = 0

    else:  # # # # # # # NIGHTTIME # # # # # # # #

        Tgwall = 0
        # CI_Tg = -999  # F_sh = []

        # Nocturnal K fluxes set to 0
        Knight = np.zeros((rows, cols))
        Kdown = np.zeros((rows, cols))
        Kwest = np.zeros((rows, cols))
        Kup = np.zeros((rows, cols))
        Keast = np.zeros((rows, cols))
        Ksouth = np.zeros((rows, cols))
        Knorth = np.zeros((rows, cols))
        KsideI = np.zeros((rows, cols))
        KsideD = np.zeros((rows, cols))
        F_sh = np.zeros((rows, cols))
        Tg = np.zeros((rows, cols))
        shadow = np.zeros((rows, cols))
        CI_Tg = deepcopy(CI)
        CI_TgG = deepcopy(CI)
        dRad = np.zeros((rows,cols))
        Kside = np.zeros((rows,cols))

        # # # # Lup # # # #
        Lup = SBC * emis_grid * ((Knight + Ta + Tg + 273.15) ** 4)
        if landcover == 1:
            Lup[lc_grid == 3] = SBC * 0.98 * (Twater + 273.15) ** 4  # nocturnal Water temp

        LupE = Lup
        LupS = Lup
        LupW = Lup
        LupN = Lup

        # # For Tg output in POIs
        TgOut = Ta + Tg

        I0 = 0
        timeadd = 0
        firstdaytime = 1

    # # # # Ldown # # # #
    Ldown = (svf + svfveg - 1) * esky * SBC * ((Ta + 273.15) ** 4) + (2 - svfveg - svfaveg) * ewall * SBC * \
                ((Ta + 273.15) ** 4) + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + \
                (2 - svf - svfveg) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)
    # Ldown = Ldown - 25 # Shown by Jonsson et al.(2006) and Duarte et al.(2006)

    if CI < 0.95:  # non - clear conditions
        c = 1 - CI
        Ldown = Ldown * (1 - c) + c * ((svf + svfveg - 1) * SBC * ((Ta + 273.15) ** 4) + (2 - svfveg - svfaveg) *
                ewall * SBC * ((Ta + 273.15) ** 4) + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) +
                (2 - svf - svfveg) * (1 - ewall) * SBC * ((Ta + 273.15) ** 4))  # NOT REALLY TESTED!!! BUT MORE CORRECT?

    # # # # Lside # # # #
    Least, Lsouth, Lwest, Lnorth = Lside_veg_v2022a(svfS, svfW, svfN, svfE, svfEveg, svfSveg, svfWveg, svfNveg,
                    svfEaveg, svfSaveg, svfWaveg, svfNaveg, azimuth, altitude, Ta, Tgwall, SBC, ewall, Ldown,
                                                      esky, t, F_sh, CI, LupE, LupS, LupW, LupN, anisotropic_sky)

    # New parameterization scheme for wall temperatures
    if wallScheme == 1:
        # albedo_g = 0.15 #TODO Change to correct
        if altitude < 0:
            wallsh_ = 0
        voxelTable = wall_surface_temperature(voxelTable, wallsh_, altitude, azimuth, timeStep, radI, radD, radG, Ldown, Lup, Ta, esky)
    # Anisotropic sky
    if (anisotropic_sky == 1):
        if 'lv' not in locals():
            # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
            skyvaultalt, skyvaultazi, _, _, _, _, _ = create_patches(patch_option)

            patch_emissivities = np.zeros(skyvaultalt.shape[0])

            x = np.transpose(np.atleast_2d(skyvaultalt))
            y = np.transpose(np.atleast_2d(skyvaultazi))
            z = np.transpose(np.atleast_2d(patch_emissivities))

            L_patches = np.append(np.append(x, y, axis=1), z, axis=1)

        else:
            L_patches = deepcopy(lv)

        # Calculate steradians for patches if it is the first model iteration
        if i == 0:
            steradians, skyalt, patch_altitude = patch_steradians(L_patches)

        # Create lv from L_patches if nighttime, i.e. lv does not exist
        if altitude < 0:
            # CI = deepcopy(CI)
            lv = deepcopy(L_patches); KupE = 0; KupS = 0; KupW = 0; KupN = 0

        # Adjust sky emissivity under semi-cloudy/hazy/cloudy/overcast conditions, i.e. CI lower than 0.95
        if CI < 0.95:
            esky_c = CI * esky + (1 - CI) * 1.
            esky = esky_c

        Ldown, Lside, Lside_sky, Lside_veg, Lside_sh, Lside_sun, Lside_ref, Least_, Lwest_, Lnorth_, Lsouth_, \
           Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside, steradians, skyalt = ani_sky(shmat, vegshmat, vbshvegshmat, altitude, azimuth, asvf, cyl, esky,
                                                                                     L_patches, wallScheme, voxelTable, voxelMaps, steradians, Ta, Tgwall, ewall, Lup, radI, radD, radG, lv, 
                                                                                     albedo_b, 0, diffsh, shadow, KupE, KupS, KupW, KupN, i)
    else:
        Lside = np.zeros((rows, cols))
        L_patches = None

    # Box and anisotropic longwave
    if cyl == 0 and anisotropic_sky == 1:
        Least += Least_
        Lwest += Lwest_
        Lnorth += Lnorth_
        Lsouth += Lsouth_

    # # # # Calculation of radiant flux density and Tmrt # # # #
    # Human body considered as a cylinder with isotropic all-sky diffuse
    if cyl == 1 and anisotropic_sky == 0: 
        Sstr = absK * (KsideI * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
                        ((Ldown + Lup) * Fup + (Lnorth + Least + Lsouth + Lwest) * Fside)
    # Human body considered as a cylinder with Perez et al. (1993) (anisotropic sky diffuse) 
    # and Martin and Berdahl (1984) (anisotropic sky longwave)
    elif cyl == 1 and anisotropic_sky == 1:
        Sstr = absK * (Kside * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
                        ((Ldown + Lup) * Fup + Lside * Fcyl + (Lnorth + Least + Lsouth + Lwest) * Fside)
    # Knorth = nan Ksouth = nan Kwest = nan Keast = nan
    else: # Human body considered as a standing cube
        Sstr = absK * ((Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
                        ((Ldown + Lup) * Fup + (Lnorth + Least + Lsouth + Lwest) * Fside)

    Tmrt = np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.2

    # Add longwave to cardinal directions for output in POI
    if (cyl == 1) and (anisotropic_sky == 1):
        Least += Least_
        Lwest += Lwest_
        Lnorth += Lnorth_
        Lsouth += Lsouth_

    return Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, \
           timeadd, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, Keast, Ksouth, Kwest, Knorth, Least, \
           Lsouth, Lwest, Lnorth, KsideI, TgOut1, TgOut, radI, radD, \
               Lside, L_patches, CI_Tg, CI_TgG, KsideD, dRad, Kside, steradians, voxelTable
