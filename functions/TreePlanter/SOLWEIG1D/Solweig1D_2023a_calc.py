from __future__ import absolute_import
import numpy as np
from ...SOLWEIGpython.daylen import daylen
from ....util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import clearnessindex_2013b
from ....util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
from ...SOLWEIGpython.cylindric_wedge import cylindric_wedge
# from .TsWaveDelay_2015a import TsWaveDelay_2015a
# from .Kup_veg_2015a import Kup_veg_2015a
from ...SOLWEIGpython.Lside_veg_v2015a import Lside_veg_v2015a
from ..SOLWEIG1D.Kside1D_veg_v2019a import Kside_veg_v2019a
#from ...SOLWEIGpython.Perez_v3_moved import Perez_v3
from ....util.SEBESOLWEIGCommonFiles.Perez_v3 import Perez_v3

from .anisotropic_sky import anisotropic_sky as ani_sky
from copy import deepcopy

def Solweig1D_2019a_calc(svf, svfveg, svfaveg, sh, vegsh,  albedo_b, absK, absL, ewall, Fside, Fup, Fcyl, altitude, azimuth, zen, jday,
                         onlyglobal, location, dectime, altmax, cyl, elvis, Ta, RH, radG, radD, radI, P,
                         Twater, TgK, Tstart, albedo_g, eground, TgK_wall, Tstart_wall, TmaxLST, TmaxLST_wall,
                         svfalfa, CI, ani, diffsh, trans, patch_option, skyp, buip, vegp, asvf, L_patches, steradians):

    # This is the core function of the SOLWEIG1D model, 2019-Jun-21
    # Fredrik Lindberg, fredrikl@gvc.gu.se, Goteborg Urban Climate Group, Gothenburg University, Sweden

    svfE = svf
    svfW = svf
    svfN = svf
    svfS = svf
    svfEveg = svfveg
    svfSveg = svfveg
    svfWveg = svfveg
    svfNveg = svfveg
    svfEaveg = svfaveg
    svfSaveg = svfaveg
    svfWaveg = svfaveg
    svfNaveg = svfaveg
    psi = trans

    # Instrument offset in degrees
    t = 0.

    # Stefan Bolzmans Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi/180

    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location['latitude'])

    shadow = sh - (1 - vegsh) * (1 - psi)

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
        if ani == 1:
            patchchoice = 1
            zenDeg = zen*(180/np.pi)
            lv, _, _ = Perez_v3(zenDeg, azimuth, radD, radI, jday, patchchoice, patch_option)   # Relative luminance

            aniLum = 0.
            for idx in range(skyp.shape[0]):
                aniLum = aniLum + skyp[idx] * lv[idx, 2]     # Total relative luminance from sky into each cell

            dRad = aniLum * radD   # Total diffuse radiation from sky into each cell

        else:
            dRad = radD * svf
            lv = 0

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
        #if landcover == 1:
        #    Tg[Tg < 0] = 0  # temporary for removing low Tg during morning 20130205

        # gvf = shadow

        # Lup = SBC * eground * ((gvf + Ta + Tg + 273.15) ** 4)
        # Ground always in shade, i.e. ground longwave radiation based on ground emissivity and air temperature
        Lup = SBC * eground * ((Ta + 273.15) ** 4)
        LupE = Lup * 0.5
        LupS = Lup * 0.5
        LupW = Lup * 0.5
        LupN = Lup * 0.5

        # Building height angle from svf
        F_sh = cylindric_wedge(zen, svfalfa, 1, 1)  # Fraction shadow on building walls based on sun alt and svf
        F_sh[np.isnan(F_sh)] = 0.5

        # # # # # # # Calculation of shortwave daytime radiative fluxes # # # # # # #
        Kdown = radI * shadow * np.sin(altitude * (np.pi / 180)) + dRad + albedo_b * (1 - svf) * \
                            (radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))

        #Kdown = radI * shadow * np.sin(altitude * (np.pi / 180)) + radD * svfbuveg + albedo_b * (1 - svfbuveg) * \
                            #(radG * (1 - F_sh) + radD * F_sh) # *sin(altitude(i) * (pi / 180))

        Kup = albedo_g * (shadow * radI * np.sin(altitude * (np.pi / 180.))) + radD * svf + \
              albedo_b * (1 - svf) * (radG * (1 - F_sh) + radD * F_sh)

        # Kup, KupE, KupS, KupW, KupN = Kup_veg_2015a(radI, radD, radG, altitude, svf, albedo_b, F_sh, gvfalb,
        #             gvfalbE, gvfalbS, gvfalbW, gvfalbN, gvfalbnosh, gvfalbnoshE, gvfalbnoshS, gvfalbnoshW, gvfalbnoshN)

        # Keast, Ksouth, Kwest, Knorth, KsideI, KsideD = Kside_veg_v2019a(radI, radD, radG, shadow, svfS, svfW, svfN, svfE,
        #             svfEveg, svfSveg, svfWveg, svfNveg, azimuth, altitude, psi, t, albedo_b, F_sh, Kup, Kup, Kup,
        #             Kup, cyl, lv, ani, diffsh, 1, 1)

        # firstdaytime = 0

    else:  # # # # # # # NIGHTTIME # # # # # # # #

        Tgwall = 0
        # CI_Tg = -999  # F_sh = []

        # Nocturnal K fluxes set to 0
        Knight = 0.
        Kdown = 0.
        Kwest = 0.
        Kup = 0.
        Keast = 0.
        Ksouth = 0.
        Knorth = 0.
        KsideI = 0.
        KsideD = 0.
        F_sh = 0.
        Tg = 0.
        shadow = 0.

        # # # # Lup # # # #
        Lup = SBC * eground * ((Knight + Ta + Tg + 273.15) ** 4)
        # if landcover == 1:
        #     Lup[lc_grid == 3] = SBC * 0.98 * (Twater + 273.15) ** 4  # nocturnal Water temp
        LupE = Lup * 0.5
        LupS = Lup * 0.5
        LupW = Lup * 0.5
        LupN = Lup * 0.5

        # # For Tg output in POIs
        # TgOut = Ta + Tg

        I0 = 0
        # timeadd = 0
        # firstdaytime = 1

    # # # # Ldown # # # #
    Ldown = svf * esky * SBC * ((Ta + 273.15) ** 4) + (1 - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) + \
            (1 - svf) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4)  # Jonsson et al.(2006)

    if CI < 0.95:  # non - clear conditions
        c = 1 - CI
        Ldown = Ldown * (1 - c) + \
                c * (svf * SBC * ((Ta + 273.15) ** 4) + (1 - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4) +
                     (1 - svf) * (1 - ewall) * esky * SBC * ((Ta + 273.15) ** 4))

    if ani == 0:
        # # # # Lside # # # #
        Least, Lsouth, Lwest, Lnorth = Lside_veg_v2015a(svfS, svfW, svfN, svfE, svfEveg, svfSveg, svfWveg, svfNveg,
                        svfEaveg, svfSaveg, svfWaveg, svfNaveg, azimuth, altitude, Ta, Tgwall, SBC, ewall, Ldown,
                                                        esky, t, F_sh, CI, LupE, LupS, LupW, LupN)

    else:

        if CI < 0.95:
            esky_c = CI * esky + (1 - CI) * 1.
            esky = esky_c

        L_patches[:, 2] = 0.

        Ldown, Lside, Lside_sky, Lside_veg, Lside_sh, Lside_sun, Lside_ref, Least_, Lwest_, Lnorth_, Lsouth_, \
            Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside, steradians, skyalt = ani_sky(skyp, buip, vegp, altitude, azimuth, asvf, cyl, esky,
                                                                                        L_patches, False, False, False, steradians, Ta, Tgwall, ewall, Lup, radI, radD, radG, lv, 
                                                                                        albedo_b, 0, diffsh, shadow, Kup)
        
        Least = 0
        Lwest = 0
        Lsouth = 0
        Lnorth = 0

    # Box and anisotropic longwave
    if cyl == 0 and ani == 1:
        Least += Least_
        Lwest += Lwest_
        Lnorth += Lnorth_
        Lsouth += Lsouth_

    Least += LupE
    Lwest += LupW
    Lnorth += LupN
    Lsouth += LupS

    # # # # Calculation of radiant flux density and Tmrt # # # #
    # if cyl == 1 and ani == 1:  # Human body considered as a cylinder with Perez et al. (1993)
    #     Sstr = absK * ((KsideI + KsideD) * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
    #                     (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)
    if cyl == 1 and ani == 1:
        Sstr = absK * (Kside * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
                        ((Ldown + Lup) * Fup + Lside * Fcyl + (Lnorth + Least + Lsouth + Lwest) * Fside)

    elif cyl == 1 and ani == 0: # Human body considered as a cylinder with isotropic all-sky diffuse
        Sstr = absK * (KsideI * Fcyl + (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) + absL * \
                        (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)

    else: # Human body considered as a standing cube
        Sstr = absK * ((Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside) +absL * \
                        (Ldown * Fup + Lup * Fup + Lnorth * Fside + Least * Fside + Lsouth * Fside + Lwest * Fside)

    Tmrt = np.sqrt(np.sqrt((Sstr / (absL * SBC)))) - 273.2

    return Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, Keast, Ksouth, Kwest, Knorth, Least, Lsouth, Lwest, \
           Lnorth, KsideI, radI, radD, shadow
