"""
@author Fredrik Lindberg, University of Gothenburg
"""

from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
from .daylen import daylen
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import (
    clearnessindex_2013b,
)
from ...util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_13 import (
    shadowingfunction_wallheight_13,
)
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_23 import (
    shadowingfunction_wallheight_23,
)
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
from .Lside_veg import Lside_veg_v2022a, Lside_veg_v2026
from .anisotropic_sky import anisotropic_sky as ani_sky
from .patch_radiation import patch_steradians
from copy import deepcopy
import time
import torch

# Wall surface temperature scheme
from .wall_surface_temperature import wall_surface_temperature

# Ground surface temperature
from .ground_surface import surfaceTemperature_calc, outgoingLongwave_calc


def Solweig_2026a_calc(
    i,
    dsm,
    scale,
    rows,
    cols,
    svf,
    svfN,
    svfW,
    svfE,
    svfS,
    svfveg,
    svfNveg,
    svfEveg,
    svfSveg,
    svfWveg,
    svfaveg,
    svfEaveg,
    svfSaveg,
    svfWaveg,
    svfNaveg,
    vegdem,
    vegdem2,
    albedo_b,
    absK,
    absL,
    ewall,
    Fside,
    Fup,
    Fcyl,
    altitude,
    azimuth,
    zen,
    jday,
    usevegdem,
    onlyglobal,
    buildings,
    location,
    psi,
    landcover,
    lc_grid,
    dectime,
    altmax,
    dirwalls,
    walls,
    cyl,
    elvis,
    Ta,
    RH,
    radG,
    radD,
    radI,
    P,
    amaxvalue,
    bush,
    Twater,
    TgK,
    Tstart,
    alb_grid,
    emis_grid,
    TgK_wall,
    Tstart_wall,
    TmaxLST,
    TmaxLST_wall,
    first,
    second,
    svfalfa,
    svfbuveg,
    firstdaytime,
    timeadd,
    timestepdec,
    Tgmap1,
    Tgmap1E,
    Tgmap1S,
    Tgmap1W,
    Tgmap1N,
    CI,
    diffsh,
    shmat,
    vegshmat,
    vbshvegshmat,
    anisotropic_sky,
    asvf,
    patch_option,
    voxelMaps,
    voxelTable,
    ws,
    wallScheme,
    timeStep,
    steradians,
    walls_scheme,
    dirwalls_scheme,
    groundScheme,
    Tg,
    Rn,
    Rn_past,
    G,
    Tm,
    cap_grid,
    diff_grid,
    a1_grid,
    a2_grid,
    a3_grid,
    shadow_past,
):
    """
    This is the core function of the SOLWEIG model
    2016-Aug-28
    Fredrik Lindberg, fredrikl@gvc.gu.se
    Goteborg Urban Climate Group
    Gothenburg University

    :Input variables:
    dsm = digital surface model
    scale = height to pixel size (2m pixel gives scale = 0.5)
    svf,svfN,svfW,svfE,svfS = SVFs for building and ground
    svfveg,svfNveg,svfEveg,svfSveg,svfWveg = Veg SVFs blocking sky
    svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg = Veg SVFs blocking buildings
    vegdem = Vegetation canopy DSM
    vegdem2 = Vegetation trunk zone DSM
    albedo_b = building wall albedo
    absK = human absorption coefficient for shortwave radiation
    absL = human absorption coefficient for longwave radiation
    ewall = Emissivity of building walls
    Fside = The angular factors between a person and the surrounding surfaces
    Fup = The angular factors between a person and the surrounding surfaces
    Fcyl = The angular factors between a culidric person and the surrounding surfaces
    altitude = Sun altitude (degree)
    azimuth = Sun azimuth (degree)
    zen = Sun zenith angle (radians)
    jday = day of year
    usevegdem = use vegetation scheme
    onlyglobal = calculate dir and diff from global shortwave (Reindl et al. 1990)
    buildings = Boolena grid to identify building pixels
    location = geographic location
    height = height of measurements point (center of gravity of human)
    psi = 1 - Transmissivity of shortwave through vegetation
    landcover = use landcover scheme !!!NEW IN 2015a!!!
    lc_grid = grid with landcoverclasses
    lc_class = table with landcover properties
    dectime = decimal time
    altmax = maximum sun altitude
    dirwalls = aspect of walls
    walls = one pixel row outside building footprint. height of building walls
    cyl = consider man as cylinder instead of cude
    elvis = dummy
    Ta = air temp
    RH
    radG = global radiation
    radD = diffuse
    radI = direct
    P = pressure
    amaxvalue = max height of buildings
    bush = grid representing bushes
    Twater = temperature of water (daily)
    TgK, Tstart, TgK_wall, Tstart_wall, TmaxLST,TmaxLST_wall,
    alb_grid, emis_grid = albedo and emmissivity on ground
    first, second = conneted to old Ts model (source area based on Smidt et al.)
    svfalfa = SVF recalculated to angle
    svfbuveg = complete SVF
    firstdaytime, timeadd, timestepdec, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N,
    CI = Clearness index
    TgOut1 = old Ts model
    diffsh, ani = Used in anisotrpic models (Wallenberg et al. 2019, 2022)
    """
    

    # # # Core program start # # #
    # Instrument offset in degrees
    t = 0.0

    # Stefan Bolzmans Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = torch.pi / 180
    device = (
        Tg.device
        if isinstance(Tg, torch.Tensor)
        else Ta.device
        if isinstance(Ta, torch.Tensor)
        else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    )

    # Find sunrise decimal hour - new from 2014a
    _, _, _, SNUP = daylen(jday, location["latitude"])

    # Vapor pressure in hPa
    ea = 6.107 * 10 ** ((7.5 * Ta) / (237.3 + Ta)) * (RH / 100.0)

    # Determination of clear-sky emissivity from Prata (1996)
    msteg = 46.5 * (ea / (Ta + 273.15))
    esky = (
        1 - (1 + msteg) * torch.exp(-((1.2 + 3.0 * msteg) ** 0.5))
    ) + elvis  # -0.04 old error from Jonsson et al.2006

    if altitude > 0:  # # # # # # DAYTIME # # # # # #
        # Clearness Index on Earth's surface after Crawford and Dunchon (1999) with a correction
        # factor for low sun elevations after Lindberg et al.(2008)
        I0, CI, Kt, I0et, CIuncorr = clearnessindex_2013b(
            zen, jday, Ta, RH / 100.0, radG, location, P
        )
        if (CI > 1) or (CI == torch.inf):
            CI = 1

        # Estimation of radD and radI if not measured after Reindl et al.(1990)
        if onlyglobal == 1:
            I0, CI, Kt, I0et, CIuncorr = clearnessindex_2013b(
                zen, jday, Ta, RH / 100.0, radG, location, P
            )
            if (CI > 1) or (CI == torch.inf):
                CI = 1

            radI, radD = diffusefraction(radG, altitude, Kt, Ta, RH)

        # Diffuse Radiation
        # Anisotropic Diffuse Radiation after Perez et al. 1993
        if anisotropic_sky == 1:
            patchchoice = 1
            zenDeg = zen * (180 / torch.pi)
            # Relative luminance
            lv, pc_, pb_ = Perez_v3(
                zenDeg, azimuth, radD, radI, jday, patchchoice, patch_option
            )
            # Total relative luminance from sky, i.e. from each patch, into each cell
            aniLum = torch.zeros((rows, cols), device=device)
            for idx in range(lv.shape[0]):
                aniLum += diffsh[:, :, idx] * lv[idx, 2]

            dRad = (
                aniLum * radD
            )  # Total diffuse radiation from sky into each cell
        else:
            dRad = radD * svfbuveg
            patchchoice = 1
            lv = None

        # Shadow  images
        if usevegdem == 1:
            vegsh, sh, _, wallsh, wallsun, wallshve, _, facesun, wallsh_ = (
                shadowingfunction_wallheight_23(
                    dsm,
                    vegdem,
                    vegdem2,
                    azimuth,
                    altitude,
                    scale,
                    amaxvalue,
                    bush,
                    walls,
                    dirwalls * torch.pi / 180.0,
                    walls_scheme,
                    dirwalls_scheme * torch.pi / 180.0,
                )
            )
            shadow = sh - (1 - vegsh) * (1 - psi)
        else:
            sh, wallsh, wallsun, facesh, facesun, wallsh_ = (
                shadowingfunction_wallheight_13(
                    dsm,
                    azimuth,
                    altitude,
                    scale,
                    walls,
                    dirwalls * torch.pi / 180.0,
                    walls_scheme,
                    dirwalls_scheme * torch.pi / 180.0,
                )
            )
            shadow = sh

        # Building height angle from svf
        F_sh = cylindric_wedge(
            zen, svfalfa, rows, cols
        )  # Fraction shadow on building walls based on sun alt and svf
        F_sh[torch.isnan(F_sh)] = 0.5

        # New estimation of Tgwall with reduction for non-clear situation based on Reindl et al. 1990
        radI0, _ = diffusefraction(I0, altitude, 1.0, Ta, RH)
        radG0 = radI0 * (torch.sin(altitude * deg2rad)) + _
        corr = (
            0.1473 * torch.log(90 - (zen / torch.pi * 180)) + 0.3454
        )  # 20070329 correction of lat, Lindberg et al. 2008
        CI_TgG = (radG / radG0) + (1 - corr)
        if (CI_TgG > 1) or (CI_TgG == torch.inf):
            CI_TgG = 1

        Tgampwall = TgK_wall * altmax + Tstart_wall
        Tgwall = Tgampwall * torch.sin(
            (
                ((dectime - torch.floor(dectime)) - SNUP / 24)
                / (TmaxLST_wall / 24 - SNUP / 24)
            )
            * torch.pi
            / 2
        )  # 2015a, based on max sun altitude
        if Tgwall < 0:  # temporary for removing low Tg during morning 20130205
            Tgwall = 0
        Tgwall = Tgwall * CI_TgG

        # # # # Kdown # # # #
        Kdown = (
            radI * shadow * torch.sin(altitude * (torch.pi / 180))
            + dRad
            + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
        )  # *sin(altitude(i) * (pi / 180))

        # # # # Ldown # # # #
        Ldown = (
            (svf + svfveg - 1) * esky * SBC * ((Ta + 273.15) ** 4)
            + (2 - svfveg - svfaveg) * ewall * SBC * ((Ta + 273.15) ** 4)
            + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4)
            + (2 - svf - svfveg)
            * (1 - ewall)
            * esky
            * SBC
            * ((Ta + 273.15) ** 4)
        )  # Jonsson et al.(2006)
        # Ldown = Ldown - 25 # Shown by Jonsson et al.(2006) and Duarte et al.(2006)
        if CI < 0.95:  # non - clear conditions
            c = 1 - CI
            Ldown = Ldown * (1 - c) + c * (
                (svf + svfveg - 1) * SBC * ((Ta + 273.15) ** 4)
                + (2 - svfveg - svfaveg) * ewall * SBC * ((Ta + 273.15) ** 4)
                + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4)
                + (2 - svf - svfveg) * (1 - ewall) * SBC * ((Ta + 273.15) ** 4)
            )  # NOT REALLY TESTED!!! BUT MORE CORRECT?

        # Surface temperature parameterization during daytime
        if groundScheme == 1:
            # calculate the ground surface temperature, and relevant heat fluxes
            Tg, Rn, Rn_past, G = surfaceTemperature_calc(
                Kdown,
                Ldown,
                Rn,
                Rn_past,
                G,
                Tg,
                Tm,
                alb_grid,
                emis_grid,
                cap_grid,
                diff_grid,
                lc_grid,
                a1_grid,
                a2_grid,
                a3_grid,
                timeStep,
                RH,
                shadow,
                shadow_past,
            )

            # # # # Lup, daytime # # # #
            (
                Lup,
                gvfalb,
                gvfalbnosh,
                LupE,
                gvfalbE,
                gvfalbnoshE,
                LupS,
                gvfalbS,
                gvfalbnoshS,
                LupW,
                gvfalbW,
                gvfalbnoshW,
                LupN,
                gvfalbN,
                gvfalbnoshN,
                gvfLsideW,
                gvfLsideS,
                gvfLsideE,
                gvfLsideN,
            ) = outgoingLongwave_calc(
                Tg,
                Tgwall,
                Ta,
                Ldown,
                emis_grid,
                alb_grid,
                buildings,
                shadow,
                wallsun,
                walls,
                rows,
                cols,
                1 / scale,
            )

        else:
            # using max sun alt instead of dfm
            Tgamp = TgK * altmax + Tstart  # Fixed 2021
            Tgdiff = Tgamp * torch.sin(
                (
                    ((dectime - torch.floor(dectime)) - SNUP / 24)
                    / (TmaxLST / 24 - SNUP / 24)
                )
                * torch.pi
                / 2
            )  # 2015 a, based on max sun altitude

            Tgdiff = Tgdiff * CI_TgG  # new estimation

            # For Tg output in POIs
            TgTemp = Tgdiff * shadow + Ta
            _, timeadd, Tg = TsWaveDelay_2015a(
                TgTemp, firstdaytime, timeadd, timestepdec, Tg
            )  # timeadd only here v2021a

            if landcover == 1:
                Tg[Tg < 0] = (
                    0  # temporary for removing low Tg during morning 20130205
                )

            ### Ground View Factors
            (
                gvfLup,
                gvfalb,
                gvfalbnosh,
                gvfLupE,
                gvfalbE,
                gvfalbnoshE,
                gvfLupS,
                gvfalbS,
                gvfalbnoshS,
                gvfLupW,
                gvfalbW,
                gvfalbnoshW,
                gvfLupN,
                gvfalbN,
                gvfalbnoshN,
                gvfSum,
                gvfNorm,
            ) = gvf_2018a(
                wallsun,
                walls,
                buildings,
                scale,
                shadow,
                first,
                second,
                dirwalls,
                Tg,
                Tgwall,
                Ta,
                emis_grid,
                ewall,
                alb_grid,
                SBC,
                albedo_b,
                rows,
                cols,
                Twater,
                lc_grid,
                landcover,
            )

            # # # # Lup, daytime # # # #
            # Surface temperature wave delay - new as from 2014a
            Lup, timeaddnotused, Tgmap1 = TsWaveDelay_2015a(
                gvfLup, firstdaytime, timeadd, timestepdec, Tgmap1
            )
            LupE, timeaddnotused, Tgmap1E = TsWaveDelay_2015a(
                gvfLupE, firstdaytime, timeadd, timestepdec, Tgmap1E
            )
            LupS, timeaddnotused, Tgmap1S = TsWaveDelay_2015a(
                gvfLupS, firstdaytime, timeadd, timestepdec, Tgmap1S
            )
            LupW, timeaddnotused, Tgmap1W = TsWaveDelay_2015a(
                gvfLupW, firstdaytime, timeadd, timestepdec, Tgmap1W
            )
            LupN, timeaddnotused, Tgmap1N = TsWaveDelay_2015a(
                gvfLupN, firstdaytime, timeadd, timestepdec, Tgmap1N
            )

        # # # # Kup # # # #
        Kup, KupE, KupS, KupW, KupN = Kup_veg_2015a(
            radI,
            radD,
            radG,
            altitude,
            svfbuveg,
            albedo_b,
            F_sh,
            gvfalb,
            gvfalbE,
            gvfalbS,
            gvfalbW,
            gvfalbN,
            gvfalbnosh,
            gvfalbnoshE,
            gvfalbnoshS,
            gvfalbnoshW,
            gvfalbnoshN,
        )
        
        Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside = Kside_veg_v2022a(
            radI,
            radD,
            radG,
            shadow,
            svfS,
            svfW,
            svfN,
            svfE,
            svfEveg,
            svfSveg,
            svfWveg,
            svfNveg,
            azimuth,
            altitude,
            psi,
            t,
            albedo_b,
            F_sh,
            KupE,
            KupS,
            KupW,
            KupN,
            cyl,
            lv,
            anisotropic_sky,
            diffsh,
            rows,
            cols,
            asvf,
            shmat,
            vegshmat,
            vbshvegshmat,
        )

        firstdaytime = 0

    else:  # # # # # # # NIGHTTIME # # # # # # # #
        # Nocturnal K fluxes set to 0
        Knight = torch.zeros((rows, cols), device=device)
        Kdown = torch.zeros((rows, cols), device=device)
        Kwest = torch.zeros((rows, cols), device=device)
        Kup = torch.zeros((rows, cols), device=device)
        Keast = torch.zeros((rows, cols), device=device)
        Ksouth = torch.zeros((rows, cols), device=device)
        Knorth = torch.zeros((rows, cols), device=device)
        KsideI = torch.zeros((rows, cols), device=device)
        KsideD = torch.zeros((rows, cols), device=device)
        F_sh = torch.zeros((rows, cols), device=device)
        shadow = torch.zeros((rows, cols), device=device)
        CI_TgG = deepcopy(CI)
        dRad = torch.zeros((rows, cols), device=device)
        Kside = torch.zeros((rows, cols), device=device)
        wallsun = torch.zeros((rows, cols), device=device)

        Tgwall = 0

        # # # # Ldown # # # #
        Ldown = (
            (svf + svfveg - 1) * esky * SBC * ((Ta + 273.15) ** 4)
            + (2 - svfveg - svfaveg) * ewall * SBC * ((Ta + 273.15) ** 4)
            + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4)
            + (2 - svf - svfveg)
            * (1 - ewall)
            * esky
            * SBC
            * ((Ta + 273.15) ** 4)
        )  # Jonsson et al.(2006)
        # Ldown = Ldown - 25 # Shown by Jonsson et al.(2006) and Duarte et al.(2006)

        if CI < 0.95:  # non - clear conditions
            c = 1 - CI
            Ldown = Ldown * (1 - c) + c * (
                (svf + svfveg - 1) * SBC * ((Ta + 273.15) ** 4)
                + (2 - svfveg - svfaveg) * ewall * SBC * ((Ta + 273.15) ** 4)
                + (svfaveg - svf) * ewall * SBC * ((Ta + 273.15 + Tgwall) ** 4)
                + (2 - svf - svfveg) * (1 - ewall) * SBC * ((Ta + 273.15) ** 4)
            )  # NOT REALLY TESTED!!! BUT MORE CORRECT?

        # Surface temperature parameterization
        if groundScheme == 1:
            # calculate the ground surface temperature, and relevant heat fluxes
            Tg, Rn, Rn_past, G = surfaceTemperature_calc(
                Kdown,
                Ldown,
                Rn,
                Rn_past,
                G,
                Tg,
                Tm,
                alb_grid,
                emis_grid,
                cap_grid,
                diff_grid,
                lc_grid,
                a1_grid,
                a2_grid,
                a3_grid,
                timeStep,
                RH,
                shadow,
                shadow_past,
            )

            # # # # Lup, daytime # # # #
            (
                Lup,
                gvfalb,
                gvfalbnosh,
                LupE,
                gvfalbE,
                gvfalbnoshE,
                LupS,
                gvfalbS,
                gvfalbnoshS,
                LupW,
                gvfalbW,
                gvfalbnoshW,
                LupN,
                gvfalbN,
                gvfalbnoshN,
                gvfLsideW,
                gvfLsideS,
                gvfLsideE,
                gvfLsideN,
            ) = outgoingLongwave_calc(
                Tg,
                Tgwall,
                Ta,
                Ldown,
                emis_grid,
                alb_grid,
                buildings,
                shadow,
                wallsun,
                walls,
                rows,
                cols,
                1 / scale,
            )

        else:
            # In the old scheme the ground surface temperature is equal to the air temperature during nighttime
            Tg = torch.ones((rows, cols), device=device) * Ta

            # # # # Lup, nighttime # # # #
            Lup = SBC * emis_grid * ((Knight + Tg + 273.15) ** 4)
            LupE = Lup
            LupS = Lup
            LupW = Lup
            LupN = Lup

        I0 = 0
        timeadd = 0
        firstdaytime = 1

    # # # # Lside # # # #
    if groundScheme == 1:
        Least = torch.clone(gvfLsideE)
        Lsouth = torch.clone(gvfLsideS)
        Lwest = torch.clone(gvfLsideW)
        Lnorth = torch.clone(gvfLsideN)
        Least_, Lsouth_, Lwest_, Lnorth_ = Lside_veg_v2026(
            svfS,
            svfW,
            svfN,
            svfE,
            svfEveg,
            svfSveg,
            svfWveg,
            svfNveg,
            svfEaveg,
            svfSaveg,
            svfWaveg,
            svfNaveg,
            azimuth,
            altitude,
            Ta,
            Tgwall,
            SBC,
            ewall,
            Ldown,
            esky,
            t,
            F_sh,
            CI,
            anisotropic_sky,
        )
    else:
        Least = torch.zeros_like(Ldown)
        Lnorth = torch.zeros_like(Ldown)
        Lwest = torch.zeros_like(Ldown)
        Lsouth = torch.zeros_like(Ldown)
        Least_, Lsouth_, Lwest_, Lnorth_ = Lside_veg_v2022a(
            svfS,
            svfW,
            svfN,
            svfE,
            svfEveg,
            svfSveg,
            svfWveg,
            svfNveg,
            svfEaveg,
            svfSaveg,
            svfWaveg,
            svfNaveg,
            azimuth,
            altitude,
            Ta,
            Tgwall,
            SBC,
            ewall,
            Ldown,
            esky,
            t,
            F_sh,
            CI,
            LupE,
            LupS,
            LupW,
            LupN,
            anisotropic_sky,
        )

    Least += Least_
    Lsouth += Lsouth_
    Lwest += Lwest_
    Lnorth += Lnorth_
    Lside = (Lsouth + Lnorth + Least + Lwest) / 4

    ### Anisotropic sky
    if anisotropic_sky == 1:
        if "lv" not in locals():
            # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
            skyvaultalt, skyvaultazi, _, _, _, _, _ = create_patches(
                patch_option
            )

            patch_emissivities = torch.zeros(skyvaultalt.shape[0], device=device)

            x = torch.transpose(torch.atleast_2d(skyvaultalt))
            y = torch.transpose(torch.atleast_2d(skyvaultazi))
            z = torch.transpose(torch.atleast_2d(patch_emissivities))

            L_patches = torch.append(torch.append(x, y, axis=1), z, axis=1)

        else:
            L_patches = deepcopy(lv)

        # Calculate steradians for patches if it is the first model iteration
        if i == 0:
            steradians, skyalt, patch_altitude = patch_steradians(L_patches)

        # Create lv from L_patches if nighttime, i.e. lv does not exist
        if altitude < 0:
            # CI = deepcopy(CI)
            lv = deepcopy(L_patches)
            KupE = 0
            KupS = 0
            KupW = 0
            KupN = 0

        # Adjust sky emissivity under semi-cloudy/hazy/cloudy/overcast conditions, i.e. CI lower than 0.95
        if CI < 0.95:
            esky_c = CI * esky + (1 - CI) * 1.0
            esky = esky_c

        (
            Ldown,
            Lside_,
            Lside_sky,
            Lside_veg,
            Lside_sh,
            Lside_sun,
            Lside_ref,
            Least_,
            Lwest_,
            Lnorth_,
            Lsouth_,
            Keast,
            Ksouth,
            Kwest,
            Knorth,
            KsideI,
            KsideD,
            Kside,
            steradians,
            skyalt,
        ) = ani_sky(
            shmat,
            vegshmat,
            vbshvegshmat,
            altitude,
            azimuth,
            asvf,
            cyl,
            esky,
            L_patches,
            wallScheme,
            voxelTable,
            voxelMaps,
            steradians,
            Ta,
            Tgwall,
            ewall,
            Lup,
            radI,
            radD,
            radG,
            lv,
            albedo_b,
            0,
            diffsh,
            shadow,
            KupE,
            KupS,
            KupW,
            KupN,
            i,
        )
        Lside += Lside_
    else:
        Lside_ = torch.zeros((rows, cols), device=device)
        L_patches = None

    # Box and anisotropic longwave
    if cyl == 0 and anisotropic_sky == 1:
        Least += Least_
        Lwest += Lwest_
        Lnorth += Lnorth_
        Lsouth += Lsouth_

    # Calculation of radiant flux density
    # Human body considered as a cylinder with isotropic all-sky diffuse
    if cyl == 1 and anisotropic_sky == 0:
        Sstr = absK * (
            KsideI * Fcyl
            + (Kdown + Kup) * Fup
            + (Knorth + Keast + Ksouth + Kwest) * Fside
        ) + absL * (
            (Ldown + Lup) * Fup + (Lnorth + Least + Lsouth + Lwest) * Fside
        )
    # Human body considered as a cylinder with Perez et al. (1993) (anisotropic sky diffuse)
    # and Martin and Berdahl (1984) (anisotropic sky longwave)
    elif cyl == 1 and anisotropic_sky == 1:
        Sstr = absK * (
            Kside * Fcyl
            + (Kdown + Kup) * Fup
            + (Knorth + Keast + Ksouth + Kwest) * Fside
        ) + absL * (
            (Ldown + Lup) * Fup
            + Lside * Fcyl
            + (Lnorth + Least + Lsouth + Lwest) * Fside
        )
    # Human body considered as a standing cube
    else:
        Sstr = absK * (
            (Kdown + Kup) * Fup + (Knorth + Keast + Ksouth + Kwest) * Fside
        ) + absL * (
            (Ldown + Lup) * Fup + (Lnorth + Least + Lsouth + Lwest) * Fside
        )

    # # # # Tmrt # # # #
    Tmrt = torch.sqrt(torch.sqrt((Sstr / (absL * SBC)))) - 273.2

    # Add longwave to cardinal directions for output in POI
    if (cyl == 1) and (anisotropic_sky == 1):
        Least += Least_
        Lwest += Lwest_
        Lnorth += Lnorth_
        Lsouth += Lsouth_

    return (
        Tmrt,
        Kdown,
        Kup,
        Ldown,
        Lup,
        Tg,
        ea,
        esky,
        I0,
        CI,
        shadow,
        firstdaytime,
        timestepdec,
        timeadd,
        Tgmap1,
        Tgmap1E,
        Tgmap1S,
        Tgmap1W,
        Tgmap1N,
        Keast,
        Ksouth,
        Kwest,
        Knorth,
        Least,
        Lsouth,
        Lwest,
        Lnorth,
        KsideI,
        radI,
        radD,
        Lside,
        L_patches,
        CI_TgG,
        KsideD,
        dRad,
        Kside,
        steradians,
        voxelTable,
        Rn,
        Rn_past,
        Tm,
        G,
    )
