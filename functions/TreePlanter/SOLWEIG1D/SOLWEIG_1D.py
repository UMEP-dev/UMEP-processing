import numpy as np
import os
from ....util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import clearnessindex_2013b
from ....util.SEBESOLWEIGCommonFiles import Solweig_v2015_metdata_noload as metload
from ..SOLWEIG1D import Solweig1D_2019a_calc as so

def tmrt_1d_fun(metfilepath,tau,lon,lat,dsm,r_range):
    # Misc
    UTC = 1
    met = np.loadtxt(metfilepath, skiprows=1, delimiter=' ')
    alt = 10

    onlyglobal = 1  # 1 if only global radiation exist
    landcovercode = 0  # according to landcovrclasses_2018a_orig.txt
    metfile = 1  # 1 if time series data is used
    useveg = 1  # 1 if vegetation should be considered
    ani = 1
    cyl = 1
    elvis = 0
    alt = 3.0

    sh = 1.  # 0 if shadowed by building
    vegsh = 0.  # 0 if shadowed by tree
    svf = 0.6
    if useveg == 1:
        svfveg = 0.8
        svfaveg = 0.9
        trans = tau
    else:
        svfveg = 1.
        svfaveg = 1.
        trans = 1.

    absK = 0.7
    absL = 0.98
    PA = 'STAND'
    sensorheight = 2.0

    # program start
    if PA == 'STAND':
        Fside = 0.22
        Fup = 0.06
        height = 1.1
        Fcyl = 0.28
    else:
        Fside = 0.166666
        Fup = 0.166667
        height = 0.75
        Fcyl = 0.20

    if metfile == 1:
        met = np.loadtxt(metfilepath, skiprows=1, delimiter=' ')
    else:
        met = np.zeros((1, 24)) - 999.
        year = 2011
        month = 6
        day = 6
        hour = 12
        minu = 30

        if (year % 4) == 0:
            if (year % 100) == 0:
                if (year % 400) == 0:
                    leapyear = 1
                else:
                    leapyear = 0
            else:
                leapyear = 1
        else:
            leapyear = 0

        if leapyear == 1:
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        doy = np.sum(dayspermonth[0:month - 1]) + day

        Ta = 25.
        RH = 50.
        radG = 880.
        radD = 150.
        radI = 950.

        met[0, 0] = year
        met[0, 1] = doy
        met[0, 2] = hour
        met[0, 3] = minu
        met[0, 11] = Ta
        met[0, 10] = RH
        met[0, 14] = radG
        met[0, 21] = radD
        met[0, 22] = radI

    location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
    YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = metload.Solweig_2015a_metdata_noload(met, location,
                                                                                                       UTC)

    svfalfa = np.arcsin(np.exp((np.log((1. - svf)) / 2.)))

    dectime_new = np.zeros((dectime.shape[0], 1))
    dec_doy = int(dectime[0])
    for i in range(dectime.shape[0]):
        dectime_new[i, 0] = dectime[i] - dec_doy

    # Creating vectors from meteorological input
    DOY = met[:, 1]
    hours = met[:, 2]
    minu = met[:, 3]
    Ta = met[:, 11]
    RH = met[:, 10]
    radG = met[:, 14]
    radD = met[:, 21]
    radI = met[:, 22]
    P = met[:, 12]
    Ws = met[:, 9]
    Twater = []
    amaxvalue = dsm.max() - dsm.min()

    # load landcover file
    sitein = os.path.dirname(os.path.abspath(__file__)) + "/landcoverclasses_2018a_orig.txt"
    f = open(sitein)
    lin = f.readlines()
    lc_class = np.zeros((lin.__len__() - 1, 6))
    for i in range(1, lin.__len__()):
        lines = lin[i].split()
        for j in np.arange(1, 7):
            lc_class[i - 1, j - 1] = float(lines[j])
    f.close()

    # ground material parameters
    ground_pos = np.where(lc_class[:, 0] == landcovercode)
    albedo_g = lc_class[ground_pos, 1]
    eground = lc_class[ground_pos, 2]
    TgK = lc_class[ground_pos, 3]
    Tstart = lc_class[ground_pos, 4]
    TmaxLST = lc_class[ground_pos, 5]

    # wall material parameters
    wall_pos = np.where(lc_class[:, 0] == 99)
    albedo_b = lc_class[wall_pos, 1]
    ewall = lc_class[wall_pos, 2]
    TgK_wall = lc_class[wall_pos, 3]
    Tstart_wall = lc_class[wall_pos, 4]
    TmaxLST_wall = lc_class[wall_pos, 5]

    # If metfile starts at night
    CI = 1.

    if ani == 1:
        skyvaultalt = np.atleast_2d([])
        skyvaultazi = np.atleast_2d([])
        skyvaultaltint = [6, 18, 30, 42, 54, 66, 78]
        skyvaultaziint = [12, 12, 15, 15, 20, 30, 60]
        for j in range(7):
            for k in range(1, int(360 / skyvaultaziint[j]) + 1):
                skyvaultalt = np.append(skyvaultalt, skyvaultaltint[j])

        skyvaultalt = np.append(skyvaultalt, 90)

        diffsh = np.zeros((145))
        svfalfadeg = svfalfa / (np.pi / 180.)
        for k in range(0, 145):
            if skyvaultalt[k] > svfalfadeg:
                diffsh[k] = 1
    else:
        diffsh = []

    tmrt_1d = np.zeros((r_range.__len__(), 2))
    i_c = 0
    for i in r_range:
        # Daily water body temperature
        if (dectime[i] - np.floor(dectime[i])) == 0 or (i == 0):
            Twater = np.mean(Ta[jday[0] == np.floor(dectime[i])])

        # Nocturnal cloudfraction from Offerle et al. 2003
        if (dectime[i] - np.floor(dectime[i])) == 0:
            daylines = np.where(np.floor(dectime) == dectime[i])
            alt = altitude[0][daylines]
            alt2 = np.where(alt > 1)
            rise = alt2[0][0]
            [_, CI, _, _, _] = clearnessindex_2013b(zen[0, i + rise + 1], jday[0, i + rise + 1],
                                                    Ta[i + rise + 1],
                                                    RH[i + rise + 1] / 100., radG[i + rise + 1], location,
                                                    P[i + rise + 1])
            if (CI > 1) or (CI == np.inf):
                CI = 1

        Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, Keast, Ksouth, Kwest, Knorth, Least, Lsouth, Lwest, \
        Lnorth, KsideI, radIo, radDo, shadow1d = so.Solweig1D_2019a_calc(svf, svfveg, svfaveg, sh, vegsh, albedo_b,
                                                                         absK, absL, ewall,
                                                                         Fside, Fup, Fcyl,
                                                                         altitude[0][i], azimuth[0][i], zen[0][i],
                                                                         jday[0][i],
                                                                         onlyglobal, location, dectime[i], altmax[0][i],
                                                                         cyl, elvis,
                                                                         Ta[i], RH[i], radG[i], radD[i], radI[i], P[i],
                                                                         Twater, TgK, Tstart, albedo_g, eground,
                                                                         TgK_wall, Tstart_wall,
                                                                         TmaxLST, TmaxLST_wall, svfalfa, CI, ani,
                                                                         diffsh, trans)

        tmrt_1d[i_c, 0] = Tmrt
        tmrt_1d[i_c, 1] = hours[i]
        i_c += 1

    return tmrt_1d, azimuth, altitude, amaxvalue
