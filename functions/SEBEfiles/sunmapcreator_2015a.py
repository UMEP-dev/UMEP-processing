from __future__ import division
from __future__ import absolute_import
from builtins import range
import numpy as np
from ...util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
from ...util.SEBESOLWEIGCommonFiles.Perez_v3 import Perez_v3
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import (
    clearnessindex_2013b,
)
from ...util.SEBESOLWEIGCommonFiles.create_patches import create_patches


def sunmapcreator_2015a(
    met, altitude, azimuth, onlyglobal, output, jday, albedo, location, zen
):
    """
    % This function creates a sun map based on hourly values of solar radiation.

    :param met:
    :param altitude: 2D array
    :param azimuth:
    :param onlyglobal:
    :param output:
    :param jday:
    :param albedo:
    :return:
    """
    np.seterr(over="raise")
    np.seterr(invalid="raise")

    # Creating skyvault of patches of constant radians (Tregeneza and
    # Sharples, 1993)
    patch_option = 1  # 145 patches
    # patch_option = 2 # 153 patches
    # patch_option = 3 # 306 patches
    # patch_option = 4 # 612 patches
    (
        skyvaultalt,
        skyvaultazi,
        annulino,
        skyvaultaltint,
        aziinterval,
        skyvaultaziint,
        azistart,
    ) = create_patches(patch_option)

    iangle2 = np.array([])
    Gyear = 0
    Dyear = 0
    Gmonth = np.zeros([1, 12])
    Dmonth = Gmonth
    for j in range(len(aziinterval)):
        iangle2 = np.append(
            iangle2, skyvaultaltint[j] * np.ones([1, aziinterval[j]])
        )

    radmatI = np.transpose(
        np.vstack((iangle2, skyvaultazi, np.zeros((13, len(iangle2)))))
    )
    radmatD = np.transpose(
        np.vstack((iangle2, skyvaultazi, np.zeros((13, len(iangle2)))))
    )
    radmatR = np.transpose(
        np.vstack((iangle2, skyvaultazi, np.zeros((13, len(iangle2)))))
    )

    iazimuth = skyvaultazi
    # Ta = met[:, 11]
    # RH = met[:, 10]

    # Main loop
    for i in range(len(met[:, 0])):
        alt = altitude[0, i]
        azi = azimuth[0, i]
        # disp(alt)
        if alt > 2:
            # Estimation of radD and radI if not measured after Reindl et al.
            # (1990)
            if onlyglobal:
                if (
                    met[i, 11] <= -999.00
                    or met[i, 10] <= -999.00
                    or np.isnan(met[i, 11])
                    or np.isnan(met[i, 10])
                ):
                    met[i, 11] = 15.0
                    met[i, 10] = 75.0

                I0, CI, Kt, I0et, CIuncorr = clearnessindex_2013b(
                    zen[0, i],
                    jday[0, i],
                    met[i, 11],
                    met[i, 10],
                    met[i, 14],
                    location,
                    -999.0,
                )
                I, D = diffusefraction(
                    met[i, 14], altitude[0, i], Kt, met[i, 11], met[i, 10]
                )
            else:
                I = met[i, 22]
                D = met[i, 21]

            G = met[i, 14]

            # Anisotropic diffuse distribution (Perez et al., 1993/Robinson &
            # Stone, 2004)
            lv, _, _ = Perez_v3(
                90 - altitude[0, i],
                azimuth[0, i],
                D,
                I,
                jday[0, i],
                1,
                patch_option,
            )

            distalt = np.abs(iangle2 - alt)
            altlevel = distalt == (np.min(np.abs(distalt)))
            distazi = np.abs(iazimuth - azi)
            azipos = distazi[altlevel] == (np.min(distazi[altlevel]))
            azipos2 = np.where(altlevel)[0][0] + np.where(azipos)[0][0]
            # azipos2 = np.where(altlevel)[0] + np.where(azipos)[0]
            # azipos2 = find(altlevel, 1)-1 + find(azipos, 1)
            radmatI[azipos2, 2] = radmatI[azipos2, 2] + I
            radmatD[:, 2] = radmatD[:, 2] + D * lv[:, 2]
            radmatR[:, 2] = radmatR[:, 2] + G * (1 / 145) * albedo
            #         Gyear=Gyear+(G*sin(altitude(i)*(pi/180)));
            #         Dyear=Dyear+D;

            if output["energymonth"] == 1:
                radmatI[azipos2, met[i, 1] + 2] = (
                    radmatI[azipos2, met[i, 1] + 2] + I
                )
                radmatD[:, met[i, 1] + 2] = (
                    radmatD[:, met[i, 1] + 2] + D * lv[:, 2]
                )
                radmatR[:, met[i, 1] + 2] = (
                    radmatR[:, met[i, 1] + 2] + G * (1 / 145) * albedo
                )

    # Adjusting the numbers if multiple years is used
    if np.shape(met)[0] > 8760:
        multiyear = np.shape(met)[0] / 8760
        radmatI[:, 2:15] = radmatI[:, 2:15] / multiyear
        radmatD[:, 2:15] = radmatD[:, 2:15] / multiyear
        radmatR[:, 2:15] = radmatR[:, 2:15] / multiyear

    return radmatI, radmatD, radmatR
