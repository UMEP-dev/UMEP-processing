from __future__ import division
from __future__ import absolute_import
from builtins import range
import torch
from ...util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
from ...util.SEBESOLWEIGCommonFiles.Perez_v3 import Perez_v3
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import (
    clearnessindex_2013b,
)
from ...util.SEBESOLWEIGCommonFiles.create_patches import create_patches


def sunmapcreator_2015a(
    met, altitude, azimuth, onlyglobal, output, jday, albedo, location, zen, device
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

    # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
    patch_option = 1  # 145 patches
    (
        skyvaultalt,
        skyvaultazi,
        annulino,
        skyvaultaltint,
        aziinterval,
        skyvaultaziint,
        azistart,
    ) = create_patches(patch_option, device)

    iangle2 = torch.tensor([])

    for j in range(len(aziinterval)):
        iangle2 = torch.append(
            iangle2, skyvaultaltint[j] * torch.ones([1, aziinterval[j]], device=device)
        )

    radmatI = torch.transpose(
        torch.vstack(
            (iangle2, skyvaultazi, torch.zeros((13, len(iangle2)), device=device))
        )
    )
    radmatD = torch.transpose(
        torch.vstack(
            (iangle2, skyvaultazi, torch.zeros((13, len(iangle2)), device=device))
        )
    )
    radmatR = torch.transpose(
        torch.vstack(
            (iangle2, skyvaultazi, torch.zeros((13, len(iangle2)), device=device))
        )
    )

    iazimuth = skyvaultazi

    # Main loop
    for i in range(len(met[:, 0])):
        alt = altitude[0, i]
        azi = azimuth[0, i]
        # disp(alt)
        if alt > 2:
            # Estimation of radD and radI if not measured after Reindl et al. (1990)
            if onlyglobal:
                if (
                    met[i, 11] <= -999.00
                    or met[i, 10] <= -999.00
                    or torch.isnan(met[i, 11])
                    or torch.isnan(met[i, 10])
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

            # Anisotropic diffuse distribution (Perez et al., 1993/Robinson & Stone, 2004)
            lv, _, _ = Perez_v3(
                90 - altitude[0, i],
                azimuth[0, i],
                D,
                I,
                jday[0, i],
                1,
                patch_option,
            )

            distalt = torch.abs(iangle2 - alt)
            altlevel = distalt == (torch.min(torch.abs(distalt)))
            distazi = torch.abs(iazimuth - azi)
            azipos = distazi[altlevel] == (torch.min(distazi[altlevel]))
            azipos2 = torch.where(altlevel)[0][0] + torch.where(azipos)[0][0]

            radmatI[azipos2, 2] = radmatI[azipos2, 2] + I
            radmatD[:, 2] = radmatD[:, 2] + D * lv[:, 2]
            radmatR[:, 2] = radmatR[:, 2] + G * (1 / 145) * albedo

            if output["energymonth"] == 1:
                radmatI[azipos2, met[i, 1] + 2] = radmatI[azipos2, met[i, 1] + 2] + I
                radmatD[:, met[i, 1] + 2] = radmatD[:, met[i, 1] + 2] + D * lv[:, 2]
                radmatR[:, met[i, 1] + 2] = (
                    radmatR[:, met[i, 1] + 2] + G * (1 / 145) * albedo
                )

    # Adjusting the numbers if multiple years is used

    if torch.shape(met)[0] > 8760:
        multiyear = torch.shape(met)[0] / 8760
        radmatI[:, 2:15] = radmatI[:, 2:15] / multiyear
        radmatD[:, 2:15] = radmatD[:, 2:15] / multiyear
        radmatR[:, 2:15] = radmatR[:, 2:15] / multiyear

    return radmatI, radmatD, radmatR
