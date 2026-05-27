from builtins import range
import numpy as np
import torch
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_13 import (
    shadowingfunction_wallheight_13,
)
from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_23 import (
    shadowingfunction_wallheight_23,
)


def SEBE_2015a_calc(
    a,
    scale,
    slope,
    aspect,
    voxelheight,
    sizey,
    sizex,
    vegdem,
    vegdem2,
    walls,
    dirwalls,
    albedo,
    psi,
    radmatI,
    radmatD,
    radmatR,
    usevegdem,
    feedback,
    wallmaxheight,
    device,
):

    # Parameters
    deg2rad = torch.pi / 180
    Knight = torch.zeros((sizex, sizey), device=device)
    Energyyearroof = torch.clone(Knight)

    if usevegdem == 1:
        # amaxvalue
        vegmax = vegdem.max()
        amaxvalue = a.max() - a.min()
        amaxvalue = torch.maximum(amaxvalue, vegmax)

        # Elevation vegdsms if buildingDEM includes ground heights
        vegdem = vegdem + a
        vegdem[vegdem == a] = 0
        vegdem2 = vegdem2 + a
        vegdem2[vegdem2 == a] = 0

        # % Bush separation
        bush = torch.logical_not((vegdem2 * vegdem)) * vegdem
    else:
        psi = 1

    # Creating wallmatrix (1 meter interval)
    wallcol, wallrow = torch.where(
        torch.transpose(walls) > 0
    )  # row and col for each wall pixel
    wallstot = torch.floor(walls * (1 / voxelheight)) * voxelheight
    # wallsections = torch.floor(torch.max(walls) * (1 / voxelheight))     # finding tallest wall
    wallsections = torch.floor(wallmaxheight * (1 / voxelheight))
    # feedback.setProgressText('torch.max(walls):' + str(torch.max(walls)))
    # feedback.setProgressText('1 / voxelheight:' + str(1 / voxelheight))
    # feedback.setProgressText('voxel:' + str(voxelheight))
    # feedback.setProgressText('wallsections:' + str(wallsections))
    # feedback.setProgressText('torch.shape(wallrow)[0]:' + str(torch.shape(wallrow)[0]))
    wallmatrix = torch.zeros(
        (torch.shape(wallrow)[0], int(wallsections)), device=device
    )
    Energyyearwall = torch.clone(wallmatrix)

    # Main loop - Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
    skyvaultaltint = torch.tensor(
        [6, 18, 30, 42, 54, 66, 78, 90], device=device
    )
    aziinterval = torch.tensor([30, 30, 24, 24, 18, 12, 6, 1], device=device)

    if usevegdem == 1:
        wallshve = torch.zeros(torch.shape(a), device=device)
        vegrow, vegcol = torch.where(
            vegdem > 0
        )  # row and col for each veg pixel
        vegdata = torch.zeros((torch.shape(vegrow)[0], 3), device=device)
        for i in range(0, vegrow.shape[0] - 1):
            vegdata[i, 0] = vegrow[i] + 1
            vegdata[i, 1] = vegcol[i] + 1
            vegdata[i, 2] = vegdem[vegrow[i], vegcol[i]]
    else:
        vegdata = 0

    index = 0
    for i in range(skyvaultaltint.size):
        for j in range(aziinterval[i]):

            feedback.setProgress(int(index * (100.0 / 145.0)))

            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break

            #################### SOLAR RADIATION POSITIONS ###################
            # Solar Incidence angle (Roofs)
            suniroof = torch.sin(slope) * torch.cos(
                radmatI[index, 0] * deg2rad
            ) * torch.cos((radmatI[index, 1] * deg2rad) - aspect) + torch.cos(
                slope
            ) * torch.sin(
                (radmatI[index, 0] * deg2rad)
            )

            suniroof[suniroof < 0] = 0

            # Solar Incidence angle (Walls)
            suniwall = torch.abs(
                torch.sin(torch.pi / 2)
                * torch.cos(radmatI[index, 0] * deg2rad)
                * torch.cos((radmatI[index, 1] * deg2rad) - dirwalls * deg2rad)
                + torch.cos(torch.pi / 2)
                * torch.sin((radmatI[index, 0] * deg2rad))
            )

            # Shadow image
            if usevegdem == 1:
                vegsh, sh, _, wallsh, wallsun, wallshve, _, facesun = (
                    shadowingfunction_wallheight_23(
                        a,
                        vegdem,
                        vegdem2,
                        radmatI[index, 1],
                        radmatI[index, 0],
                        scale,
                        amaxvalue,
                        bush,
                        walls,
                        dirwalls * deg2rad,
                    )
                )
                shadow = torch.clone(sh - (1.0 - vegsh) * (1.0 - psi))
            else:
                sh, wallsh, wallsun, facesh, facesun = (
                    shadowingfunction_wallheight_13(
                        a,
                        radmatI[index, 1],
                        radmatI[index, 0],
                        scale,
                        walls,
                        dirwalls * deg2rad,
                    )
                )
                shadow = torch.clone(sh)

            # roof irradiance calculation
            # direct radiation
            if radmatI[index, 2] > 0:
                I = shadow * radmatI[index, 2] * suniroof
            else:
                I = torch.clone(Knight)

            # roof diffuse and reflected radiation
            D = radmatD[index, 2] * shadow
            R = radmatR[index, 2] * (shadow * -1 + 1)

            Energyyearroof = torch.clone(Energyyearroof + D + R + I)

            # WALL IRRADIANCE
            # direct radiation
            if radmatI[index, 2] > 0:
                Iw = radmatI[index, 2] * suniwall  # wall
            else:
                Iw = torch.clone(Knight)

            # wall diffuse and reflected radiation
            Dw = radmatD[index, 2] * facesun
            Rw = radmatR[index, 2] * facesun

            # for each wall level (voxelheight interval)
            wallsun = torch.floor(wallsun * (1 / voxelheight)) * voxelheight
            wallsh = torch.floor(wallsh * (1 / voxelheight)) * voxelheight
            if usevegdem == 1:
                wallshve = (
                    torch.floor(wallshve * (1 / voxelheight)) * voxelheight
                )

            wallmatrix = wallmatrix * 0

            for p in range(torch.shape(wallmatrix)[0]):
                if wallsun[wallrow[p], wallcol[p]] > 0:  # Sections in sun
                    if (
                        wallsun[int(wallrow[p]), int(wallcol[p])]
                        == wallstot[int(wallrow[p]), int(wallcol[p])]
                    ):  # All sections in sun
                        wallmatrix[
                            p,
                            0 : int(
                                wallstot[int(wallrow[p]), int(wallcol[p])]
                                / voxelheight
                            ),
                        ] = (
                            Iw[wallrow[p], wallcol[p]]
                            + Dw[wallrow[p], wallcol[p]]
                            + Rw[wallrow[p], wallcol[p]]
                        )
                    else:
                        wallmatrix[
                            p,
                            int(
                                (
                                    wallstot[wallrow[p], wallcol[p]]
                                    - wallsun[wallrow[p], wallcol[p]]
                                )
                                / voxelheight
                            )
                            - 1 : int(
                                wallstot[wallrow[p], wallcol[p]] / voxelheight
                            ),
                        ] = (
                            Iw[wallrow[p], wallcol[p]]
                            + Dw[wallrow[p], wallcol[p]]
                            + Rw[wallrow[p], wallcol[p]]
                        )

                if (
                    usevegdem == 1 and wallshve[wallrow[p], wallcol[p]] > 0
                ):  # sections in vegetation shade
                    wallmatrix[
                        p,
                        0 : int(
                            (
                                wallshve[int(wallrow[p]), int(wallcol[p])]
                                + wallsh[int(wallrow[p]), int(wallcol[p])]
                            )
                            / voxelheight
                        ),
                    ] = (
                        Iw[wallrow[p], wallcol[p]] + Dw[wallrow[p], wallcol[p]]
                    ) * psi

                if (
                    wallsh[wallrow[p], wallcol[p]] > 0
                ):  # sections in building shade
                    wallmatrix[
                        p,
                        0 : int(wallsh[wallrow[p], wallcol[p]] / voxelheight),
                    ] = Rw[wallrow[p], wallcol[p]]

            Energyyearwall = Energyyearwall + torch.clone(wallmatrix)

            index = index + 1

    # Including radiation from ground on walls as well as removing pixels high than walls
    # fix_print_with_import
    wallmatrixbol = (Energyyearwall > 0).astype(float)
    Energyyearwall = (
        Energyyearwall + (torch.sum(radmatR[:, 2]) * albedo) / 2
    ) * wallmatrixbol

    Energyyearroof /= 1000
    Energyyearwall /= 1000
    Energyyearwall = torch.transpose(
        torch.vstack(
            (wallrow + 1, wallcol + 1, torch.transpose(Energyyearwall))
        )
    )  # adding 1 to wallrow and wallcol so that the tests pass

    seberesult = {
        "Energyyearroof": Energyyearroof,
        "Energyyearwall": Energyyearwall,
        "vegdata": vegdata,
    }

    return seberesult
