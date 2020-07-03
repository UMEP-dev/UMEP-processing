from builtins import range
import numpy as np
from ..util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_13 import shadowingfunction_wallheight_13
from ..util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_23 import shadowingfunction_wallheight_23
import linecache
import sys

def SEBE_2015a_calc(a, scale, slope, aspect, voxelheight, sizey, sizex, vegdem, vegdem2, walls, dirwalls, albedo, psi, 
                radmatI, radmatD, radmatR, usevegdem, feedback):

    # Parameters
    deg2rad = np.pi/180
    Knight = np.zeros((sizex, sizey))
    Energyyearroof = np.copy(Knight)

    if usevegdem == 1:
        # amaxvalue
        vegmax = vegdem.max()
        amaxvalue = a.max() - a.min()
        amaxvalue = np.maximum(amaxvalue, vegmax)

        # Elevation vegdsms if buildingDEM includes ground heights
        vegdem = vegdem+a
        vegdem[vegdem == a] = 0
        vegdem2 = vegdem2+a
        vegdem2[vegdem2 == a] = 0

        #% Bush separation
        bush = np.logical_not((vegdem2*vegdem))*vegdem
    else:
        psi = 1

    # Creating wallmatrix (1 meter interval)
    wallcol, wallrow = np.where(np.transpose(walls) > 0)    # row and col for each wall pixel
    wallstot = np.floor(walls * (1 / voxelheight)) * voxelheight
    wallsections = np.floor(np.max(walls) * (1 / voxelheight))     # finding tallest wall
    wallmatrix = np.zeros((np.shape(wallrow)[0], int(wallsections)))
    Energyyearwall = np.copy(wallmatrix)

    # Main loop - Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
    skyvaultaltint = np.array([6, 18, 30, 42, 54, 66, 78, 90])
    aziinterval = np.array([30, 30, 24, 24, 18, 12, 6, 1])

    if usevegdem == 1:
        wallshve = np.zeros(np.shape(a))
        vegrow, vegcol = np.where(vegdem > 0)  # row and col for each veg pixel
        vegdata = np.zeros((np.shape(vegrow)[0], 3))
        for i in range(0, vegrow.shape[0] - 1):
            vegdata[i, 0] = vegrow[i] + 1
            vegdata[i, 1] = vegcol[i] + 1
            vegdata[i, 2] = vegdem[vegrow[i], vegcol[i]]
    else:
        vegdata = 0

    index = 0
    for i in range(skyvaultaltint.size):
        for j in range(aziinterval[i]):

            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break

            #################### SOLAR RADIATION POSITIONS ###################
            #Solar Incidence angle (Roofs)
            suniroof = np.sin(slope) * np.cos(radmatI[index, 0] * deg2rad) * \
                        np.cos((radmatI[index, 1]*deg2rad)-aspect) + \
                        np.cos(slope) * np.sin((radmatI[index, 0] * deg2rad))

            suniroof[suniroof < 0] = 0

            # Solar Incidence angle (Walls)
            suniwall = np.abs(np.sin(np.pi/2) * np.cos(radmatI[index, 0] * deg2rad) *
                                np.cos((radmatI[index, 1] * deg2rad) - dirwalls*deg2rad) + np.cos(np.pi/2) *
                                np.sin((radmatI[index, 0] * deg2rad)))

            # Shadow image
            if usevegdem == 1:
                vegsh, sh, _, wallsh, wallsun, wallshve, _, facesun = shadowingfunction_wallheight_23(a,
                                    vegdem, vegdem2, radmatI[index, 1], radmatI[index, 0], scale, amaxvalue,
                                                                                bush, walls, dirwalls * deg2rad)
                shadow = np.copy(sh-(1.-vegsh)*(1.-psi))
            else:
                sh, wallsh, wallsun, facesh, facesun = shadowingfunction_wallheight_13(a, radmatI[index, 1],
                                                            radmatI[index, 0], scale, walls, dirwalls * deg2rad)
                shadow = np.copy(sh)

            # roof irradiance calculation
            # direct radiation
            if radmatI[index, 2] > 0:
                I = shadow * radmatI[index, 2] * suniroof
            else:
                I = np.copy(Knight)

            # roof diffuse and reflected radiation
            D = radmatD[index, 2] * shadow
            R = radmatR[index, 2] * (shadow*-1 + 1)

            Energyyearroof = np.copy(Energyyearroof+D+R+I)

            # WALL IRRADIANCE
            # direct radiation
            if radmatI[index, 2] > 0:
                Iw = radmatI[index, 2] * suniwall    # wall
            else:
                Iw = np.copy(Knight)

            # wall diffuse and reflected radiation
            Dw = radmatD[index, 2] * facesun
            Rw = radmatR[index, 2] * facesun

            # for each wall level (voxelheight interval)
            wallsun = np.floor(wallsun*(1/voxelheight)) * voxelheight
            wallsh = np.floor(wallsh*(1/voxelheight)) * voxelheight
            if usevegdem == 1:
                wallshve = np.floor(wallshve*(1/voxelheight)) * voxelheight

            wallmatrix = wallmatrix * 0

            for p in range(np.shape(wallmatrix)[0]):
                if wallsun[wallrow[p], wallcol[p]] > 0:    # Sections in sun
                    if wallsun[int(wallrow[p]), int(wallcol[p])] == wallstot[int(wallrow[p]), int(wallcol[p])]:  # All sections in sun
                        wallmatrix[p, 0:int(wallstot[int(wallrow[p]), int(wallcol[p])] / voxelheight)] = Iw[wallrow[p], wallcol[p]] + Dw[wallrow[p], wallcol[p]] + Rw[wallrow[p], wallcol[p]]
                    else:
                        wallmatrix[p, int((wallstot[wallrow[p], wallcol[p]] - wallsun[wallrow[p], wallcol[p]]) / voxelheight) - 1:int(wallstot[wallrow[p], wallcol[p]] / voxelheight)] = Iw[wallrow[p], wallcol[p]] + Dw[wallrow[p], wallcol[p]] + Rw[wallrow[p], wallcol[p]]

                if usevegdem == 1 and wallshve[wallrow[p], wallcol[p]] > 0:    # sections in vegetation shade
                    wallmatrix[p, 0:int((wallshve[int(wallrow[p]), int(wallcol[p])] + wallsh[int(wallrow[p]), int(wallcol[p])]) / voxelheight)] = (Iw[wallrow[p], wallcol[p]] + Dw[wallrow[p], wallcol[p]])*psi

                if wallsh[wallrow[p], wallcol[p]] > 0:    # sections in building shade
                    wallmatrix[p, 0:int(wallsh[wallrow[p], wallcol[p]] / voxelheight)] = Rw[wallrow[p], wallcol[p]]

            Energyyearwall = Energyyearwall + np.copy(wallmatrix)

            index = index + 1

            feedback.setProgress(int(index * (100. / 145.)))


    # Including radiation from ground on walls as well as removing pixels high than walls
    # fix_print_with_import
    wallmatrixbol = (Energyyearwall > 0).astype(float)
    Energyyearwall = (Energyyearwall + (np.sum(radmatR[:, 2]) * albedo)/2) * wallmatrixbol

    Energyyearroof /= 1000
    Energyyearwall /= 1000
    Energyyearwall = np.transpose(np.vstack((wallrow + 1, wallcol + 1, np.transpose(Energyyearwall))))    # adding 1 to wallrow and wallcol so that the tests pass

    seberesult = {'Energyyearroof': Energyyearroof, 'Energyyearwall': Energyyearwall, 'vegdata': vegdata}

    return seberesult