try:
    import torch
except:
    pass

from ..util import shadowingfunctions_torch as shadow
from ..util.SEBESOLWEIGCommonFiles.create_patches_torch import create_patches


def _to_tensor(x, device, dtype=torch.float32):
    if isinstance(x, torch.Tensor):
        return x.to(device)
    if x is None:
        return None
    return torch.tensor(x, dtype=dtype, device=device)


def _sync_if_cuda(device):
    if device is not None and device.type == "cuda":
        torch.cuda.synchronize()


# from ..functions.wallalgorithms import findwalls
from . import wallalgorithms_torch as wa
from . import svf_for_voxels_torch as svfv
from ..util.SEBESOLWEIGCommonFiles import (
    shadowingfunction_wallheight_13_torch as shb,
)
from ..util.SEBESOLWEIGCommonFiles import (
    shadowingfunction_wallheight_23_torch as shbv,
)

# remove
from ..util.misc import saveraster
from osgeo.gdalconst import *
from osgeo import gdal, osr


def annulus_weight(altitude, aziinterval, device):
    n = torch.tensor(90.0, device=device)
    steprad = (360.0 / aziinterval) * (torch.pi / 180.0)
    annulus = 91.0 - altitude
    w = (
        (1.0 / (2.0 * torch.pi))
        * torch.sin(torch.pi / (2.0 * n))
        * torch.sin((torch.pi * (2.0 * annulus - 1.0)) / (2.0 * n))
    )
    weight = steprad * w

    return weight


def svf_angles_100121(device):
    azi1 = torch.arange(1.0, 360.0, 360.0 / 16.0)  # %22.5
    azi2 = torch.arange(12.0, 360.0, 360.0 / 16.0)  # %22.5
    azi3 = torch.arange(5.0, 360.0, 360.0 / 32.0)  # %11.25
    azi4 = torch.arange(2.0, 360.0, 360.0 / 32.0)  # %11.25
    azi5 = torch.arange(4.0, 360.0, 360.0 / 40.0)  # %9
    azi6 = torch.arange(7.0, 360.0, 360.0 / 48.0)  # %7.50
    azi7 = torch.arange(6.0, 360.0, 360.0 / 48.0)  # %7.50
    azi8 = torch.arange(1.0, 360.0, 360.0 / 48.0)  # %7.50
    azi9 = torch.arange(4.0, 359.0, 360.0 / 52.0)  # %6.9231
    azi10 = torch.arange(5.0, 360.0, 360.0 / 52.0)  # %6.9231
    azi11 = torch.arange(1.0, 360.0, 360.0 / 48.0)  # %7.50
    azi12 = torch.arange(0.0, 359.0, 360.0 / 44.0)  # %8.1818
    azi13 = torch.arange(3.0, 360.0, 360.0 / 44.0)  # %8.1818
    azi14 = torch.arange(2.0, 360.0, 360.0 / 40.0)  # %9
    azi15 = torch.arange(7.0, 360.0, 360.0 / 32.0)  # %10
    azi16 = torch.arange(3.0, 360.0, 360.0 / 24.0)  # %11.25
    azi17 = torch.arange(10.0, 360.0, 360.0 / 16.0)  # %15
    azi18 = torch.arange(19.0, 360.0, 360.0 / 12.0)  # %22.5
    azi19 = torch.arange(17.0, 360.0, 360.0 / 8.0)  # %45
    azi20 = 0.0  # %360
    iazimuth = torch.tensor(
        torch.hstack(
            (
                azi1,
                azi2,
                azi3,
                azi4,
                azi5,
                azi6,
                azi7,
                azi8,
                azi9,
                azi10,
                azi11,
                azi12,
                azi13,
                azi14,
                azi15,
                azi16,
                azi17,
                azi18,
                azi19,
                azi20,
            )
        ),
        device=device,
    )
    aziinterval = torch.tensor(
        torch.hstack(
            (
                16.0,
                16.0,
                32.0,
                32.0,
                40.0,
                48.0,
                48.0,
                48.0,
                52.0,
                52.0,
                48.0,
                44.0,
                44.0,
                40.0,
                32.0,
                24.0,
                16.0,
                12.0,
                8.0,
                1.0,
            )
        ),
        device=device,
    )
    angleresult = {"iazimuth": iazimuth, "aziinterval": aziinterval}

    return angleresult


def svfForProcessing153(
    dsm,
    vegdem,
    vegdem2,
    scale,
    usevegdem,
    pixel_resolution,
    wallScheme,
    demlayer,
    feedback,
    device=torch.device("cpu"),
):

    with torch.no_grad():

        dsm = _to_tensor(dsm, device)
        vegdem = _to_tensor(vegdem, device)
        vegdem2 = _to_tensor(vegdem2, device)
        demlayer = _to_tensor(demlayer, device)

        rows = dsm.shape[0]
        cols = dsm.shape[1]
        svf = torch.zeros([rows, cols], device=device)
        svfE = torch.zeros([rows, cols], device=device)
        svfS = torch.zeros([rows, cols], device=device)
        svfW = torch.zeros([rows, cols], device=device)
        svfN = torch.zeros([rows, cols], device=device)
        svfveg = torch.zeros((rows, cols), device=device)
        svfEveg = torch.zeros((rows, cols), device=device)
        svfSveg = torch.zeros((rows, cols), device=device)
        svfWveg = torch.zeros((rows, cols), device=device)
        svfNveg = torch.zeros((rows, cols), device=device)
        svfaveg = torch.zeros((rows, cols), device=device)
        svfEaveg = torch.zeros((rows, cols), device=device)
        svfSaveg = torch.zeros((rows, cols), device=device)
        svfWaveg = torch.zeros((rows, cols), device=device)
        svfNaveg = torch.zeros((rows, cols), device=device)

        # % amaxvalue
        vegmax = vegdem.max()
        amaxvalue = dsm.max()
        amaxvalue = torch.maximum(amaxvalue, vegmax)

        # % Elevation vegdems if buildingDSM inclused ground heights
        vegdem = vegdem + dsm
        vegdem[vegdem == dsm] = 0
        vegdem2 = vegdem2 + dsm
        vegdem2[vegdem2 == dsm] = 0
        # % Bush separation
        bush = torch.logical_not((vegdem2 * vegdem)) * vegdem

        # index = int(0)

        # patch_option = 1 # 145 patches
        patch_option = 2  # 153 patches
        # patch_option = 3 # 306 patches
        # patch_option = 4 # 612 patches

        # Create patches based on patch_option
        (
            skyvaultalt,
            skyvaultazi,
            annulino,
            skyvaultaltint,
            aziinterval,
            skyvaultaziint,
            azistart,
        ) = create_patches(patch_option, device)

        skyvaultalt = _to_tensor(skyvaultalt, device)
        skyvaultazi = _to_tensor(skyvaultazi, device)
        annulino = _to_tensor(annulino, device)
        skyvaultaltint = _to_tensor(skyvaultaltint, device)
        aziinterval = _to_tensor(aziinterval, device)
        skyvaultaziint = _to_tensor(skyvaultaziint, device)
        azistart = _to_tensor(azistart, device)

        skyvaultaziint = torch.tensor(
            [360 / patches for patches in aziinterval], device=device
        )
        iazimuth = torch.zeros(
            int(torch.sum(aziinterval).item()), device=device
        )

        shmat = torch.zeros(
            (rows, cols, int(torch.sum(aziinterval).item())), device=device
        )
        vegshmat = torch.zeros(
            (rows, cols, int(torch.sum(aziinterval).item())), device=device
        )
        vbshvegshmat = torch.zeros(
            (rows, cols, int(torch.sum(aziinterval).item())), device=device
        )

        # Preparations for wall temperature scheme
        if wallScheme:
            feedback.setProgressText("Estimating view factors for wall voxels")
            (
                voxelTable,
                voxelId_list,
                wall_dict,
                walls,
                aspect,
                uniqueWallIDs,
                wall2d_id,
                voxel_height,
            ) = svfv.wallscheme_prepare(
                dsm, scale, pixel_resolution, feedback, device=device
            )

            # Rasters to fill with values in loop
            all_buildIDSeen = torch.zeros(
                (rows, cols, skyvaultalt.shape[0]), device=device
            )
            all_voxelHeight = torch.zeros(
                (rows, cols, skyvaultalt.shape[0]), device=device
            )
            all_voxelId = torch.zeros(
                (rows, cols, skyvaultalt.shape[0]), device=device
            )
        else:
            voxelTable = 0
            allbuildIDSeen = 0
            allvoxelHeight = 0
            all_voxelId = 0
            walls = 0

        index = 0
        for j in range(0, skyvaultaltint.shape[0]):
            for k in range(0, int(360 / skyvaultaziint[j].item())):
                iazimuth[index] = k * skyvaultaziint[j] + azistart[j]
                if iazimuth[index] > 360.0:
                    iazimuth[index] = iazimuth[index] - 360.0
                index = index + 1
        aziintervalaniso = torch.ceil(aziinterval / 2.0)
        index = 0
        for i in range(0, skyvaultaltint.shape[0]):
            for j in range(0, int(aziinterval[int(i)].item())):
                if feedback.isCanceled():
                    feedback.setProgressText("Calculation cancelled")
                    break
                altitude = torch.tensor(
                    float(skyvaultaltint[int(i)].item()), device=device
                )
                azimuth = torch.tensor(
                    float(iazimuth[int(index)].item()), device=device
                )

                # Casting shadow
                if wallScheme:
                    if usevegdem == 1:
                        (
                            vegsh,
                            sh,
                            vbshvegsh,
                            wallsh,
                            wallsun,
                            wallshve,
                            facesh,
                            facesun,
                        ) = shbv.shadowingfunction_wallheight_23(
                            dsm,
                            vegdem,
                            vegdem2,
                            azimuth,
                            altitude,
                            scale,
                            amaxvalue,
                            bush,
                            walls,
                            aspect * torch.pi / 180,
                        )
                        vegshmat[:, :, index] = vegsh
                        vbshvegshmat[:, :, index] = vbshvegsh
                    else:
                        sh, wallsh, wallsun, facesh, facesun = (
                            shb.shadowingfunction_wallheight_13(
                                dsm,
                                azimuth,
                                altitude,
                                scale,
                                walls,
                                aspect * torch.pi / 180.0,
                            )
                        )
                        vegsh = torch.ones(
                            (sh.shape[0], sh.shape[1]),
                            device=device,
                            dtype=torch.float32,
                        )
                        vbshvegsh = torch.ones(
                            (sh.shape[0], sh.shape[1]),
                            device=device,
                            dtype=torch.float32,
                        )
                        vegshmat[:, :, index] = vegsh
                        vbshvegshmat[:, :, index] = vbshvegsh
                else:
                    if usevegdem == 1:
                        shadowresult = shadow.shadowingfunction_20(
                            dsm,
                            vegdem,
                            vegdem2,
                            azimuth,
                            altitude,
                            scale,
                            amaxvalue,
                            bush,
                            feedback,
                            1,
                            device,
                        )

                        vegsh = torch.tensor(
                            shadowresult["vegsh"], device=device
                        )
                        vbshvegsh = torch.tensor(
                            shadowresult["vbshvegsh"], device=device
                        )
                        vegshmat[:, :, index] = vegsh
                        vbshvegshmat[:, :, index] = vbshvegsh
                        sh = torch.tensor(shadowresult["sh"], device=device)
                    else:
                        sh = shadow.shadowingfunctionglobalradiation(
                            dsm, azimuth, altitude, scale, feedback, 1, device
                        )

                shmat[:, :, index] = sh

                # Wall temperature scheme, i.e. finding out which voxel is seen from each pixel, where direction is patch azimuth and altitude
                if wallScheme:
                    (
                        all_buildIDSeen[:, :, index],
                        all_voxelHeight[:, :, index],
                        all_voxelId[:, :, index],
                    ) = shadow.shadowingfunction_findwallID(
                        dsm,
                        azimuth,
                        altitude,
                        scale,
                        walls,
                        uniqueWallIDs,
                        demlayer,
                        wall2d_id,
                        voxel_height,
                        voxelId_list,
                        facesh,
                        wall_dict,
                        sh,
                        device,
                    )

                # Calculate svfs
                for k in range(
                    int(annulino[int(i)].item()) + 1,
                    int(annulino[int(i + 1.0)].item()) + 1,
                ):

                    weight = annulus_weight(k, aziinterval[i], device) * sh
                    svf = svf + weight
                    weight = (
                        annulus_weight(k, aziintervalaniso[i], device) * sh
                    )
                    if (azimuth >= 0) and (azimuth < 180):
                        svfE = svfE + weight
                    if (azimuth >= 90) and (azimuth < 270):
                        svfS = svfS + weight
                    if (azimuth >= 180) and (azimuth < 360):
                        svfW = svfW + weight
                    if (azimuth >= 270) or (azimuth < 90):
                        svfN = svfN + weight

                if usevegdem == 1:
                    for k in torch.arange(
                        annulino[int(i)] + 1, (annulino[int(i + 1.0)]) + 1
                    ):
                        # % changed to include 90
                        weight = annulus_weight(k, aziinterval[i], device)
                        svfveg = svfveg + weight * vegsh
                        svfaveg = svfaveg + weight * vbshvegsh
                        weight = annulus_weight(k, aziintervalaniso[i], device)
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
                feedback.setProgress(
                    int(index * (100.0 / torch.sum(aziinterval)))
                )

        svfS = svfS + 3.0459e-004
        svfW = svfW + 3.0459e-004
        # % Last azimuth is 90. Hence, manual add of last annuli for svfS and SVFW
        # %Forcing svf not be greater than 1 (some MATLAB crazyness)
        svf[(svf > 1.0)] = 1.0
        svfE[(svfE > 1.0)] = 1.0
        svfS[(svfS > 1.0)] = 1.0
        svfW[(svfW > 1.0)] = 1.0
        svfN[(svfN > 1.0)] = 1.0

        if usevegdem == 1:
            last = torch.zeros((rows, cols), device=device)
            last[(vegdem2 == 0.0)] = 3.0459e-004
            svfSveg = svfSveg + last
            svfWveg = svfWveg + last
            svfSaveg = svfSaveg + last
            svfWaveg = svfWaveg + last
            # %Forcing svf not be greater than 1 (some MATLAB crazyness)
            svfveg[(svfveg > 1.0)] = 1.0
            svfEveg[(svfEveg > 1.0)] = 1.0
            svfSveg[(svfSveg > 1.0)] = 1.0
            svfWveg[(svfWveg > 1.0)] = 1.0
            svfNveg[(svfNveg > 1.0)] = 1.0
            svfaveg[(svfaveg > 1.0)] = 1.0
            svfEaveg[(svfEaveg > 1.0)] = 1.0
            svfSaveg[(svfSaveg > 1.0)] = 1.0
            svfWaveg[(svfWaveg > 1.0)] = 1.0
            svfNaveg[(svfNaveg > 1.0)] = 1.0

        svfresult = {
            "svf": svf,
            "svfE": svfE,
            "svfS": svfS,
            "svfW": svfW,
            "svfN": svfN,
            "svfveg": svfveg,
            "svfEveg": svfEveg,
            "svfSveg": svfSveg,
            "svfWveg": svfWveg,
            "svfNveg": svfNveg,
            "svfaveg": svfaveg,
            "svfEaveg": svfEaveg,
            "svfSaveg": svfSaveg,
            "svfWaveg": svfWaveg,
            "svfNaveg": svfNaveg,
            "shmat": shmat,
            "vegshmat": vegshmat,
            "vbshvegshmat": vbshvegshmat,
            "voxelIds": all_voxelId,
            "voxelTable": voxelTable,
            "walls": walls,
        }

        # ,
        # 'vbshvegshmat': vbshvegshmat, 'wallshmat': wallshmat, 'wallsunmat': wallsunmat,
        # 'wallshvemat': wallshvemat, 'facesunmat': facesunmat}
        return svfresult


def svfForProcessing655(
    dsm,
    vegdem,
    vegdem2,
    scale,
    usevegdem,
    feedback,
    device=torch.device("cpu"),
):

    with torch.no_grad():

        dsm = _to_tensor(dsm, device)
        vegdem = _to_tensor(vegdem, device)
        vegdem2 = _to_tensor(vegdem2, device)

        rows = dsm.shape[0]
        cols = dsm.shape[1]
        svf = torch.zeros([rows, cols], device=device)
        svfE = torch.zeros([rows, cols], device=device)
        svfS = torch.zeros([rows, cols], device=device)
        svfW = torch.zeros([rows, cols], device=device)
        svfN = torch.zeros([rows, cols], device=device)
        svfveg = torch.zeros((rows, cols), device=device)
        svfEveg = torch.zeros((rows, cols), device=device)
        svfSveg = torch.zeros((rows, cols), device=device)
        svfWveg = torch.zeros((rows, cols), device=device)
        svfNveg = torch.zeros((rows, cols), device=device)
        svfaveg = torch.zeros((rows, cols), device=device)
        svfEaveg = torch.zeros((rows, cols), device=device)
        svfSaveg = torch.zeros((rows, cols), device=device)
        svfWaveg = torch.zeros((rows, cols), device=device)
        svfNaveg = torch.zeros((rows, cols), device=device)

        # % amaxvalue
        vegmax = vegdem.max()
        amaxvalue = dsm.max()
        amaxvalue = torch.maximum(amaxvalue, vegmax)

        # % Elevation vegdems if buildingDSM inclused ground heights
        vegdem = vegdem + dsm
        vegdem[vegdem == dsm] = 0
        vegdem2 = vegdem2 + dsm
        vegdem2[vegdem2 == dsm] = 0
        # % Bush separation
        bush = torch.logical_not((vegdem2 * vegdem)) * vegdem

        noa = torch.tensor(19.0, device=device)
        # % No. of anglesteps minus 1
        step = 89.0 / noa
        iangle = torch.tensor(
            torch.hstack((torch.arange(step / 2.0, 89.0, step), 90.0)),
            device=device,
        )
        annulino = torch.tensor(
            torch.hstack((torch.round(torch.arange(0.0, 89.0, step)), 90.0)),
            device=device,
        )
        angleresult = svf_angles_100121(device)
        aziinterval = angleresult["aziinterval"].to(device)
        iazimuth = angleresult["iazimuth"].to(device)
        aziintervalaniso = torch.ceil((aziinterval / 2.0))
        index = 1

        for i in range(0, iangle.shape[0] - 1):
            for j in range(0, int(aziinterval[int(i)].item())):
                if feedback.isCanceled():
                    feedback.setProgressText("Calculation cancelled")
                    break
                altitude = float(iangle[int(i)].item())
                azimuth = float(iazimuth[int(index) - 1].item())

                # Casting shadow
                if usevegdem == 1:
                    shadowresult = shadow.shadowingfunction_20(
                        dsm,
                        vegdem,
                        vegdem2,
                        azimuth,
                        altitude,
                        scale,
                        amaxvalue,
                        bush,
                        feedback,
                        1,
                        device,
                    )

                    vegsh = torch.tensor(shadowresult["vegsh"], device=device)
                    vbshvegsh = torch.tensor(
                        shadowresult["vbshvegsh"], device=device
                    )
                    sh = torch.tensor(shadowresult["sh"], device=device)
                else:
                    sh = shadow.shadowingfunctionglobalradiation(
                        dsm, azimuth, altitude, scale, feedback, 1, device
                    )

                # Calculate svfs
                for k in range(
                    int(annulino[int(i)].item()) + 1,
                    int(annulino[int(i + 1.0)].item()) + 1,
                ):
                    weight = annulus_weight(k, aziinterval[i], device) * sh
                    svf = svf + weight
                    weight = (
                        annulus_weight(k, aziintervalaniso[i], device) * sh
                    )
                    if (azimuth >= 0) and (azimuth < 180):
                        svfE = svfE + weight
                    if (azimuth >= 90) and (azimuth < 270):
                        svfS = svfS + weight
                    if (azimuth >= 180) and (azimuth < 360):
                        svfW = svfW + weight
                    if (azimuth >= 270) or (azimuth < 90):
                        svfN = svfN + weight

                if usevegdem == 1:
                    for k in torch.arange(
                        annulino[int(i)] + 1, (annulino[int(i + 1.0)]) + 1
                    ):
                        # % changed to include 90
                        weight = annulus_weight(k, aziinterval[i], device)
                        svfveg = svfveg + weight * vegsh
                        svfaveg = svfaveg + weight * vbshvegsh
                        weight = annulus_weight(k, aziintervalaniso[i], device)
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
                feedback.setProgress(int(index * (100.0 / 655.0)))

        svfS = svfS + 3.0459e-004
        svfW = svfW + 3.0459e-004
        # % Last azimuth is 90. Hence, manual add of last annuli for svfS and SVFW
        # %Forcing svf not be greater than 1 (some MATLAB crazyness)
        svf[(svf > 1.0)] = 1.0
        svfE[(svfE > 1.0)] = 1.0
        svfS[(svfS > 1.0)] = 1.0
        svfW[(svfW > 1.0)] = 1.0
        svfN[(svfN > 1.0)] = 1.0

        if usevegdem == 1:
            last = torch.zeros((rows, cols), device=device)
            last[(vegdem2 == 0.0)] = 3.0459e-004
            svfSveg = svfSveg + last
            svfWveg = svfWveg + last
            svfSaveg = svfSaveg + last
            svfWaveg = svfWaveg + last
            # %Forcing svf not be greater than 1 (some MATLAB crazyness)
            svfveg[(svfveg > 1.0)] = 1.0
            svfEveg[(svfEveg > 1.0)] = 1.0
            svfSveg[(svfSveg > 1.0)] = 1.0
            svfWveg[(svfWveg > 1.0)] = 1.0
            svfNveg[(svfNveg > 1.0)] = 1.0
            svfaveg[(svfaveg > 1.0)] = 1.0
            svfEaveg[(svfEaveg > 1.0)] = 1.0
            svfSaveg[(svfSaveg > 1.0)] = 1.0
            svfWaveg[(svfWaveg > 1.0)] = 1.0
            svfNaveg[(svfNaveg > 1.0)] = 1.0

        svfresult = {
            "svf": svf,
            "svfE": svfE,
            "svfS": svfS,
            "svfW": svfW,
            "svfN": svfN,
            "svfveg": svfveg,
            "svfEveg": svfEveg,
            "svfSveg": svfSveg,
            "svfWveg": svfWveg,
            "svfNveg": svfNveg,
            "svfaveg": svfaveg,
            "svfEaveg": svfEaveg,
            "svfSaveg": svfSaveg,
            "svfWaveg": svfWaveg,
            "svfNaveg": svfNaveg,
        }
        return svfresult
