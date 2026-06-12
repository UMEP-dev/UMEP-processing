from __future__ import absolute_import
from .Kvikt_veg import Kvikt_veg

try:
    import torch

except:
    pass


def Kside_veg_v2022a(
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
    albedo,
    F_sh,
    KupE,
    KupS,
    KupW,
    KupN,
    cyl,
    lv,
    anisotropic_diffuse,
    diffsh,
    rows,
    cols,
    asvf,
    shmat,
    vegshmat,
    vbshvegshmat,
):
    device = (
        altitude.device
        if isinstance(altitude, torch.Tensor)
        else (
            azimuth.device
            if isinstance(azimuth, torch.Tensor)
            else torch.device(
                "cuda"
                if torch.cuda.is_available()
                else (
                    "xpu"
                    if (hasattr(torch, "xpu") and torch.xpu.is_available())
                    else "cpu"
                )
            )
        )
    )

    vikttot = 4.4897
    aziE = azimuth + t
    aziS = azimuth - 90 + t
    aziW = azimuth - 180 + t
    aziN = azimuth - 270 + t
    deg2rad = torch.pi / 180.0
    deg2rad = torch.tensor(deg2rad)

    KsideD = torch.zeros((rows, cols), device=device)
    Kref_sun = torch.zeros((rows, cols), device=device)
    Kref_sh = torch.zeros((rows, cols), device=device)
    Kref_veg = torch.zeros((rows, cols), device=device)
    Kside = torch.zeros((rows, cols), device=device)

    Kref_veg_n = torch.zeros((rows, cols), device=device)
    Kref_veg_s = torch.zeros((rows, cols), device=device)
    Kref_veg_e = torch.zeros((rows, cols), device=device)
    Kref_veg_w = torch.zeros((rows, cols), device=device)

    Kref_sh_n = torch.zeros((rows, cols), device=device)
    Kref_sh_s = torch.zeros((rows, cols), device=device)
    Kref_sh_e = torch.zeros((rows, cols), device=device)
    Kref_sh_w = torch.zeros((rows, cols), device=device)

    Kref_sun_n = torch.zeros((rows, cols), device=device)
    Kref_sun_s = torch.zeros((rows, cols), device=device)
    Kref_sun_e = torch.zeros((rows, cols), device=device)
    Kref_sun_w = torch.zeros((rows, cols), device=device)

    KeastRef = torch.zeros((rows, cols), device=device)
    KwestRef = torch.zeros((rows, cols), device=device)
    KnorthRef = torch.zeros((rows, cols), device=device)
    KsouthRef = torch.zeros((rows, cols), device=device)
    diffRadE = torch.zeros((rows, cols), device=device)
    diffRadS = torch.zeros((rows, cols), device=device)
    diffRadW = torch.zeros((rows, cols), device=device)
    diffRadN = torch.zeros((rows, cols), device=device)

    if cyl == 1:
        KsideI = shadow * radI * torch.cos(altitude * deg2rad)
        KeastI = torch.zeros((rows, cols), device=device)
        KsouthI = torch.zeros((rows, cols), device=device)
        KwestI = torch.zeros((rows, cols), device=device)
        KnorthI = torch.zeros((rows, cols), device=device)
    else:
        KeastI = torch.where(
            (azimuth > (360 - t)) | (azimuth <= (180 - t)),
            radI
            * shadow
            * torch.cos(altitude * deg2rad)
            * torch.sin(aziE * deg2rad),
            torch.zeros((rows, cols), device=device),
        )
        KsouthI = torch.where(
            (azimuth > (90 - t)) & (azimuth <= (270 - t)),
            radI
            * shadow
            * torch.cos(altitude * deg2rad)
            * torch.sin(aziS * deg2rad),
            torch.zeros((rows, cols), device=device),
        )
        KwestI = torch.where(
            (azimuth > (180 - t)) & (azimuth <= (360 - t)),
            radI
            * shadow
            * torch.cos(altitude * deg2rad)
            * torch.sin(aziW * deg2rad),
            torch.zeros((rows, cols), device=device),
        )
        KnorthI = torch.where(
            (azimuth <= (90 - t)) | (azimuth > (270 - t)),
            radI
            * shadow
            * torch.cos(altitude * deg2rad)
            * torch.sin(aziN * deg2rad),
            torch.zeros((rows, cols), device=device),
        )
        KsideI = shadow * 0

    viktveg, viktwall = Kvikt_veg(svfE, svfEveg, vikttot)
    svfviktbuvegE = viktwall + (viktveg * (1 - psi))

    viktveg, viktwall = Kvikt_veg(svfS, svfSveg, vikttot)
    svfviktbuvegS = viktwall + (viktveg * (1 - psi))

    viktveg, viktwall = Kvikt_veg(svfW, svfWveg, vikttot)
    svfviktbuvegW = viktwall + (viktveg * (1 - psi))

    viktveg, viktwall = Kvikt_veg(svfN, svfNveg, vikttot)
    svfviktbuvegN = viktwall + (viktveg * (1 - psi))

    if anisotropic_diffuse == 1:
        anisotropic_sky = True

        patch_altitude = lv[:, 0]
        patch_azimuth = lv[:, 1]
        if anisotropic_sky:
            patch_luminance = lv[:, 2]
        else:
            patch_luminance = (
                torch.ones((patch_altitude.shape[0]), device=device)
                / patch_altitude.shape[0]
            )

        skyalt, skyalt_c = torch.unique(patch_altitude, return_counts=True)
        radTot = torch.zeros(1, device=device)
        steradian = torch.zeros((patch_altitude.shape[0]), device=device)
        for i in range(patch_altitude.shape[0]):
            if skyalt_c[skyalt == patch_altitude[i]] > 1:
                steradian[i] = (
                    (360 / skyalt_c[skyalt == patch_altitude[i]]) * deg2rad
                ) * (
                    torch.sin(
                        (patch_altitude[i] + patch_altitude[0]) * deg2rad
                    )
                    - torch.sin(
                        (patch_altitude[i] - patch_altitude[0]) * deg2rad
                    )
                )
            else:
                steradian[i] = (
                    (360 / skyalt_c[skyalt == patch_altitude[i]]) * deg2rad
                ) * (
                    torch.sin((patch_altitude[i]) * deg2rad)
                    - torch.sin(
                        (patch_altitude[i - 1] + patch_altitude[0]) * deg2rad
                    )
                )

            radTot += (
                patch_luminance[i]
                * steradian[i]
                * torch.sin(patch_altitude[i] * deg2rad)
            )

        lumChi = (patch_luminance * radD) / radTot

        if cyl == 1:
            for idx in range(patch_azimuth.shape[0]):
                anglIncC = torch.cos(
                    patch_altitude[idx] * deg2rad
                ) * torch.cos(torch.tensor(0.0))
                KsideD += (
                    diffsh[:, :, idx] * lumChi[idx] * anglIncC * steradian[idx]
                )

                sunlit_surface = (
                    albedo * (radI * torch.cos(altitude * deg2rad))
                    + (radD * 0.5)
                ) / torch.pi
                shaded_surface = (albedo * radD * 0.5) / torch.pi

                temp_vegsh = (vegshmat[:, :, idx] == 0) | (
                    vbshvegshmat[:, :, idx] == 0
                )
                Kref_veg += (
                    shaded_surface
                    * temp_vegsh
                    * steradian[idx]
                    * torch.cos(patch_altitude[idx] * deg2rad)
                )

                temp_vbsh = (1 - shmat[:, :, idx]) * vbshvegshmat[:, :, idx]
                temp_sh = temp_vbsh == 1

                sunlit_patches, shaded_patches = shaded_or_sunlit(
                    altitude,
                    azimuth,
                    patch_altitude[idx],
                    patch_azimuth[idx],
                    asvf,
                )
                Kref_sun += (
                    sunlit_surface
                    * sunlit_patches
                    * temp_sh
                    * steradian[idx]
                    * torch.cos(patch_altitude[idx] * deg2rad)
                )
                Kref_sh += (
                    shaded_surface
                    * shaded_patches
                    * temp_sh
                    * steradian[idx]
                    * torch.cos(patch_altitude[idx] * deg2rad)
                )

            Kside = KsideI + KsideD + Kref_sun + Kref_sh + Kref_veg

            Keast = KupE * 0.5
            Kwest = KupW * 0.5
            Knorth = KupN * 0.5
            Ksouth = KupS * 0.5

        else:
            for idx in range(patch_azimuth.shape[0]):
                if (patch_azimuth[idx] > 360) or (patch_azimuth[idx] <= 180):
                    anglIncE = torch.cos(
                        patch_altitude[idx] * deg2rad
                    ) * torch.cos((90 - patch_azimuth[idx] + t) * deg2rad)
                    diffRadE += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncE
                        * steradian[idx]
                    )

                if (patch_azimuth[idx] > 90) and (patch_azimuth[idx] <= 270):
                    anglIncS = torch.cos(
                        patch_altitude[idx] * deg2rad
                    ) * torch.cos((180 - patch_azimuth[idx] + t) * deg2rad)
                    diffRadS += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncS
                        * steradian[idx]
                    )

                if (patch_azimuth[idx] > 180) and (patch_azimuth[idx] <= 360):
                    anglIncW = torch.cos(
                        patch_altitude[idx] * deg2rad
                    ) * torch.cos((270 - patch_azimuth[idx] + t) * deg2rad)
                    diffRadW += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncW
                        * steradian[idx]
                    )

                if (patch_azimuth[idx] > 270) or (patch_azimuth[idx] <= 90):
                    anglIncN = torch.cos(
                        patch_altitude[idx] * deg2rad
                    ) * torch.cos((0 - patch_azimuth[idx] + t) * deg2rad)
                    diffRadN += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncN
                        * steradian[idx]
                    )

                sunlit_surface = (
                    albedo * (radI * torch.cos(altitude * deg2rad))
                    + (radD * 0.5)
                ) / torch.pi
                shaded_surface = (albedo * radD * 0.5) / torch.pi

                temp_vegsh = (vegshmat[:, :, idx] == 0) | (
                    vbshvegshmat[:, :, idx] == 0
                )
                Kref_veg += (
                    shaded_surface
                    * temp_vegsh
                    * steradian[idx]
                    * torch.cos(patch_altitude[idx] * deg2rad)
                )

                if (patch_azimuth[idx] > 360) or (patch_azimuth[idx] < 180):
                    Kref_veg_e += (
                        shaded_surface
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                        * temp_vegsh
                        * torch.cos((90 - patch_azimuth[idx] + t) * deg2rad)
                    )
                if (patch_azimuth[idx] > 90) and (patch_azimuth[idx] < 270):
                    Kref_veg_s += (
                        shaded_surface
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                        * temp_vegsh
                        * torch.cos((180 - patch_azimuth[idx] + t) * deg2rad)
                    )
                if (patch_azimuth[idx] > 180) and (patch_azimuth[idx] < 360):
                    Kref_veg_w += (
                        shaded_surface
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                        * temp_vegsh
                        * torch.cos((270 - patch_azimuth[idx] + t) * deg2rad)
                    )
                if (patch_azimuth[idx] > 270) or (patch_azimuth[idx] < 90):
                    Kref_veg_n += (
                        shaded_surface
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                        * temp_vegsh
                        * torch.cos((0 - patch_azimuth[idx] + t) * deg2rad)
                    )

                temp_vbsh = (1 - shmat[:, :, idx]) * vbshvegshmat[:, :, idx]
                temp_sh = temp_vbsh == 1
                azimuth_difference = torch.abs(azimuth - patch_azimuth[idx])

                if (azimuth_difference > 90) and (azimuth_difference < 270):
                    sunlit_patches, shaded_patches = shaded_or_sunlit(
                        altitude,
                        azimuth,
                        patch_altitude[idx],
                        patch_azimuth[idx],
                        asvf,
                    )
                    Kref_sun += (
                        sunlit_surface
                        * sunlit_patches
                        * temp_sh
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                    )
                    Kref_sh += (
                        shaded_surface
                        * shaded_patches
                        * temp_sh
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                    )

                    if (patch_azimuth[idx] > 360) or (
                        patch_azimuth[idx] < 180
                    ):
                        Kref_sun_e += (
                            sunlit_surface
                            * sunlit_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (90 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                        Kref_sh_e += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (90 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 90) and (
                        patch_azimuth[idx] < 270
                    ):
                        Kref_sun_s += (
                            sunlit_surface
                            * sunlit_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (180 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                        Kref_sh_s += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (180 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 180) and (
                        patch_azimuth[idx] < 360
                    ):
                        Kref_sun_w += (
                            sunlit_surface
                            * sunlit_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (270 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                        Kref_sh_w += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (270 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 270) or (patch_azimuth[idx] < 90):
                        Kref_sun_n += (
                            sunlit_surface
                            * sunlit_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((0 - patch_azimuth[idx] + t) * deg2rad)
                        )
                        Kref_sh_n += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((0 - patch_azimuth[idx] + t) * deg2rad)
                        )
                else:
                    Kref_sh += (
                        shaded_surface
                        * temp_sh
                        * steradian[idx]
                        * torch.cos(patch_altitude[idx] * deg2rad)
                    )

                    if (patch_azimuth[idx] > 360) or (
                        patch_azimuth[idx] < 180
                    ):
                        Kref_sh_e += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (90 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 90) and (
                        patch_azimuth[idx] < 270
                    ):
                        Kref_sh_s += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (180 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 180) and (
                        patch_azimuth[idx] < 360
                    ):
                        Kref_sh_w += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos(
                                (270 - patch_azimuth[idx] + t) * deg2rad
                            )
                        )
                    if (patch_azimuth[idx] > 270) or (patch_azimuth[idx] < 90):
                        Kref_sh_n += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((0 - patch_azimuth[idx] + t) * deg2rad)
                        )

            Keast = (
                KeastI
                + diffRadE
                + Kref_sun_e
                + Kref_sh_e
                + Kref_veg_e
                + KupE * 0.5
            )
            Kwest = (
                KwestI
                + diffRadW
                + Kref_sun_w
                + Kref_sh_w
                + Kref_veg_w
                + KupW * 0.5
            )
            Knorth = (
                KnorthI
                + diffRadN
                + Kref_sun_n
                + Kref_sh_n
                + Kref_veg_n
                + KupN * 0.5
            )
            Ksouth = (
                KsouthI
                + diffRadS
                + Kref_sun_s
                + Kref_sh_s
                + Kref_veg_s
                + KupS * 0.5
            )

    else:
        KeastDG = (
            radD * (1 - svfviktbuvegE)
            + albedo * (svfviktbuvegE * (radG * (1 - F_sh) + radD * F_sh))
            + KupE
        ) * 0.5
        Keast = KeastI + KeastDG

        KsouthDG = (
            radD * (1 - svfviktbuvegS)
            + albedo * (svfviktbuvegS * (radG * (1 - F_sh) + radD * F_sh))
            + KupS
        ) * 0.5
        Ksouth = KsouthI + KsouthDG

        KwestDG = (
            radD * (1 - svfviktbuvegW)
            + albedo * (svfviktbuvegW * (radG * (1 - F_sh) + radD * F_sh))
            + KupW
        ) * 0.5
        Kwest = KwestI + KwestDG

        KnorthDG = (
            radD * (1 - svfviktbuvegN)
            + albedo * (svfviktbuvegN * (radG * (1 - F_sh) + radD * F_sh))
            + KupN
        ) * 0.5
        Knorth = KnorthI + KnorthDG

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside
