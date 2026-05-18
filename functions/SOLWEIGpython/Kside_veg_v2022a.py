from __future__ import absolute_import
import numpy as np
from .Kvikt_veg import Kvikt_veg
from . import sunlit_shaded_patches
import torch


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
        else azimuth.device
        if isinstance(azimuth, torch.Tensor)
        else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    )

    # New reflection equation 2012-05-25
    vikttot = 4.4897
    aziE = azimuth + t
    aziS = azimuth - 90 + t
    aziW = azimuth - 180 + t
    aziN = azimuth - 270 + t
    deg2rad = torch.pi / 180
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

    ### Direct radiation ###
    if cyl == 1:  ### Kside with cylinder ###
        KsideI = shadow * radI * torch.cos(altitude * deg2rad)
        KeastI = 0
        KsouthI = 0
        KwestI = 0
        KnorthI = 0
    else:  ### Kside with weights ###
        if azimuth > (360 - t) or azimuth <= (180 - t):
            KeastI = (
                radI
                * shadow
                * torch.cos(altitude * deg2rad)
                * torch.sin(aziE * deg2rad)
            )
        else:
            KeastI = 0
        if azimuth > (90 - t) and azimuth <= (270 - t):
            KsouthI = (
                radI
                * shadow
                * torch.cos(altitude * deg2rad)
                * torch.sin(aziS * deg2rad)
            )
        else:
            KsouthI = 0
        if azimuth > (180 - t) and azimuth <= (360 - t):
            KwestI = (
                radI
                * shadow
                * torch.cos(altitude * deg2rad)
                * torch.sin(aziW * deg2rad)
            )
        else:
            KwestI = 0
        if azimuth <= (90 - t) or azimuth > (270 - t):
            KnorthI = (
                radI
                * shadow
                * torch.cos(altitude * deg2rad)
                * torch.sin(aziN * deg2rad)
            )
        else:
            KnorthI = 0

        KsideI = shadow * 0

    ### Diffuse and reflected radiation ###
    [viktveg, viktwall] = Kvikt_veg(svfE, svfEveg, vikttot)
    svfviktbuvegE = viktwall + (viktveg) * (1 - psi)

    [viktveg, viktwall] = Kvikt_veg(svfS, svfSveg, vikttot)
    svfviktbuvegS = viktwall + (viktveg) * (1 - psi)

    [viktveg, viktwall] = Kvikt_veg(svfW, svfWveg, vikttot)
    svfviktbuvegW = viktwall + (viktveg) * (1 - psi)

    [viktveg, viktwall] = Kvikt_veg(svfN, svfNveg, vikttot)
    svfviktbuvegN = viktwall + (viktveg) * (1 - psi)

    ### Anisotropic Diffuse Radiation after Perez et al. 1993 ###
    if anisotropic_diffuse == 1:

        anisotropic_sky = True

        patch_altitude = lv[:, 0]
        patch_azimuth = lv[:, 1]
        if anisotropic_sky:
            patch_luminance = lv[:, 2]
        else:
            patch_luminance = torch.zeros((patch_altitude.shape[0]), device=device)
            patch_luminance[:] = 1.0 / patch_luminance.shape[0]

        # Unique altitudes for patches
        skyalt, skyalt_c = torch.unique(patch_altitude, return_counts=True)

        radTot = torch.zeros(1, device=device)

        # Calculation of steradian for each patch
        steradian = torch.zeros((patch_altitude.shape[0]), device=device)
        for i in range(patch_altitude.shape[0]):
            # If there are more than one patch in a band
            if skyalt_c[skyalt == patch_altitude[i]] > 1:
                steradian[i] = (
                    (360 / skyalt_c[skyalt == patch_altitude[i]]) * deg2rad
                ) * (
                    torch.sin((patch_altitude[i] + patch_altitude[0]) * deg2rad)
                    - torch.sin((patch_altitude[i] - patch_altitude[0]) * deg2rad)
                )
            # If there is only one patch in band, i.e. 90 degrees
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
            )  # Radiance fraction normalization

        lumChi = (
            patch_luminance * radD
        ) / radTot  # Radiance fraction normalization

        if cyl == 1:
            for idx in range(patch_azimuth.shape[0]):
                # Angle of incidence, torch.cos(0) because cylinder - always perpendicular
                anglIncC = torch.cos(patch_altitude[idx] * deg2rad) * torch.cos(
                    0
                )  # * torch.sin(torch.pi / 2) \
                # + torch.sin(patch_altitude[idx] * deg2rad) * torch.cos(torch.pi / 2)
                # Diffuse vertical radiation
                KsideD += (
                    diffsh[:, :, idx] * lumChi[idx] * anglIncC * steradian[idx]
                )

                # Shortwave reflected on sunlit surfaces
                # sunlit_surface = ((albedo * radG) / torch.pi)
                sunlit_surface = (
                    albedo * (radI * torch.cos(altitude * deg2rad)) + (radD * 0.5)
                ) / torch.pi
                # Shortwave reflected on shaded surfaces and vegetation
                shaded_surface = (albedo * radD * 0.5) / torch.pi

                # Shortwave radiation reflected on vegetation - based on diffuse shortwave radiation
                temp_vegsh = (vegshmat[:, :, idx] == 0) | (
                    vbshvegshmat[:, :, idx] == 0
                )
                Kref_veg += (
                    shaded_surface
                    * temp_vegsh
                    * steradian[idx]
                    * torch.cos(patch_altitude[idx] * deg2rad)
                )

                # Shortwave radiation reflected on buildings (shaded and sunlit) - based on global and diffuse shortwave radiation
                temp_vbsh = (1 - shmat[:, :, idx]) * vbshvegshmat[:, :, idx]
                temp_sh = temp_vbsh == 1  # & (vbshvegshmat[:,:,idx] == 1)

                sunlit_patches, shaded_patches = (
                    sunlit_shaded_patches.shaded_or_sunlit(
                        altitude,
                        azimuth,
                        patch_altitude[idx],
                        patch_azimuth[idx],
                        asvf,
                    )
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

            # Keast  = (albedo * (svfviktbuvegE * (radG * (1 - F_sh) + radD * F_sh)) + KupE) * 0.5
            # Ksouth = (albedo * (svfviktbuvegS * (radG * (1 - F_sh) + radD * F_sh)) + KupS) * 0.5
            # Kwest  = (albedo * (svfviktbuvegW * (radG * (1 - F_sh) + radD * F_sh)) + KupW) * 0.5
            # Knorth = (albedo * (svfviktbuvegN * (radG * (1 - F_sh) + radD * F_sh)) + KupN) * 0.5
        else:  # Box
            diffRadE = torch.zeros((rows, cols), device=device)
            diffRadS = torch.zeros((rows, cols), device=device)
            diffRadW = torch.zeros((rows, cols), device=device)
            diffRadN = torch.zeros((rows, cols), device=device)

            for idx in range(patch_azimuth.shape[0]):
                if (patch_azimuth[idx] > 360) or (patch_azimuth[idx] <= 180):
                    anglIncE = torch.cos(patch_altitude[idx] * deg2rad) * torch.cos(
                        (90 - patch_azimuth[idx] + t) * deg2rad
                    )  # * torch.sin(torch.pi / 2) \
                    # + torch.sin(patch_altitude[idx] * deg2rad) * torch.cos(torch.pi / 2)
                    diffRadE += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncE
                        * steradian[idx]
                    )  # * 0.5

                if (patch_azimuth[idx] > 90) and (patch_azimuth[idx] <= 270):
                    anglIncS = torch.cos(patch_altitude[idx] * deg2rad) * torch.cos(
                        (180 - patch_azimuth[idx] + t) * deg2rad
                    )  # * torch.sin(torch.pi / 2) \
                    # + torch.sin(patch_altitude[idx] * deg2rad) * torch.cos(torch.pi / 2)
                    diffRadS += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncS
                        * steradian[idx]
                    )  # * 0.5

                if (patch_azimuth[idx] > 180) and (patch_azimuth[idx] <= 360):
                    anglIncW = torch.cos(patch_altitude[idx] * deg2rad) * torch.cos(
                        (270 - patch_azimuth[idx] + t) * deg2rad
                    )  # * torch.sin(torch.pi / 2) \
                    # + torch.sin(patch_altitude[idx] * deg2rad) * torch.cos(torch.pi / 2)
                    diffRadW += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncW
                        * steradian[idx]
                    )  # * 0.5

                if (patch_azimuth[idx] > 270) or (patch_azimuth[idx] <= 90):
                    anglIncN = torch.cos(patch_altitude[idx] * deg2rad) * torch.cos(
                        (0 - patch_azimuth[idx] + t) * deg2rad
                    )  # * torch.sin(torch.pi / 2) \
                    # + torch.sin(patch_altitude[idx] * deg2rad) * torch.cos(torch.pi / 2)
                    diffRadN += (
                        diffsh[:, :, idx]
                        * lumChi[idx]
                        * anglIncN
                        * steradian[idx]
                    )  # * 0.5

                # Shortwave reflected on sunlit surfaces
                # sunlit_surface = ((albedo * radG) / torch.pi)
                sunlit_surface = (
                    albedo * (radI * torch.cos(altitude * deg2rad)) + (radD * 0.5)
                ) / torch.pi
                # Shortwave reflected on shaded surfaces and vegetation
                shaded_surface = (albedo * radD * 0.5) / torch.pi

                # Shortwave radiation reflected on vegetation - based on diffuse shortwave radiation
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

                # Shortwave radiation reflected on buildings (shaded and sunlit) - based on global and diffuse shortwave radiation
                temp_vbsh = (1 - shmat[:, :, idx]) * vbshvegshmat[:, :, idx]
                temp_sh = temp_vbsh == 1  # & (vbshvegshmat[:,:,idx] == 1)
                azimuth_difference = torch.abs(azimuth - patch_azimuth[idx])

                if (azimuth_difference > 90) and (azimuth_difference < 270):
                    sunlit_patches, shaded_patches = (
                        sunlit_shaded_patches.shaded_or_sunlit(
                            altitude,
                            azimuth,
                            patch_altitude[idx],
                            patch_azimuth[idx],
                            asvf,
                        )
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
                            * torch.cos((90 - patch_azimuth[idx] + t) * deg2rad)
                        )
                        Kref_sh_e += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((90 - patch_azimuth[idx] + t) * deg2rad)
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
                            * torch.cos((180 - patch_azimuth[idx] + t) * deg2rad)
                        )
                        Kref_sh_s += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((180 - patch_azimuth[idx] + t) * deg2rad)
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
                            * torch.cos((270 - patch_azimuth[idx] + t) * deg2rad)
                        )
                        Kref_sh_w += (
                            shaded_surface
                            * shaded_patches
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((270 - patch_azimuth[idx] + t) * deg2rad)
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
                            * torch.cos((90 - patch_azimuth[idx] + t) * deg2rad)
                        )
                    if (patch_azimuth[idx] > 90) and (
                        patch_azimuth[idx] < 270
                    ):
                        Kref_sh_s += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((180 - patch_azimuth[idx] + t) * deg2rad)
                        )
                    if (patch_azimuth[idx] > 180) and (
                        patch_azimuth[idx] < 360
                    ):
                        Kref_sh_w += (
                            shaded_surface
                            * steradian[idx]
                            * torch.cos(patch_altitude[idx] * deg2rad)
                            * temp_sh
                            * torch.cos((270 - patch_azimuth[idx] + t) * deg2rad)
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

    return Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside
    # return Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kref_sh, Kref_sun, Kref_veg, Kside
