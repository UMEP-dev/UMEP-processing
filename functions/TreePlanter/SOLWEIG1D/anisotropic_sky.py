import numpy as np
import pandas as pd

from copy import deepcopy

from . import emissivity_models
from . import patch_radiation
from . import sunlit_shaded_patches
from . import patch_radiation

def anisotropic_sky(skyp, buip, vegp, solar_altitude, solar_azimuth, asvf, cyl, 
                    esky, L_patches, wallScheme, voxelTable, voxelMaps, steradians, Ta, Tgwall, ewall, Lup,
                    radI, radD, radG, lv, albedo, anisotropic_diffuse, diffsh, shadow, Kup):
    
    '''This function estimates short- and longwave radiation from an anisotropic sky....'''
    # Stefan-Boltzmann's Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi / 180

    # Define shortwave variables
    KsideD = 0
    KsideI = 0
    Kside = 0
    Kref_sun = 0
    Kref_sh = 0
    Kref_veg = 0

    Keast = 0
    Kwest = 0
    Knorth = 0
    Ksouth = 0

    # Define longwave variables
    Ldown_sky = 0
    Ldown_veg = 0
    Ldown_sun = 0
    Ldown_sh = 0
    Ldown_ref = 0

    Lside_sky = 0
    Lside_veg = 0
    Lside_sun = 0
    Lside_sh = 0
    Lside_ref = 0

    Least = 0
    Lwest = 0
    Lnorth = 0
    Lsouth = 0

    # Unique altitudes for patches
    skyalt, skyalt_c = np.unique(L_patches[:, 0], return_counts=True)

    # Patch altitudes and azimuths
    patch_altitude = L_patches[:, 0]
    patch_azimuth = L_patches[:, 1]

    # Longwave radiation
    # Current model for anisotropic sky emissivity
    emis_m = 2

    # Unsworth & Monteith (1975)
    if emis_m == 1:
        patch_emissivity_normalized, esky_band = emissivity_models.model1(L_patches, esky, Ta)
    # Martin & Berdahl (1984)
    elif emis_m == 2:
        patch_emissivity_normalized, esky_band = emissivity_models.model2(L_patches, esky, Ta)
    # Bliss (1961)
    elif emis_m == 3:
        patch_emissivity_normalized, esky_band = emissivity_models.model3(L_patches, esky, Ta)

    # True = anisotropic sky, False = isotropic sky
    anisotropic_sky_ = True

    # Longwave based on spectral flux density (divide by pi)
    Ldown = np.zeros((patch_altitude.shape[0]))
    Lside = np.zeros((patch_altitude.shape[0]))
    Lnormal = np.zeros((patch_altitude.shape[0]))
    for temp_altitude in skyalt:
        # Anisotropic sky
        if anisotropic_sky_:
            temp_emissivity = esky_band[skyalt == temp_altitude]
        # Isotropic sky but with patches (need to switch anisotropic_sky to False)
        else:
            temp_emissivity = esky
        # Estimate longwave radiation on a horizontal surface (Ldown), vertical surface (Lside) and perpendicular (Lnormal)
        Ldown[patch_altitude == temp_altitude] = ((temp_emissivity * SBC * ((Ta + 273.15) ** 4)) / np.pi) * steradians[patch_altitude == temp_altitude] * np.sin(temp_altitude * deg2rad)
        Lside[patch_altitude == temp_altitude] = ((temp_emissivity * SBC * ((Ta + 273.15) ** 4)) / np.pi) * steradians[patch_altitude == temp_altitude] * np.cos(temp_altitude * deg2rad)
        Lnormal[patch_altitude == temp_altitude] = ((temp_emissivity * SBC * ((Ta + 273.15) ** 4)) / np.pi) * steradians[patch_altitude == temp_altitude]

    Lsky_normal = deepcopy(L_patches)
    Lsky_down = deepcopy(L_patches)
    Lsky_side = deepcopy(L_patches)

    Lsky_normal[:,2] = Lnormal
    Lsky_down[:,2] = Ldown
    Lsky_side[:,2] = Lside

    # Shortwave radiation
    if solar_altitude > 0:
        # Patch luminance
        patch_luminance = lv[:, 2]
        # Variable to fill with total sky diffuse shortwave radiation
        radTot = np.zeros(1)
        # Radiance fraction normalization
        for i in np.arange(patch_altitude.shape[0]):
            radTot += (patch_luminance[i] * steradians[i] * np.sin(patch_altitude[i] * deg2rad)) # Radiance fraction normalization
        lumChi = (patch_luminance * radD) / radTot

    for i in np.arange(patch_altitude.shape[0]):
        # Calculations for patches on sky, shmat = 1 = sky is visible
        # temp_sky = ((shmat[:, :, i] == 1) & (vegshmat[:, :, i] == 1))
        
        # Calculations for patches that are vegetation, vegshmat = 0 = shade from vegetation
        # temp_vegsh = ((vegshmat[:, :, i] == 0) | (vbshvegshmat[:, :, i] == 0))

        # Calculations for patches that are buildings, shmat = 0 = shade from buildings
        # temp_vbsh = (1 - shmat[:, :, i]) * vbshvegshmat[:, :, i]
        # temp_sh = (temp_vbsh == 1)
        # if wallScheme == 1:
        #     # temp_vbsh = (voxelMaps[:, :, i] > 0) * vbshvegshmat[:, :, i]
        #     # temp_sh_w = (temp_vbsh == 1) * voxelMaps[:, :, i]
        #     temp_sh_w = temp_sh * voxelMaps[:, :, i]
        #     temp_sh_roof = temp_sh * (voxelMaps[:, :, i] == 0)

        # Estimate sunlit and shaded patches
        sunlit_patches, shaded_patches = sunlit_shaded_patches.shaded_or_sunlit(solar_altitude, solar_azimuth, patch_altitude[i], patch_azimuth[i], asvf)

        if cyl == 1:
            # Angle of incidence, np.cos(0) because cylinder - always perpendicular
            angle_of_incidence = np.cos(patch_altitude[i] * deg2rad) * np.cos(0) # * np.sin(np.pi / 2) \
            # Angle of incidence to horizontal surface
            angle_of_incidence_h = np.sin(patch_altitude[i] * deg2rad)

            ### CALCULATIONS FOR LONGWAVE RADIATION ###
            # Longwave radiation from sky
            Lside_sky_temp, Ldown_sky_temp, Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_sky(skyp[i], Lsky_side[i, 2], Lsky_down[i, 2], patch_azimuth[i])          
            
            Lside_sky += Lside_sky_temp
            Ldown_sky += Ldown_sky_temp
            Least += Least_temp
            Lsouth += Lsouth_temp
            Lwest += Lwest_temp
            Lnorth += Lnorth_temp

            # Longwave radiation from vegetation
            Lside_veg_temp, Ldown_veg_temp, \
                Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_veg(vegp[i], steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                    patch_altitude[i], patch_azimuth[i], ewall, Ta)                   
            
            Lside_veg += Lside_veg_temp
            Ldown_veg += Ldown_veg_temp
            Least += Least_temp
            Lsouth += Lsouth_temp
            Lwest += Lwest_temp
            Lnorth += Lnorth_temp            

            # Longwave radiation from buildings
            # if wallScheme == 0:
            azimuth_difference = np.abs(solar_azimuth - patch_azimuth[i])

            Lside_sun_temp, Lside_sh_temp, \
            Ldown_sun_temp, Ldown_sh_temp, \
            Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_buildings(buip[i], steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                        patch_azimuth[i], sunlit_patches, shaded_patches, 
                                                                                                        azimuth_difference, solar_altitude, ewall, Ta, Tgwall)

            # else:
            #     azimuth_difference = np.abs(solar_azimuth - patch_azimuth[i])
            #     # print('Building pixels = ' + str(temp_sh.sum()))
            #     # print('Building pixels wall scheme = ' + str(np.sum(temp_sh_w > 0)))
            #     Lside_sun_temp, Lside_sh_temp, \
            #     Ldown_sun_temp, Ldown_sh_temp, \
            #     Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_buildings_wallScheme(temp_sh_w, voxelTable, steradians[i], angle_of_incidence, angle_of_incidence_h, 
            #                                                                                              patch_azimuth[i])                

            #     Lside_sun_r_temp, Lside_sh_r_temp, \
            #     Ldown_sun_r_temp, Ldown_sh_r_temp, \
            #     Least_r_temp, Lsouth_r_temp, Lwest_r_temp, Lnorth_r_temp = patch_radiation.longwave_from_buildings(temp_sh_roof, steradians[i], angle_of_incidence, angle_of_incidence_h, 
            #                                                                                              patch_azimuth[i], sunlit_patches, shaded_patches, 
            #                                                                                              azimuth_difference, solar_altitude, ewall, Ta, Tgwall)

            #     Lside_sun_temp += Lside_sun_r_temp
            #     Lside_sh_temp += Lside_sh_r_temp
            #     Ldown_sun_temp += Ldown_sun_r_temp
            #     Ldown_sh_temp += Ldown_sh_r_temp
            #     Least_temp += Least_r_temp
            #     Lsouth_temp += Lsouth_r_temp
            #     Lwest_temp += Lwest_r_temp
            #     Lnorth_temp += Lnorth_r_temp

            Lside_sun += Lside_sun_temp
            Lside_sh += Lside_sh_temp
            Ldown_sun += Ldown_sun_temp
            Ldown_sh += Ldown_sh_temp
            Least += Least_temp
            Lsouth += Lsouth_temp
            Lwest += Lwest_temp
            Lnorth += Lnorth_temp 

            ### CALCULATIONS FOR SHORTWAVE RADIATION ###
            if solar_altitude > 0:
                # Shortwave radiation from sky
                KsideD += skyp[i] * lumChi[i] * angle_of_incidence * steradians[i]
                
                if skyp[i]:
                    Keast_temp, Ksouth_temp, Kwest_temp, Knorth_temp = patch_radiation.cardinal_shortwave(lumChi[i], angle_of_incidence, skyp[i], steradians[i], 1, patch_azimuth[i], patch_altitude[i])

                # Shortwave reflected on sunlit surfaces
                # sunlit_surface = ((albedo * radG) / np.pi)
                sunlit_surface = ((albedo * (radI * np.cos(solar_altitude * deg2rad)) + (radD * 0.5)) / np.pi)
                # Shortwave reflected on shaded surfaces and vegetation
                shaded_surface = ((albedo * radD * 0.5) / np.pi)
                
                # Shortwave radiation from vegetation
                Kref_veg += shaded_surface * vegp[i] * steradians[i] * angle_of_incidence
                
                if vegp[i]:
                    Keast_temp, Ksouth_temp, Kwest_temp, Knorth_temp = patch_radiation.cardinal_shortwave(shaded_surface, angle_of_incidence, vegp[i], steradians[i], 1, patch_azimuth[i], patch_altitude[i])                

                # Shortwave radiation from buildings
                Kref_sun += sunlit_surface * sunlit_patches * buip[i] * steradians[i] * angle_of_incidence
                Kref_sh += shaded_surface * shaded_patches * buip[i] * steradians[i] * angle_of_incidence

                if buip[i]:
                    if sunlit_patches == 1:
                        Keast_temp, Ksouth_temp, Kwest_temp, Knorth_temp = patch_radiation.cardinal_shortwave(sunlit_surface, angle_of_incidence, buip[i], steradians[i], sunlit_patches, patch_azimuth[i], patch_altitude[i])
                    elif shaded_patches == 1:
                        Keast_temp, Ksouth_temp, Kwest_temp, Knorth_temp = patch_radiation.cardinal_shortwave(shaded_surface, angle_of_incidence, buip[i], steradians[i], shaded_patches, patch_azimuth[i], patch_altitude[i])

                Keast += Keast_temp
                Kwest += Kwest_temp
                Knorth += Knorth_temp
                Ksouth += Ksouth_temp

    # Calculate reflected longwave in each patch
    for idx in np.arange(patch_altitude.shape[0]):
        # Angle of incidence, np.cos(0) because cylinder - always perpendicular
        angle_of_incidence = np.cos(patch_altitude[idx] * deg2rad) * np.cos(0) # * np.sin(np.pi / 2) \

        # Angle of incidence to horizontal surface
        angle_of_incidence_h = np.sin(patch_altitude[idx] * deg2rad)

        # Patches with reflected surfaces
        temp_sh = ((buip[idx] == 1) | (vegp[idx] == 1))
        
        Lside_ref_temp, Ldown_ref_temp, \
        Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.reflected_longwave(temp_sh, steradians[idx], angle_of_incidence, angle_of_incidence_h, patch_azimuth[idx], Ldown_sky, Lup, ewall)

        Lside_ref += Lside_ref_temp
        Ldown_ref += Ldown_ref_temp
        Least += Least_temp
        Lsouth += Lsouth_temp
        Lwest += Lwest_temp
        Lnorth += Lnorth_temp

    # Sum of all Lside components (sky, vegetation, sunlit and shaded buildings, reflected)
    Lside = Lside_sky + Lside_veg + Lside_sh + Lside_sun + Lside_ref

    # Sum of all Lside components (sky, vegetation, sunlit and shaded buildings, reflected)
    Ldown = Ldown_sky + Ldown_veg + Ldown_sh + Ldown_sun + Ldown_ref

    ### Direct radiation ###
    if cyl == 1: ### Kside with cylinder ###
        KsideI = shadow * radI * np.cos(solar_altitude * deg2rad)

    if solar_altitude > 0:

        Kside = KsideI + KsideD + Kref_sun + Kref_sh + Kref_veg

        Keast_ = (Kup * 0.5)
        Kwest_ = (Kup * 0.5)
        Knorth_ = (Kup * 0.5)
        Ksouth_ = (Kup * 0.5)

        if cyl == 0:
            Keast_ += Keast
            Kwest_ += Kwest
            Knorth_ += Knorth
            Ksouth_ += Ksouth

    return Ldown, Lside, Lside_sky, Lside_veg, Lside_sh, Lside_sun, Lside_ref, Least, Lwest, Lnorth, Lsouth, \
           Keast_, Ksouth_, Kwest_, Knorth_, KsideI, KsideD, Kside, steradians, skyalt

    