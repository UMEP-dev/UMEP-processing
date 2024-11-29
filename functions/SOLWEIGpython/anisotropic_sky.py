import numpy as np
import pandas as pd

from copy import deepcopy

from . import emissivity_models
from . import patch_radiation
from . import sunlit_shaded_patches
from . import patch_radiation

def anisotropic_sky(shmat, vegshmat, vbshvegshmat, solar_altitude, solar_azimuth, asvf, cyl, 
                    esky, L_patches, wallScheme, voxelTable, voxelMaps, steradians, Ta, Tgwall, ewall, Lup,
                    radI, radD, radG, lv, albedo, anisotropic_diffuse, diffsh, shadow, KupE, KupS, KupW, KupN, current_step):
    
    '''This function estimates short- and longwave radiation from an anisotropic sky....'''
    # Stefan-Boltzmann's Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi / 180

    # Shape of rasters
    rows = shmat.shape[0]
    cols = shmat.shape[1]

    # Define shortwave variables
    KsideD = np.zeros((rows, cols))
    KsideI = np.zeros((rows, cols))
    Kside = np.zeros((rows, cols))
    Kref_sun = np.zeros((rows, cols))
    Kref_sh = np.zeros((rows, cols))
    Kref_veg = np.zeros((rows, cols))

    Keast = np.zeros((rows, cols))
    Kwest = np.zeros((rows, cols))
    Knorth = np.zeros((rows, cols))
    Ksouth = np.zeros((rows, cols))

    # Define longwave variables
    Ldown_sky = np.zeros((rows, cols))
    Ldown_veg = np.zeros((rows, cols))
    Ldown_sun = np.zeros((rows, cols))
    Ldown_sh = np.zeros((rows, cols))
    Ldown_ref = np.zeros((rows, cols))

    Lside_sky = np.zeros((rows, cols))
    Lside_veg = np.zeros((rows, cols))
    Lside_sun = np.zeros((rows, cols))
    Lside_sh = np.zeros((rows, cols))
    Lside_ref = np.zeros((rows, cols))

    Least = np.zeros((rows, cols))
    Lwest = np.zeros((rows, cols))
    Lnorth = np.zeros((rows, cols))
    Lsouth = np.zeros((rows, cols))

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
        temp_sky = ((shmat[:, :, i] == 1) & (vegshmat[:, :, i] == 1))
        
        # Calculations for patches that are vegetation, vegshmat = 0 = shade from vegetation
        temp_vegsh = ((vegshmat[:, :, i] == 0) | (vbshvegshmat[:, :, i] == 0))

        # Calculations for patches that are buildings, shmat = 0 = shade from buildings
        temp_vbsh = (1 - shmat[:, :, i]) * vbshvegshmat[:, :, i]
        temp_sh = (temp_vbsh == 1)
        if wallScheme == 1:
            # temp_vbsh = (voxelMaps[:, :, i] > 0) * vbshvegshmat[:, :, i]
            # temp_sh_w = (temp_vbsh == 1) * voxelMaps[:, :, i]
            temp_sh_w = temp_sh * voxelMaps[:, :, i]
            temp_sh_roof = temp_sh * (voxelMaps[:, :, i] == 0)

        # Estimate sunlit and shaded patches
        sunlit_patches, shaded_patches = sunlit_shaded_patches.shaded_or_sunlit(solar_altitude, solar_azimuth, patch_altitude[i], patch_azimuth[i], asvf)

        if cyl == 1:
            # Angle of incidence, np.cos(0) because cylinder - always perpendicular
            angle_of_incidence = np.cos(patch_altitude[i] * deg2rad) * np.cos(0) # * np.sin(np.pi / 2) \
            # Angle of incidence to horizontal surface
            angle_of_incidence_h = np.sin(patch_altitude[i] * deg2rad)

            ### CALCULATIONS FOR LONGWAVE RADIATION ###
            # Longwave radiation from sky
            Lside_sky_temp, Ldown_sky_temp, Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_sky(temp_sky, Lsky_side[i, 2], Lsky_down[i, 2], patch_azimuth[i])          
            
            Lside_sky += Lside_sky_temp
            Ldown_sky += Ldown_sky_temp
            Least += Least_temp
            Lsouth += Lsouth_temp
            Lwest += Lwest_temp
            Lnorth += Lnorth_temp

            # Longwave radiation from vegetation
            Lside_veg_temp, Ldown_veg_temp, \
                Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_veg(temp_vegsh, steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                    patch_altitude[i], patch_azimuth[i], ewall, Ta)                   
            
            Lside_veg += Lside_veg_temp
            Ldown_veg += Ldown_veg_temp
            Least += Least_temp
            Lsouth += Lsouth_temp
            Lwest += Lwest_temp
            Lnorth += Lnorth_temp            

            # Longwave radiation from buildings
            if wallScheme == 0:
                azimuth_difference = np.abs(solar_azimuth - patch_azimuth[i])

                Lside_sun_temp, Lside_sh_temp, \
                Ldown_sun_temp, Ldown_sh_temp, \
                Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_buildings(temp_sh, steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                         patch_azimuth[i], sunlit_patches, shaded_patches, 
                                                                                                         azimuth_difference, solar_altitude, ewall, Ta, Tgwall)

            else:
                azimuth_difference = np.abs(solar_azimuth - patch_azimuth[i])
                # print('Building pixels = ' + str(temp_sh.sum()))
                # print('Building pixels wall scheme = ' + str(np.sum(temp_sh_w > 0)))
                Lside_sun_temp, Lside_sh_temp, \
                Ldown_sun_temp, Ldown_sh_temp, \
                Least_temp, Lsouth_temp, Lwest_temp, Lnorth_temp = patch_radiation.longwave_from_buildings_wallScheme(temp_sh_w, voxelTable, steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                         patch_azimuth[i])                

                Lside_sun_r_temp, Lside_sh_r_temp, \
                Ldown_sun_r_temp, Ldown_sh_r_temp, \
                Least_r_temp, Lsouth_r_temp, Lwest_r_temp, Lnorth_r_temp = patch_radiation.longwave_from_buildings(temp_sh_roof, steradians[i], angle_of_incidence, angle_of_incidence_h, 
                                                                                                         patch_azimuth[i], sunlit_patches, shaded_patches, 
                                                                                                         azimuth_difference, solar_altitude, ewall, Ta, Tgwall)

                Lside_sun_temp += Lside_sun_r_temp
                Lside_sh_temp += Lside_sh_r_temp
                Ldown_sun_temp += Ldown_sun_r_temp
                Ldown_sh_temp += Ldown_sh_r_temp
                Least_temp += Least_r_temp
                Lsouth_temp += Lsouth_r_temp
                Lwest_temp += Lwest_r_temp
                Lnorth_temp += Lnorth_r_temp

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
                KsideD += temp_sky * lumChi[i] * angle_of_incidence * steradians[i]
                
                # Shortwave reflected on sunlit surfaces
                # sunlit_surface = ((albedo * radG) / np.pi)
                sunlit_surface = ((albedo * (radI * np.cos(solar_altitude * deg2rad)) + (radD * 0.5)) / np.pi)
                # Shortwave reflected on shaded surfaces and vegetation
                shaded_surface = ((albedo * radD * 0.5) / np.pi)
                
                # Shortwave radiation from vegetation
                Kref_veg += shaded_surface * temp_vegsh * steradians[i] * angle_of_incidence
                
                # Shortwave radiation from buildings
                Kref_sun += sunlit_surface * sunlit_patches * temp_sh * steradians[i] * angle_of_incidence
                Kref_sh += shaded_surface * shaded_patches * temp_sh * steradians[i] * angle_of_incidence

    # Calculate reflected longwave in each patch
    for idx in np.arange(patch_altitude.shape[0]):
        # Angle of incidence, np.cos(0) because cylinder - always perpendicular
        angle_of_incidence = np.cos(patch_altitude[idx] * deg2rad) * np.cos(0) # * np.sin(np.pi / 2) \

        # Angle of incidence to horizontal surface
        angle_of_incidence_h = np.sin(patch_altitude[idx] * deg2rad)

        # Patches with reflected surfaces
        temp_sh = ((shmat[:,:,idx] == 0) | (vegshmat[:,:,idx] == 0) | (vbshvegshmat[:,:,idx] == 0))
        
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

        Keast = (KupE * 0.5)
        Kwest = (KupW * 0.5)
        Knorth = (KupN * 0.5)
        Ksouth = (KupS * 0.5)

    return Ldown, Lside, Lside_sky, Lside_veg, Lside_sh, Lside_sun, Lside_ref, Least, Lwest, Lnorth, Lsouth, \
           Keast, Ksouth, Kwest, Knorth, KsideI, KsideD, Kside, steradians, skyalt

    