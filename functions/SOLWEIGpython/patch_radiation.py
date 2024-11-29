import numpy as np
from . import sunlit_shaded_patches

def shortwave_from_sky(sky, angle_of_incidence, lumChi, steradian, patch_azimuth, cyl):
    '''Calculates the amount of diffuse shortwave radiation from the sky for a patch with:
    angle of incidence = angle_of_incidence
    luminance = lumChi
    steradian = steradian'''

    # Diffuse vertical radiation
    diffuse_shortwave_radiation = sky * lumChi * angle_of_incidence * steradian

    return diffuse_shortwave_radiation

def longwave_from_sky(sky, Lsky_side, Lsky_down, patch_azimuth):

    # Degrees to radians
    deg2rad = np.pi / 180    

    # Longwave radiation from sky to vertical surface
    Ldown_sky = sky * Lsky_down

    # Longwave radiation from sky to horizontal surface
    Lside_sky = sky * Lsky_side

    #
    Least = np.zeros((sky.shape[0], sky.shape[1])) 
    Lsouth = np.zeros((sky.shape[0], sky.shape[1]))
    Lwest = np.zeros((sky.shape[0], sky.shape[1]))
    Lnorth = np.zeros((sky.shape[0], sky.shape[1]))

    # Portion into cardinal directions to be used for standing box or POI output
    if (patch_azimuth > 360) or (patch_azimuth < 180):
        Least = sky * Lsky_side * np.cos((90 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 90) and (patch_azimuth < 270):
        Lsouth = sky * Lsky_side * np.cos((180 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 180) and (patch_azimuth < 360):
        Lwest = sky * Lsky_side * np.cos((270 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 270) or (patch_azimuth < 90):
        Lnorth = sky * Lsky_side * np.cos((0 - patch_azimuth) * deg2rad)

    return Lside_sky, Ldown_sky, Least, Lsouth, Lwest, Lnorth

def longwave_from_veg(vegetation, steradian, angle_of_incidence, angle_of_incidence_h, patch_altitude, patch_azimuth, ewall, Ta):
    '''Calculates the amount of longwave radiation from vegetation for a patch with:
    angle of incidence = angle_of_incidence
    steradian = steradian
    if a patch is vegetation = vegetation
    amount of radiation from vegetated patch = vegetation_surface'''

    # Stefan-Boltzmann's Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi / 180    
    
    # Longwave radiation from vegetation surface (considered vertical)
    vegetation_surface = ((ewall * SBC * ((Ta + 273.15) ** 4)) / np.pi)

    # Longwave radiation reaching a vertical surface
    Lside_veg = vegetation_surface * steradian * angle_of_incidence * vegetation

    # Longwave radiation reaching a horizontal surface
    Ldown_veg = vegetation_surface * steradian * angle_of_incidence_h * vegetation

    #
    Least = np.zeros((vegetation.shape[0], vegetation.shape[1])) 
    Lsouth = np.zeros((vegetation.shape[0], vegetation.shape[1]))
    Lwest = np.zeros((vegetation.shape[0], vegetation.shape[1]))
    Lnorth = np.zeros((vegetation.shape[0], vegetation.shape[1]))

    # Portion into cardinal directions to be used for standing box or POI output
    if (patch_azimuth > 360) or (patch_azimuth < 180):
        Least = vegetation_surface * steradian * np.cos(patch_altitude * deg2rad) * vegetation * np.cos((90 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 90) and (patch_azimuth < 270):
        Lsouth = vegetation_surface * steradian * np.cos(patch_altitude * deg2rad) * vegetation * np.cos((180 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 180) and (patch_azimuth < 360):
        Lwest = vegetation_surface * steradian * np.cos(patch_altitude * deg2rad) * vegetation * np.cos((270 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 270) or (patch_azimuth < 90):
        Lnorth = vegetation_surface * steradian * np.cos(patch_altitude * deg2rad) * vegetation * np.cos((0 - patch_azimuth) * deg2rad)    

    return Lside_veg, Ldown_veg, Least, Lsouth, Lwest, Lnorth

def longwave_from_buildings(building, steradian, angle_of_incidence, angle_of_incidence_h, patch_azimuth, sunlit_patches, shaded_patches, azimuth_difference, solar_altitude, ewall, Ta, Tgwall):

    # Stefan-Boltzmann's Constant
    SBC = 5.67051e-8

    # Degrees to radians
    deg2rad = np.pi / 180    

    #
    Least = np.zeros((building.shape[0], building.shape[1])) 
    Lsouth = np.zeros((building.shape[0], building.shape[1]))
    Lwest = np.zeros((building.shape[0], building.shape[1]))
    Lnorth = np.zeros((building.shape[0], building.shape[1]))

    # Longwave radiation from sunlit surfaces
    sunlit_surface = ((ewall * SBC * ((Ta + Tgwall + 273.15) ** 4)) / np.pi)
    # Longwave radiation from shaded surfaces
    shaded_surface = ((ewall * SBC * ((Ta + 273.15) ** 4)) / np.pi)
    if ((azimuth_difference > 90) and (azimuth_difference < 270) and (solar_altitude > 0)):
        # Calculate which patches defined as buildings that are sunlit or shaded
        # sunlit_patches, shaded_patches = sunlit_shaded_patches.shaded_or_sunlit(solar_altitude, solar_azimuth, patch_altitude, patch_azimuth, asvf)
        
        # Calculate longwave radiation from sunlit walls to vertical surface
        Lside_sun = sunlit_surface * sunlit_patches * steradian * angle_of_incidence * building
        # Calculate longwave radiation from shaded walls to vertical surface
        Lside_sh = shaded_surface * shaded_patches * steradian * angle_of_incidence * building
        
        # Calculate longwave radiation from sunlit walls to horizontal surface
        Ldown_sun = sunlit_surface * sunlit_patches * steradian * angle_of_incidence_h * building
        # Calculate longwave radiation from shaded walls to horizontal surface
        Ldown_sh = shaded_surface * shaded_patches * steradian * angle_of_incidence_h * building
        
        # Portion into cardinal directions to be used for standing box or POI output
        if (patch_azimuth > 360) or (patch_azimuth < 180):
            Least = sunlit_surface * sunlit_patches * steradian * angle_of_incidence * building * np.cos((90 - patch_azimuth) * deg2rad)
            Least += shaded_surface * shaded_patches * steradian * angle_of_incidence * building * np.cos((90 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 90) and (patch_azimuth < 270):
            Lsouth = sunlit_surface * sunlit_patches * steradian * angle_of_incidence * building * np.cos((180 - patch_azimuth) * deg2rad)
            Lsouth += shaded_surface * shaded_patches * steradian * angle_of_incidence * building * np.cos((180 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 180) and (patch_azimuth < 360):
            Lwest = sunlit_surface * sunlit_patches * steradian * angle_of_incidence * building * np.cos((270 - patch_azimuth) * deg2rad)
            Lwest += shaded_surface * shaded_patches * steradian * angle_of_incidence * building * np.cos((270 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 270) or (patch_azimuth < 90):
            Lnorth = sunlit_surface * sunlit_patches * steradian * angle_of_incidence * building * np.cos((0 - patch_azimuth) * deg2rad)
            Lnorth += shaded_surface * shaded_patches * steradian * angle_of_incidence * building * np.cos((0 - patch_azimuth) * deg2rad)

    else:
        # Calculate longwave radiation from shaded walls reaching a vertical surface
        Lside_sh = shaded_surface * steradian * angle_of_incidence * building
        Lside_sun = np.zeros((Lside_sh.shape[0], Lside_sh.shape[1]))
        
        # Calculate longwave radiation from shaded walls reaching a horizontal surface
        Ldown_sh = shaded_surface * steradian * angle_of_incidence_h * building
        Ldown_sun = np.zeros((Lside_sh.shape[0], Lside_sh.shape[1]))

        # Portion into cardinal directions to be used for standing box or POI output
        if (patch_azimuth > 360) or (patch_azimuth < 180):
            Least = shaded_surface * steradian * angle_of_incidence * building * np.cos((90 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 90) and (patch_azimuth < 270):
            Lsouth = shaded_surface * steradian * angle_of_incidence * building * np.cos((180 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 180) and (patch_azimuth < 360):
            Lwest = shaded_surface * steradian * angle_of_incidence * building * np.cos((270 - patch_azimuth) * deg2rad)
        if (patch_azimuth > 270) or (patch_azimuth < 90):
            Lnorth = shaded_surface * steradian * angle_of_incidence * building * np.cos((0 - patch_azimuth) * deg2rad)

    return Lside_sun, Lside_sh, Ldown_sun, Ldown_sh, Least, Lsouth, Lwest, Lnorth

def longwave_from_buildings_wallScheme(voxelMaps, voxelTable, steradian, angle_of_incidence, angle_of_incidence_h, patch_azimuth):

    # Degrees to radians
    deg2rad = np.pi / 180    

    #
    Lside = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Lside_sh = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Ldown = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Ldown_sh = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Least = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1])) 
    Lsouth = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Lwest = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))
    Lnorth = np.zeros((voxelMaps.shape[0], voxelMaps.shape[1]))    
    
    # print(voxelMaps)
    #print(voxelTable.head())
    l = list(np.unique(voxelMaps)[1:])
    # print(l)
    a = dict(voxelTable.loc[l, 'LongwaveRadiation'])
    # print(a)
    patch_radiation = np.vectorize(a.get)(voxelMaps).astype(float)
    patch_radiation[np.isnan(patch_radiation)] = 0
    Lside += patch_radiation * steradian * angle_of_incidence
    Ldown += patch_radiation * steradian * angle_of_incidence_h

    # Portion into cardinal directions to be used for standing box or POI output
    if (patch_azimuth > 360) or (patch_azimuth < 180):
        Least = patch_radiation * steradian * angle_of_incidence * np.cos((90 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 90) and (patch_azimuth < 270):
        Lsouth = patch_radiation * steradian * angle_of_incidence * np.cos((180 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 180) and (patch_azimuth < 360):
        Lwest = patch_radiation * steradian * angle_of_incidence * np.cos((270 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 270) or (patch_azimuth < 90):
        Lnorth = patch_radiation * steradian * angle_of_incidence * np.cos((0 - patch_azimuth) * deg2rad)    

    return Lside, Lside_sh, Ldown, Ldown_sh, Least, Lsouth, Lwest, Lnorth

def reflected_longwave(reflecting_surface, steradian, angle_of_incidence, angle_of_incidence_h, patch_azimuth, Ldown_sky, Lup, ewall):

    # Degrees to radians
    deg2rad = np.pi / 180  

    # Calculate reflected longwave in each patch
    reflected_radiation = (((Ldown_sky+Lup) * (1-ewall)*0.5) / np.pi)

    # Reflected longwave radiation reaching vertical surfaces
    Lside_ref = reflected_radiation * steradian * angle_of_incidence * reflecting_surface
    
    # Reflected longwave radiation reaching horizontal surfaces
    Ldown_ref = reflected_radiation * steradian * angle_of_incidence_h * reflecting_surface

    #
    Least = np.zeros((reflecting_surface.shape[0], reflecting_surface.shape[1])) 
    Lsouth = np.zeros((reflecting_surface.shape[0], reflecting_surface.shape[1]))
    Lwest = np.zeros((reflecting_surface.shape[0], reflecting_surface.shape[1]))
    Lnorth = np.zeros((reflecting_surface.shape[0], reflecting_surface.shape[1]))

    # Portion into cardinal directions to be used for standing box or POI output
    if (patch_azimuth > 360) or (patch_azimuth < 180):
        Least = reflected_radiation * steradian * angle_of_incidence * reflecting_surface * np.cos((90 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 90) and (patch_azimuth < 270):
        Lsouth = reflected_radiation * steradian * angle_of_incidence * reflecting_surface * np.cos((180 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 180) and (patch_azimuth < 360):
        Lwest = reflected_radiation * steradian * angle_of_incidence * reflecting_surface * np.cos((270 - patch_azimuth) * deg2rad)
    if (patch_azimuth > 270) or (patch_azimuth < 90):
        Lnorth = reflected_radiation * steradian * angle_of_incidence * reflecting_surface * np.cos((0 - patch_azimuth) * deg2rad)

    return Lside_ref, Ldown_ref, Least, Lsouth, Lwest, Lnorth

def patch_steradians(L_patches):
    ''''This function calculates the steradians of the patches'''

    # Degrees to radians
    deg2rad = np.pi / 180  

    # Unique altitudes for patches
    skyalt, skyalt_c = np.unique(L_patches[:, 0], return_counts=True)

    # Altitudes of the Robinson & Stone patches
    patch_altitude = L_patches[:, 0]

    # Calculation of steradian for each patch
    steradian = np.zeros((patch_altitude.shape[0]))
    for i in range(patch_altitude.shape[0]):
        # If there are more than one patch in a band
        if skyalt_c[skyalt == patch_altitude[i]] > 1:
            steradian[i] = ((360 / skyalt_c[skyalt == patch_altitude[i]]) * deg2rad) * (np.sin((patch_altitude[i] + patch_altitude[0]) * deg2rad) \
            - np.sin((patch_altitude[i] - patch_altitude[0]) * deg2rad))
        # If there is only one patch in band, i.e. 90 degrees
        else:
            steradian[i] = ((360 / skyalt_c[skyalt == patch_altitude[i]]) * deg2rad) * (np.sin((patch_altitude[i]) * deg2rad) \
                - np.sin((patch_altitude[i-1] + patch_altitude[0]) * deg2rad))     
            
    return steradian, skyalt, patch_altitude