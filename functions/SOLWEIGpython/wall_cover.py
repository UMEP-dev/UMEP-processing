import numpy as np

def get_wall_cover(voxelTable, lcgrid, dsm, lc_params):
    '''Function to set thermal properties of wall used in surface temperature scheme of walls'''

    # Y-position of wall pixel in raster
    ypos = voxelTable['ypos'].to_numpy().astype(int)
    # X-position of wall pixel in raster
    xpos = voxelTable['xpos'].to_numpy().astype(int)
    # Empty array to store wall code
    wallCode = np.zeros((voxelTable.shape[0]))
    # Empty array to store thermal effusivity
    wallTu = np.zeros((voxelTable.shape[0]))
    # Empty array to store thermal conductivity
    wallTd = np.zeros((voxelTable.shape[0]))
    # Empty array to store albedo
    wallAlbedo = np.zeros((voxelTable.shape[0]))
    # Empty array to store emissivity
    wallEmissivity = np.zeros((voxelTable.shape[0]))
    # Empty array to store wall thickness
    wallThickness = np.zeros((voxelTable.shape[0]))    
    # Search kernel used to find wall code
    domain = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    # Loop through all wall pixels in voxel table
    for i in range(voxelTable.shape[0]):
        # Temporary lc_grid based on kernel and wall y and x position
        temp_lc = lcgrid[ypos[i]-1:ypos[i]+2, xpos[i]-1:xpos[i]+2] * domain
        # Temporary dsm based on kernel and wall y and x position
        temp_dsm = dsm[ypos[i]-1:ypos[i]+2, xpos[i]-1:xpos[i]+2] * domain
        # Temporary code based on highest pixel in temp_dsm where temp_lc is a building
        temp_code = temp_lc[((temp_lc > 99) & (temp_dsm == temp_dsm.max()))]

        # If more than one option in temp_code, use the first one #TODO CHANGE TO MOST COMMON
        if len(temp_code) > 1:
            temp_code = temp_code[0].astype(int)
        elif len(temp_code) == 1:
            temp_code = temp_code[0].astype(int)
        # If no wall type specified in land cover for specific wall pixel, set to concrete
        elif temp_code.size == 0:
            temp_code = 101

        # Save wall code in array
        wallCode[i] = temp_code
        # Volumetric heat capacity
        temp_Tc = lc_params['Specific_heat']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Thermal conductivity
        temp_Tk = lc_params['Thermal_conductivity']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Material density
        temp_D = lc_params['Density']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Calculate thermal effusivity
        temp_Tu = np.sqrt(temp_Tc * temp_D * temp_Tk)
        # Calculate thermal diffusivity
        temp_Td = temp_Tk/(temp_Tc*temp_D)
        # Save thermal effusivity of wall pixel in array
        wallTu[i] = temp_Tu
        # Save thermal diffusivity of wall pixel in array
        wallTd[i] = temp_Td
        # Get wall albedo
        temp_a = lc_params['Albedo']['Material']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Save wall albedo
        wallAlbedo[i] = temp_a
        # Get wall emissivity
        temp_e = lc_params['Emissivity']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Save wall emissivity
        wallEmissivity[i] = temp_e
        # Get wall thickness
        temp_wt = lc_params['Wall_thickness']['Value'][lc_params['Names']['Value'][str(temp_code)]]
        # Save wall thickness
        wallThickness[i] = temp_wt

    voxelTable['thermalEffusivity'] = wallTu
    voxelTable['thermalDiffusivity'] = wallTd
    voxelTable['wallAlbedo'] = wallAlbedo
    voxelTable['wallEmissivity'] = wallEmissivity
    voxelTable['wallThickness'] = wallThickness

    return voxelTable
