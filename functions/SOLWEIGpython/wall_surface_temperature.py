import numpy as np
import pandas as pd
import math
from ...functions.SOLWEIGpython.wall_cover import get_wall_cover
from .cylindric_wedge import cylindric_wedge_voxel
from .Lvikt_veg import Lvikt_veg

# Stefan Boltzmans Constant
SBC = 5.67051e-8

def load_walls(voxelTable, solweig_parameters, wall_type, wallaspect, Ta, timeStep, albedo_b, emissivity_b, albedo_grid, landcover, lcgrid, dsm):
    '''This function loads the voxel data created in sky view factor calculator and turns it into a Pandas DataFrame.'''
    
    # Load data as Pandas DataFrame and add column names
    voxelTable = pd.DataFrame(voxelTable, columns = ['voxelId', 'voxelHeight', 'wallHeight', 'wallHeight_exact', 'wallId', 'ypos', 'xpos', 
                                                        'SVF_height', 'SVF', 'SVF_fix', 'svfbu', 'svfveg', 'svfaveg'])
    
    # Initiate/declare new columns used by SOLWEIG and parameterization scheme for wall surface temperatures
    voxelTable['wallTemperature'] = Ta # Initial wall temperature is set to air temperature
    voxelTable['timeStep'] = timeStep
    
    # Add columns filled later
    columns_to_add = ['SVF_ground', 'svfbu_at_ground', 'svfaveg_at_ground', 'wallAspect', 'wallEmissivity', 'wallThickness', 
                      'wallAlbedo', 'thermalEffusivity', 'thermalDiffusivity', 'groundAlbedo', 'wallShade', 'wallShadeHeight', 
                      'LongwaveRadiation', 'K_in', 'L_in', 'Lwallsun', 'Lwallsh', 'Lrefl', 'Lveg', 'Lground', 'Lsky', 'esky', 'voxelHeightMasl']
    for col in columns_to_add:
        voxelTable[col] = 0
    
    # tmp = svf + svfveg - 1.
    # tmp[tmp < 0.] = 0.
    # %matlab crazyness around 0
    # svfalfa = np.arcsin(np.exp((np.log((1. - tmp)) / 2.)))
    tmp = voxelTable['SVF_fix'].to_numpy() + voxelTable['svfveg'].to_numpy() - 1.
    tmp[tmp < 0.] = 0.
    voxelTable['svfalfa'] = np.arcsin(np.exp((np.log((1. - tmp)) / 2.)))

    # Add wall aspect of each voxel
    # temp_table = np.column_stack([voxelTable['voxelId'].to_numpy(), voxelTable['ypos'].to_numpy(), voxelTable['xpos'].to_numpy()])
    temp_table = np.column_stack([voxelTable['wallId'].to_numpy(), voxelTable['ypos'].to_numpy(), voxelTable['xpos'].to_numpy()])
    temp_table = np.unique(temp_table, axis = 0)
    # Add wall aspect
    for i in np.arange(temp_table.shape[0]):
        temp_aspect = wallaspect[temp_table[i,1].astype(int), temp_table[i,2].astype(int)]
        voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'wallAspect'] = temp_aspect
        temp_building = voxelTable.loc[((voxelTable['wallId'] == temp_table[i,0]) & (voxelTable['voxelHeight'] == voxelTable['voxelHeight'].min())), 'svfbu'].copy().to_numpy()[0]
        temp_veg = voxelTable.loc[((voxelTable['wallId'] == temp_table[i,0]) & (voxelTable['voxelHeight'] == voxelTable['voxelHeight'].min())), 'svfaveg'].copy().to_numpy()[0]
        temp_albedo = albedo_grid[temp_table[i,1].astype(int), temp_table[i,2].astype(int)]
        temp_dsm = dsm[temp_table[i,1].astype(int), temp_table[i,2].astype(int)]
        
        voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'svfbu_at_ground'] = temp_building
        voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'svfaveg_at_ground'] = temp_veg
        voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'groundAlbedo'] = temp_albedo
        voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'voxelHeightMasl'] = voxelTable.loc[voxelTable['wallId'] == temp_table[i,0], 'voxelHeight'].to_numpy() + temp_dsm
    
    # Set voxelId to index
    voxelTable = voxelTable.set_index('voxelId')

    # Calculate fractions
    # Non-sky fraction = (1 - svf at ground level) - 0.5 (ground surface seen from wall)
    building_fraction = 1 - voxelTable['svfbu_at_ground'].to_numpy() - 0.5
    building_fraction[building_fraction < 0] = 0.
    veg_fraction = 1 - voxelTable['svfaveg_at_ground'].to_numpy() - 0.5
    veg_fraction[veg_fraction < 0] = 0.
    voxelTable['building_fraction'] = building_fraction
    voxelTable['veg_fraction'] = veg_fraction
    sky_fraction = voxelTable['SVF_fix'].to_numpy()
    ground_fraction = 1 - sky_fraction - building_fraction - veg_fraction
    voxelTable['ground_fraction'] = ground_fraction

    voxelTable['total_fraction'] = building_fraction + sky_fraction + ground_fraction + veg_fraction

    # Set thermal effusivity according to wall type (either from GUI or from landcover)
    if landcover == 1:
        unique_landcover = np.unique(lcgrid)
        # Unique wall codes for unique wall types, e.g. brick, concrete, wood, etc.
        unique_walls = unique_landcover[unique_landcover > 99].astype(int)
        # Get wall properties for wall surface temperature scheme
        if unique_walls.size > 1:
            voxelTable = get_wall_cover(voxelTable, lcgrid, dsm, solweig_parameters)
        elif unique_walls.size == 1:
            # Specific heat capacity
            wallTc = solweig_parameters['Specific_heat']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]
            # Thermal conductivity
            wallTk = solweig_parameters['Thermal_conductivity']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]
            # Material density
            wallD = solweig_parameters['Density']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]
            # Thermal effusivity
            wallTu = np.sqrt(wallTc * wallD * wallTk)
            # Thermal diffusivity
            wallTd = wallTk/(wallTc*wallD)
            # Set thermal effusivity
            voxelTable['thermalEffusivity'] = wallTu      
            # Set thermal diffusivity
            voxelTable['thermalDiffusivity'] = wallTd
            # Set wall albedo
            # voxelTable['wallAlbedo'] = solweig_parameters['Albedo']['Material']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]   
            voxelTable['wallAlbedo'] = albedo_b
            # Set wall emissivity
            # voxelTable['wallEmissivity'] = solweig_parameters['Emissivity']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]
            voxelTable['wallEmissivity'] = emissivity_b
            # Set thickness
            voxelTable['wallThickness'] = solweig_parameters['Wall_thickness']['Value'][solweig_parameters['Names']['Value'][str(unique_walls[0])]]
        else:
            landcover = 0
    
    if landcover == 0:
        # Specific heat capacity
        wallTc = solweig_parameters['Specific_heat']['Value'][solweig_parameters['Names']['Value'][wall_type]]
        # Thermal conductivity
        wallTk = solweig_parameters['Thermal_conductivity']['Value'][solweig_parameters['Names']['Value'][wall_type]]
        # Material density
        wallD = solweig_parameters['Density']['Value'][solweig_parameters['Names']['Value'][wall_type]]
        # Thermal effusivity
        wallTu = np.sqrt(wallTc * wallD * wallTk)    
        # Set thermal effusivity
        voxelTable['thermalEffusivity'] = wallTu
        # Calculate thermal diffusivity
        wallTd = wallTk/(wallTc*wallD)
        # Set thermal diffusivity
        voxelTable['thermalDiffusivity'] = wallTd
        # Get wall albedo
        # voxelTable['wallAlbedo'] = solweig_parameters['Albedo']['Material']['Value'][solweig_parameters['Names']['Value'][wall_type]]
        voxelTable['wallAlbedo'] = albedo_b
        # Get wall emissivity
        # voxelTable['wallEmissivity'] = solweig_parameters['Emissivity']['Value'][solweig_parameters['Names']['Value'][wall_type]]
        voxelTable['wallEmissivity'] = emissivity_b
        # Get wall thickness
        voxelTable['wallThickness'] = solweig_parameters['Wall_thickness']['Value'][solweig_parameters['Names']['Value'][wall_type]]

    eqTime = True

    ### REMEMEBER TO TURN OFF FOR KOLUMBUS ###
    if eqTime:
        voxelTable['timeStep'] = voxelTable['wallThickness'].to_numpy()**2/(np.pi**2 * voxelTable['thermalDiffusivity'].to_numpy())

    return voxelTable, wallaspect

def step_heating(q, e, t):
    '''Function to calculate delta surface temperature based on heat flux (q), thermal effusivity (e) and time (t)'''
    return (2*q)/e*np.sqrt(t/np.pi)

def surface_temperature_calc(effusivity, t, Kin, Lin, Ta, wall_emissivity, Ts_previous):
    '''Function to get surface temperature'''
       
    dT = np.zeros((effusivity.shape[0]))
    Ts = np.zeros((effusivity.shape[0]))
    # Calculate heat flux based on Ts in previous time step and Kin, Lin and Ta from current time step
    Lout_temp = wall_emissivity * SBC * (Ts_previous + 273.15)**4
    # sensible_temp = 20 * (Ts_previous - Ta)
    energy_in_temp = Kin + Lin - Lout_temp #- sensible_temp
    # Calculate dT (Ts - Ta) with step wise heating
    dT = step_heating(energy_in_temp, effusivity, t)

    # Surface temperature
    Ts_current = Ta + dT
    # Estimate surface temperature with Ts of current timestep, i.e. not Ts_previous
    Lout_temp = wall_emissivity * SBC * (Ts_current + 273.15)**4
    energy_in_temp = Kin + Lin - Lout_temp #- sensible_temp
    dT = step_heating(energy_in_temp, effusivity, t)
    # Surface temperature after second iteration
    Ts = Ta + dT
        
    return Ts, dT

def wall_surface_temperature(voxelTable, wallsh, altitude, azimuth, timeStep, K_direct, K_diff, K_down, Ldown, Lup, Ta, esky):
    '''Wall surface temperature parameterization
    
    This parameterization scheme estimates a wall temperature based on the 
    incoming energy (short- and longwave radiation), outgoing longwave radiation and sensible heat.'''

    # Stefan Boltzmans Constant
    SBC = 5.67051e-8

    # Degrees fo radians conversion factor
    deg2rad = np.pi/180

    # Estimate shadow on wall, i.e. how far up the wall that the shade stretches (or how far down the wall that the sun reaches)
    voxelTable['wallShade'] = 0 # Starting value is zero, i.e. the entire wall is shaded
    # If sun is above horizon, get wall shadow height on wall from wallsh calculated in shadowingfunction_wallheight_23/13 in Solweig_2022a_calc_forprocessing.py
    if altitude > 0:
        # Starting value of wallShadeHeight
        voxelTable['wallShadeHeight'] = 0
        # Find number of unique walls in model domain
        unique_walls = np.unique(voxelTable['wallId'])
        # Find wall shade height of each unique wall
        for unique_wall in unique_walls:
            # temp_sh is the shadow height of a wall
            temp_sh = wallsh[voxelTable.loc[voxelTable.wallId == unique_wall, 'ypos'].to_numpy().astype(int)[0], voxelTable.loc[voxelTable.wallId == unique_wall, 'xpos'].to_numpy().astype(int)[0]]
            
            # Change to sunlit (1) for voxels that are above wall shade height
            voxelTable.loc[(voxelTable['wallId'] == unique_wall) & (voxelTable['voxelHeight'] >= temp_sh) & (voxelTable['wallHeight_exact'] > temp_sh), 'wallShade'] = 1
            # Add wall shade height to pandas dataframe
            voxelTable.loc[(voxelTable['wallId'] == unique_wall), 'wallShadeHeight'] = temp_sh
    # If sun is below horizon, everything is in "shade"
    else:
        voxelTable['wallShadeHeight'] = voxelTable['wallHeight_exact']

    # Ldown and Lup for wall pixels used when estimating the amount of longwave received from surrounding surfaces (ground and reflected)
    Ldown_array = np.zeros((voxelTable.shape[0]))
    Lup_array = np.zeros((voxelTable.shape[0]))
    for idx in np.arange(voxelTable.shape[0]):
        temp_y = voxelTable.iloc[idx]['ypos'].astype(int)
        temp_x = voxelTable.iloc[idx]['xpos'].astype(int)
        Ldown_array[idx] = Ldown[temp_y, temp_x]
        Lup_array[idx] = Lup[temp_y, temp_x]

    # Estimate longwave in based on old methods by Bogren et al.
    # If sun above horizon
    if altitude > 0:
        # Cylindric wedge based on building height angle from svf and solar zenith angle
        F_sh = cylindric_wedge_voxel((90 - altitude) * deg2rad, voxelTable['svfalfa'].to_numpy())  # Fraction shadow on building walls based on sun alt and svf
        F_sh[np.isnan(F_sh)] = 0.5

        # Only half of hemisphere visible from a wall, therefore the entire seen area can be sunlit, compared to an open space where only half can be sunlit
        F_sh = 2. * F_sh - 1.  #(cylindric_wedge scaled 0-1)

        # Fraction of sunlit and shaded building surfaces seen by a building wall based on its aspect and the azimuth of the sun
        wallSun = np.abs(voxelTable['wallAspect'].to_numpy() - azimuth)
        wallSun = wallSun/180.
        wallSun[wallSun > 1.] = 2 - (wallSun[wallSun > 1.])
        wallSun = 0.2 + wallSun * 0.6 # Scaling wall shadows to avoid totally shaded and totally sunlit walls

        # Temperature in shade and sun based on mean values for all voxels, i.e. shade = mean of all shaded and sun = mean of all sunlit
        ts_shade = voxelTable.loc[voxelTable['wallShade'] == 0, 'wallTemperature'].to_numpy().mean()
        ts_sun = voxelTable.loc[voxelTable['wallShade'] == 1, 'wallTemperature'].to_numpy().mean()

        # Calculation of longwave radiation received by the voxels based on sunlit and shaded surroundings
        Lwallsun = (SBC * voxelTable['wallEmissivity'].to_numpy() * ((ts_sun + 273.15)**4) * voxelTable['building_fraction'].to_numpy() * (1. - F_sh)) * wallSun # (1 - svf - 0.5) istället för viktwall?
        Lwallsh = (SBC * voxelTable['wallEmissivity'].to_numpy() * ((ts_shade + 273.15)**4) * voxelTable['building_fraction'].to_numpy() * (1. - F_sh)) * (1 - wallSun) # (1 - svf - 0.5) istället för viktwall?
        
        Lwallsh += SBC * voxelTable['wallEmissivity'].to_numpy() * ((ts_shade + 273.15)**4) * voxelTable['building_fraction'].to_numpy() * F_sh # * 0.5 # (1 - svf - 0.5) istället för viktwall?

    # If sun below horizon
    else:
        Lwallsun = 0.
        ts_shade = voxelTable.loc[voxelTable['wallShade'] == 0, 'wallTemperature'].to_numpy().mean()
        Lwallsh = SBC * voxelTable['wallEmissivity'].to_numpy() * (ts_shade + 273.15)**4 * voxelTable['building_fraction'].to_numpy() # * 0.5 # (1 - svf - 0.5) istället för viktwall?

    # Received longwave radiation from vegetation
    Lveg = SBC * voxelTable['wallEmissivity'].to_numpy() * (Ta + 273.15)**4 * voxelTable['veg_fraction'].to_numpy() # * 0.5 # (1 - svf - 0.5) istället för viktwall?
    # Received longwave radiation from the sky
    Lsky = SBC * esky * ((Ta + 273.15)**4) * voxelTable['SVF_fix'].to_numpy() # * 0.5 # svf istället för viktwall?
    # Received reflected longwave radiation
    Lrefl = (1. - voxelTable['wallEmissivity'].to_numpy()) * (Ldown_array + Lup_array) * voxelTable['building_fraction'].to_numpy() # * 0.5 # (1 - svf - 0.5) istället för viktwall?
    # Received longwave radiation from ground
    Lground = Lup_array * voxelTable['ground_fraction'].to_numpy()
    
    # Total amount of longwave radiation received by a wall surface
    L_in = Lwallsun + Lwallsh + Lrefl + Lveg + Lground + Lsky #TODO ground??
    
    voxelTable['Lwallsun'] = Lwallsun
    voxelTable['Lwallsh'] = Lwallsh
    voxelTable['Lrefl'] = Lrefl
    voxelTable['Lveg'] = Lveg
    voxelTable['Lground'] = Lground
    voxelTable['Lsky'] = Lsky
    voxelTable['esky'] = esky
    if altitude > 0:
        voxelTable['F_sh'] = F_sh
        voxelTable['wallSun'] = wallSun
    else:
        voxelTable['F_sh'] = 0.
        voxelTable['wallSun'] = 0.

    # Angle of incidence (wall aspect vs solar position (altitude and azimuth))
    sun_x = 1 * np.cos(math.radians(altitude)) * np.cos(math.radians(azimuth)) * np.cos(deg2rad * voxelTable['wallAspect'].to_numpy())
    sun_y = 1 * np.cos(math.radians(altitude)) * np.sin(math.radians(azimuth)) * np.sin(deg2rad * voxelTable['wallAspect'].to_numpy())
    cf = sun_x + sun_y
    cf[cf <= 0] = 0.

    # If sun above horizon
    if altitude > 0:
        K_in = (1 - voxelTable['wallAlbedo'].to_numpy()) * (K_direct * cf * voxelTable['wallShade'].to_numpy() + 
                                                      K_diff * voxelTable['SVF_fix'].to_numpy() + 
                                                       K_down * voxelTable['wallAlbedo'].to_numpy() * voxelTable['building_fraction'] + 
                                                       (K_down * voxelTable['groundAlbedo']) * voxelTable['ground_fraction'])
    # If sun below horizon
    else:
        K_in = 0.    

    voxelTable['K_in'] = K_in
    voxelTable['L_in'] = L_in

    # Calculate wall surface temperature with step-wise heating
    Ts_previous = voxelTable['wallTemperature'].to_numpy()
    if voxelTable['timeStep'].unique().size == 1:
        timeStep = voxelTable.iloc[0]['timeStep']
    else:
        timeStep = voxelTable['timeStep'].to_numpy()
    Ts, dT = surface_temperature_calc(voxelTable['thermalEffusivity'].to_numpy(), timeStep, K_in, L_in, Ta, voxelTable['wallEmissivity'].to_numpy(), Ts_previous)
    voxelTable['wallTemperature'] = Ts
        
    # Convert wall temperature to longwave radiation
    voxelTable['LongwaveRadiation'] = ((voxelTable['wallEmissivity'] * SBC * ((voxelTable['wallTemperature'] + 273.15) ** 4)) / np.pi)

    # Position of the sun for current timestep (can be removed?)
    voxelTable['sunAltitude'] = altitude
    voxelTable['sunAzimuth'] = azimuth

    return voxelTable