# This is the main function of the SOLWEIG model
# 2025-Mar-21
# Fredrik Lindberg, fredrikl@gvc.gu.se
# Goteborg Urban Climate Group
# Gothenburg University

#sommon imports
from __future__ import absolute_import
from ...util.umep_solweig_export_component import read_solweig_config
from ...util.SEBESOLWEIGCommonFiles.Solweig_v2015_metdata_noload import Solweig_2015a_metdata_noload
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import clearnessindex_2013b

from ...functions.SOLWEIGpython import Solweig_2025a_calc_forprocessing as so
#from ...functions.SOLWEIGpython import WriteMetadataSOLWEIG # Not needed anymore?
from ...functions.SOLWEIGpython import PET_calculations as p
from ...functions.SOLWEIGpython import UTCI_calculations as utci
from ...functions.SOLWEIGpython.CirclePlotBar import PolarBarPlot
from ...functions.SOLWEIGpython.wall_surface_temperature import load_walls
from ...functions.SOLWEIGpython.wallOfInterest import wallOfInterest
from ...functions.SOLWEIGpython.wallsAsNetCDF import walls_as_netcdf
from ...functions.SOLWEIGpython.Tgmaps_v1 import Tgmaps_v1

import numpy as np
import json
import zipfile
import pandas as pd
import matplotlib.pylab as plt
from shutil import copyfile

# imports for osgeo/qgis dependency
try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    from ...util.misc import saveraster
    from qgis.core import QgsVectorLayer
except:
    pass

# imports for standalone
try:
    from umep import common
    from rasterio.transform import xy
    import pyproj
    from tqdm import tqdm
except:
    pass


# import numpy as np
# from .daylen import daylen
# from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b import clearnessindex_2013b
# from ...util.SEBESOLWEIGCommonFiles.diffusefraction import diffusefraction
# from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_13 import shadowingfunction_wallheight_13
# from ...util.SEBESOLWEIGCommonFiles.shadowingfunction_wallheight_23 import shadowingfunction_wallheight_23
# from .gvf_2018a import gvf_2018a
# from .cylindric_wedge import cylindric_wedge
# from .TsWaveDelay_2015a import TsWaveDelay_2015a
# from .Kup_veg_2015a import Kup_veg_2015a
# # from .Lside_veg_v2015a import Lside_veg_v2015a
# # from .Kside_veg_v2019a import Kside_veg_v2019a
# from .Kside_veg_v2022a import Kside_veg_v2022a
# from ...util.SEBESOLWEIGCommonFiles.Perez_v3 import Perez_v3
# from ...util.SEBESOLWEIGCommonFiles.create_patches import create_patches

def solweig_run(configPath, feedback):

    """
    Input:
    configPath : config file including geodata paths and settings.
    feedback : To communicate with qgis gui. Set to None if standalone
    """
    
    # Load config file
    configDict = read_solweig_config(configPath)

    # Load parameters settings for SOLWEIG
    with open(configDict['para_json_path'], "r") as jsn:
        param = json.load(jsn)

    standAlone = int(configDict['standalone'])

    #Load DSM
    if standAlone == 1: #TODO How to do this?
        dsm, dsm_transf, dsm_crs = common.load_raster(configDict['filepath_dsm'], bbox=None) 
        scale = 1 / dsm_transf.a
        # dsm_height, dsm_width = dsm.shape  # y rows by x cols
        # y is flipped - so return max for lower row
        minx, miny = xy(dsm_transf, dsm.shape[0], 0)
        # Define the source and target CRS
        source_crs = pyproj.CRS(dsm_crs)
        target_crs = pyproj.CRS(4326)  # WGS 84
        # Create a transformer object
        transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
        # Perform the transformation
        lon, lat = transformer.transform(minx, miny) 
        nd = -9999 #TODO: extract nodatavalue from rasterio
    else:
        gdal_dsm = gdal.Open(configDict['filepath_dsm'])
        dsm = gdal_dsm.ReadAsArray().astype(float)
        lat = float(configDict['lat'])
        lon = float(configDict['lon']) 
        scale = float(configDict['scale'])
        nd = gdal_dsm.GetRasterBand(1).GetNoDataValue()

    rows = dsm.shape[0]
    cols = dsm.shape[1]
    
    # response to issue #85
    dsm[dsm == nd] = 0.
    if dsm.min() < 0:
        dsmraise = np.abs(dsm.min())
        dsm = dsm + dsmraise
    else:
        dsmraise = 0

    alt = np.median(dsm)
    if alt < 0:
        alt = 3

    #Vegetation
    transVeg = param["Tree_settings"]["Value"]["Transmissivity"] 
    trunkratio = param["Tree_settings"]["Value"]["Trunk_ratio"] 
    usevegdem = int(configDict['usevegdem'])
    if usevegdem == 1:
        if standAlone == 0:
            vegdsm = gdal.Open(configDict['filepath_cdsm']).ReadAsArray().astype(float)
        else:
            vegdsm, _ , _ = common.load_raster(configDict['filepath_cdsm'], bbox=None)
        if configDict['filepath_tdsm'] is None:
            if standAlone == 1:
                vegdsm2 = gdal.Open(configDict['filepath_tdsm']).ReadAsArray().astype(float)
            else:
                vegdsm2, _ , _ = common.load_raster(configDict['filepath_tdsm'], bbox=None)
        else:
            vegdsm2 = vegdsm * trunkratio

    #Land cover
    landcover = int(configDict['landcover'])
    if landcover == 1:
        if standAlone == 0:
            lcgrid = gdal.Open(configDict['filepath_lc']).ReadAsArray().astype(float)
        else:
            lcgrid, _ , _ = common.load_raster(configDict['filepath_lc'], bbox=None)
    else:
        lcgrid = 0

    #DEM for buildings #TODO: fix nodata in standalone
    demforbuild = int(configDict['demforbuild'])
    if demforbuild == 1:
        if standAlone == 0:
            gdal_dem = gdal.Open(configDict['filepath_dem']).ReadAsArray().astype(float)
            dem = gdal_dem.ReadAsArray().astype(float)
            nd = gdal_dem.GetRasterBand(1).GetNoDataValue()
        else:
            dem, _ , _ = common.load_raster(configDict['filepath_dem'], bbox=None)
            nd = -9999 #TODO: standAlone nd exposure

        # response to issue and #230
        dem[dem == nd] = 0.
        if dem.min() < 0:
            demraise = np.abs(dem.min())
            dem = dem + demraise
        else:
            demraise = 0

    #SVF
    zip = zipfile.ZipFile(configDict['input_svf'], 'r')
    zip.extractall(configDict['working_dir'])
    zip.close()

    if standAlone == 0:
        svf = gdal.Open(configDict['working_dir'] + "/svf.tif").ReadAsArray().astype(float)
        svfN = gdal.Open(configDict['working_dir'] + "/svfN.tif").ReadAsArray().astype(float)
        svfS = gdal.Open(configDict['working_dir'] + "/svfS.tif").ReadAsArray().astype(float)
        svfE = gdal.Open(configDict['working_dir'] + "/svfE.tif").ReadAsArray().astype(float)
        svfW = gdal.Open(configDict['working_dir'] + "/svfW.tif").ReadAsArray().astype(float)
    else:
        svf, _ , _ = common.load_raster(configDict['working_dir'] + "/svf.tif", bbox=None)
        svfN, _ , _ = common.load_raster(configDict['working_dir'] + "/svfN.tif", bbox=None)
        svfS, _ , _ = common.load_raster(configDict['working_dir'] + "/svfS.tif", bbox=None)
        svfE, _ , _ = common.load_raster(configDict['working_dir'] + "/svfE.tif", bbox=None)
        svfW, _ , _ = common.load_raster(configDict['working_dir'] + "/svfW.tif", bbox=None)    

    if usevegdem == 1:
        if standAlone == 0:
            svfveg = gdal.Open(configDict['working_dir'] + "/svfveg.tif").ReadAsArray().astype(float)
            svfNveg = gdal.Open(configDict['working_dir'] + "/svfNveg.tif").ReadAsArray().astype(float)
            svfSveg = gdal.Open(configDict['working_dir'] + "/svfSveg.tif").ReadAsArray().astype(float)
            svfEveg = gdal.Open(configDict['working_dir'] + "/svfEveg.tif").ReadAsArray().astype(float)
            svfWveg = gdal.Open(configDict['working_dir'] + "/svfWveg.tif").ReadAsArray().astype(float)

            svfaveg = gdal.Open(configDict['working_dir'] + "/svfaveg.tif").ReadAsArray().astype(float)
            svfNaveg = gdal.Open(configDict['working_dir'] + "/svfNaveg.tif").ReadAsArray().astype(float)
            svfSaveg = gdal.Open(configDict['working_dir'] + "/svfSaveg.tif").ReadAsArray().astype(float)
            svfEaveg = gdal.Open(configDict['working_dir'] + "/svfEaveg.tif").ReadAsArray().astype(float)
            svfWaveg = gdal.Open(configDict['working_dir'] + "/svfWaveg.tif").ReadAsArray().astype(float)
        else:
            svfveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfveg.tif", bbox=None)
            svfNveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfNveg.tif", bbox=None)
            svfSveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfSveg.tif", bbox=None)
            svfEveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfEveg.tif", bbox=None)
            svfWveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfWveg.tif", bbox=None)

            svfaveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfaveg.tif", bbox=None)
            svfNaveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfNaveg.tif", bbox=None)
            svfSaveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfSaveg.tif", bbox=None)
            svfEaveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfEaveg.tif", bbox=None)
            svfWaveg, _ , _ = common.load_raster(configDict['working_dir'] + "/svfWaveg.tif", bbox=None)
    else:
        svfveg = np.ones((rows, cols))
        svfNveg = np.ones((rows, cols))
        svfSveg = np.ones((rows, cols))
        svfEveg = np.ones((rows, cols))
        svfWveg = np.ones((rows, cols))
        svfaveg = np.ones((rows, cols))
        svfNaveg = np.ones((rows, cols))
        svfSaveg = np.ones((rows, cols))
        svfEaveg = np.ones((rows, cols))
        svfWaveg = np.ones((rows, cols))

    tmp = svf + svfveg - 1.
    tmp[tmp < 0.] = 0.
    # %matlab crazyness around 0
    svfalfa = np.arcsin(np.exp((np.log((1. - tmp)) / 2.)))

    if standAlone == 0:
        wallheight = gdal.Open(configDict['filepath_wh']).ReadAsArray().astype(float)
        wallaspect = gdal.Open(configDict['filepath_wa']).ReadAsArray().astype(float)
    else:
        wallheight, _ , _ = common.load_raster(configDict['filepath_wh'], bbox=None)
        wallaspect, _ , _ = common.load_raster(configDict['filepath_wa'], bbox=None)

    # Metdata
    headernum = 1
    delim = ' '
    Twater = []

    metdata = np.loadtxt(configDict['input_met'],skiprows=headernum, delimiter=delim) 

    location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
    YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = Solweig_2015a_metdata_noload(metdata,location, int(configDict['utc']))

    DOY = metdata[:, 1]
    hours = metdata[:, 2]
    minu = metdata[:, 3]
    Ta = metdata[:, 11]
    RH = metdata[:, 10]
    radG = metdata[:, 14]
    radD = metdata[:, 21]
    radI = metdata[:, 22]
    P = metdata[:, 12]
    Ws = metdata[:, 9]

    # POIs check TODO: fix for standalone
    if configDict['poi_file'] is not '': # usePOI:
        header = 'yyyy id   it imin dectime altitude azimuth kdir kdiff kglobal kdown   kup    keast ksouth ' \
                    'kwest knorth ldown   lup    least lsouth lwest  lnorth   Ta      Tg     RH    Esky   Tmrt    ' \
                    'I0     CI   Shadow  SVF_b  SVF_bv KsideI PET UTCI  CI_Tg   CI_TgG  KsideD  Lside   diffDown    Kside'
        poi_field = configDict['poi_field'] #self.parameterAsString(parameters, self.POI_FIELD, context)

        vlayer = QgsVectorLayer(configDict['poi_file'], 'point', 'ogr')
        idx = vlayer.fields().indexFromName(poi_field)
        numfeat = vlayer.featureCount()
        poiname = []
        poisxy = np.zeros((numfeat, 3)) - 999
        ind = 0
        for f in vlayer.getFeatures():  # looping through each POI
            y = f.geometry().centroid().asPoint().y()
            x = f.geometry().centroid().asPoint().x()

            poiname.append(f.attributes()[idx])
            poisxy[ind, 0] = ind
            poisxy[ind, 1] = np.round((x - minx) * scale)
            if miny >= 0:
                poisxy[ind, 2] = np.round((miny + rows * (1. / scale) - y) * scale)
            else:
                poisxy[ind, 2] = np.round((miny + rows * (1. / scale) - y) * scale)

            ind += 1

        for k in range(0, poisxy.shape[0]):
            poi_save = []  # np.zeros((1, 33))
            data_out = configDict['output_dir'] + '/POI_' + str(poiname[k]) + '.txt'
            np.savetxt(data_out, poi_save,  delimiter=' ', header=header, comments='')
        
        # Num format for POI output
        numformat = '%d %d %d %d %.5f ' + '%.2f ' * 36

        # Other PET variables
        sensorheight = param['Wind_Height']['Value']['magl']
        age = param['PET_settings']['Value']['Age']
        mbody = param['PET_settings']['Value']['Weight']
        ht = param['PET_settings']['Value']['Height']
        clo = param['PET_settings']['Value']['clo']
        activity = param['PET_settings']['Value']['Activity']
        sex = param['PET_settings']['Value']['Sex']
    else:
        poisxy = None
    
    #Posture settings
    if  param['Tmrt_params']['Value']['posture'] == "Standing":
        Fside = param["Posture"]["Standing"]["Value"]["Fside"]
        Fup = param["Posture"]["Standing"]["Value"]["Fup"]
        height = param["Posture"]["Standing"]["Value"]["height"]
        Fcyl = param["Posture"]["Standing"]["Value"]["Fcyl"] 
        pos = 1
    else:
        Fside = param["Posture"]["Sitting"]["Value"]["Fside"]
        Fup = param["Posture"]["Sitting"]["Value"]["Fup"]
        height = param["Posture"]["Sitting"]["Value"]["height"]
        Fcyl = param["Posture"]["Sitting"]["Value"]["Fcyl"]
        pos = 0

    # Radiative surface influence, Rule of thumb by Schmid et al. (1990).
    first = np.round(height)
    if first == 0.:
        first = 1.
    second = np.round((height * 20.))

    if usevegdem == 1:
        # Conifer or deciduous
        if configDict['conifer_bool']:
            leafon = np.ones((1, DOY.shape[0]))
        else:
            leafon = np.zeros((1, DOY.shape[0]))
            if param["Tree_settings"]["Value"]["First_day_leaf"] > param["Tree_settings"]["Value"]["Last_day_leaf"]:
                leaf_bool = ((DOY > param["Tree_settings"]["Value"]["First_day_leaf"]) | (DOY < param["Tree_settings"]["Value"]["Last_day_leaf"]))
            else:
                leaf_bool = ((DOY > param["Tree_settings"]["Value"]["First_day_leaf"]) & (DOY < param["Tree_settings"]["Value"]["Last_day_leaf"]))
            leafon[0, leaf_bool] = 1

        # % Vegetation transmittivity of shortwave radiation
        psi = leafon * transVeg
        psi[leafon == 0] = 0.5
        # amaxvalue
        vegmax = vegdsm.max()
        amaxvalue = dsm.max() - dsm.min()
        amaxvalue = np.maximum(amaxvalue, vegmax)

        # Elevation vegdsms if buildingDEM includes ground heights
        vegdsm = vegdsm + dsm
        vegdsm[vegdsm == dsm] = 0
        vegdsm2 = vegdsm2 + dsm
        vegdsm2[vegdsm2 == dsm] = 0

        # % Bush separation
        bush = np.logical_not((vegdsm2 * vegdsm)) * vegdsm

        svfbuveg = (svf - (1. - svfveg) * (1. - transVeg))  # % major bug fixed 20141203
    else:
        psi = leafon * 0. + 1.
        svfbuveg = svf
        bush = np.zeros([rows, cols])
        amaxvalue = 0

    # Initialization of maps
    Knight = np.zeros((rows, cols))
    Tgmap1 = np.zeros((rows, cols))
    Tgmap1E = np.zeros((rows, cols))
    Tgmap1S = np.zeros((rows, cols))
    Tgmap1W = np.zeros((rows, cols))
    Tgmap1N = np.zeros((rows, cols))

    # Create building boolean raster from either land cover or height rasters
    if demforbuild == 0:
        buildings = np.copy(lcgrid)
        buildings[buildings == 7] = 1
        buildings[buildings == 6] = 1
        buildings[buildings == 5] = 1
        buildings[buildings == 4] = 1
        buildings[buildings == 3] = 1
        buildings[buildings == 2] = 0
    else:
        buildings = dsm - dem
        buildings[buildings < 2.] = 1.
        buildings[buildings >= 2.] = 0.

    if int(configDict['savebuild']) == 1:
        if standAlone == 0:
            saveraster(gdal_dsm, configDict['output_dir'] + '/buildings.tif', buildings)
        else:
            common.save_raster(configDict['output_dir'] + '/buildings.tif', buildings, dsm_transf, dsm_crs)

    # Import shadow matrices (Anisotropic sky)
    anisotropic_sky = int(configDict['aniso'])
    if anisotropic_sky == 1:  #UseAniso
        data = np.load(configDict['input_aniso'])
        shmat = data['shadowmat']
        vegshmat = data['vegshadowmat']
        vbshvegshmat = data['vbshmat']
        if usevegdem == 1:
            diffsh = np.zeros((rows, cols, shmat.shape[2]))
            for i in range(0, shmat.shape[2]):
                diffsh[:, :, i] = shmat[:, :, i] - (1 - vegshmat[:, :, i]) * (1 - transVeg) # changes in psi not implemented yet
        else:
            diffsh = shmat

        # Estimate number of patches based on shadow matrices
        if shmat.shape[2] == 145:
            patch_option = 1 # patch_option = 1 # 145 patches
        elif shmat.shape[2] == 153:
            patch_option = 2 # patch_option = 2 # 153 patches
        elif shmat.shape[2] == 306:
            patch_option = 3 # patch_option = 3 # 306 patches
        elif shmat.shape[2] == 612:
            patch_option = 4 # patch_option = 4 # 612 patches

        # asvf to calculate sunlit and shaded patches
        asvf = np.arccos(np.sqrt(svf))

        # Empty array for steradians
        steradians = np.zeros((shmat.shape[2]))
    else:
        # anisotropic_sky = 0
        diffsh = None
        shmat = None
        vegshmat = None
        vbshvegshmat = None
        asvf = None
        patch_option = 0
        steradians = 0
    
    # % Ts parameterisation maps
    if landcover == 1.:
        # Get land cover properties for Tg wave (land cover scheme based on Bogren et al. 2000, explained in Lindberg et al., 2008 and Lindberg, Onomura & Grimmond, 2016)
        [TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST, TmaxLST_wall] = Tgmaps_v1(lcgrid.copy(), param)
    else:
        TgK = Knight + param['Ts_deg']['Value']['Cobble_stone_2014a']
        Tstart = Knight - param['Tstart']['Value']['Cobble_stone_2014a']
        TmaxLST = param['TmaxLST']['Value']['Cobble_stone_2014a']
        alb_grid = Knight + param['Albedo']['Effective']['Value']['Cobble_stone_2014a']
        emis_grid = Knight + param['Emissivity']['Value']['Cobble_stone_2014a']
        TgK_wall = param['Ts_deg']['Value']['Walls']
        Tstart_wall = param['Tstart']['Value']['Walls']
        TmaxLST_wall = param['TmaxLST']['Value']['Walls']

    # Import data for wall temperature parameterization TODO: fix for standalone
    wallScheme = int(configDict['wallscheme'])
    if wallScheme == 1:
        wallData = np.load(configDict['input_wall'])
        voxelMaps = wallData['voxelId']
        voxelTable = wallData['voxelTable']
        # Get wall type set in GUI
        wall_type = configDict['walltype'] # str(100 + int(self.parameterAsString(parameters, self.WALL_TYPE, context))) #TODO

        from ...functions import wallalgorithms as wa
        walls_scheme = wa.findwalls_sp(dsm, 2, np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]]))

        dirwalls_scheme = wa.filter1Goodwin_as_aspect_v3(walls_scheme.copy(), scale, dsm, None, 100./180.)
        
        # Used in wall temperature parameterization scheme
        first_timestep = (pd.to_datetime(YYYY[0][0], format='%Y') + pd.to_timedelta(DOY[0] - 1, unit='d') + 
            pd.to_timedelta(hours[0], unit='h') + pd.to_timedelta(minu[0], unit='m'))
        second_timestep = (pd.to_datetime(YYYY[0][1], format='%Y') + pd.to_timedelta(DOY[1] - 1, unit='d') + 
            pd.to_timedelta(hours[1], unit='h') + pd.to_timedelta(minu[1], unit='m'))
        
        timeStep = (second_timestep - first_timestep).seconds

        # Load voxelTable as Pandas DataFrame
        voxelTable, dirwalls_scheme = load_walls(voxelTable, param, wall_type, dirwalls_scheme, Ta[0], timeStep, alb_grid, landcover, lcgrid, dsm)
        # Unique wall types
        #thermal_effusivity = voxelTable['thermalEffusivity'].unique()

        # Empty wall temperature matrix for parameterization scheme
        # wallScheme_feedback = "Running with wall parameterization scheme. Walls divided into " + str(int(voxelTable.shape[0])) + " unique voxels."
        # feedback.setProgressText(wallScheme_feedback)

        # Use wall of interest
        if configDict['woi_file'] is not None:
            (dsm_minx, dsm_x_size, dsm_x_rotation, dsm_miny, dsm_y_rotation, dsm_y_size) = gdal_dsm.GetGeoTransform() #TODO: fix for standalone

            woi_field = configDict['woi_field'] #self.parameterAsStrings(parameters, self.WOI_FIELD, context)
            woisxy, woiname = wallOfInterest(configDict['woi_file'], woi_field, minx, miny, scale, rows, configDict['output_dir'], dsm_minx, dsm_x_size, dsm_miny, dsm_y_size)

        # Create pandas datetime object to be used when createing an xarray DataSet where wall temperatures/radiation is stored and eventually saved as a NetCDf
        if configDict["wallnetcdf"] == 1:
            met_for_xarray = (pd.to_datetime(YYYY[0][:], format='%Y') + pd.to_timedelta(DOY - 1, unit='d') + 
                            pd.to_timedelta(hours, unit='h') + pd.to_timedelta(minu, unit='m'))
    else:
        wallScheme = 0
        voxelMaps = 0
        voxelTable = 0
        timeStep = 0
        #thermal_effusivity = 0
        walls_scheme = np.ones((rows, cols)) * 10.
        dirwalls_scheme = np.ones((rows, cols)) * 10.

    # Initialisation of time related variables
    if Ta.__len__() == 1:
        timestepdec = 0
    else:
        timestepdec = dectime[1] - dectime[0]
    timeadd = 0.
    # timeaddE = 0.
    # timeaddS = 0.
    # timeaddW = 0.
    # timeaddN = 0.
    firstdaytime = 1.


    # Save svf #TODO is this for the patch images?
    if anisotropic_sky == 1:
        if not poisxy is None:
            patch_characteristics = np.zeros((shmat.shape[2], poisxy.shape[0]))
            for idx in range(poisxy.shape[0]):
                for idy in range(shmat.shape[2]):
                    # Calculations for patches on sky, shmat = 1 = sky is visible
                    temp_sky = ((shmat[:,:,idy] == 1) & (vegshmat[:,:,idy] == 1))
                    # Calculations for patches that are vegetation, vegshmat = 0 = shade from vegetation
                    temp_vegsh = ((vegshmat[:,:,idy] == 0) | (vbshvegshmat[:,:,idy] == 0))
                    # Calculations for patches that are buildings, shmat = 0 = shade from buildings
                    temp_vbsh = (1 - shmat[:,:,idy]) * vbshvegshmat[:,:,idy]
                    temp_sh = (temp_vbsh == 1)
                    if wallScheme:
                        temp_sh_w = temp_sh * voxelMaps[:, :, idy]
                        temp_sh_roof = temp_sh * (voxelMaps[:, :, idy] == 0)
                    else:
                        temp_sh_w = 0
                        temp_sh_roof = 0
                    # Sky patch
                    if temp_sky[int(poisxy[idx, 2]), int(poisxy[idx, 1])]:
                        patch_characteristics[idy,idx] = 1.8
                    # Vegetation patch
                    elif (temp_vegsh[int(poisxy[idx, 2]), int(poisxy[idx, 1])]):
                        patch_characteristics[idy,idx] = 2.5
                    # Building patch
                    elif (temp_sh[int(poisxy[idx, 2]), int(poisxy[idx, 1])]):
                        if wallScheme:
                            if (temp_sh_w[int(poisxy[idx, 2]), int(poisxy[idx, 1])]):
                                patch_characteristics[idy, idx] = 4.5
                            elif (temp_sh_roof[int(poisxy[idx, 2]), int(poisxy[idx, 1])]):
                                patch_characteristics[idy, idx] = 6.0
                        else:
                            patch_characteristics[idy,idx] = 4.5
                    # Roof patch
                    #elif

    # If metfile starts at night
    CI = 1.

    # reading variables from config and parameters that is not yet presented
    cyl = int(configDict["cyl"])
    albedo_b = param['Albedo']['Effective']['Value']['Walls']
    ewall = param['Emissivity']['Value']['Walls']
    onlyglobal = int(configDict['onlyglobal'])
    elvis = 0.0
    absK = param['Tmrt_params']['Value']['absK'] 
    absL = param['Tmrt_params']['Value']['absL'] 
    
    # Main loop
    tmrtplot = np.zeros((rows, cols))
    TgOut1 = np.zeros((rows, cols))

    # Initiate array for I0 values
    if np.unique(DOY).shape[0] > 1:
        unique_days = np.unique(DOY)
        first_unique_day = DOY[DOY == unique_days[0]]
        I0_array = np.zeros((first_unique_day.shape[0]))
    else:
        first_unique_day = DOY.copy()
        I0_array = np.zeros((DOY.shape[0]))

    if standAlone == 1:
        progress = tqdm(total=Ta.__len__())

    for i in np.arange(0, Ta.__len__()):
        if feedback is not None:
            feedback.setProgress(int(i * (100. / Ta.__len__()))) # move progressbar forward
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break
        else:
            progress.update(1)

        # Daily water body temperature
        if landcover == 1:
            if ((dectime[i] - np.floor(dectime[i]))) == 0 or (i == 0):
                Twater = np.mean(Ta[jday[0] == np.floor(dectime[i])])
        # Nocturnal cloudfraction from Offerle et al. 2003
        if (dectime[i] - np.floor(dectime[i])) == 0:
            daylines = np.where(np.floor(dectime) == dectime[i])
            if daylines.__len__() > 1:
                alt = altitude[0][daylines]
                alt2 = np.where(alt > 1)
                rise = alt2[0][0]
                [_, CI, _, _, _] = clearnessindex_2013b(zen[0, i + rise + 1], jday[0, i + rise + 1],
                                                        Ta[i + rise + 1],
                                                        RH[i + rise + 1] / 100., radG[i + rise + 1], location,
                                                        P[i + rise + 1])
                if (CI > 1.) or (CI == np.inf):
                    CI = 1.
            else:
                CI = 1.

        if altitude[0][i] > 0:
            radI[i] = radI[i]/np.sin(altitude[0][i] * np.pi/180)
        else:
            radG[i] = 0.
            radD[i] = 0.
            radI[i] = 0.


        Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, timeadd, \
                Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, Keast, Ksouth, Kwest, Knorth, Least, \
                Lsouth, Lwest, Lnorth, KsideI, TgOut1, TgOut, radIout, radDout, \
                Lside, Lsky_patch_characteristics, CI_Tg, CI_TgG, KsideD, \
                    dRad, Kside, steradians, voxelTable = so.Solweig_2025a_calc(
                    i, dsm, scale, rows, cols, svf, svfN, svfW, svfE, svfS, svfveg,
                    svfNveg, svfEveg, svfSveg, svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg, \
                    vegdsm, vegdsm2, albedo_b, absK, absL, ewall, Fside, Fup, Fcyl, altitude[0][i],
                    azimuth[0][i], zen[0][i], jday[0][i], usevegdem, onlyglobal, buildings, location,
                    psi[0][i], landcover, lcgrid, dectime[i], altmax[0][i], wallaspect,
                    wallheight, cyl, elvis, Ta[i], RH[i], radG[i], radD[i], radI[i], P[i], amaxvalue,
                    bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
                    TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timestepdec, 
                    Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI, TgOut1, diffsh, shmat, vegshmat, vbshvegshmat, 
                    anisotropic_sky, asvf, patch_option, voxelMaps, voxelTable, Ws[i], wallScheme, timeStep, steradians, walls_scheme, dirwalls_scheme)

        # Save I0 for I0 vs. Kdown output plot to check if UTC is off
        if i < first_unique_day.shape[0]:
            I0_array[i] = I0
        elif i == first_unique_day.shape[0]:
            # Output I0 vs. Kglobal plot
            radG_for_plot = radG[DOY == first_unique_day[0]]
            # hours_for_plot = hours[DOY == first_unique_day[0]]
            dectime_for_plot = dectime[DOY == first_unique_day[0]]
            fig, ax = plt.subplots()
            ax.plot(dectime_for_plot, I0_array, label='I0')
            ax.plot(dectime_for_plot, radG_for_plot, label='Kglobal')
            ax.set_ylabel('Shortwave radiation [$Wm^{-2}$]')
            ax.set_xlabel('Decimal time')
            ax.set_title('UTC' + str(configDict['utc']))
            ax.legend()
            fig.savefig(configDict['output_dir'] + '/metCheck.png', dpi=150)                

        tmrtplot = tmrtplot + Tmrt

        if altitude[0][i] > 0:
            w = 'D'
        else:
            w = 'N'

        # Write to POIs
        if not poisxy is None:              
            for k in range(0, poisxy.shape[0]):
                poi_save = np.zeros((1, 41))
                poi_save[0, 0] = YYYY[0][i]
                poi_save[0, 1] = jday[0][i]
                poi_save[0, 2] = hours[i]
                poi_save[0, 3] = minu[i]
                poi_save[0, 4] = dectime[i]
                poi_save[0, 5] = altitude[0][i]
                poi_save[0, 6] = azimuth[0][i]
                poi_save[0, 7] = radIout
                poi_save[0, 8] = radDout
                poi_save[0, 9] = radG[i]
                poi_save[0, 10] = Kdown[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 11] = Kup[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 12] = Keast[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 13] = Ksouth[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 14] = Kwest[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 15] = Knorth[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 16] = Ldown[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 17] = Lup[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 18] = Least[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 19] = Lsouth[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 20] = Lwest[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 21] = Lnorth[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 22] = Ta[i]
                poi_save[0, 23] = TgOut[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 24] = RH[i]
                poi_save[0, 25] = esky
                poi_save[0, 26] = Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 27] = I0
                poi_save[0, 28] = CI
                poi_save[0, 29] = shadow[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 30] = svf[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 31] = svfbuveg[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 32] = KsideI[int(poisxy[k, 2]), int(poisxy[k, 1])]
                # Recalculating wind speed based on powerlaw
                WsPET = (1.1 / sensorheight) ** 0.2 * Ws[i]
                WsUTCI = (10. / sensorheight) ** 0.2 * Ws[i]
                resultPET = p._PET(Ta[i], RH[i], Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])], WsPET,
                                    mbody, age, ht, activity, clo, sex)
                poi_save[0, 33] = resultPET
                resultUTCI = utci.utci_calculator(Ta[i], RH[i], Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])],
                                                    WsUTCI)
                poi_save[0, 34] = resultUTCI
                poi_save[0, 35] = CI_Tg
                poi_save[0, 36] = CI_TgG
                poi_save[0, 37] = KsideD[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 38] = Lside[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 39] = dRad[int(poisxy[k, 2]), int(poisxy[k, 1])]
                poi_save[0, 40] = Kside[int(poisxy[k, 2]), int(poisxy[k, 1])]
                data_out = configDict['output_dir'] + '/POI_' + str(poiname[k]) + '.txt'
                # f_handle = file(data_out, 'a')
                f_handle = open(data_out, 'ab')
                np.savetxt(f_handle, poi_save, fmt=numformat)
                f_handle.close()

        # If wall temperature parameterization scheme is in use
        if configDict['wallscheme'] == 1: # folderWallScheme: TODO: Fix for standalone
            # Store wall data for output
            if not woisxy is None:
                for k in range(0, woisxy.shape[0]):
                    temp_wall = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallTemperature'].to_numpy()
                    K_in = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'K_in'].to_numpy()
                    L_in = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'L_in'].to_numpy()
                    wallShade = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallShade'].to_numpy()
                    temp_all = np.concatenate([temp_wall, K_in, L_in, wallShade])
                    # temp_all = np.concatenate([temp_wall])
                    # wall_data = np.zeros((1, 7 + temp_wall.shape[0]))
                    wall_data = np.zeros((1, 7 + temp_all.shape[0]))
                    # Part of file name (wallid), i.e. WOI_wallid.txt
                    data_out = configDict['output_dir'] + '/WOI_' + str(woiname[k]) + '.txt'                    
                    if i == 0:
                        # Output file header
                        #header = 'yyyy id   it imin dectime Ta  SVF Ts'
                        header = 'yyyy id   it imin dectime Ta  SVF' + ' Ts' * temp_wall.shape[0] + ' Kin' * K_in.shape[0] + ' Lin' * L_in.shape[0] + ' shade' * wallShade.shape[0]
                        # Part of file name (wallid), i.e. WOI_wallid.txt
                        # woiname = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallId'].to_numpy()[0]
                        woi_save = []  # 
                        np.savetxt(data_out, woi_save,  delimiter=' ', header=header, comments='')                        
                    # Fill wall_data with variables
                    wall_data[0, 0] = YYYY[0][i] 
                    wall_data[0, 1] = jday[0][i]
                    wall_data[0, 2] = hours[i]
                    wall_data[0, 3] = minu[i]
                    wall_data[0, 4] = dectime[i]
                    wall_data[0, 5] = Ta[i]
                    wall_data[0, 6] = svf[int(woisxy[k, 2]), int(woisxy[k, 1])]
                    wall_data[0, 7:] = temp_all

                    # Num format for output file data
                    woi_numformat = '%d %d %d %d %.5f %.2f %.2f' + ' %.2f' * temp_all.shape[0]
                    # Open file, add data, save
                    f_handle = open(data_out, 'ab')
                    np.savetxt(f_handle, wall_data, fmt=woi_numformat)
                    f_handle.close()                

            # Save wall temperature/radiation as NetCDF TODO: fix for standAlone?
            if configDict['wallnetcdf'] == '1': # wallNetCDF:
                netcdf_output = configDict['outputDir'] + '/walls.nc'
                walls_as_netcdf(voxelTable, rows, cols, met_for_xarray, i, dsm, configDict['filepath_dsm'], netcdf_output)

        if hours[i] < 10:
            XH = '0'
        else:
            XH = ''
        if minu[i] < 10:
            XM = '0'
        else:
            XM = ''

        time_code = (str(int(YYYY[0, i])) + "_" + str(int(DOY[i])) + "_" + XH + str(int(hours[i])) + XM + str(int(minu[i])) + w)

        if configDict['outputtmrt'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Tmrt_' + time_code + '.tif', Tmrt)
            else:
                common.save_raster(configDict['output_dir'] + '/Tmrt_' + time_code + '.tif', Tmrt, dsm_transf, dsm_crs)
        if configDict['outputkup'] == '1':
            if standAlone == 0:           
                saveraster(gdal_dsm, configDict['output_dir'] + '/Kup_' + time_code + '.tif', Kup)
            else:
                common.save_raster(configDict['output_dir'] + '/Kup_' + time_code + '.tif', Kup, dsm_transf, dsm_crs)
        if configDict['outputkdown'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Kdown_' + time_code + '.tif', Kdown)
            else:
                common.save_raster(configDict['output_dir'] + '/Kdown_' + time_code + '.tif', Kdown, dsm_transf, dsm_crs)
        if configDict['outputlup'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Lup_' + time_code + '.tif', Lup)
            else:
                common.save_raster(configDict['output_dir'] + '/Lup_' + time_code + '.tif', Lup, dsm_transf, dsm_crs)
        if configDict['outputldown'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Ldown_' + time_code + '.tif', Ldown)
            else:
                common.save_raster(configDict['output_dir'] + '/Ldown_' + time_code + '.tif', Ldown, dsm_transf, dsm_crs)
        if configDict['outputsh'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Shadow_' + time_code + '.tif', shadow)
            else:
                common.save_raster(configDict['output_dir'] + '/Shadow_' + time_code + '.tif', shadow, dsm_transf, dsm_crs)
            
        if configDict['outputkdiff'] == '1':
            if standAlone == 0:
                saveraster(gdal_dsm, configDict['output_dir'] + '/Kdiff_' + time_code + '.tif', dRad)  
            else:
                common.save_raster(configDict['output_dir'] + '/Kdiff_' + time_code + '.tif', dRad, dsm_transf, dsm_crs)         

        # Sky view image of patches
        if ((anisotropic_sky == 1) & (i == 0) & (not poisxy is None)):
                for k in range(poisxy.shape[0]):
                    Lsky_patch_characteristics[:,2] = patch_characteristics[:,k]
                    skyviewimage_out = configDict['output_dir'] + '/POI_' + str(poiname[k]) + '.png'
                    PolarBarPlot(Lsky_patch_characteristics, altitude[0][i], azimuth[0][i], 'Hemisphere partitioning', skyviewimage_out, 0, 5, 0)

    # Save files for Tree Planter
    if configDict['outputtreeplanter'] == '1': # outputTreeplanter:
        if feedback is not None:
            feedback.setProgressText("Saving files for Tree Planter tool")
        # Save DSM
        copyfile(configDict['filepath_dsm'], configDict['output_dir'] + '/DSM.tif')

        # Save CDSM
        if usevegdem == 1:
            copyfile(configDict['filepath_cdsm'],  configDict['output_dir'] + '/CDSM.tif')

        albedo_g = param['Albedo']['Effective']['Value']['Cobble_stone_2014a']
        eground = param['Emissivity']['Value']['Cobble_stone_2014a']

        # Saving settings from SOLWEIG for SOLWEIG1D in TreePlanter
        settingsHeader = 'UTC, posture, onlyglobal, landcover, anisotropic, cylinder, albedo_walls, albedo_ground, emissivity_walls, emissivity_ground, absK, absL, elevation, patch_option'
        settingsFmt = '%i', '%i', '%i', '%i', '%i', '%i', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%1.2f', '%i'
        settingsData = np.array([[int(configDict['utc']), pos, onlyglobal, landcover, anisotropic_sky, cyl, albedo_b, albedo_g, ewall, eground, absK, absL, alt, patch_option]])
        # print(settingsData)
        np.savetxt(configDict['output_dir'] + '/treeplantersettings.txt', settingsData, fmt=settingsFmt, header=settingsHeader, delimiter=' ')

    # Copying met file for SpatialTC
    copyfile(configDict['input_met'], configDict['output_dir'] + '/metforcing.txt')
    
    tmrtplot = tmrtplot / Ta.__len__()  # fix average Tmrt instead of sum, 20191022
    if standAlone == 0:
        saveraster(gdal_dsm, configDict['output_dir'] + '/Tmrt_average.tif', tmrtplot)
    else:
        common.save_raster(configDict['output_dir'] + '/Tmrt_average.tif', tmrtplot, dsm_transf, dsm_crs)
