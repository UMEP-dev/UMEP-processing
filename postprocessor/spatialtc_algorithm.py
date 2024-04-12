# -*- coding: utf-8 -*-


__author__ = 'Fredrik Lindberg'
__date__ = '2021-02-05'
__copyright__ = '(C) 2021 by Fredrik Lindberg'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterEnum,
                       QgsProcessingException,
                       QgsFeature,
                       QgsVectorFileWriter,
                       QgsVectorDataProvider,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterDefinition)
from qgis.PyQt.QtGui import QIcon
from osgeo import gdal, osr, ogr
from osgeo.gdalconst import *
import os
import numpy as np
import pandas as pd
import inspect
from pathlib import Path
import sys
from ..util.misc import saveraster
from ..functions.SOLWEIGpython import PET_calculations as pet
from ..functions.SOLWEIGpython import UTCI_calculations as utci
from ..functions.SOLWEIGpython.COMFA.radiationfunctionsCOMFA import COMFA_RAD_SPATIAL_TC as COMFA_rad
from ..functions.SOLWEIGpython.COMFA.COMFA_BUDGET import COMFA_Mact
from ..functions.SOLWEIGpython.COMFA.COMFA_BUDGET import COMFA_BUDGET

def load_grid(filepath, feedback):
    # Function that loads QGIS input raster into numpy array
    # provider = data.dataProvider()  # Provider for raster data
    # temp_filename = str(provider.dataSourceUri()) # File path for data
    try:
        temp_grid = gdal.Open(filepath) # Open raster with GDAL
        feedback.setProgressText("Successfully loaded " + filepath)
    except:
        raise QgsProcessingException("Error: Could not load " + filepath + ". File does not exist. Check path!")
    
    # Return gdal raster layer as numpy array, number of rows and columns in raster
    return temp_grid.ReadAsArray().astype(float), temp_grid.ReadAsArray().astype(float).shape[0], temp_grid.ReadAsArray().astype(float).shape[1] 

def get_latlon(raster, gdal_raster):
    old_cs = osr.SpatialReference()
    raster_ref = raster.crs().toWkt()
    old_cs.ImportFromWkt(raster_ref)

    wgs84_wkt = """
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]]"""

    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(wgs84_wkt)

    transform = osr.CoordinateTransformation(old_cs, new_cs)
    widthx = gdal_raster.RasterXSize
    heightx = gdal_raster.RasterYSize
    geotransform = gdal_raster.GetGeoTransform()
    minx = geotransform[0]
    miny = geotransform[3] + widthx * geotransform[4] + heightx * geotransform[5]
    lonlat = transform.TransformPoint(minx, miny)
    gdalver = float(gdal.__version__[0])
    if gdalver == 3.:
        lon = lonlat[1] #changed to gdal 3
        lat = lonlat[0] #changed to gdal 3
    else:
        lon = lonlat[0] #changed to gdal 2
        lat = lonlat[1] #changed to gdal 2

    return lat, lon

class ProcessingSpatialTCAlgorithm(QgsProcessingAlgorithm):
    """
    This class is a postprocessing algoritm to calculate Spatial thermal comfort maps
    """
    # SOLWEIG_DIR = 'SOLWEIG_DIR'
    TMRT_MAP = 'TMRT_MAP'
    BUILDING_MAP = 'BUILDING_MAP'
    UROCK_MAP = 'UROCK_MAP'
    TC_TYPE = 'TC_TYPE'
    METDATA = 'METDATA'

    # KUP_MAP = 'KUP_MAP'
    # KDOWN_MAP = 'KDOWN_MAP'
    # KDIFF_MAP = 'KDIFF_MAP'
    # LUP_MAP = 'LUP_MAP'
    # LDOWN_MAP = 'LDOWN_MAP'

    #PET parameters
    AGE = 'AGE'
    ACTIVITY = 'ACTIVITY'
    CLO = 'CLO'
    WEIGHT = 'WEIGHT'
    HEIGHT = 'HEIGHT'
    SEX = 'SEX'
    # SENSOR_HEIGHT = 'SENSOR_HEIGHT'

    COMFA = 'COMFA'

    # Output
    TC_OUT = 'TC_OUT'


    def initAlgorithm(self, config):
        # self.addParameter(QgsProcessingParameterFile(self.SOLWEIG_DIR,
        #                                              self.tr('Path to SOLWEIG output folder'),
        #                                              QgsProcessingParameterFile.Folder,
        #                                              optional=False))
        self.addParameter(QgsProcessingParameterRasterLayer(self.TMRT_MAP,
                                                            self.tr('Mean radiant temperature raster (SOLWEIG)'),
                                                             '', 
                                                             optional=False))
        self.addParameter(QgsProcessingParameterRasterLayer(self.UROCK_MAP,
                                                            self.tr('Pedestrian wind speed raster (URock)'),
                                                             '', 
                                                             optional=False))                   
        # self.addParameter(QgsProcessingParameterRasterLayer(self.BUILDING_MAP,
        #                                                     self.tr('Raster to exclude building pixels from analysis (SOLWEIG)'),
        #                                                      '', 
        #                                                      optional=False))
        # self.addParameter(QgsProcessingParameterFile(self.METDATA,
        #                                             self.tr('Input meteorological file (.txt). Same as used for Tmrt map.'), extension='txt',
        #                                             optional=False))
        self.varType = ((self.tr('Physiological Equivalent Temperature (PET)'), '0'),
                         (self.tr('Universal Thermal Climate Index (UTCI)'), '1'),
                         (self.tr('COMfort FormulA (COMFA)'), '2'))
        self.addParameter(QgsProcessingParameterEnum(self.TC_TYPE,
                                                     self.tr('Thermal Index to map'),
                                                     options=[i[0] for i in self.varType],
                                                     defaultValue=0))
        
        #PET parameters
        age = QgsProcessingParameterNumber(self.AGE, self.tr('Age (yy)'),
                QgsProcessingParameterNumber.Integer,
                QVariant(35), optional=True, minValue=0, maxValue=120)
        age.setFlags(age.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(age)
        act = QgsProcessingParameterNumber(self.ACTIVITY, self.tr('Activity (W)'),
                QgsProcessingParameterNumber.Double,
                QVariant(80), optional=True, minValue=0, maxValue=1000)
        act.setFlags(act.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(act)
        clo = QgsProcessingParameterNumber(self.CLO, self.tr('Clothing (clo)'),
                QgsProcessingParameterNumber.Double,
                QVariant(0.9), optional=True, minValue=0, maxValue=10)
        clo.setFlags(clo.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(clo)
        wei = QgsProcessingParameterNumber(self.WEIGHT, self.tr('Weight (kg)'),
                QgsProcessingParameterNumber.Integer,
                QVariant(75), optional=True, minValue=0, maxValue=500) 
        wei.setFlags(wei.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(wei)
        hei = QgsProcessingParameterNumber(self.HEIGHT, self.tr('Height (cm)'),
                QgsProcessingParameterNumber.Integer,
                QVariant(180), optional=True, minValue=0, maxValue=250) 
        hei.setFlags(hei.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(hei)
        sex = QgsProcessingParameterEnum(
            self.SEX, self.tr('Sex'), ['Male', 'Female'], optional=True, defaultValue=0)
        sex.setFlags(sex.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(sex)

        # COMFA or COMFA-kid
        comfa = QgsProcessingParameterBoolean(self.COMFA,
            self.tr("COMFA-kid (Cheng and Brown, 2020)"), defaultValue=False)
        comfa.setFlags(comfa.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(comfa)        

        # Output
        self.addParameter(QgsProcessingParameterRasterDestination(self.TC_OUT,
                                                                  self.tr("Output thermal comfort raster"),
                                                                  None,
                                                                  optional=False))


    def processAlgorithm(self, parameters, context, feedback):
        
        # InputParameters
        # solweigDir = self.parameterAsString(parameters, self.SOLWEIG_DIR, context)
        
        tmrt = self.parameterAsRasterLayer(parameters, self.TMRT_MAP, context)
        # Kup = self.parameterAsRasterLayer(parameters, self.KUP_MAP, context)
        # Kdown = self.parameterAsRasterLayer(parameters, self.KDOWN_MAP, context)
        # Kdiff = self.parameterAsRasterLayer(parameters, self.KDIFF_MAP, context)
        # Lup = self.parameterAsRasterLayer(parameters, self.LUP_MAP, context)
        # Ldown = self.parameterAsRasterLayer(parameters, self.LDOWN_MAP, context)
        ws = self.parameterAsRasterLayer(parameters, self.UROCK_MAP, context)
        #sensorheight = self.parameterAsDouble(parameters, self.SENSOR_HEIGHT, context)
        sensorheight = 1.5 # Should be retrieved from URock RunInfo
        # buildings = self.parameterAsRasterLayer(parameters, self.BUILDING_MAP, context)
        tcTypeStr = self.parameterAsString(parameters, self.TC_TYPE, context)
        tcType = int(tcTypeStr)
        if tcType == 0:
            thermal_index = 'Physiological Equivalent Temperature (PET)'
        elif tcType == 1:
            thermal_index = 'Universal Thermal Climate Index (UTCI)'
        elif tcType == 2:
            thermal_index = 'COMfort FormulA (COMFA)'
            comfa_kid = self.parameterAsBool(parameters, self.COMFA, context)
        # metdata = self.parameterAsString(parameters, self.METDATA, context)

        outputRaster = self.parameterAsOutputLayer(parameters, self.TC_OUT, context)

        mbody = None
        ht = None
        clo = None
        age = None
        activity = None
        sex = None

        feedback.setProgressText("Initializing...")
        # Get SOLWEIG output folder path from Tmrt raster path
        provider = tmrt.dataProvider()
        filepath_tmrt = str(provider.dataSourceUri()) # Path for Tmrt raster
        solweig_path = filepath_tmrt.split('Tmrt')[0] # Path to SOLWEIG output folder, i.e. where Tmrt raster is located
        _, solweigfile = os.path.split(filepath_tmrt)
        # solweig_path = os.path.dirname(filepath_tmrt) # issue 31 

        # LOAD raster data
        try: #response to issue #51
            # Load buildings raster (should be in SOLWEIG output folder)
            filepath_build = solweig_path + '/buildings.tif' 
            gdal_buildings = gdal.Open(filepath_build)
            build = gdal_buildings.ReadAsArray().astype(float)
        except:
            raise QgsProcessingException("Error: No building raster found. It should be located in the same folder as the output"
                                         " folder specified when SOLWEIG was executed. Make sure you ticked in 'Save necessary rasters for" 
                                         " Tree planter and Spatial TC tools' when you ran SOLWEIG.")
        
        # Get latlon from grid coordinate system
        lat, lon = get_latlon(tmrt, gdal_buildings)

        # LOAD Metdata
        try:
            metdata = np.loadtxt(solweig_path + '/metforcing.txt', skiprows=1, delimiter=' ')
        except:
            raise QgsProcessingException("Error: Make sure format of meteorological file is correct. You can "
                                                        "prepare your data by using 'Prepare Existing Data' in "
                                                        "the Pre-processor")

        if metdata.shape[1] == 24:
            feedback.setProgressText("Meteorological data successfully loaded")
        else:
            raise QgsProcessingException("Error: Wrong number of columns in meteorological data. Use the forcing file "
                                                        "that was used when Tmrt was calculated")

        # Creating vectors from meteorological input
        # yyyy = metdata[:,0]
        # DOY = metdata[:, 1]
        # hours = metdata[:, 2]
        # minu = metdata[:, 3]
        # Ta = metdata[:, 11]
        # RH = metdata[:, 10]
        # radG = self.metdata[:, 14]
        # radD = self.metdata[:, 21]
        # radI = self.metdata[:, 22]
        # P = self.metdata[:, 12]
        # Ws = self.metdata[:, 9]

        # Derive metdata from Trmt raster name
        yyyyTmrt = int(solweigfile.split('_')[-3]) #int(filepath_tmrt[-18:-14]) #issue 571
        doyTmrt = int(solweigfile.split('_')[-2]) #int(filepath_tmrt[-13:-10])
        hoursTmrt = int(filepath_tmrt[-9:-7])
        minuTmrt = int(filepath_tmrt[-7:-5])

        if np.where(metdata[:,2] == hoursTmrt).__len__() == 1:
            posMet = np.where(metdata[:,0] == yyyyTmrt) and np.where(metdata[:,1] == doyTmrt) and np.where(metdata[:,2] == hoursTmrt)  
        else:
            posMet = np.where(metdata[:,0] == yyyyTmrt) and np.where(metdata[:,1] == doyTmrt) and np.where(metdata[:,2] == hoursTmrt) and np.where(metdata[:,3] == minuTmrt) 
        
        Ta = metdata[posMet, 11][0][0]
        RH = metdata[posMet, 10][0][0]
        YYYY = metdata[posMet, 0][0][0]
        jday = metdata[posMet, 1][0][0]
        hours = metdata[posMet, 2][0][0]
        minu = metdata[posMet, 3][0][0]

        # Timestamp
        current_time = str((pd.to_datetime(YYYY, format='%Y') + pd.to_timedelta(jday - 1, unit='d') + 
                 pd.to_timedelta(hours, unit='h') + pd.to_timedelta(minu, unit='m')))
        
        feedback.setProgressText('Estimating ' + thermal_index + ' on ' + current_time)
        feedback.setProgressText('Location: ' + str(np.around(lat, decimals=2)) + ' latitude,  ' + str(np.around(lon, decimals=2)) + ' longitude')
        feedback.setProgressText("Air temperature derived from meteorological data is: " + str(Ta) + ' ' + u'\N{DEGREE SIGN}C')
        feedback.setProgressText("Relative Humidity derived from meteorological data is: " + str(RH) + '%')
        feedback.setProgressText("Incoming shortwave radiation derived from meteorological data is: " + str(metdata[posMet[0], 14][0]) + u' Wm\u00b2')
        # feedback.setProgressText("Wind speed derived from meteorological data is: " + str(metdata[posMet[0], 9][0]) + ' m/s')

        # Loading Kup, Kdown, Kdiff, Lup and Ldown if calculating COMFA
        if tcType == 2:
            filepath = filepath_tmrt.split('Tmrt')
            # filepath = os.path.dirname(filepath_tmrt) # issue 31  filepath_tmrt.split('Tmrt')
            # Load Kup, Kdown, Lup, Ldown grids
            Kup, rows, cols = load_grid(filepath[0] + 'Kup' + filepath[1], feedback)
            Kdown, _, __ = load_grid(filepath[0] + 'Kdown' + filepath[1], feedback)
            Kdiff, _, __ = load_grid(filepath[0] + 'Kdiff' + filepath[1], feedback)
            Lup, _, __ = load_grid(filepath[0] + 'Lup' + filepath[1], feedback)
            Ldown, _, __ = load_grid(filepath[0] + 'Ldown' + filepath[1], feedback)
            settingsSolweig = np.loadtxt(filepath[0] + '/treeplantersettings.txt', skiprows=1, delimiter=' ')
            UTC = int(settingsSolweig[0])
            alt = settingsSolweig[12]
            location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
        else:
            # Load Tmrt grid for PET and UTCI
            tmrtGrid, rows, cols = load_grid(filepath_tmrt, feedback)

        # Load wind speed grid (URock)
        provider = ws.dataProvider()
        filename_ws = str(provider.dataSourceUri())
        wsGrid, rows2, cols2 = load_grid(filename_ws, feedback)

        rows3 = build.shape[0] # Number of rows in building grid
        cols3 = build.shape[1] # Number of columns in building grid

        if not (rows == rows2) & (cols == cols2):
            raise QgsProcessingException("Error: Wind speed raster not same domain as Tmrt raster: All rasters must be of same extent and resolution")

        if not (rows == rows3) & (cols == cols3):
            raise QgsProcessingException("Error: Buildings raster not same domain as Tmrt raster: All rasters must be of same extent and resolution")

        # Physiological variables for PET and COMFA
        mbody = self.parameterAsDouble(parameters, self.WEIGHT, context) # Body weight in kg
        clo = self.parameterAsDouble(parameters, self.CLO, context) # Clothing in clo
        age = self.parameterAsDouble(parameters, self.AGE, context) # Age in years
        activity = self.parameterAsDouble(parameters, self.WEIGHT, context) # Activity in watt
        sex = self.parameterAsInt(parameters, self.SEX, context) + 1 # Sex, #TODO CHECK SO SAME FOR PET AND COMFA

        if tcType == 0:
            feedback.setProgressText("Calculating PET for all ground level pixels")
            # Other PET variables
            ht = self.parameterAsDouble(parameters, self.HEIGHT, context) / 100. # Body height in meters

            pet.mbody = mbody
            pet.age = age
            pet.height = ht
            pet.activity = activity
            pet.sex = sex
            pet.clo = clo

            # Wind speed
            #WsPET = (10. / sensorheight) ** 0.2 * wsGrid

            result = pet.calculate_PET_grid(Ta, RH, tmrtGrid, wsGrid, pet, feedback)

        elif tcType == 1:
            feedback.setProgressText("Calculating UTCI for all ground level pixels")
            # Recalculating wind speed based on power law
            WsUTCI = (10. / sensorheight) ** 0.2 * wsGrid
            result = utci.utci_calculator_grid(Ta, RH, tmrtGrid, WsUTCI, feedback)

        elif tcType == 2:
            # If True = COMFA-kid (Cheng & Brown, 2020), if False = regular COMFA
            if comfa_kid:
                feedback.setProgressText("Calculating energy balance (COMFA-kid (Cheng and Brown, 2020)) for all ground level pixels")                
            else:
                feedback.setProgressText("Calculating energy balance (COMFA) for all ground level pixels")
            # COMFA parameters
            # Settings
            # Atr = 0.7 # atmospheric transmittance
            alpha = 0.37 # Albedo of the cylinder
            emis = 0.95 # Emissivity of the cylinder
            Aeff = 0.78 # Effective area of body. 0.78 for standing from Campbell and Normal (1998) (0.70 for sitting)                        
            L = 0.1 # length of cylinder (cm)
            D = 0.01 # Diameter of cylinder (cm)
            ht = self.parameterAsDouble(parameters, self.HEIGHT, context) # Height of person in cm
            # Calculate COMFA radiation using SOLWEIG output L, D, Lin, Lup, Kin, Kup, emis, alpha, Aeff, Kd, metdata, location, utc
            Rabs = COMFA_rad(L, D, Ldown, Lup, Kdown, Kup, emis, alpha, Aeff, Kdiff, metdata[posMet[0],:], location, UTC)
            # Recalculating wind speed based on powerlaw
            WsCOMFA = (10 / sensorheight) ** 0.2 * wsGrid
            # Mact (Wm-2) based on Cheng & Brown (2020)
            Mact, Mact_PET = COMFA_Mact(mbody, ht, sex, age, activity, 'W')
            # Activity speed
            va = 0.
            # Rco: rco = Icl * 186.6. Or convert from conductivity based on 1 Clo = 0.1555 (m 2K W-1).
            rco = clo * 186.6
            # Rcvo: Can convert from a conductivity value of Icl, where Icl is found in tables for specific clothing ensembles from ISO (2007). 
            # 1) convert from Icl to Re,cl (m2 kPa W-1):  Re,cl = Icl (m 2K W-1)*0.18 (constant from ISO 9920 pg 12 based on 1 or 2-layer clothing ensembles). 
            # 2) Convert Re,cl to rcvo:  rcvo = Re,cl*18,400, where 18400 is a conversion factor from Re,cl (vapour resistance, in m2kPaW-1) 		
            # to rcvo, using Lv = 2.5*106 J kg-1, rho = 1.16 kg m-3, and Pa = 98kPa.
            rcvo = clo * 0.18 * 18400

            MET = np.zeros((Rabs.shape[0], Rabs.shape[1]))
            CONV = np.zeros((Rabs.shape[0], Rabs.shape[1]))
            EVAP = np.zeros((Rabs.shape[0], Rabs.shape[1]))
            TREMIT = np.zeros((Rabs.shape[0], Rabs.shape[1]))
            for y in np.arange(Rabs.shape[0]):
                for x in np.arange(Rabs.shape[1]):
                    MET[y,x], CONV[y,x], EVAP[y,x], TREMIT[y,x] = COMFA_BUDGET(Mact, Ta, RH, WsCOMFA[y,x], va, rco, rcvo, mbody, ht, age, comfa_kid)
            result = MET + Rabs - CONV - EVAP - TREMIT     
        
        result[build == 0] = -9999

        saveraster(gdal_buildings, outputRaster, result)

        feedback.setProgressText("Processing finished.")

        return {self.TC_OUT: outputRaster}
    
    def name(self):
        return 'Outdoor Thermal Comfort: Spatial Thermal Comfort'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Post-Processor'

    def shortHelpString(self):
        return self.tr('The <b>Spatial TC (Thermal Comfort)</b> plugin can be used to create maps of spatial variation of thermal comfort indices by exploiting output from the SOLWEIG model together with the URock model in UMEP. All tools need to be executed via the Processing toolbox in QGIS and no filenames and folder structures should be changed. <br>'
        '\n'
        'The <b>Tmrt raster</b> should be located in the SOLWEIG output folder. Remember to tick "Save necessary raster(s) for the TreePlanter and Spatial TC tools" when running SOLWEIG.'
        '\n'
        '--------------\n'
        'Full manual available via the <b>Help</b>-button.')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/post_processor/Outdoor%20Thermal%20Comfort%20SOLWEIG%20Analyzer.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Post-Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/stc.png")
        return icon

    def createInstance(self):
        return ProcessingSpatialTCAlgorithm()
