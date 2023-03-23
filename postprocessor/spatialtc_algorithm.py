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
import inspect
from pathlib import Path
import sys
from ..util.misc import saveraster
from ..functions.SOLWEIGpython import PET_calculations as pet

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

    #PET parameters
    AGE = 'AGE'
    ACTIVITY = 'ACTIVITY'
    CLO = 'CLO'
    WEIGHT = 'WEIGHT'
    HEIGHT = 'HEIGHT'
    SEX = 'SEX'
    SENSOR_HEIGHT = 'SENSOR_HEIGHT'

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
                                                            self.tr('Pedestrain wind speed raster (URock)'),
                                                             '', 
                                                             optional=False))
        self.addParameter(QgsProcessingParameterRasterLayer(self.BUILDING_MAP,
                                                            self.tr('Raster to exclude building pixels from analysis (SOLWEIG)'),
                                                             '', 
                                                             optional=False))
        self.addParameter(QgsProcessingParameterFile(self.METDATA,
                                                    self.tr('Input meteorological file (.txt). Same as used for Tmrt map.'), extension='txt',
                                                    optional=False))
        self.varType = ((self.tr('Physilogical Equivalent Temperature (PET)'), '0'),
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

        # Output
        self.addParameter(QgsProcessingParameterRasterDestination(self.TC_OUT,
                                                                  self.tr("Output thermal comfort raster"),
                                                                  None,
                                                                  optional=False))


    def processAlgorithm(self, parameters, context, feedback):
        
        # InputParameters
        # solweigDir = self.parameterAsString(parameters, self.SOLWEIG_DIR, context)
        
        tmrt = self.parameterAsRasterLayer(parameters, self.TMRT_MAP, context)
        ws = self.parameterAsRasterLayer(parameters, self.UROCK_MAP, context)
        buildings = self.parameterAsRasterLayer(parameters, self.BUILDING_MAP, context)
        tcTypeStr = self.parameterAsString(parameters, self.TC_TYPE, context) 
        metdata = self.parameterAsString(parameters, self.METDATA, context)

        outputRaster = self.parameterAsOutputLayer(parameters, self.TC_OUT, context)

        mbody = None
        ht = None
        clo = None
        age = None
        activity = None
        sex = None

        feedback.setProgressText("Initializing...")
        # LOAD raster data
        provider = buildings.dataProvider()
        filepath_dsm = str(provider.dataSourceUri())
        gdal_dsm = gdal.Open(filepath_dsm)
        build = gdal_dsm.ReadAsArray().astype(float)

        # LOAD Metdata
        try:
            metdata = np.loadtxt(metdata,skiprows=1, delimiter=' ')
        except:
            raise QgsProcessingException("Error: Make sure format of meteorological file is correct. You can"
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
        provider = tmrt.dataProvider()
        filepath_tmrt = str(provider.dataSourceUri())
        yyyyTmrt = int(filepath_tmrt[-18:-14])
        doyTmrt = int(filepath_tmrt[-13:-10])
        hoursTmrt = int(filepath_tmrt[-9:-7])
        minuTmrt = int(filepath_tmrt[-7:-5])

        if np.where(metdata[:,2] == hoursTmrt).__len__() == 1:
            posMet = np.where(metdata[:,0] == yyyyTmrt) and np.where(metdata[:,1] == doyTmrt) and np.where(metdata[:,2] == hoursTmrt)  
        else:
            posMet = np.where(metdata[:,0] == yyyyTmrt) and np.where(metdata[:,1] == doyTmrt) and np.where(metdata[:,2] == hoursTmrt) and np.where(metdata[:,3] == minuTmrt) 
        
        Ta = metdata[posMet, 11][0][0]
        RH = metdata[posMet, 10][0][0]

        feedback.setProgressText("Air temperature derived from meteorological data is: " + str(Ta))
        feedback.setProgressText("Relative Humidity derived from meteorological data is: " + str(RH))

        # load grids
        provider = tmrt.dataProvider()
        filename_trmt = str(provider.dataSourceUri())
        gdal_tmrt = gdal.Open(filename_trmt)
        tmrtGrid = gdal_tmrt.ReadAsArray().astype(float)
        rows = tmrtGrid.shape[0]
        cols = tmrtGrid.shape[1]

        provider = ws.dataProvider()
        filename_ws = str(provider.dataSourceUri())
        gdal_ws = gdal.Open(filename_ws)
        wsGrid = gdal_ws.ReadAsArray().astype(float)
        rows2 = tmrtGrid.shape[0]
        cols2 = tmrtGrid.shape[1]

        rows3 = build.shape[0]
        cols3 = build.shape[1]

        if not (rows == rows2) & (cols == cols2):
            raise QgsProcessingException("Error: Wind speed raster not same domain as Tmrt raster: All rasters must be of same extent and resolution")

        if not (rows == rows3) & (cols == cols3):
            raise QgsProcessingException("Error: Buildings raster not same domain as Tmrt raster: All rasters must be of same extent and resolution")

        tcType = int(tcTypeStr)
        if tcType == 0:
            feedback.setProgressText("Calculating PET for all ground level pixels")
            # Other PET variables
            mbody = self.parameterAsDouble(parameters, self.WEIGHT, context)
            ht = self.parameterAsDouble(parameters, self.HEIGHT, context) / 100.
            clo = self.parameterAsDouble(parameters, self.CLO, context)
            age = self.parameterAsDouble(parameters, self.AGE, context)
            activity = self.parameterAsDouble(parameters, self.WEIGHT, context)
            sex = self.parameterAsInt(parameters, self.SEX, context) + 1

            pet.mbody = mbody
            pet.age = age
            pet.height = ht
            pet.activity = activity
            pet.sex = sex
            pet.clo = clo

            result = pet.calculate_PET_grid(Ta, RH, tmrtGrid, wsGrid, pet, feedback)

        elif tcType == 1:
            raise QgsProcessingException("UTCI not available yet. Coming soon...")
        elif tcType == 2:
            raise QgsProcessingException("COMFA not available yet. Coming soon...")        
        
        result[build == 0] = -9999

        saveraster(gdal_dsm, outputRaster, result)

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
