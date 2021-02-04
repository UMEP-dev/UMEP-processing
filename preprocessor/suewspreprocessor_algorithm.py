# -*- coding: utf-8 -*-

__author__ = 'Fredrik Lindberg'
__date__ = '2021-02-04'
__copyright__ = '(C) 2021 by Fredrik Lindberg'

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterString,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterDefinition,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterField,
                       QgsProcessingException,
                       QgsVectorLayer,
                       QgsFeature,
                       QgsVectorFileWriter,
                       QgsVectorDataProvider,
                       QgsField)

from qgis.PyQt.QtGui import QIcon
from osgeo import gdal, osr, ogr
from osgeo.gdalconst import *
import os
import numpy as np
import inspect
from pathlib import Path
import sys
from ..util import misc


class ProcessingSUEWSPreprocessorAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm is an altered version of SUEWS Prepare for the processing toolbox
    """

    INPUT_POLYGONLAYER = 'INPUT_POLYGONLAYER'
    ID_FIELD = 'ID_FIELD'
    LOD0 = 'LOD0'
    LOD1 = 'LOD1'
    LOD2_GRASS = 'LOD2_GRASS'
    LOD2_SOIL = 'LOD2"_SOIL'
    INPUT_BUILD = 'INPUT_BUILD'
    INPUT_LC = 'INPUT_LC'
    INPUT_DEC = 'INPUT_DEC'
    METFILE = 'METFILE'
    POP_NIGHT = 'POP_NIGHT'
    POP_DAY = 'POP_DAY'
    DAYLIGHT_START = 'DAYLIGHT_START'
    DAYLIGHT_END = 'DAYLIGHT_END'
    FILE_CODE = 'FILE_CODE'
    UTC = 'UTC'
    LEAF_CYCLE = 'LEAF_CYCLE'
    SOIL_STATE = 'SOIL_STATE'
    OUTPUT_DIR = 'OUTPUT_DIR'
    
    def initAlgorithm(self, config):
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT_POLYGONLAYER,
                                                            self.tr('Vector polygon grid'),
                                                            [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterField(self.ID_FIELD,
                                                      self.tr('ID field'),
                                                      '',
                                                      self.INPUT_POLYGONLAYER,
                                                      QgsProcessingParameterField.Numeric))
        self.lod0Opt = ((self.tr('Urban/dense'), '0'),
                        (self.tr('Urban/medium'), '1'),
                        (self.tr('Urban/low'), '2'),
                        (self.tr('Non-urban'), '3'))
        self.addParameter(QgsProcessingParameterEnum(self.LOD0,
                                                    self.tr('General (first guess) area description'),
                                                    options=[i[0] for i in self.lod0Opt],
                                                    defaultValue=0))
        self.lod1Opt = ((self.tr('London/UK'), '0'),
                        (self.tr('Helsinki/Finland'), '1'),
                        (self.tr('Vancouver/Cananda'), '2'),
                        (self.tr('Sacramento/US'), '3'),
                        (self.tr('Beijing/China'), '4'))
        self.addParameter(QgsProcessingParameterEnum(self.LOD1,
                                                    self.tr('Improved regional parameterisation'),
                                                    options=[i[0] for i in self.lod1Opt],
                                                    optional=True))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_BUILD,
                                                     self.tr('Building morphology file'), 
                                                     extension='txt'))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_LC,
                                                     self.tr('Land cover fraction file'),
                                                     extension='txt'))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_DEC,
                                                     self.tr('Tree morphology file'),
                                                     extension='txt'))
        # self.addParameter(QgsProcessingParameterFile(self.INPUT_CON,
        #     self.tr('Tree morphology file'), extension = 'txt'))
        self.addParameter(QgsProcessingParameterFile(self.METFILE,
                                                     self.tr('Meteorological forcing file'),
                                                     extension='txt'))
        self.addParameter(QgsProcessingParameterField(self.POP_NIGHT,
                                                      self.tr('Population density (Night-time, census data)'),
                                                      '',
                                                      self.INPUT_POLYGONLAYER,
                                                      QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterField(self.POP_DAY,
                                                      self.tr('Population density (daytime)'),
                                                      '',
                                                      self.INPUT_POLYGONLAYER,
                                                      QgsProcessingParameterField.Numeric,
                                                      optional=True))
        self.addParameter(QgsProcessingParameterNumber(self.DAYLIGHT_START,
                                                      self.tr('Start of the day light savings'),
                                                      QgsProcessingParameterNumber.Integer,
                                                      QVariant(85),
                                                      False,
                                                      minValue=1,
                                                      maxValue=365))
        self.addParameter(QgsProcessingParameterNumber(self.DAYLIGHT_END,
                                                       self.tr('End of the day light savings'),
                                                       QgsProcessingParameterNumber.Integer,
                                                       QVariant(302),
                                                       False,
                                                       minValue=1,
                                                       maxValue=365))
        self.addParameter(QgsProcessingParameterNumber(self.UTC,
                                                       self.tr('UTC offset'),
                                                       QgsProcessingParameterNumber.Integer,
                                                       QVariant(0),
                                                       False,
                                                       minValue=-12,
                                                       maxValue=12))
        self.addParameter(QgsProcessingParameterString(self.FILE_CODE,
                                                       self.tr('File code')))

        # Advanced parameters
        self.leafnum = ((self.tr('Winter (0%)'), '0'),
                        (self.tr('Early Spring (25%)'), '1'),
                        (self.tr('Spring (50%)'), '2'),
                        (self.tr('Late Spring (75%)'), '3'),
                        (self.tr('Summer (100%)'), '4'),
                        (self.tr('Early Autumn (75%)'), '5'),
                        (self.tr('Autumn (50%)'), '6'),
                        (self.tr('Late Autumn (25%)'), '7'))
        leafstate = QgsProcessingParameterEnum(self.LEAF_CYCLE,
                                                self.tr("Leaf cycle at model start"),
                                                defaultValue=0,
                                                options=[i[0] for i in self.leafnum])
        leafstate.setFlags(leafstate.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(leafstate)
        soilstate = QgsProcessingParameterNumber(self.SOIL_STATE,
                                                self.tr('Soil moisture state (%)'),
                                                      QgsProcessingParameterNumber.Integer,
                                                      QVariant(100),
                                                      False,
                                                      minValue=0,
                                                      maxValue=100)
        soilstate.setFlags(soilstate.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(soilstate)
        self.lod2_grass = ((self.tr('Long'), '0'),
                          (self.tr('Short'), '1'),
                          (self.tr('Wheat'), '2'),
                          (self.tr('Rice'), '3'))
        lod2_grass = QgsProcessingParameterEnum(self.LOD2_GRASS,
                                                self.tr("Type of dominating grass land cover"),
                                                optional=True,
                                                options=[i[0] for i in self.lod2_grass])
        lod2_grass.setFlags(lod2_grass.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(lod2_grass)
        self.lod2_soil = ((self.tr('Clay'), '0'),
                          (self.tr('Loam'), '1'),
                          (self.tr('Silt'), '2'),
                          (self.tr('Sand'), '3'))
        lod2_soil = QgsProcessingParameterEnum(self.LOD2_SOIL,
                                                self.tr("Type of dominating soil land cover"),
                                                optional=True,
                                                options=[i[0] for i in self.lod2_soil])
        lod2_soil.setFlags(lod2_soil.flags() | QgsProcessingParameterDefinition.FlagAdvanced)
        self.addParameter(lod2_soil)

        # output
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_DIR,
                                                                  self.tr('Output folder')))


    def processAlgorithm(self, parameters, context, feedback):
        # InputParameters
        inputPolygonlayer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGONLAYER, context)
        idField = self.parameterAsFields(parameters, self.ID_FIELD, context)
        # lod1 = self.parameterAsString(parameters, self.LOD0, context)
        # inputInterval = self.parameterAsDouble(parameters, self.INPUT_INTERVAL, context)
        # useDsmBuild = self.parameterAsBool(parameters, self.USE_DSMBUILD, context)
        # outputDir = self.parameterAsString(parameters, self.OUTPUT_DIR, context)
        
        if parameters['OUTPUT_DIR'] == 'TEMPORARY_OUTPUT':
            if not (os.path.isdir(outputDir)):
                os.mkdir(outputDir)

        poly = inputPolygonlayer
        poly_field = idField
        vlayer = inputPolygonlayer
        prov = vlayer.dataProvider()
        fields = prov.fields()
        idx = vlayer.fields().indexFromName(poly_field[0])
        dir_poly = self.plugin_dir + '/data/poly_temp.shp'
        nGrids = vlayer.featureCount()
        index = 1
        
        feedback.setProgressText("Number of grids to analyse: " + str(nGrids))

        feedback.setProgressText("Adding result to layer attribute table") 


        return {self.OUTPUT_DIR: outputDir}


    
    def name(self):
        return 'Urban Energy Balance: SUEWS Pre-processor'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Pre-Processor'

    def shortHelpString(self):
        return self.tr('UNDER CONSTRUCTION.\n'
        '-------------\n')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/pre-processor/SUEWS%20Prepare.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/SuewsLogo.png")
        return icon

    def createInstance(self):
        return ProcessingSUEWSPreprocessorAlgorithm()