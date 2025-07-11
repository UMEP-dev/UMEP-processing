# -*- coding: utf-8 -*-

"""
/***************************************************************************
 ProcessingUMEP
                                 A QGIS plugin
 UMEP for processing toolbox
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2020-04-02
        copyright            : (C) 2020 by Fredrik Lindberg
        email                : fredrikl@gvc.gu.se
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Fredrik Lindberg'
__date__ = '2020-04-02'
__copyright__ = '(C) 2020 by Fredrik Lindberg'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication, QDate, QTime, Qt, QVariant
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingParameterFileDestination,
                       QgsProcessingParameterFile,
                       QgsProcessingException,
                       QgsProcessingParameterRasterLayer)
from processing.gui.wrappers import WidgetWrapper
from qgis.PyQt.QtWidgets import QDateEdit, QTimeEdit
import numpy as np
from osgeo import gdal, osr
from osgeo.gdalconst import *
import os
from qgis.PyQt.QtGui import QIcon
import inspect
from pathlib import Path
from ..functions.SEBEfiles import SEBE_2015a_calc_forprocessing as sebe
from ..functions.SEBEfiles.sunmapcreator_2015a import sunmapcreator_2015a
from ..functions.SEBEfiles import WriteMetaDataSEBE
from ..util.SEBESOLWEIGCommonFiles.Solweig_v2015_metdata_noload import Solweig_2015a_metdata_noload
from ..util.misc import get_ders, saveraster, createTSlist


class ProcessingSEBEAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm is a processing version of SEBE
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.
    
    INPUT_DSM = 'INPUT_DSM'
    INPUT_CDSM = 'INPUT_CDSM'
    INPUT_TDSM = 'INPUT_TDSM'
    INPUT_HEIGHT = 'INPUT_HEIGHT'
    INPUT_ASPECT = 'INPUT_ASPECT'
    TRANS_VEG = 'TRANS_VEG'
    INPUT_THEIGHT = 'INPUT_THEIGHT'
    UTC = 'UTC'
    ALBEDO = 'ALBEDO'
    ONLYGLOBAL = 'ONLYGLOBAL'
    INPUT_MET = 'INPUTMET'
    SAVESKYIRR = 'SAVESKYIRR'
    IRR_FILE = 'IRR_FILE'
    OUTPUT_DIR = 'OUTPUT_DIR'
    OUTPUT_ROOF = 'OUTPUT_ROOF'
    

    def initAlgorithm(self, config):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DSM,
            self.tr('Input building and ground DSM'), None, False))
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_CDSM,
            self.tr('Vegetation Canopy DSM'), '', True))
        self.addParameter(QgsProcessingParameterNumber(self.TRANS_VEG, 
            self.tr('Transmissivity of light through vegetation (%):'), 
            QgsProcessingParameterNumber.Integer,
            QVariant(3), True, minValue=0, maxValue=100))
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_TDSM,
            self.tr('Vegetation Trunk zone DSM'), '', True))
        self.addParameter(QgsProcessingParameterNumber(self.INPUT_THEIGHT, 
            self.tr("Trunk zone height (percent of Canopy Height)"), 
            QgsProcessingParameterNumber.Double,
            QVariant(25.0), True, minValue=0.1, maxValue=99.9))
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_HEIGHT,
            self.tr('Wall height raster'), '', False))
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_ASPECT,
            self.tr('Wall aspect raster'), '', False))
        self.addParameter(QgsProcessingParameterNumber(self.ALBEDO, 
            self.tr('Albedo:'), QgsProcessingParameterNumber.Double,
            QVariant(0.15), False, minValue=0, maxValue=1))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_MET,
            self.tr('Input meteorological file'), extension = 'txt'))
        self.addParameter(QgsProcessingParameterBoolean(self.ONLYGLOBAL,
            self.tr("Estimate diffuse and direct shortwave radiation from global radiation"), defaultValue=False))
        self.sorted_utclist, _ = createTSlist() #response to #104
        lista = []
        for i in self.sorted_utclist:
            lista.append((str(i['utc_offset_str']), str(i['utc_offset'])))
        self.addParameter(QgsProcessingParameterEnum(self.UTC,
                                                self.tr('Coordinated Universal Time (UTC) '),
                                                options=[i[0] for i in lista],
                                                defaultValue=14))
        self.addParameter(QgsProcessingParameterBoolean(self.SAVESKYIRR,
            self.tr("Save sky irradiance distribution"), defaultValue=False))
        self.addParameter(QgsProcessingParameterFileDestination(self.IRR_FILE,
             self.tr('Sky irradiance distribution'), self.tr('txt files (*.txt)')))
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_DIR,
                                                     'Output folder'))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_ROOF,
            self.tr("Roof irradiance raster (kWh)"), optional=True,
            createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        # InputParameters
        outputDir = self.parameterAsString(parameters, self.OUTPUT_DIR, context)
        dsmlayer = self.parameterAsRasterLayer(parameters, self.INPUT_DSM, context) 
        transVeg = self.parameterAsDouble(parameters, self.TRANS_VEG, context) 
        vegdsm = self.parameterAsRasterLayer(parameters, self.INPUT_CDSM, context) 
        vegdsm2 = self.parameterAsRasterLayer(parameters, self.INPUT_TDSM, context) 
        whlayer = self.parameterAsRasterLayer(parameters, self.INPUT_HEIGHT, context) 
        walayer = self.parameterAsRasterLayer(parameters, self.INPUT_ASPECT, context) 
        trunkr = self.parameterAsDouble(parameters, self.INPUT_THEIGHT, context) 
        onlyglobal = self.parameterAsBool(parameters, self.ONLYGLOBAL, context)
        utcpos = self.parameterAsString(parameters, self.UTC, context)
        albedo = self.parameterAsDouble(parameters, self.ALBEDO, context)
        inputMet = self.parameterAsString(parameters, self.INPUT_MET, context)
        saveskyirr = self.parameterAsBool(parameters, self.SAVESKYIRR, context)
        irrFile = self.parameterAsFileOutput(parameters, self.IRR_FILE, context)
        outputRoof = self.parameterAsOutputLayer(parameters, self.OUTPUT_ROOF, context)

        if parameters['OUTPUT_DIR'] == 'TEMPORARY_OUTPUT':
            if not (os.path.isdir(outputDir)):
                os.mkdir(outputDir)

        provider = dsmlayer.dataProvider()
        filepath_dsm = str(provider.dataSourceUri())
        self.gdal_dsm = gdal.Open(filepath_dsm)
        self.dsm = self.gdal_dsm.ReadAsArray().astype(float)
        sizex = self.dsm.shape[0]
        sizey = self.dsm.shape[1]

        # response to issue #85
        nd = self.gdal_dsm.GetRasterBand(1).GetNoDataValue()
        self.dsm[self.dsm == nd] = 0.
        if self.dsm.min() < 0:
            self.dsm = self.dsm + np.abs(self.dsm.min())

        # response to issue #104
        self.sorted_utclist
        utc = self.sorted_utclist[int(utcpos)]['utc_offset']

        # Get latlon from grid coordinate system
        old_cs = osr.SpatialReference()
        dsm_ref = dsmlayer.crs().toWkt()
        old_cs.ImportFromWkt(dsm_ref)

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
        width = self.gdal_dsm.RasterXSize
        height = self.gdal_dsm.RasterYSize
        geotransform = self.gdal_dsm.GetGeoTransform()
        minx = geotransform[0]
        miny = geotransform[3] + width*geotransform[4] + height*geotransform[5]
        lonlat = transform.TransformPoint(minx, miny)
        gdalver = float(gdal.__version__[0])
        if gdalver == 3.:
            lon = lonlat[1] #changed to gdal 3
            lat = lonlat[0] #changed to gdal 3
        else:
            lon = lonlat[0] #changed to gdal 2
            lat = lonlat[1] #changed to gdal 2
        self.scale = 1 / geotransform[1]

        feedback.setProgressText('Longitude derived from DSM: ' + str(lon))
        feedback.setProgressText('Latitude derived from DSM: ' + str(lat))

        trunkfile = 0
        trunkratio = 0
        psi = transVeg / 100.0

        if vegdsm:
            usevegdem = 1
            feedback.setProgressText('Vegetation scheme activated')
            # vegdsm = self.parameterAsRasterLayer(parameters, self.INPUT_CDSM, context) 
            # if vegdsm is None:
            #     raise QgsProcessingException("Error: No valid vegetation DSM selected")

            # load raster
            gdal.AllRegister()
            provider = vegdsm.dataProvider()
            filePathOld = str(provider.dataSourceUri())
            dataSet = gdal.Open(filePathOld)
            vegdsm = dataSet.ReadAsArray().astype(float)
            filePath_cdsm = filePathOld
            vegsizex = vegdsm.shape[0]
            vegsizey = vegdsm.shape[1]

            if not (vegsizex == sizex) & (vegsizey == sizey):
                raise QgsProcessingException("Error in Vegetation Canopy DSM: All rasters must be of same extent and resolution")

            if vegdsm2:
                # vegdsm2 = self.parameterAsRasterLayer(parameters, self.INPUT_TDSM, context) 

                # if vegdsm2 is None:
                #     raise QgsProcessingException("Error: No valid Trunk zone DSM selected")

                # load raster
                gdal.AllRegister()
                provider = vegdsm2.dataProvider()
                filePathOld = str(provider.dataSourceUri())
                filePath_tdsm = filePathOld
                dataSet = gdal.Open(filePathOld)
                vegdsm2 = dataSet.ReadAsArray().astype(float)
            else:
                trunkratio = trunkr / 100.0
                vegdsm2 = vegdsm * trunkratio
                filePath_tdsm = None

            vegsizex = vegdsm2.shape[0]
            vegsizey = vegdsm2.shape[1]

            if not (vegsizex == sizex) & (vegsizey == sizey):  # &
                raise QgsProcessingException("Error in Trunk Zone DSM: All rasters must be of same extent and resolution")
        else:
            vegdsm = 0
            vegdsm2 = 0
            usevegdem = 0
            filePath_cdsm = None
            filePath_tdsm = None

        # wall height layer
        # if whlayer is None:
        #     raise QgsProcessingException("Error: No valid wall height raster layer is selected")
        provider = whlayer.dataProvider()
        filepath_wh = str(provider.dataSourceUri())
        self.gdal_wh = gdal.Open(filepath_wh)
        wheight = self.gdal_wh.ReadAsArray().astype(float)
        vhsizex = wheight.shape[0]
        vhsizey = wheight.shape[1]
        if not (vhsizex == sizex) & (vhsizey == sizey):
            raise QgsProcessingException("Error in Wall height raster: All rasters must be of same extent and resolution")

        wallmaxheight = self.gdal_wh.GetRasterBand(1).GetStatistics(True,True)[1]

        # wall aspectlayer
        # if walayer is None:
        #     raise QgsProcessingException("Error: No valid wall aspect raster layer is selected")
        provider = walayer.dataProvider()
        filepath_wa = str(provider.dataSourceUri())
        self.gdal_wa = gdal.Open(filepath_wa)
        waspect = self.gdal_wa.ReadAsArray().astype(float)
        vasizex = waspect.shape[0]
        vasizey = waspect.shape[1]
        if not (vasizex == sizex) & (vasizey == sizey):
            raise QgsProcessingException("Error in Wall aspect raster: All rasters must be of same extent and resolution")

        voxelheight = geotransform[1]  # float

        # Metdata
        headernum = 1
        delim = ' '

        try:
            self.metdata = np.loadtxt(inputMet, skiprows=headernum, delimiter=delim)
        except:
            QgsProcessingException("Error: Make sure format of meteorological file is correct. You can"
                                                        "prepare your data by using 'Prepare Existing Data' in "
                                                        "the Pre-processor")

        testwhere = np.where((self.metdata[:, 14] < 0.0) | (self.metdata[:, 14] > 1300.0))
        if testwhere[0].__len__() > 0:
             QgsProcessingException("Error: Kdown - beyond what is expected at line: " + str(testwhere[0] + 1))

        if self.metdata.shape[1] == 24:
            feedback.setProgressText("Meteorological data succefully loaded")
        else:
            QgsProcessingException("Error: Wrong number of columns in meteorological data. You can "
                                                        "prepare your data by using 'Prepare Existing Data' in "
                                                        "the Pre-processor")

        alt = np.median(self.dsm)
        if alt < 0:
            alt = 3

        feedback.setProgressText("Calculating sun positions for each time step")
        location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
        YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = \
            Solweig_2015a_metdata_noload(self.metdata, location, utc)

        feedback.setProgressText("Distributing irradiance on sky vault")
        output = {'energymonth': 0, 'energyyear': 1, 'suitmap': 0}
        radmatI, radmatD, radmatR = sunmapcreator_2015a(self.metdata, altitude, azimuth,
                                                        onlyglobal, output, jday, albedo, location, zen)

        if saveskyirr:
            metout = np.zeros((145, 4))
            metout[:, 0] = radmatI[:, 0]
            metout[:, 1] = radmatI[:, 1]
            metout[:, 2] = radmatI[:, 2]
            metout[:, 3] = radmatD[:, 2]
            header = '%altitude azimuth radI radD'
            numformat = '%6.2f %6.2f %6.2f %6.2f'
            np.savetxt(irrFile, metout, fmt=numformat, header=header, comments='')

        building_slope, building_aspect = get_ders(self.dsm, self.scale)

        WriteMetaDataSEBE.writeRunInfo(outputDir, filepath_dsm, self.gdal_dsm, usevegdem,
                                        filePath_cdsm, trunkfile, filePath_tdsm, lat, lon, utc,
                                        inputMet, albedo, onlyglobal, trunkratio, psi, sizex, sizey)

        # Main function
        feedback.setProgressText("Executing main model")
        seberesult = sebe.SEBE_2015a_calc(self.dsm, self.scale, building_slope,
                    building_aspect, voxelheight, sizey, sizex, vegdsm, vegdsm2, wheight,
                    waspect, albedo, psi, radmatI, radmatD, radmatR, usevegdem, feedback, wallmaxheight)

        Energyyearroof = seberesult["Energyyearroof"]
        Energyyearwall = seberesult["Energyyearwall"]
        vegdata = seberesult["vegdata"]

        feedback.setProgressText("SEBE: Model calculation finished. Saving to disk")

        
        if outputRoof:
            saveraster(self.gdal_dsm, outputRoof, Energyyearroof)

        saveraster(self.gdal_dsm, outputDir + '/dsm.tif', self.dsm)
        filenameroof = outputDir + '/Energyyearroof.tif'
        saveraster(self.gdal_dsm, filenameroof, Energyyearroof)
        filenamewall = outputDir + '/Energyyearwall.txt'
        header = '%row col irradiance'
        numformat = '%4d %4d ' + '%6.2f ' * (Energyyearwall.shape[1] - 2)
        np.savetxt(filenamewall, Energyyearwall, fmt=numformat, header=header, comments='')
        if usevegdem == 1:
            filenamewall = outputDir + '/Vegetationdata.txt'
            header = '%row col height'
            numformat = '%4d %4d %6.2f'
            np.savetxt(filenamewall, vegdata, fmt=numformat, header=header, comments='')

        return {self.OUTPUT_DIR: outputDir, self.IRR_FILE: irrFile, self.OUTPUT_ROOF: outputRoof}
    
    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Solar Radiation: Solar Energy of Builing Envelopes (SEBE)'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Processor'

    def shortHelpString(self):
        return self.tr('The SEBE plugin (Solar Energy on Building Envelopes) can be used to calculate pixel wise potential solar energy using ground and building digital surface models (DSM). SEBE is also able to estimate irradiance on building walls. Optionally, vegetation DSMs could also be used.<br>'
                        '--------------\n'
                        'Full manual available via the <b>Help</b>-button.')
    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/processor/Solar%20Radiation%20Solar%20Energy%20on%20Building%20Envelopes%20(SEBE).html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/sebeicon.png")
        return icon

    def createInstance(self):
        return ProcessingSEBEAlgorithm()
