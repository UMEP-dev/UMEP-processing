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
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterString,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterFolderDestination,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingException,
                       QgsProcessingParameterDateTime,               
                       QgsProcessingParameterRasterLayer)

from processing.gui.wrappers import WidgetWrapper
from qgis.PyQt.QtWidgets import QDateEdit, QTimeEdit
import numpy as np
from osgeo import gdal, osr
from osgeo.gdalconst import *
import os
from ..functions import dailyshading as dsh
from qgis.PyQt.QtGui import QIcon
import inspect
from pathlib import Path
import datetime
from ..util.misc import createTSlist


class ProcessingShadowGeneratorAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm is a processing version of ShadowGenerator
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
    ONE_SHADOW = 'ONE_SHADOW'
    ITERTIME = 'ITERTIME'
    DATEINI = 'DATEINI'
    TIMEINI = 'TIMEINI'
    UTC = 'UTC'
    DST = 'DST'
    OUTPUT_DIR = 'OUTPUT_DIR'
    OUTPUT_FILE = 'OUTPUT_FILE'

    def initAlgorithm(self, config):
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_DSM,
                self.tr('Input building and ground DSM'),
                None,
                False))
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_CDSM,
                self.tr('Vegetation Canopy DSM'),
                '',
                True))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.TRANS_VEG,
                self.tr('Transmissivity of light through vegetation (%):'), 
                QgsProcessingParameterNumber.Integer,
                QVariant(3), 
                True, 
                minValue=0, 
                maxValue=100))
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_TDSM,
                self.tr('Vegetation Trunk zone DSM'), '', True))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.INPUT_THEIGHT,
                self.tr("Trunk zone height (percent of Canopy Height)"),
                QgsProcessingParameterNumber.Double,
                QVariant(25.0),
                True, 
                minValue=0.1, 
                maxValue=99.9))
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_HEIGHT,
                self.tr('Wall height raster (required if facade shadow should be claculated)'),
                '', 
                True))
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_ASPECT,
                self.tr('Wall aspect raster (required if facade shadow should be claculated)'),
                '', 
                True))
        self.sorted_utclist, _ = createTSlist() #Inculde all UTC times
        lista = []
        for i in self.sorted_utclist:
            lista.append((str(i['utc_offset_str']), str(i['utc_offset'])))
        self.addParameter(QgsProcessingParameterEnum(self.UTC,
                                                self.tr('Coordinated Universal Time (UTC) '),
                                                options=[i[0] for i in lista],
                                                defaultValue=15))        
        
        # self.addParameter(
        #     QgsProcessingParameterNumber(
        #         self.UTC,
        #         self.tr('Coordinated Universal Time (UTC) '),
        #         QgsProcessingParameterNumber.Integer,
        #         QVariant(0),
        #         True, 
        #         minValue=-12, 
        #         maxValue=12)) 
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.DST,
                self.tr("Add Daylight savings time"), 
                defaultValue=False))       
        self.addParameter(QgsProcessingParameterDateTime(self.DATEINI,
            self.tr('Date'),
            QgsProcessingParameterDateTime.Date))
        self.addParameter(
            QgsProcessingParameterNumber(
                self.ITERTIME,
                self.tr('Time interval between casting of each shadow (minutes)'),
                QgsProcessingParameterNumber.Integer,
                QVariant(30),
                True, 
                minValue=0.1,
                maxValue=360))
        self.addParameter(
            QgsProcessingParameterBoolean(
                self.ONE_SHADOW,
                self.tr("Cast only one shadow"), 
                defaultValue=False))
        self.addParameter(QgsProcessingParameterDateTime(self.TIMEINI,
            self.tr('Time for single shadow'),
            QgsProcessingParameterDateTime.Time))
        self.addParameter(
            QgsProcessingParameterFolderDestination(
                self.OUTPUT_DIR,
                'Output folder'))
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT_FILE,
                self.tr("Aggregated (or single) shadow raster"), 
                optional=True,
                createByDefault=False))


    def processAlgorithm(self, parameters, context, feedback):
        # InputParameters
        outputDir = self.parameterAsString(parameters, self.OUTPUT_DIR, context)
        outputFile = self.parameterAsOutputLayer(parameters, self.OUTPUT_FILE, context)
        dsmlayer = self.parameterAsRasterLayer(parameters, self.INPUT_DSM, context)
        transVeg = self.parameterAsDouble(parameters, self.TRANS_VEG, context)
        vegdsm = self.parameterAsRasterLayer(parameters, self.INPUT_CDSM, context)
        vegdsm2 = self.parameterAsRasterLayer(parameters, self.INPUT_TDSM, context)
        whlayer = self.parameterAsRasterLayer(parameters, self.INPUT_HEIGHT, context) 
        walayer = self.parameterAsRasterLayer(parameters, self.INPUT_ASPECT, context) 
        trunkr = self.parameterAsDouble(parameters, self.INPUT_THEIGHT, context) 
        utcpos = self.parameterAsString(parameters, self.UTC, context) 
        dst = self.parameterAsBool(parameters, self.DST, context)
        myDate = self.parameterAsString(parameters, self.DATEINI, context)
        oneShadow = self.parameterAsDouble(parameters, self.ONE_SHADOW, context) 
        myTime = self.parameterAsString(parameters, self.TIMEINI, context)
        iterShadow = self.parameterAsDouble(parameters, self.ITERTIME, context)

        if parameters['OUTPUT_DIR'] == 'TEMPORARY_OUTPUT':
            if not os.path.isdir(outputDir):
                os.mkdir(outputDir)

        if dst:
            dst = 1
        else:
            dst = 0

        provider = dsmlayer.dataProvider()
        filepath_dsm = str(provider.dataSourceUri())
        gdal_dsm = gdal.Open(filepath_dsm)
        dsm = gdal_dsm.ReadAsArray().astype(float)

        # response to issue #85
        nd = gdal_dsm.GetRasterBand(1).GetNoDataValue()
        dsm[dsm == nd] = 0.
        if dsm.min() < 0:
            dsm = dsm + np.abs(dsm.min())

        self.sorted_utclist
        utc = self.sorted_utclist[int(utcpos)]['utc_offset']

        sizex = dsm.shape[0]
        sizey = dsm.shape[1]

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

        width = gdal_dsm.RasterXSize
        height = gdal_dsm.RasterYSize
        gt = gdal_dsm.GetGeoTransform()
        minx = gt[0]
        miny = gt[3] + width*gt[4] + height*gt[5]
        lonlat = transform.TransformPoint(minx, miny)
        geotransform = gdal_dsm.GetGeoTransform()
        scale = 1 / geotransform[1]

        gdalver = float(gdal.__version__[0])
        if gdalver >= 3.:
            lon = lonlat[1] #changed to gdal 3
            lat = lonlat[0] #changed to gdal 3
        else:
            lon = lonlat[0] #changed to gdal 2
            lat = lonlat[1] #changed to gdal 2

        feedback.setProgressText('Longitude derived from DSM: ' + str(lon))
        feedback.setProgressText('Latitude derived from DSM: ' + str(lat))

        trans = transVeg / 100.0

        if vegdsm:
            usevegdem = 1
            feedback.setProgressText('Vegetation scheme activated')

            # load raster
            gdal.AllRegister()
            provider = vegdsm.dataProvider()
            filePathOld = str(provider.dataSourceUri())
            dataSet = gdal.Open(filePathOld)
            vegdsm = dataSet.ReadAsArray().astype(float)

            vegsizex = vegdsm.shape[0]
            vegsizey = vegdsm.shape[1]

            if not (vegsizex == sizex) & (vegsizey == sizey):
                raise QgsProcessingException("Error in Vegetation Canopy DSM: All rasters must be of same extent and resolution")

            if vegdsm2:
                gdal.AllRegister()
                provider = vegdsm2.dataProvider()
                filePathOld = str(provider.dataSourceUri())
                dataSet = gdal.Open(filePathOld)
                vegdsm2 = dataSet.ReadAsArray().astype(float)
            else:
                trunkratio = trunkr / 100.0
                vegdsm2 = vegdsm * trunkratio

            vegsizex = vegdsm2.shape[0]
            vegsizey = vegdsm2.shape[1]

            if not (vegsizex == sizex) & (vegsizey == sizey):  # &
                raise QgsProcessingException("Error in Trunk Zone DSM: All rasters must be of same extent and resolution")
        else:
            vegdsm = 0
            vegdsm2 = 0
            usevegdem = 0

        if whlayer and walayer:
            feedback.setProgressText('Facade shadow scheme activated')
            wallsh = 1
            provider = whlayer.dataProvider()
            filepath_wh = str(provider.dataSourceUri())
            self.gdal_wh = gdal.Open(filepath_wh)
            wheight = self.gdal_wh.ReadAsArray().astype(float)
            vhsizex = wheight.shape[0]
            vhsizey = wheight.shape[1]
            if not (vhsizex == sizex) & (vhsizey == sizey):  # &
                raise QgsProcessingException("Error in Wall height raster: All rasters must be of same extent and resolution")

            provider = walayer.dataProvider()
            filepath_wa = str(provider.dataSourceUri())
            self.gdal_wa = gdal.Open(filepath_wa)
            waspect = self.gdal_wa.ReadAsArray().astype(float)
            vasizex = waspect.shape[0]
            vasizey = waspect.shape[1]
            if not (vasizex == sizex) & (vasizey == sizey):
                raise QgsProcessingException("Error in Wall aspect raster: All rasters must be of same extent and resolution")
        else:
            wallsh = 0
            wheight = 0
            waspect = 0

        if outputDir is 'None':
            raise QgsProcessingException("Error: No selected folder")
        else:
            startDate = datetime.datetime.strptime(myDate, '%Y-%m-%d')
            year = startDate.year
            month = startDate.month
            day = startDate.day
            # UTC = utc #self.dlg.spinBoxUTC.value()
            if oneShadow: #self.dlg.shadowCheckBox.isChecked():
                onetime = 1
                onetimetime = datetime.datetime.strptime(myTime, '%H:%M:%S.%f')
                hour = onetimetime.hour
                minu = onetimetime.minute
                sec = onetimetime.second
            else:
                onetime = 0
                hour = 0
                minu = 0
                sec = 0

            tv = [year, month, day, hour, minu, sec]

            timeInterval = iterShadow # self.dlg.intervalTimeEdit.time()
            shadowresult = dsh.dailyshading(dsm, vegdsm, vegdsm2, scale, lon, lat, sizex, sizey, tv, utc, usevegdem,
                                            timeInterval, onetime, feedback, outputDir, gdal_dsm, trans,
                                            dst, wallsh, wheight, waspect)
            
            shfinal = shadowresult["shfinal"]

        if outputFile:
            dsh.saveraster(gdal_dsm, outputFile, shfinal)

        feedback.setProgressText("ShadowGenerator: Shadow grid(s) successfully generated")

        return {self.OUTPUT_DIR: outputDir, self.OUTPUT_FILE: outputFile}
    
    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Solar Radiation: Shadow Generator'

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
        return self.tr('The Shadow generator plugin can be used to generate pixel wise shadow analysis using ground and '
               'building digital surface models (DSM). Optionally, vegetation DSMs could also be used. '
               'The methodology that is used to generate shadows originates from Ratti and Richens (1990) '
               'and is further developed and described in Lindberg and Grimmond (2011).<br>'
               '\n'
               '------------------<br>'
               'Lindberg, F., Grimmond, C.S.B., 2011a. The influence of vegetation and building morphology on shadow patterns and mean radiant temperatures in urban areas: model development and evaluation. Theoret. Appl. Climatol. 105, 311–323 <br>'
               '\n'
               'Ratti CF, Richens P (1999) Urban texture analysis with image processing techniques. In: Proceedings of the CAADFutures99, Atalanta, GA'
               '\n'
               'Full manual available via the <b>Help</b>-button.')


    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/processor/Solar%20Radiation%20Daily%20Shadow%20Pattern.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/ShadowIcon.png")
        return icon

    def createInstance(self):
        return ProcessingShadowGeneratorAlgorithm()


class DateWidget(WidgetWrapper):
    """
    QDateEdit widget with calendar pop up
    """
    def createWidget(self):
        self._combo = QDateEdit() #QDateTimeEdit()
        self._combo.setCalendarPopup(True)

        today = QDate.currentDate()
        self._combo.setDate(today)

        return self._combo

    def value(self):
        date_chosen = self._combo.dateTime()
        return date_chosen.toString(Qt.ISODate)

class TimeWidget(WidgetWrapper):
    """
    QTimeEdit widget with calendar pop up
    """
    def createWidget(self):
        self._combo = QTimeEdit() #QDateTimeEdit()
        self._combo.setCalendarPopup(True)

        today = QTime.currentTime()
        self._combo.setTime(today)

        return self._combo

    def value(self):
        time_chosen = self._combo.dateTime()
        return time_chosen.toString(Qt.ISODate)
