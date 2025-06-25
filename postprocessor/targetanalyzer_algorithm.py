# -*- coding: utf-8 -*-


__author__ = 'Fredrik Lindberg'
__date__ = '2021-02-05'
__copyright__ = '(C) 2021 by Fredrik Lindberg'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication, QVariant, QDate, Qt
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterBoolean,
                       QgsProcessingParameterNumber,
                       QgsProcessingParameterRasterDestination,
                       QgsProcessingParameterString,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterField,
                       QgsProcessingParameterEnum,
                       QgsProcessingException,
                       QgsProcessingParameterDateTime,
                       QgsFeature,
                       QgsVectorFileWriter,
                       QgsVectorDataProvider,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterDefinition,
                       QgsVectorLayer,
                       QgsField)
from qgis.PyQt.QtGui import QIcon

from processing.gui.wrappers import WidgetWrapper

from qgis.PyQt.QtWidgets import QDateEdit
from osgeo import gdal, osr, ogr
from osgeo.gdalconst import GDT_Float32
import os
import numpy as np
import inspect
from pathlib import Path
import shutil
import datetime
from ..util.misc import saveraster
#from ..util.umep_uwg_export_component import read_uwg_file

try:
    from target_py import Target
    from target_py.ui.utils import read_config
except:
    pass


class ProcessingTARGETAnalyzerAlgorithm(QgsProcessingAlgorithm):
    """
    This class is a processing version of UWGAnalyzer but only for generating aggregated grids
    """
    INPUT_FOLDER = 'INPUT_FOLDER'
    OUTPUT_FOLDER = 'OUTPUT_FOLDER'
    INPUT_POLYGONLAYER = 'INPUT_POLYGONLAYER'
    ID_FIELD = 'ID_FIELD'
    SINGLE_DAY_BOOL = 'SINGLE_DAY_BOOL'
    SINGLE_DAY = 'SINGLE_DAY'
    IRREGULAR = 'IRREGULAR'
    PIXELSIZE = 'PIXELSIZE'
    STAT_TYPE = 'STAT_TYPE'
    ADD_ATTRIBUTES ='ADD_ATTRIBUTES'
    TIME_OF_DAY ='TIME_OF_DAY'
 
    # Output
    UWG_GRID_OUT = 'UWG_GRID_OUT'


    def initAlgorithm(self, config):

        self.plugin_dir = os.path.dirname(__file__)

        self.addParameter(QgsProcessingParameterFile(self.INPUT_FOLDER,
            self.tr('Path to TARGET Run name folder'),
            QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterBoolean(self.SINGLE_DAY_BOOL,
            self.tr("Examine single night"), defaultValue=False, optional=True))
        self.addParameter(QgsProcessingParameterDateTime(self.SINGLE_DAY,
            self.tr('Month and day when single night begins (year is irrelevant)'),
            QgsProcessingParameterDateTime.Date))
        self.statType = ((self.tr('Mean'), '0'),
                         (self.tr('Maximun'), '1'),
                         (self.tr('Median'), '2'),
                         (self.tr('75% percentile'), '3'),
                         (self.tr('95% percentile'), '4'))
        self.addParameter(QgsProcessingParameterEnum(self.STAT_TYPE,
                                                     self.tr('Statistic measure'),
                                                     options=[i[0] for i in self.statType],
                                                     defaultValue=0))
        self.dayType =  ((self.tr('Diurnal'), '0'),
                    (self.tr('Daytime'), '1'),
                    (self.tr('Nighttime'), '2'))
        self.addParameter(QgsProcessingParameterEnum(self.TIME_OF_DAY,
                                                     self.tr('Time of day'),
                                                     options=[i[0] for i in self.dayType],
                                                     defaultValue=0))
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT_POLYGONLAYER,
                                                              self.tr('Vector polygon grid'), 
                                                              [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterField(self.ID_FIELD,
                                                      self.tr('ID field'),
                                                      '', 
                                                      self.INPUT_POLYGONLAYER, 
                                                      QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterBoolean(self.IRREGULAR,
                                                        self.tr("Polygon grid irregular (not squared)"), 
                                                        defaultValue=False))
        self.addParameter(QgsProcessingParameterNumber(self.PIXELSIZE,
                                                       self.tr('Pixelsize if irregular grid is used (meter)'),
                                                       QgsProcessingParameterNumber.Integer,
                                                       QVariant(10), False, minValue=1))

        # Output
        self.addParameter(QgsProcessingParameterBoolean(self.ADD_ATTRIBUTES,
                                                        self.tr("Add results to vector polygon grid attribute table"), 
                                                        defaultValue=False))
        self.addParameter(QgsProcessingParameterRasterDestination(self.UWG_GRID_OUT,
                                                                  self.tr("Output raster from statistical analysis"),
                                                                  None,
                                                                  optional=True,
                                                                  createByDefault=False))

    def processAlgorithm(self, parameters, context, feedback):
        
        # InputParameters
        targetIn = self.parameterAsString(parameters, self.INPUT_FOLDER, context)
        #uwgOut = self.parameterAsString(parameters, self.OUTPUT_FOLDER, context)
        # variaIn = self.parameterAsString(parameters, self.VARIA_IN, context)
        # startday = self.parameterAsString(parameters, self.DATEINISTART, context)
        singleNight = self.parameterAsBool(parameters, self.SINGLE_DAY_BOOL, context)
        # endday = self.parameterAsString(parameters, self.DATEINIEND, context)
        inputPolygonlayer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGONLAYER, context)
        idField = self.parameterAsFields(parameters, self.ID_FIELD, context)
        irreg = self.parameterAsBool(parameters, self.IRREGULAR, context)
        statTypeStr = self.parameterAsString(parameters, self.STAT_TYPE, context)
        # dayTypeStr = self.parameterAsString(parameters, self.TIME_OF_DAY, context)
        pixelsize = self.parameterAsDouble(parameters, self.PIXELSIZE, context)
        addAttributes = self.parameterAsBool(parameters, self.ADD_ATTRIBUTES, context)
        outputStat = self.parameterAsOutputLayer(parameters, self.UWG_GRID_OUT, context)
        dayTypeStr = self.parameterAsString(parameters, self.TIME_OF_DAY, context)

        feedback.setProgressText("Initializing...")

        # statType = int(statTypeStr)

        feedback.setProgressText("Model input directory: " + targetIn + "/input")
        feedback.setProgressText("Model output directory: " + targetIn + "/output")

        cfM = read_config(targetIn + '/config.ini')

        runName = cfM['run_name']
        feedback.setProgressText("Run name: " + runName)


        #fileList = os.listdir(uwgIn)
        # a = fileList[0].find("_")
        # prefix = fileList[0][0:a]

        # feedback.setProgressText("Prefix: " + prefix)

        # uwgDict = read_uwg_file(uwgIn, fileList[0][:-4])
        # mm = uwgDict['Month']
        # dd = uwgDict['Day']
        # nDays = uwgDict['nDay']

        # Load rural data
        #inputDir + "/output/" + runName + 'metdata_UMEP.txt'
        sitein = targetIn + "/output/" + runName + '_metdata_UMEP.txt'
        dataref = np.genfromtxt(sitein, skip_header=1)
        yyyy = dataref[0,0]

        # start = datetime.date(int(yyyy), int(mm), int(dd))
        # end = start + datetime.timedelta(days=int(nDays))
        start = datetime.datetime.strptime(cfM['date1'], "%Y,%m,%d,%H").strftime("%Y-%m-%d") #string
        end = datetime.datetime.strptime(cfM['date2'], "%Y,%m,%d,%H").strftime("%Y-%m-%d") #string
        startDate = datetime.datetime.strptime(cfM['date1'], "%Y,%m,%d,%H")
        endDate = datetime.datetime.strptime(cfM['date2'], "%Y,%m,%d,%H")

        feedback.setProgressText('Days indentified as modelled by TARGET: ' + start + ' to ' + end)

        if singleNight:
            singledate = self.parameterAsString(parameters, self.SINGLE_DAY, context)
            startDate = datetime.datetime.strptime(singledate, '%Y-%m-%d')
            mm = startDate.month
            dd = startDate.day
            startDate = datetime.date(int(yyyy), int(mm), int(dd))
            nDays = 1
            endDate = startDate + datetime.timedelta(days=int(nDays))
            if startDate >= end or startDate < start:
                raise QgsProcessingException('Single date has to be within the modelled days OR the same as the last modelled day.')
            else:
                feedback.setProgressText('Single day and following night to examine: ' + startDate.strftime('%d %b'))
        else:
            #uwgDict = read_uwg_file(uwgIn, fileList[0][:-4])
            #mm = uwgDict['Month']
            #dd = uwgDict['Day']
            #nDays = uwgDict['nDay']
            #startDate = datetime.date(int(cfM['date1'][0]), int(cfM['date1'][1]), int(cfM['date1'][2]))
            #endDate = datetime.date(int(cfM['date2'][0]), int(cfM['date2'][1]), int(cfM['date2'][2]))
            #endDate = startDate + datetime.timedelta(days=int(nDays))
            feedback.setProgressText('Dates to be analyzed: ' + startDate.strftime('%d %b') + ' to ' + endDate.strftime('%d %b'))

        # poly = inputPolygonlayer
        poly_field = idField
        vlayer = inputPolygonlayer
        # prov = vlayer.dataProvider()

        path=vlayer.dataProvider().dataSourceUri()
        # polygonpath = path [:path.rfind('|')] # work around. Probably other solution exists
        if path.rfind('|') > 0:
            polygonpath = path [:path.rfind('|')] # work around. Probably other solution exists
        else:
            polygonpath = path

        # fields = prov.fields()
        idx = vlayer.fields().indexFromName(poly_field[0])

        # load, cut data and calculate statistics
        statvectemp = [0]
        statresult = [0]
        idvec = [0]

        startD = int(startDate.strftime('%j'))
        endD = int(endDate.strftime('%j'))

       
        # for i in range(0, self.idgrid.shape[0]): # loop over vector grid instead
        index = 1
        nGrids = vlayer.featureCount()
        for f in vlayer.getFeatures():
            feedback.setProgress(int((index * 100) / nGrids))
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break
            index += 1

            gid = str(int(f.attributes()[idx]))

            feedback.setProgressText("Processing grid: " + str(gid))

            datawhole = np.genfromtxt(targetIn + '/output/' + runName + '_' + gid + '_UMEP_TARGET.txt', skip_header=1)

            # cut TARGET data
            start = np.min(np.where(datawhole[:, 1] == startD))
            if endD > np.max(datawhole[:, 1]):
                ending = np.max(np.where(datawhole[:, 1] == endD - 1))
            else:
                ending = np.min(np.where(datawhole[:, 1] == endD))



            data1 = datawhole[start:int(ending + 12), :] # + 12 to include whole final night 


            # Select depending of time of day for modelled data
            if dayTypeStr == '1':
                data1 = data1[np.where(data1[:, 14] > 1.), :]
                data1 = data1[0][:]
            if dayTypeStr == '2':
                data1 = data1[np.where(data1[:, 14] < 1.), :]
                data1 = data1[0][:]
            
            #data1 = data1[np.where(data1[:, 14] < 1.), :] # include only nighttime. 14 is position for global radiation
            #data1 = data1[0][:]

            # cut ref data
            if endD > np.max(dataref[:, 1]):
                ending = np.max(np.where(dataref[:, 1] == endD - 1))
            else:
                ending = np.min(np.where(dataref[:, 1] == endD))
            data2 = dataref[start:int(ending + 12), :] # + 12 to include whole final night 

            # Select depending of time of day for ref data
            if dayTypeStr == '1':
                data2 = data2[np.where(data2[:, 14] > 1.), :]
                data2 = data2[0][:]
            if dayTypeStr == '2':
                data2 = data2[np.where(data2[:, 14] < 1.), :]
                data2 = data2[0][:]

            # data2 = data2[np.where(data2[:, 14] < 1.), :] # include only nighttime. 14 is position for global radiation
            # data2 = data2[0][:]

            vardatauwg = data1[:, 11] # 11 is temperature column
            vardataref = data2[:, 11] 
            vardata = vardatauwg - vardataref

            if statTypeStr == '0':
                statresult = np.nanmean(vardata)
                header = 'mean'
            if statTypeStr == '1':
                statresult = np.nanmax(vardata)
                header = 'max'
            if statTypeStr == '2':
                statresult = np.nanpercentile(vardata, 50)
                header = 'median'
            if statTypeStr == '3':
                statresult = np.nanpercentile(vardata, 75)
                header = '75precentile'
            if statTypeStr == '4':
                statresult = np.nanpercentile(vardata, 95)
                header = '95precentile'

            statvectemp = np.vstack((statvectemp, statresult))
            idvec = np.vstack((idvec, int(gid)))

        statvector = statvectemp[1:, :]
        # fix_print_with_import
        statmat = np.hstack((idvec[1:, :], statvector))
        statmat[statmat < -500] = -9999 #Response to #107

        if addAttributes:
            self.addattributes(vlayer, statmat, header)

        if irreg:
            resx = pixelsize
        else:
            for f in vlayer.getFeatures():  # Taking first polygon. Could probably be done nicer
                # geom = f.geometry().asPolygon()
                geom = f.geometry().asMultiPolygon()
                break
            resx = np.abs(geom[0][0][0][0] - geom[0][0][2][0])  # x
            resy = np.abs(geom[0][0][0][1] - geom[0][0][2][1])  # y

            if not resx == resy:
                raise QgsProcessingException("Polygons not squared in current CRS")

        if os.path.isfile(self.plugin_dir + '/tempgrid.tif'): # response to issue 103
            try:
                shutil.rmtree(self.plugin_dir + '/tempgrid.tif')
            except OSError:
                os.remove(self.plugin_dir + '/tempgrid.tif')
        
        extent = vlayer.extent()
        crs = vlayer.crs().toWkt()
        self.rasterize(polygonpath, str(self.plugin_dir + '/tempgrid.tif'), str(poly_field[0]), resx, crs, extent)

        dataset = gdal.Open(self.plugin_dir + '/tempgrid.tif')
        idgrid_array = dataset.ReadAsArray().astype(float)

        gridout = np.zeros((idgrid_array.shape[0], idgrid_array.shape[1]))

        for i in range(0, statmat.shape[0]):
            gridout[idgrid_array == statmat[i, 0]] = statmat[i, 1]

        if outputStat:
            saveraster(dataset, outputStat, gridout)

        feedback.setProgressText("Processing finished.")

        return {self.UWG_GRID_OUT: outputStat}

    def rasterize(self, src, dst, attribute, resolution, crs, extent, all_touch=False, na=-9999):

        # Open shapefile, retrieve the layer
        src_data = ogr.Open(src)
        layer = src_data.GetLayer()

        # Use transform to derive coordinates and dimensions
        xmax = extent.xMaximum()
        xmin = extent.xMinimum()
        ymax = extent.yMaximum()
        ymin = extent.yMinimum()

        # Create the target raster layer
        cols = int((xmax - xmin)/resolution)
        # rows = int((ymax - ymin)/resolution) + 1
        rows = int((ymax - ymin)/resolution)  # issue 164
        trgt = gdal.GetDriverByName("GTiff").Create(dst, cols, rows, 1, GDT_Float32)
        trgt.SetGeoTransform((xmin, resolution, 0, ymax, 0, -resolution))

        # Add crs
        refs = osr.SpatialReference()
        refs.ImportFromWkt(crs)
        trgt.SetProjection(refs.ExportToWkt())

        # Set no value
        band = trgt.GetRasterBand(1)
        band.SetNoDataValue(na)

        # Set options
        if all_touch is True:
            ops = ["-at", "ATTRIBUTE=" + attribute]
        else:
            ops = ["ATTRIBUTE=" + attribute]

        # Finally rasterize
        gdal.RasterizeLayer(trgt, [1], layer, options=ops)

        # Close target an source rasters
        del trgt
        del src_data 

    
    def addattributes(self, vlayer, matdata, header):
        current_index_length = len(vlayer.dataProvider().attributeIndexes())
        caps = vlayer.dataProvider().capabilities()

        if caps & QgsVectorDataProvider.AddAttributes:
            vlayer.dataProvider().addAttributes([QgsField(header, QVariant.Double)])
            attr_dict = {}
            for y in range(0, matdata.shape[0]):
                attr_dict.clear()
                idx = int(matdata[y, 0])
                attr_dict[current_index_length] = float(matdata[y, 1])
                vlayer.dataProvider().changeAttributeValues({y: attr_dict})

            vlayer.updateFields()
        else:
            raise QgsProcessingException("Vector Layer does not support adding attributes")
    
    def name(self):
        return 'Urban Heat Island: TARGET Analyzer'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Post-Processor'

    def shortHelpString(self):
        return self.tr('The <b>TARGET Analyzer</b> plugin can be used to make basic grid analysis of model results generated by the TARGET model.<br>'
        '\n'
        '--------------\n'
        'Full manual available via the <b>Help</b>-button.')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/post_processor/Urban%20Heat%20Island%20TARGET%20Analyser.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Post-Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/icon_uwg.png")
        return icon

    def createInstance(self):
        return ProcessingTARGETAnalyzerAlgorithm()

class DateWidgetStart(WidgetWrapper):
    def createWidget(self):
        self._combo = QDateEdit()
        self._combo.setCalendarPopup(True)

        today = QDate(2010, 1, 1)
        self._combo.setDate(today)

        return self._combo

    def value(self):
        date_chosen = self._combo.dateTime()
        return date_chosen.toString(Qt.ISODate)

class DateWidgetEnd(WidgetWrapper):
    def createWidget(self):
        self._combo = QDateEdit()
        self._combo.setCalendarPopup(True)

        today = QDate(2010, 12, 31)
        self._combo.setDate(today)

        return self._combo

    def value(self):
        date_chosen = self._combo.dateTime()
        return date_chosen.toString(Qt.ISODate)
