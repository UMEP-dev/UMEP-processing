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
from osgeo.gdalconst import *
import os
import numpy as np
import inspect
from pathlib import Path
from ..util import f90nml
import shutil
from ..util.misc import saveraster, get_resolution_from_file, SUEWS_txt_to_df, extract_suews_years
from .params_dict import *
import yaml
import datetime

class ProcessingSuewsAnalyzerAlgorithm(QgsProcessingAlgorithm):
    """
    This class is a processing version of SuewsAnalyzer but only for generating aggregated grids
    """
    SUEWS_NL = 'SUEWS_NL'
    VARIA_IN = 'VARIA_IN'
    INPUT_POLYGONLAYER = 'INPUT_POLYGONLAYER'
    ID_FIELD = 'ID_FIELD'
    YEAR = 'YEAR'
    DATEINISTART = 'DATEINISTART'
    DATEINIEND = 'DATEINIEND'
    IRREGULAR = 'IRREGULAR'
    PIXELSIZE = 'PIXELSIZE'
    TIME_OF_DAY = 'TIME_OF_DAY'
    STAT_TYPE = 'STAT_TYPE'
    ADD_ATTRIBUTES ='ADD_ATTRIBUTES'
 
    # Output
    SUEWS_GRID_OUT = 'SUEWS_GRID_OUT'


    def initAlgorithm(self, config):

        self.plugin_dir = os.path.dirname(__file__)

    
        var_list = []
        idx = 0
        for item in list(sorted(list(params_dict.keys()))):
            item_desc = item + f" ({params_dict[item]['description']})"
            item = tuple([item_desc, str(idx)])
            var_list.append(item)
            idx = idx+1

        self.varType = tuple(var_list)

        self.addParameter(QgsProcessingParameterFile(self.SUEWS_NL,
                                                     self.tr('Yaml-file'),
                                                     extension='yml',
                                                     optional=False))
        self.addParameter(QgsProcessingParameterEnum(self.VARIA_IN,
                                                     self.tr('Variable to post-process'),
                                                     options=[i[0] for i in self.varType],
                                                     defaultValue=26))
        self.dayType =  ((self.tr('Diurnal'), '0'),
                         (self.tr('Daytime'), '1'),
                         (self.tr('Nighttime'), '2'))
        self.addParameter(QgsProcessingParameterEnum(self.TIME_OF_DAY,
                                                     self.tr('Time of day'),
                                                     options=[i[0] for i in self.dayType],
                                                     defaultValue=0))
        self.statType = ((self.tr('Mean'), '0'),
                         (self.tr('Minimum'), '1'),
                         (self.tr('Maximum'), '2'),
                         (self.tr('Median'), '3'),
                         (self.tr('IQR'), '4'))
        self.addParameter(QgsProcessingParameterEnum(self.STAT_TYPE,
                                                     self.tr('Statistic measure'),
                                                     options=[i[0] for i in self.statType],
                                                     defaultValue=0))
        paramS = QgsProcessingParameterString(self.DATEINISTART, 'Start date')
        paramS.setMetadata({'widget_wrapper': {'class': DateWidgetStart}})
        self.addParameter(paramS)
        paramE = QgsProcessingParameterString(self.DATEINIEND, 'End date')
        paramE.setMetadata({'widget_wrapper': {'class': DateWidgetEnd}})
        self.addParameter(paramE)
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
                                                        self.tr("Add results to Vector polygon grid attribute table"), 
                                                        defaultValue=False))
        self.addParameter(QgsProcessingParameterRasterDestination(self.SUEWS_GRID_OUT,
                                                                  self.tr("Output raster from statistical analysis"),
                                                                  None,
                                                                  optional=False))

    def processAlgorithm(self, parameters, context, feedback):
        
        # InputParameters
        suewsNL = self.parameterAsString(parameters, self.SUEWS_NL, context)
        variaIn = self.parameterAsString(parameters, self.VARIA_IN, context)
        startday = self.parameterAsString(parameters, self.DATEINISTART, context)
        endday = self.parameterAsString(parameters, self.DATEINIEND, context)
        inputPolygonlayer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGONLAYER, context)
        idField = self.parameterAsFields(parameters, self.ID_FIELD, context)
        irreg = self.parameterAsBool(parameters, self.IRREGULAR, context)
        statTypeStr = self.parameterAsString(parameters, self.STAT_TYPE, context)
        dayTypeStr = self.parameterAsString(parameters, self.TIME_OF_DAY, context)
        pixelsize = self.parameterAsDouble(parameters, self.PIXELSIZE, context)
        addAttributes = self.parameterAsBool(parameters, self.ADD_ATTRIBUTES, context)
        outputStat = self.parameterAsOutputLayer(parameters, self.SUEWS_GRID_OUT, context)

        feedback.setProgressText("Initializing...")

        with open(suewsNL, 'r') as f:
                yaml_dict = yaml.load(f, Loader=yaml.SafeLoader)

        self.fileoutputpath = yaml_dict['model']['control']['output_file']
        
        if self.fileoutputpath.startswith("."):
            yamlfolder = self.yamlPath[0][:-15]
            self.fileoutputpath = yamlfolder + self.fileoutputpath[1:]

        resolutionFilesOut = get_resolution_from_file(yaml_dict['model']['control']['output_file'])
        self.resout = int(float(resolutionFilesOut) / 60)

        years = extract_suews_years(self.fileoutputpath)
        self.YYYY = str(startday).split('-')[0]

        if self.YYYY not in years:
            raise QgsProcessingException(f"Selected year '{self.YYYY}' is not present in the data.\n Availible years are {str(years)}")


        self.id = int(variaIn) #self.dlg.comboBox_SpatialVariable.currentIndex() - 1

        poly_field = idField


        if startday >= endday:
            raise QgsProcessingException('Start date is greater or equal than end date')

 
        # load, cut data and calculate statistics
        statvectemp = [0]
        statresult = [0]
        idvec = [0]
        vlayer = inputPolygonlayer #QgsVectorLayer(poly.source(), "polygon", "ogr")
        prov = vlayer.dataProvider()
        fields = prov.fields()
        path=vlayer.dataProvider().dataSourceUri()
        if path.rfind('|') > 0:
            polygonpath = path [:path.rfind('|')] # work around. Probably other solution exists
        else:
            polygonpath = path

        grid_list = [feature[poly_field[0]] for feature in vlayer.getFeatures()]
        
        for grid in grid_list:

            datawhole = SUEWS_txt_to_df(self.fileoutputpath + '/' + str(grid) + '_' +
                                      str(self.YYYY) + '_SUEWS_' + str(self.resout) + '.txt')
            
            data1 = datawhole.loc[startday:endday]
            
            if dayTypeStr == '0': # Diurnal
                pass
            elif dayTypeStr == '1': # Daytime:
                    data1 = data1[data1.loc[:, 'Zenith'] < 90]

            elif dayTypeStr == '2': # Nighttime:
                    data1 = data1[data1.loc[:, 'Zenith'] > 90.]
            var = self.varType[self.id][0].split(' ')[0] # get position and slice variable
            vardata = data1.loc[:, var]

            if statTypeStr == '0':
                statresult = np.nanmean(vardata)
                suffix = '_mean'
            if statTypeStr == '1':
                statresult = np.nanmin(vardata)
                suffix = '_min'
            if statTypeStr == '2':
                statresult = np.nanmax(vardata)
                suffix = '_max'
            if statTypeStr == '3':
                statresult = np.nanmedian(vardata)
                suffix = '_median'
            if statTypeStr == '4':
                statresult = np.nanpercentile(vardata, 75) - np.percentile(vardata, 25)
                suffix = '_IQR'
            statvectemp = np.vstack((statvectemp, statresult))
            idvec = np.vstack((idvec, int(grid)))

        statvector = statvectemp[1:, :]
        statmat = np.hstack((idvec[1:, :], statvector))
 
        header = var + suffix
        if addAttributes:
            self.addattributes(vlayer, statmat, header)

        # if self.dlg.addResultToGeotiff.isChecked():
        extent = vlayer.extent()
        xmax = extent.xMaximum()
        xmin = extent.xMinimum()
        ymax = extent.yMaximum()
        ymin = extent.yMinimum()

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
                return

       
        if os.path.isfile(self.plugin_dir + '/tempgrid.tif'): # response to issue 103
            try:
                shutil.rmtree(self.plugin_dir + '/tempgrid.tif')
            except OSError:
                os.remove(self.plugin_dir + '/tempgrid.tif')

        crs = vlayer.crs().toWkt()
        self.rasterize(polygonpath, str(self.plugin_dir + '/tempgrid.tif'), str(poly_field[0]), resx, crs, extent)

        dataset = gdal.Open(self.plugin_dir + '/tempgrid.tif')
        idgrid_array = dataset.ReadAsArray().astype(float)

        gridout = np.zeros((idgrid_array.shape[0], idgrid_array.shape[1]))

        for i in range(0, statmat.shape[0]):
            gridout[idgrid_array == statmat[i, 0]] = statmat[i, 1]

        saveraster(dataset, outputStat, gridout)

        feedback.setProgressText("Processing finished.")

        return {self.SUEWS_GRID_OUT: outputStat}

    def rasterize(self, src, dst, attribute, resolution, crs, extent, all_touch=False, na=-9999):

        # Open shapefile, retrieve the layer
        # print(src)
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
        trgt = gdal.GetDriverByName("GTiff").Create(dst, cols, rows, 1, gdal.GDT_Float32)
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
            raise QgsProcessingException("Error", "Vector Layer does not support adding attributes")
    
    def name(self):
        return 'Urban Energy Balance: SUEWS Analyzer'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Post-Processor'

    def shortHelpString(self):
        return self.tr('The <b>SUEWS Analyzer</b> plugin can be used to make basic grid analysis of model results generated by the SUEWS model.<br>'
        '\n'
        '--------------\n'
        'Full manual available via the <b>Help</b>-button.')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/post_processor/Urban%20Energy%20Balance%20SUEWS%20Analyser.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Post-Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/SuewsLogo.png")
        return icon

    def createInstance(self):
        return ProcessingSuewsAnalyzerAlgorithm()

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
