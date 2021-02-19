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
import sys
from ..util import f90nml
import shutil
from ..util.misc import saverasternd, saveraster

# def saverasternd(gdal_data, filename, raster):
#     rows = gdal_data.RasterYSize
#     cols = gdal_data.RasterXSize

#     outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
#     outBand = outDs.GetRasterBand(1)

#     # write the data
#     outBand.WriteArray(raster, 0, 0)
#     # flush data to disk, set the NoData value and calculate stats
#     outBand.FlushCache()
#     # outBand.SetNoDataValue(-9999)

#     # georeference the image and set the projection
#     outDs.SetGeoTransform(gdal_data.GetGeoTransform())
#     outDs.SetProjection(gdal_data.GetProjection())

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
        writeoutoption = 2 # hardcoded to simple output
        dataunit = self.plugin_dir + '/SUEWS_OutputFormatOption' + str(int(writeoutoption)) + '.txt'  # File moved to plugin directory
        f = open(dataunit)
        lin = f.readlines()
        self.lineunit = lin[3].split(";")
        self.linevar = lin[1].split(";")
        self.linevarlong = lin[2].split(";")
        f.close()

        listA = [tuple([a, str(i)]) for i, a in enumerate(self.linevarlong)]
        self.varType = tuple(listA)

        self.addParameter(QgsProcessingParameterFile(self.SUEWS_NL,
                                                     self.tr('SUEWS RunControl namelist'),
                                                     extension='nml',
                                                     optional=False))
        self.addParameter(QgsProcessingParameterEnum(self.VARIA_IN,
                                                     self.tr('Variable to post-process'),
                                                     options=[i[0] for i in self.varType],
                                                     defaultValue=13))
        self.dayType =  ((self.tr('Diurnal'), '0'),
                         (self.tr('Daytime'), '1'),
                         (self.tr('Nighttime'), '2'))
        self.addParameter(QgsProcessingParameterEnum(self.TIME_OF_DAY,
                                                     self.tr('Time of day'),
                                                     options=[i[0] for i in self.dayType],
                                                     defaultValue=0))
        self.statType = ((self.tr('Mean'), '0'),
                         (self.tr('Minimun'), '1'),
                         (self.tr('Maximun'), '2'),
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

        statType = int(statTypeStr)

        # read nml
        # self.fileDialognml.open()
        # result = self.fileDialognml.exec_()
        # if result == 1:
        # self.nmlPath = self.fileDialognml.selectedFiles()
        # self.dlg.textModelFolder.setText(self.nmlPath[0])
        nml = f90nml.read(suewsNL)

        self.fileinputpath = nml['runcontrol']['fileinputpath']
        if self.fileinputpath.startswith("."):
            nmlfolder = self.nmlPath[0][:-15]
            self.fileinputpath = nmlfolder + self.fileinputpath[1:]

        self.fileoutputpath = nml['runcontrol']['fileoutputpath']
        if self.fileoutputpath.startswith("."):
            nmlfolder = self.nmlPath[0][:-15]
            self.fileoutputpath = nmlfolder + self.fileoutputpath[1:]

        resolutionFilesOut = nml['runcontrol']['resolutionfilesout']
        self.resout = int(float(resolutionFilesOut) / 60)
        self.fileCode = nml['runcontrol']['filecode']
        self.multiplemetfiles = nml['runcontrol']['multiplemetfiles']
        resolutionFilesIn = nml['runcontrol']['resolutionFilesIn']
        self.resin = int(resolutionFilesIn / 60)

        tstep = nml['runcontrol']['tstep']
        self.tstep = int(float(tstep) / 60)
        writeoutoption = nml['runcontrol']['writeoutoption']
        if writeoutoption != 2:
            raise QgsProcessingException('WriteOutOption=2 must have been specified in RunControl.nml to use the current version of SUEWS Analyzer.')

        feedback.setProgressText("Model input directory: " + self.fileinputpath)
        feedback.setProgressText("Model output directory: " + self.fileoutputpath)

        # mm = 0 #This doesn't work when hourly file starts with e.g.15
        # while mm < 60:
        #     self.dlg.comboBox_mm.addItem(str(mm))
        #     mm += self.resout

        sitein = self.fileinputpath + 'SUEWS_SiteSelect.txt'
        f = open(sitein)
        lin = f.readlines()
        self.YYYY = -99
        self.gridcodemetID = -99
        index = 2
        loop_out = ''
        gridcodemetmat = [0]
        yeartest = []
        while loop_out != '-9':
            lines = lin[index].split()
            if not int(lines[1]) == self.YYYY:
                self.YYYY = int(lines[1])
                feedback.setProgressText("Year identified for possible analysis: " + str(self.YYYY))
                yeartest.append(self.YYYY)
                # self.dlg.comboBox_POIYYYY.addItem(str(self.YYYY))
                # self.dlg.comboBox_SpatialYYYY.addItem(str(self.YYYY))

            # if not np.any(int(lines[0]) == gridcodemetmat):
            #     self.gridcodemetID = int(lines[0])
                # self.dlg.comboBox_POIField.addItem(str(self.gridcodemetID))
                # self.dlg.comboBox_POIField_2.addItem(str(self.gridcodemetID))
                # gridtemp = [self.gridcodemetID]
                # gridcodemetmat = np.vstack((gridcodemetmat, gridtemp))

            if index == 2:
                if self.multiplemetfiles == 0:
                    self.gridcodemet = ''
                else:
                    self.gridcodemet = lines[0]
                data_in = self.fileinputpath + self.fileCode + self.gridcodemet + '_' + str(self.YYYY) + '_data_' + str(self.resin) + '.txt'
                self.met_data = np.genfromtxt(data_in, skip_header=1, missing_values='**********', filling_values=-9999)  # , skip_footer=2

            lines = lin[index + 1].split()
            loop_out = lines[0]
            index += 1

        f.close()
        # self.idgrid = gridcodemetmat[1:, :]

        # dataunit = self.plugin_dir + '/SUEWS_OutputFormatOption' + str(int(writeoutoption)) + '.txt'  # File moved to plugin directory
        # f = open(dataunit)
        # lin = f.readlines()
        # self.lineunit = lin[3].split(";")
        # self.linevar = lin[1].split(";")
        # self.linevarlong = lin[2].split(";")
        # f.close()

        # for i in range(0, self.linevarlong.__len__()):
        #     self.dlg.comboBox_POIVariable.addItem(self.linevarlong[i])
        #     self.dlg.comboBox_POIVariable_2.addItem(self.linevarlong[i])
        #     self.dlg.comboBox_SpatialVariable.addItem(self.linevarlong[i])

        # self.dlg.runButtonPlot.setEnabled(1)
        # self.dlg.runButtonSpatial.setEnabled(1)
        
        
        # if self.dlg.comboBox_SpatialVariable.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No analyzing variable is selected")
        #     return
        # else:
        self.id = int(variaIn) #self.dlg.comboBox_SpatialVariable.currentIndex() - 1

        # if self.dlg.comboBox_SpatialYYYY.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No Year is selected")
        #     return

        # if self.dlg.comboBox_SpatialDOYMin.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No Minimum DOY is selected")
        #     return

        # if self.dlg.comboBox_SpatialDOYMax.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No Maximum DOY is selected")
        #     return

        poly = inputPolygonlayer
        # if poly is None:
        #     QMessageBox.critical(self.dlg, "Error", "No valid Polygon layer is selected")
        #     return
        # if not poly.geometryType() == 2:
        #     QMessageBox.critical(self.dlg, "Error", "No valid Polygon layer is selected")
        #     return

        poly_field = idField
        # if poly_field is None:
        #     QMessageBox.critical(self.dlg, "Error", "An attribute with unique fields/records must be selected (same as used in the model run to analyze)")
        #     return

        # if not (self.dlg.addResultToGrid.isChecked() or self.dlg.addResultToGeotiff.isChecked()):
        #     QMessageBox.critical(self.dlg, "Error", "No output method has been selected (Add results to polygon grid OR Save as GeoTIFF)")
        #     return

        # if self.dlg.comboBox_SpatialDOYMin.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No Minimum DOY is selected")
        #     return
        # else:
        # startday = int(self.dlg.comboBox_SpatialDOYMin.currentText())

        # if self.dlg.comboBox_SpatialDOYMax.currentText() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No Maximum DOY is selected")
        #     return
        # else:
        # endday = int(self.dlg.comboBox_SpatialDOYMax.currentText())

        if startday >= endday:
            raise QgsProcessingException('Start date is greater or equal than end date')

        # if startday > endday:
        #     QMessageBox.critical(self.dlg, "Error", "Start day happens after end day")
        #     return

        # if startday == endday:
        #     QMessageBox.critical(self.dlg, "Error", "End day must be higher than start day")
        #     return

        # if self.dlg.checkBox_TOD.isChecked():
        #     if self.dlg.comboBox_HH.currentText() == ' ':
        #         QMessageBox.critical(self.dlg, "Error", "No Hour specified")
        #         return
        #     if self.dlg.comboBox_mm.currentText() == ' ':
        #         QMessageBox.critical(self.dlg, "Error", "No Minute specified")
        #         return

        # load, cut data and calculate statistics
        statvectemp = [0]
        statresult = [0]
        idvec = [0]
        vlayer = inputPolygonlayer #QgsVectorLayer(poly.source(), "polygon", "ogr")
        prov = vlayer.dataProvider()
        fields = prov.fields()
        # idx = vlayer.fieldNameIndex(poly_field)
        idx = vlayer.fields().indexFromName(poly_field[0])
        typetest = fields.at(idx).type()
        # if typetest == 10:
        #     raise QgsProcessingException("ID field must be either integer or float")
        starty = int(startday[0:4]) #date.year()
        startm = int(startday[5:7]) #date.month()
        startd = int(startday[8:10]) #date.day()
        if (starty % 4) == 0:
            if (starty % 100) == 0:
                if (starty % 400) == 0:
                    leapyear = 1
                else:
                    leapyear = 0
            else:
                leapyear = 1
        else:
            leapyear = 0
        if leapyear == 1:
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        startD = sum(dayspermonth[0:int(startm - 1)]) + startd
        endy = int(endday[0:4]) #date.year()
        endm = int(endday[5:7]) #date.month()
        endd = int(endday[8:10]) #date.day()
        if (endy % 4) == 0:
            if (endy % 100) == 0:
                if (endy % 400) == 0:
                    leapyear = 1
                else:
                    leapyear = 0
            else:
                leapyear = 1
        else:
            leapyear = 0
        if leapyear == 1:
            dayspermonth = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        endD = sum(dayspermonth[0:int(endm - 1)]) + endd

        if not starty == endy:
            raise QgsProcessingException('Startdate and enddate must happen within the same year. Multiple years will be possible in future versions')

        if not starty in yeartest:
            raise QgsProcessingException('Selected timeperiod not present in output data. Choose a period within the year(s): ' + str(yeartest[:]))

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
            datawhole = np.genfromtxt(self.fileoutputpath + '/' + self.fileCode + gid + '_'
                                     + str(self.YYYY) + '_SUEWS_' + str(self.resout) + '.txt', skip_header=1,
                                     missing_values='**********', filling_values=-9999)

            writeoption = datawhole.shape[1]
            if writeoption > 38:
                altpos = 52
            else:
                altpos = 25
            feedback.setProgressText("Processing grid: " + str(gid))
            
            start = np.min(np.where(datawhole[:, 1] == startD))
            if endD > np.max(datawhole[:, 1]):
                ending = np.max(np.where(datawhole[:, 1] == endD - 1))
            else:
                ending = np.min(np.where(datawhole[:, 1] == endD))
            data1 = datawhole[start:ending, :]

            # if self.dlg.checkBox_TOD.isChecked():
            #     hh = self.dlg.comboBox_HH.currentText()
            #     hhdata = np.where(data1[:, 2] == int(hh))
            #     data1 = data1[hhdata, :]
            #     minute = self.dlg.comboBox_mm.currentText()
            #     mmdata = np.where(data1[0][:, 3] == int(minute))
            #     data1 = data1[0][mmdata, :]
            #     data1 = data1[0][:]
            # else:
            
            if dayTypeStr == '1':
                data1 = data1[np.where(data1[:, altpos] < 90.), :]
                data1 = data1[0][:]
            if dayTypeStr == '2':
                data1 = data1[np.where(data1[:, altpos] > 90.), :]
                data1 = data1[0][:]

            vardata = data1[:, self.id]

            if statTypeStr == '0':
                statresult = np.nanmean(vardata)
            if statTypeStr == '1':
                statresult = np.nanmin(vardata)
            if statTypeStr == '2':
                statresult = np.nanmax(vardata)
            if statTypeStr == '3':
                statresult = np.nanmedian(vardata)
            if statTypeStr == '4':
                statresult = np.nanpercentile(vardata, 75) - np.percentile(vardata, 25)

            statvectemp = np.vstack((statvectemp, statresult))
            idvec = np.vstack((idvec, int(gid)))

        statvector = statvectemp[1:, :]
        # fix_print_with_import
        statmat = np.hstack((idvec[1:, :], statvector))
        
        # numformat2 = '%8d %5.3f'
        # header2 = 'id value'
        # np.savetxt(self.plugin_dir + 'test.txt', statmat,
        #            fmt=numformat2, delimiter=' ', header=header2, comments='')

        header = self.linevar[self.id]
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

        # polyname = self.dlg.comboBox_Polygrid.currentText()
        # polyname = self.layerComboManagerPolygrid.currentText()
        # if self.dlg.textOutput.text() == 'Not Specified':
        #     QMessageBox.critical(self.dlg, "Error", "No output filename for GeoTIFF is added")
        #     return
        # else:
        #     filename = self.dlg.textOutput.text()

        if os.path.isfile(self.plugin_dir + '/tempgrid.tif'): # response to issue 103
            try:
                shutil.rmtree(self.plugin_dir + '/tempgrid.tif')
            except OSError:
                os.remove(self.plugin_dir + '/tempgrid.tif')

        crs = vlayer.crs().toWkt()
        self.rasterize(str(poly.source()), str(self.plugin_dir + '/tempgrid.tif'), str(poly_field[0]), resx, crs, extent)

        dataset = gdal.Open(self.plugin_dir + '/tempgrid.tif')
        idgrid_array = dataset.ReadAsArray().astype(np.float)

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

    # def saveraster(self, gdal_data, filename, raster):
    #     rows = gdal_data.RasterYSize
    #     cols = gdal_data.RasterXSize

    #     outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
    #     outBand = outDs.GetRasterBand(1)

    #     # write the data
    #     outBand.WriteArray(raster, 0, 0)
    #     # flush data to disk, set the NoData value and calculate stats
    #     outBand.FlushCache()
    #     outBand.SetNoDataValue(-9999)

    #     # georeference the image and set the projection
    #     outDs.SetGeoTransform(gdal_data.GetGeoTransform())
    #     outDs.SetProjection(gdal_data.GetProjection())

    #     del outDs, outBand  
    
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
            QMessageBox.critical(None, "Error", "Vector Layer does not support adding attributes")
    
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
