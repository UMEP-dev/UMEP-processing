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
                       QgsProcessingParameterDateTime,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterField,
                       QgsProcessingException,
                       QgsVectorLayer,
                       QgsFeature,
                       QgsVectorFileWriter,
                       QgsVectorDataProvider,
                       QgsField)

from qgis.PyQt.QtGui import QIcon
# from osgeo import gdal, osr, ogr
from osgeo.gdalconst import *
import os
import sys
import numpy as np
import inspect
from pathlib import Path
import shutil
import traceback
import math
from ..util.umep_uwg_export_component import get_uwg_file, read_uwg_file
try:
    from uwg import UWG
except:
    pass


class ProcessingUWGProcessorAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm make use of UWG for the processing toolbox
    """
    
    INPUT_FOLDER = 'INPUT_FOLDER'
    INPUT_POLYGONLAYER = 'INPUT_POLYGONLAYER'
    ID_FIELD = 'ID_FIELD'
    START_DATE = 'START_DATE'
    NDAYS = 'NDAYS'
    INPUT_MET = 'INPUT_MET'
    UMEP_OUTPUT = 'UMEP_OUTPUT'
    OUTPUT_DIR = 'OUTPUT_DIR'
    OUTPUT_FORMAT = 'OUTPUT_FORMAT'
    DTSIM = 'DTSIM'
    EXCLUDE_RURAL = 'EXCLUDE_RURAL'


    def initAlgorithm(self, config):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_FOLDER,
            self.tr('Path to folder where UWG input files are located'),
            QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT_POLYGONLAYER,
            self.tr('Vector data including location(s) to model'),
            [QgsProcessing.TypeVector]))
        self.addParameter(QgsProcessingParameterField(self.ID_FIELD,
            self.tr('ID field'),
            '',
            self.INPUT_POLYGONLAYER,
            QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterDateTime(self.START_DATE,
            self.tr('Start date of simulation'),
            QgsProcessingParameterDateTime.Date))
        self.addParameter(QgsProcessingParameterNumber(self.NDAYS,
            self.tr('Number of days to run simulation'),
            QgsProcessingParameterNumber.Integer,
            QVariant(5), False, minValue=1, maxValue=365))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_MET,
            self.tr('Input meteorological file (*.epw)'),
            extension = 'epw'))
        self.addParameter(QgsProcessingParameterNumber(self.DTSIM,
            self.tr('Simulation time step in seconds'),
            QgsProcessingParameterNumber.Integer,
            QVariant(300), False, minValue=1, maxValue=1440))
        self.addParameter(QgsProcessingParameterBoolean(self.EXCLUDE_RURAL,
            self.tr('Exculde grids with very small building fractions (< 0.5%)'), defaultValue=False))
        # output
        self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_DIR,
            self.tr('Output folder')))
        self.addParameter(QgsProcessingParameterBoolean(self.OUTPUT_FORMAT,
            self.tr('Save output in UMEP specific format. Leave ticked off to store in epw-format.')))


    def processAlgorithm(self, parameters, context, feedback):

        try:
            import uwg
        except:
            pass
            raise QgsProcessingException("uwg python library not found: Instructions on how to install missing python libraries using the pip command: https://umep-docs.readthedocs.io/en/latest/Getting_Started.html")

        # InputParameters
        inputDir = self.parameterAsString(parameters, self.INPUT_FOLDER, context)
        inputPolygonlayer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGONLAYER, context)
        idField = self.parameterAsFields(parameters, self.ID_FIELD, context)
        startDate = self.parameterAsString(parameters, self.START_DATE, context)
        nDays =  self.parameterAsDouble(parameters, self.NDAYS, context)
        inputMet = self.parameterAsString(parameters, self.INPUT_MET, context)
        outputDir = self.parameterAsString(parameters, self.OUTPUT_DIR, context)
        umepformat = self.parameterAsBoolean(parameters, self.OUTPUT_FORMAT, context)
        dtSim = self.parameterAsDouble(parameters, self.DTSIM, context)
        excludeRural = self.parameterAsBoolean(parameters, self.EXCLUDE_RURAL, context)
        
        if parameters['OUTPUT_DIR'] == 'TEMPORARY_OUTPUT':
            if not (os.path.isdir(outputDir)):
                os.mkdir(outputDir)

        fileList = os.listdir(inputDir)
        a = fileList[0].find("_")
        prefix = fileList[0][0:a]

        # poly = inputPolygonlayer
        poly_field = idField
        vlayer = inputPolygonlayer
        # prov = vlayer.dataProvider()
        # fields = prov.fields()
        idx = vlayer.fields().indexFromName(poly_field[0])
        nGrids = vlayer.featureCount()

        feedback.setProgressText("Number of grids to calculate: " + str(nGrids))

        mm = startDate[5:7]
        dd = startDate[8:10]

        index = 1

        if umepformat:
            header = '%iy  id  it imin   Q*      QH      QE      Qs      Qf    Wind    RH     Td     press   rain ' \
                        '   Kdn    snow    ldown   fcld    wuh     xsmd    lai_hr  Kdiff   Kdir    Wd'
            numformat = '%d %d %d %d %.2f %.2f %.2f %.2f %.2f %.5f %.2f %.2f %.2f %.2f %.2f %.2f %.2f ' \
                            '%.2f %.2f %.2f %.2f %.2f %.2f %.2f' 

        for f in vlayer.getFeatures():  # looping through each vector object
            feedback.setProgress(int((index * 100) / nGrids))
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break

            numtype = math.modf(f.attributes()[idx])
            if numtype[0] == 0.0:
                attr = int(f.attributes()[idx])
            
            ## generate input files for UWG
            uwgDict = read_uwg_file(inputDir, prefix + '_' + str(attr))
            uwgDict['Month'] = mm
            uwgDict['Day'] = dd
            uwgDict['nDay'] = nDays
            uwgDict['dtSim'] = dtSim
            get_uwg_file(uwgDict, inputDir, prefix + '_' + str(attr))
            
            # inputUWGfile = inputDir + '/' + prefix + '_' + str(f.attributes()[idx])  + '.uwg'
            epw_path = inputDir + '/' + prefix + '_' + str(attr)  + '.epw'
            uwg_path = inputDir + '/' + prefix + '_' + str(attr)  + '_UWG.epw'
            param_path = inputDir + '/' + prefix + '_' + str(attr)  + '.uwg'

            shutil.copy(inputMet, epw_path)

            if index == 1:
                if umepformat:
                    epwdata = np.genfromtxt(inputMet, skip_header=8, delimiter=',', filling_values=99999)
                    umep_forcing = self.epw2UMEP(epwdata)
                    np.savetxt(outputDir + '/metdata_UMEP.txt', umep_forcing, fmt=numformat, header=header, comments='')
                # else:
                #     shutil.copy(inputMet, )

            # run model
            if excludeRural:
                if uwgDict['bldDensity'] < 0.005:
                    feedback.setProgressText("Grid: " + str(attr) + ' not calculated. Less than 0.005 in building fraction.')

                    if umepformat:
                        epwdata_uwg = np.genfromtxt(epw_path, skip_header=8, delimiter=',', filling_values=99999)
                        umep_uwg = self.epw2UMEP(epwdata_uwg)
                        np.savetxt(outputDir + '/' + prefix + '_' + str(attr) +  '_UMEP_UWG.txt', umep_uwg, fmt=numformat, header=header, comments='')
                        os.remove(epw_path)
                    else:
                        shutil.move(epw_path, Path(outputDir + '/' + prefix + '_' + str(attr)  + '_UWG.epw'))
                else:
                    feedback.setProgressText("UWG calculating grid: " + str(attr))
                    try:
                        model = UWG.from_param_file(param_path, epw_path=epw_path)
                        model.generate()
                        model.simulate()
                        model.write_epw()

                        if umepformat:
                            epwdata_uwg = np.genfromtxt(uwg_path, skip_header=8, delimiter=',', filling_values=99999)
                            umep_uwg = self.epw2UMEP(epwdata_uwg)
                            np.savetxt(outputDir + '/' + prefix + '_' + str(attr) +  '_UMEP_UWG.txt', umep_uwg, fmt=numformat, header=header, comments='')
                            os.remove(uwg_path)
                        else:
                            shutil.move(uwg_path, Path(outputDir + '/' + prefix + '_' + str(attr)  + '_UWG.epw'))
                        
                        os.remove(epw_path)

                    except Exception as e:
                        feedback.pushWarning("Calculating grid " + str(attr) + ' failed: ' + str(e))
                        feedback.pushWarning('To get the full traceback error message, open the Python console in QGIS and re-run the simulation.')
                        feedback.pushWarning('If you cannot solve the error yourself, report an issue to our code reporitory (see UMEP-Manual for details).')
                        print('Traceback error message while caclulation grid: ' + str(attr))
                        print(traceback.format_exc())
                
            else:
                feedback.setProgressText("UWG calculating grid: " + str(attr))
                try:
                    model = UWG.from_param_file(param_path, epw_path=epw_path)
                    model.generate()
                    model.simulate()
                    model.write_epw()

                    if umepformat:
                        epwdata_uwg = np.genfromtxt(uwg_path, skip_header=8, delimiter=',', filling_values=99999)
                        umep_uwg = self.epw2UMEP(epwdata_uwg)
                        np.savetxt(outputDir + '/' + prefix + '_' + str(attr) +  '_UMEP_UWG.txt', umep_uwg, fmt=numformat, header=header, comments='')
                        os.remove(uwg_path)
                    else:
                        shutil.move(uwg_path, Path(outputDir + '/' + prefix + '_' + str(attr)  + '_UWG.epw'))
                    
                    os.remove(epw_path)

                except Exception as e:
                    feedback.pushWarning("Calculating grid " + str(attr) + ' failed: ' + str(e))
                    feedback.pushWarning('To get the full traceback error message, open the Python console in QGIS and re-run the simulation.')
                    feedback.pushWarning('If you cannot solve the error yourself, report an issue to our code reporitory (see UMEP-Manual for details).')
                    print('Traceback error message while caclulation grid: ' + str(attr))
                    print(traceback.format_exc())

            index += 1

        return {self.OUTPUT_DIR: outputDir}

    def epw2UMEP(self, met_old):
        met_new = np.zeros((met_old.shape[0], 24)) - 999

        # yyyy
        met_new[:, 0] = 1985
        met_new[met_old.shape[0] - 1, 0] = 1986

        # hour
        met_new[:, 2] = met_old[:, 3]
        test = met_new[:, 2] == 24
        met_new[test, 2] = 0

        # day of year
        mm = met_old[:, 1]
        dd = met_old[:, 2]
        rownum = met_old.shape[0]
        for i in range(0, rownum):
            yy = int(met_new[i, 0])
            if (yy % 4) == 0:
                if (yy % 100) == 0:
                    if (yy % 400) == 0:
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
            met_new[i, 1] = sum(dayspermonth[0:int(mm[i] - 1)]) + dd[i]

        met_new[np.where(met_new[:, 2] == 0), 1] = met_new[np.where(met_new[:, 2] == 0), 1] + 1
        met_new[met_old.shape[0] - 1, 1] = 1

        # minute
        met_new[:, 3] = 0

        # met variables
        met_new[:, 11] = met_old[:, 6]  # Ta
        met_new[:, 10] = met_old[:, 8]  # Rh
        met_new[:, 12] = met_old[:, 9] / 1000.  # P
        met_new[:, 16] = met_old[:, 12]  # Ldown
        met_new[:, 14] = met_old[:, 13]  # Kdown
        met_new[:, 22] = met_old[:, 14]  # Kdir
        met_new[:, 21] = met_old[:, 15]  # Kdiff
        met_new[:, 23] = met_old[:, 20]  # Wdir
        met_new[:, 9] = met_old[:, 21]  # Ws
        met_new[:, 13] = met_old[:, 33]  # Rain
        met_new[np.where(met_new[:, 13] == 999), 13] = 0

        return met_new
    
    def name(self):
        return 'Urban Heat Island: Urban Weather Generator'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Processor'

    def shortHelpString(self):
        return self.tr('<b>THIS PLUGIN IS EXPERIMENTAL</b>'
        '\n'
        'The <b>Urban Weather Generator</b> plugin can be used to model the urban heat island effect. Possibilities to model mutiple grids or a single location is available.<br>'
        '\n'
        'For more detailed information during execution, open the QGIS Python console (Plugins>Python Console).'
        '\n'
        '<b>NOTE</b>: This plugin requires the uwg python library. Instructions on how to install missing python libraries using the <b>pip</b> command can be found here: '
        '<a href="https://umep-docs.readthedocs.io/en/latest/Getting_Started.html">https://umep-docs.readthedocs.io/en/latest/Getting_Started.html</a>'
        '\n'
        'If you are having issues that certain grids fails to be calculated you can try to reduce the simulation time step, preferably to 150 or 100 seconds. '
        'This will increase computation time but make the model more stable.'
        '\n'
        'You can also increase stability by ticking in the box to exclude grids with very low building fraction.'
        '\n'
        '----------------------\n'
        'Full manual is available via the <b>Help</b>-button.')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/processor/Urban%20Heat%20Island%20Urban%20Weather%20Generator.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/icon_uwg.png")
        return icon

    def createInstance(self):
        return ProcessingUWGProcessorAlgorithm()