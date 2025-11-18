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
from osgeo import osr
from osgeo.gdalconst import *
import os
import sys
import numpy as np
import pandas as pd
import inspect
from pathlib import Path
import shutil
import traceback
import math
import datetime
from ..util.umep_target_export_component import write_config_file
from ..util.misc import xy2latlon

#from ..functions.target import Target
# from target import Target
# from .target.scripts.toolkit import Target
try:
    from target_py import Target
    from target_py.ui.utils import read_config
except:
    pass


class ProcessingTargetProcessorAlgorithm(QgsProcessingAlgorithm):
    """
    This algorithm make use of UWG for the processing toolbox
    """
    
    INPUT_FOLDER = 'INPUT_FOLDER'
    INPUT_POLYGONLAYER = 'INPUT_POLYGONLAYER'
    # ID_FIELD = 'ID_FIELD'
    START_DATE = 'START_DATE'
    START_DATE_INTEREST = 'START_DATE_INTEREST'
    STOP_DATE_INTEREST = 'STOP_DATE_INTEREST'
    INPUT_MET = 'INPUT_MET'
    OUTPUT_DIR = 'OUTPUT_DIR'
    OUTPUT_CSV = 'OUTPUT_CSV'
    OUTPUT_UMEP = 'OUTPUT_UMEP'
    # DTSIM = 'DTSIM'
    MOD_LDOWN = 'MOD_LDOWN'
    RUN_NAME = 'RUN_NAME'


    def initAlgorithm(self, config):
        self.addParameter(QgsProcessingParameterFile(self.INPUT_FOLDER,
            self.tr('Path to folder where TARGET input files are located (Site name folder)'),
            QgsProcessingParameterFile.Folder))
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT_POLYGONLAYER,
            self.tr('Vector polygon grid'),
            [QgsProcessing.TypeVector]))
        # self.addParameter(QgsProcessingParameterField(self.ID_FIELD,
        #     self.tr('ID field'), '', self.INPUT_POLYGONLAYER,
        #     QgsProcessingParameterField.Numeric))
        self.addParameter(QgsProcessingParameterString(self.RUN_NAME, 
            self.tr('Run name')))
        self.addParameter(QgsProcessingParameterDateTime(self.START_DATE,
            self.tr('Start date of simulation (should be a minimum of 24 hours prior to date for period of interest)'),
            QgsProcessingParameterDateTime.Date))
        self.addParameter(QgsProcessingParameterDateTime(self.START_DATE_INTEREST,
            self.tr('Start date for period of interest'),
            QgsProcessingParameterDateTime.Date))
        self.addParameter(QgsProcessingParameterDateTime(self.STOP_DATE_INTEREST,
            self.tr('End date for period of interest'),
            QgsProcessingParameterDateTime.Date))
        # self.addParameter(QgsProcessingParameterNumber(self.NDAYS,
        #     self.tr('Number of days to run from period of interest'),
        #     QgsProcessingParameterNumber.Integer,
        #     QVariant(5), False, minValue=1, maxValue=365))
        self.addParameter(QgsProcessingParameterFile(self.INPUT_MET,
            self.tr('Input meteorological file (UMEP-formatted textfile)'),
            extension = 'txt'))
        # self.addParameter(QgsProcessingParameterNumber(self.DTSIM,
        #     self.tr('Simulation time step in minutes'),
        #     QgsProcessingParameterNumber.Integer,
        #     QVariant(300), False, minValue=1, maxValue=1440))
        self.addParameter(QgsProcessingParameterBoolean(self.MOD_LDOWN,
            self.tr('Estimate incoming longwave radiation from air temperture and relative humidity'), 
            defaultValue=False))
        # output
        # self.addParameter(QgsProcessingParameterFolderDestination(self.OUTPUT_DIR,
            # self.tr('Output folder')))
        self.addParameter(QgsProcessingParameterBoolean(self.OUTPUT_CSV,
            self.tr('Save output as .csv text files')))
        self.addParameter(QgsProcessingParameterBoolean(self.OUTPUT_UMEP,
            self.tr('Save output in UMEP-specific format (required for the TARGET Analyzer)')))
        
        self.plugin_dir = os.path.dirname(__file__)


    def processAlgorithm(self, parameters, context, feedback):

        # InputParameters
        inputDir = self.parameterAsString(parameters, self.INPUT_FOLDER, context)
        inputPolygonlayer = self.parameterAsVectorLayer(parameters, self.INPUT_POLYGONLAYER, context)
        startDate = self.parameterAsString(parameters, self.START_DATE, context)
        startDateInterest = self.parameterAsString(parameters, self.START_DATE_INTEREST, context)
        stopDateInterest = self.parameterAsString(parameters, self.STOP_DATE_INTEREST, context)
        inputMet = self.parameterAsString(parameters, self.INPUT_MET, context)
        # outputDir = self.parameterAsString(parameters, self.OUTPUT_DIR, context)
        outputCSV = self.parameterAsBoolean(parameters, self.OUTPUT_CSV, context)
        # dtSim = self.parameterAsDouble(parameters, self.DTSIM, context)
        mod_Ldown = self.parameterAsBoolean(parameters, self.MOD_LDOWN, context)
        runName = self.parameterAsString(parameters, self.RUN_NAME, context)
        umepformat = self.parameterAsBoolean(parameters, self.OUTPUT_UMEP,context)
        # getting extent, lat lon, and number of x and y grids
        vlayer = inputPolygonlayer
        ext = vlayer.extent()
        xExt = ext.xMaximum() - ext.xMinimum()
        yExt = ext.yMaximum() - ext.yMinimum()

        grid_crs = osr.SpatialReference()
        grid_crs.ImportFromWkt(vlayer.crs().toWkt())
        grid_unit = grid_crs.GetAttrValue('UNIT')
        possible_units_metre = ['metre', 'Metre', 'metres', 'Metres', 'meter', 'Meter', 'meters', 'Meters', 'm'] 
        if not grid_unit in possible_units_metre:
            raise QgsProcessingException('Error! Raster projection is not in meters or feet. Please reproject.')

        features = vlayer.getFeatures()
        for feature in features:
            geom = feature.geometry()
            gridsize = np.round(geom.length() / 4)
            nGridX = int(xExt / gridsize)
            nGridY = int(yExt / gridsize)
            break

        latmin, lonmin = xy2latlon(vlayer.crs().toWkt(), ext.xMinimum(), ext.yMinimum())
        latmax, lonmax = xy2latlon(vlayer.crs().toWkt(), ext.xMaximum(), ext.yMaximum())

        # Converting UMEP met-file to target met-file
        try:
            metfile= pd.read_csv(inputMet, sep='\s+')
        except:
            raise QgsProcessingException("Error: Make sure format of meteorological file is correct. You can"
                                                "prepare your data by using 'Prepare Existing Data' in "
                                                "the Pre-processor")
        
        metfile['datetime'] = pd.to_datetime(metfile[['%iy','id','it','imin']].astype(str).agg('-'.join, axis=1), format='%Y-%j-%H-%M')  
        metfile=metfile[['datetime', 'Td', 'RH', 'Wind', 'press','Kdn','ldown']]
        metfile.columns = ['datetime', 'Ta', 'RH', 'WS', 'P','Kd','Ld']
        startmetfile = metfile['datetime'].min()
        endmetfile = metfile['datetime'].max()
        t = metfile['datetime'].iloc[1]-metfile['datetime'].iloc[0]
        tdiff = int(t.total_seconds()/60)
        metfile.set_index('datetime', inplace=True)
        metfile.to_csv(inputDir + '/input/MET/' + runName + '_metdata.txt',sep=',',header='datetime,Ta,RH,WS,P,Kd,Ld')

        # creating config.ini
        cfM = read_config(self.plugin_dir + '/configtarget.ini')

        cfM['work_dir'] = os.path.dirname(inputDir)
        cfM['para_json_path'] = inputDir + '/parameters.json'
        cfM['site_name'] = os.path.basename(inputDir)
        cfM['run_name'] = runName
        cfM['inpt_met_file'] = runName + '_metdata.txt'
        cfM['inpt_lc_file'] = 'lc_target.txt'
        cfM['date_fmt'] = '%Y-%m-%d %H:%M:%S'
        cfM['timestep'] = str(tdiff) # timestep is set from input forcingfile
        cfM['include roofs'] = 'Y'

        if mod_Ldown:
            cfM['mod_ldown'] = 'Y'
        else:
            cfM['mod_ldown'] = 'N'
            if -999 in metfile['Ld'].values:
                cfM['mod_ldown'] = 'Y'
                feedback.pushWarning("-999 found in Ldown. Ldown will be modelled.")

        cfM['domaindim'] = str(nGridX) + ',' + str(nGridY)
        cfM['latedge'] = str(latmin)
        cfM['lonedge'] = str(lonmin)
        cfM['latresolution'] = str(abs(latmin - latmax))
        cfM['lonresolution'] = str(abs(lonmin - lonmax))
        cfM['date1a'] = datetime.datetime.strptime(startDate, "%Y-%m-%d").strftime("%Y,%m,%d") + ',0'
        cfM['date1'] = datetime.datetime.strptime(startDateInterest, "%Y-%m-%d").strftime("%Y,%m,%d") + ',0'
        cfM['date2'] = datetime.datetime.strptime(stopDateInterest, "%Y-%m-%d").strftime("%Y,%m,%d") + ',0'

        #TODO: Check so that dates are within forcing data
        date1a = datetime.datetime.strptime(startDate, "%Y-%m-%d")
        date1 = datetime.datetime.strptime(startDateInterest, "%Y-%m-%d")
        date2 = datetime.datetime.strptime(stopDateInterest, "%Y-%m-%d")

        if date1a >= date1:
            raise QgsProcessingException("Error: Start date should be at least 24h before Start date for period of interest")
        if date1 >= date2:
            raise QgsProcessingException("Error: Start date of Interest should be later than End date for period of interest")

        if startmetfile <= date1a <= endmetfile:
            test =4
        else:
            raise QgsProcessingException("Error: Start date selected is not within dates available in the meteorological focring file")
        if startmetfile <= date2 <= endmetfile:
            test =4
        else:
            raise QgsProcessingException("Error: End date selected is not within dates available in the meteorological focring file")

        write_config_file(cfM, inputDir)

        ### run model ###
        feedback.setProgressText("Model calculation started.")
        start = datetime.datetime.now()
        tar = Target(
            os.path.join(inputDir, "config.ini"),            # passing the simulation's config file
            progress=True                                    # show progress bars in console
        )
        tar.load_config()

        # run simulation
        if outputCSV:
            tar.run_simulation(save_csv=True)      # save model results in csv format
        else:
            tar.run_simulation(save_csv=False)
        
        # save parameters and config used for simulation
        tar.save_simulation_parameters()
        end = datetime.datetime.now() - start
        feedback.setProgressText("Model calculation finished. Output is found in " + inputDir + "/output")
        feedback.setProgressText("Model calculation time: " + str(end.total_seconds()) + " seconds")

        #saving output as umep formatted metfile
        if umepformat:
            feedback.setProgressText('Saving data in UMEP formatted text-files.')
            header = '%iy  id  it imin   Q*      QH      QE      Qs      Qf    Wind    RH     Td     press   rain ' \
                        '   Kdn    snow    ldown   fcld    wuh     xsmd    lai_hr  Kdiff   Kdir    Wd'
            numformat = '%d %d %d %d %.2f %.2f %.2f %.2f %.2f %.5f %.2f %.2f %.2f %.2f %.2f %.2f %.2f ' \
                            '%.2f %.2f %.2f %.2f %.2f %.2f %.2f'

            dataout = np.load(inputDir + "/output/" + runName + '.npy', allow_pickle=True)
            dfin = pd.read_csv(inputMet, sep='\s+')
            dfin['datetime'] = pd.to_datetime(dfin[['%iy', 'id', 'it', 'imin']].astype(str).agg('-'.join, axis=1), format='%Y-%j-%H-%M')
            dfin.set_index('datetime', inplace=True)
            np.savetxt(inputDir + "/output/" + runName + '_metdata_UMEP.txt', dfin, fmt=numformat, header=header, comments='')
            index = 1

            for f in range(0, dataout.shape[1]):  # looping through each grid saving data
                feedback.setProgress(int((index * 100) / dataout.shape[1]))
                if feedback.isCanceled():
                    feedback.setProgressText("Calculation cancelled")
                    break
                gridID = dataout[:,f,0][0][0]
                dfout001 = pd.DataFrame(dataout[:,f,0])
                dfout001['date'] = pd.to_datetime(dfout001['date'])
                dfout001.set_index('date', inplace=True)
                dfoutmod = dfout001['Ta'].loc[date1:date2]
                dfin.loc[date1:date2,['Td']] = dfoutmod
                np.savetxt(inputDir + "/output/" + runName + '_' + str(int(gridID)) + '_UMEP_TARGET.txt', dfin, fmt=numformat, header=header, comments='')
                
        feedback.setProgressText("Process finished")

        return {self.INPUT_FOLDER: inputDir}
    

    def name(self):
        return 'Urban Heat Island: TARGET'

    def displayName(self):
        return self.tr(self.name())

    def group(self):
        return self.tr(self.groupId())

    def groupId(self):
        return 'Processor'

    def shortHelpString(self):
        return self.tr('<b>THIS TOOL IS EXPERIMENTAL</b>'
        '\n'
        '<b>TARGET</b> (The Air-temperature Response to Green blue-infrastructure Evaluaition Tool) is a simple modelling framework used to examine '
        'intra urban climate. It has specifically been developed as an efficient, easy-to-use model albe to investigate '
        'heat mitigation effects of green and blue infrastructure within urban areas but can also be used to model the canopy urban heat island '
        '<a href="https://gmd.copernicus.org/articles/12/785/2019/">(Broadbent et al. 2019)</a>. '
        'Possibilities to model mutiple grids or a single location is available.\n'
        '\n'
        'For more detailed information during execution, open the QGIS Python console (Plugins>Python Console).'
        '\n'
        '<b>NOTE</b>: This plugin requires the <b>target-py</b> python library. If UMEP was installed without issues, target-py should be installed on your system. If missing on your system, instructions on how to install missing python libraries using the pip command can be found here: '
        '<a href="https://umep-docs.readthedocs.io/en/latest/Getting_Started.html">https://umep-docs.readthedocs.io/en/latest/Getting_Started.html</a>'
        '\n'
        '----------------------\n'
        'Full manual is available via the <b>Help</b>-button.')

    def helpUrl(self):
        url = "https://umep-docs.readthedocs.io/en/latest/processor/Urban%20Heat%20Island%20TARGET.html"
        return url

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def icon(self):
        cmd_folder = Path(os.path.split(inspect.getfile(inspect.currentframe()))[0]).parent
        icon = QIcon(str(cmd_folder) + "/icons/icon_uwg.png")
        return icon

    def createInstance(self):
        return ProcessingTargetProcessorAlgorithm()