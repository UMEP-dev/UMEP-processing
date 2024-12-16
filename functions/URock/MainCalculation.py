#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 11:39:39 2021

@author: Jérémy Bernard, University of Gothenburg
"""

from .GlobalVariables import * 

from . import H2gisConnection
from . import loadData
from . import saveData
from . import Obstacles
from . import Zones
from . import CalculatesIndicators
from . import InitWindField
from . import DataUtil
from . import WindSolver
import time
try :
    from numba import jit
except ImportError:
    exit("'numba' Python package is missing")
#import copy as cp
from pathlib import Path
from qgis.core import QgsProcessingException

import os

def main(javaEnvironmentPath,
         pluginDirectory,
         outputFilePath,
         buildingFilePath,
         srid,
         outputFilename = OUTPUT_FILENAME,
         vegetationFilePath = "",
         z_ref = Z_REF,
         v_ref = V_REF,
         windDirection = WIND_DIRECTION,
         prefix = PREFIX_NAME,
         meshSize = MESH_SIZE,
         dz = DZ,
         alongWindZoneExtend = ALONG_WIND_ZONE_EXTEND,
         crossWindZoneExtend = CROSS_WIND_ZONE_EXTEND,
         verticalExtend = VERTICAL_EXTEND,
         tempoDirectory = TEMPO_DIRECTORY,
         inputDirectory = INPUT_DIRECTORY,
         outputDirectory = OUTPUT_DIRECTORY,
         cadTriangles = CAD_TRIANGLE_FILENAME,
         cadTreesIntersection = CAD_VEG_INTERSECTION_FILENAME,
         onlyInitialization = ONLY_INITIALIZATION,
         maxIterations = MAX_ITERATIONS,
         thresholdIterations = THRESHOLD_ITERATIONS,
         idFieldBuild = ID_FIELD_BUILD,
         buildingHeightField = HEIGHT_FIELD,
         vegetationBaseHeight = VEGETATION_CROWN_BASE_HEIGHT,
         vegetationTopHeight = VEGETATION_CROWN_TOP_HEIGHT,
         idVegetation = ID_VEGETATION,
         vegetationAttenuationFactor = VEGETATION_ATTENUATION_FACTOR,
         saveRockleZones = SAVE_ROCKLE_ZONES,
         z_out = Z_OUT,
         outputRaster = None,
         feedback = None,
         saveRaster = True,
         saveVector = True,
         saveNetcdf = True,
         debug = DEBUG,
         profileType = PROFILE_TYPE,
         verticalProfileFile = None):
    # If the function is called within QGIS, a feedback is sent into the QGIS interface
    if feedback:
        feedback.setProgressText('Initiating algorithm')
    
    ################################ INIT OUTPUT VARIABLES ############################
    # Define dictionaries of input and output relative directories
    outputDataRel = {}

    # Blocks and stacked blocks
    outputDataRel["blocks"] = os.path.join(tempoDirectory, "blocks.geojson")
    outputDataRel["stacked_blocks"] = os.path.join(tempoDirectory, "stackedBlocks.geojson")
    outputDataRel["vegetation"] = os.path.join(tempoDirectory, "vegetation.geojson")

    # Rotated geometries
    outputDataRel["rotated_stacked_blocks"] = os.path.join(tempoDirectory, "rotated_stacked_blocks.geojson")
    outputDataRel["rotated_vegetation"] = os.path.join(tempoDirectory, "vegetationRotated.geojson")
    outputDataRel["upwind_facades"] = os.path.join(tempoDirectory, "upwind_facades.geojson")
    outputDataRel["downwind_facades"] = os.path.join(tempoDirectory, "downwind_facades.geojson")
    
    # Created zones
    outputDataRel["displacement"] = os.path.join(tempoDirectory, "displacementZones.geojson")
    outputDataRel["displacement_vortex"] = os.path.join(tempoDirectory, "displacementVortexZones.geojson")
    outputDataRel["cavity"] = os.path.join(tempoDirectory, "cavity.geojson")
    outputDataRel["wake"] = os.path.join(tempoDirectory, "wake.geojson")
    outputDataRel["street_canyon"] = os.path.join(tempoDirectory, "streetCanyon.geojson")
    outputDataRel["rooftop_perpendicular"] = os.path.join(tempoDirectory, "rooftopPerp.geojson")
    outputDataRel["rooftop_corner"] = os.path.join(tempoDirectory, "rooftopCorner.geojson")
    outputDataRel["vegetation_built"] = os.path.join(tempoDirectory, "vegetationBuilt.geojson")
    outputDataRel["vegetation_open"] = os.path.join(tempoDirectory, "vegetationOpen.geojson")
    
    # Grid points
    outputDataRel["point3D_BuildZone"] = os.path.join(tempoDirectory, "point3D_BuildZone")
    outputDataRel["point3D_VegZone"] = os.path.join(tempoDirectory, "point3D_VegZone")
    outputDataRel["point3D_All"] = os.path.join(tempoDirectory, "point3D_All")
    
    # Put 2D grid points in the output directory
    outputDataRel["point_2DRockleZone"] = os.path.join(outputFilePath, "Rockle_zones")
    
    # Convert relative to absolute paths
    outputDataAbs = {i : os.path.abspath(outputDataRel[i]) for i in outputDataRel}
    
    ############################################################################
    ################################ SCRIPT ####################################
    ############################################################################
    # ----------------------------------------------------------------------
    # 1. SET H2GIS DATABASE ENVIRONMENT AND LOAD DATA
    # ----------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Creates an H2GIS Instance and load data')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    #Download H2GIS
    #H2gisConnection.downloadH2gis(dbDirectory = pluginDirectory)
    #Initialize a H2GIS database connection
    dBDir = os.path.join(Path(pluginDirectory).parent, 'functions','URock')
    #print(dBDir)
    if DEBUG:
        db_suffix = ""
    else:
        db_suffix = str(time.time()).replace(".", "_")
    cursor, conn, localH2InstanceDir = \
        H2gisConnection.startH2gisInstance(dbDirectory = dBDir,
                                           dbInstanceDir = tempoDirectory,
                                           suffix = db_suffix)
        
    # Load data
    loadData.loadData(fromCad = False, 
                      prefix = prefix,
                      idFieldBuild = idFieldBuild,
                      buildingHeightField = buildingHeightField,
                      vegetationBaseHeight = vegetationBaseHeight,
                      vegetationTopHeight = vegetationTopHeight,
                      idVegetation = idVegetation,
                      vegetationAttenuationFactor = vegetationAttenuationFactor,
                      cursor = cursor,
                      buildingFilePath = buildingFilePath,
                      vegetationFilePath = vegetationFilePath,
                      srid = srid)
    
    timeStartCalculation = time.time()
    
    # -----------------------------------------------------------------------------------
    # 2. CREATES OBSTACLE GEOMETRIES ----------------------------------------------------
    # -----------------------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Creates the stacked blocks used as obstacles')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Create the stacked blocks
    blockTable, stackedBlockTable = \
        Obstacles.createsBlocks(cursor = cursor, 
                                inputBuildings = BUILDING_TABLE_NAME,
                                prefix = prefix)
    
    # Save the blocks, stacked blocks and vegetation as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                          , tableName = blockTable,
                           filedir = outputDataAbs["blocks"]        , delete = True)
        saveData.saveTable(cursor = cursor                          , tableName = VEGETATION_TABLE_NAME,
                           filedir = outputDataAbs["vegetation"]    , delete = True)
    
    # -----------------------------------------------------------------------------------
    # 3. ROTATES OBSTACLES TO THE RIGHT DIRECTION AND CALCULATES GEOMETRY PROPERTIES ----
    # -----------------------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Rotates obstacles to the right direction and calculates geometry properties')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Define a set of obstacles in a dictionary before the rotation
    dicOfObstacles = {BUILDING_TABLE_NAME       : stackedBlockTable,
                      VEGETATION_TABLE_NAME     : VEGETATION_TABLE_NAME}
    
    # Rotate obstacles
    dicRotatedTables, rotationCenterCoordinates = \
        Obstacles.windRotation(cursor = cursor,
                               dicOfInputTables = dicOfObstacles,
                               rotateAngle = windDirection,
                               rotationCenterCoordinates = None,
                               prefix = prefix)
    
    # Get the rotated block and vegetation table names
    rotatedStackedBlocks = dicRotatedTables[BUILDING_TABLE_NAME]
    rotatedVegetation = dicRotatedTables[VEGETATION_TABLE_NAME]
    
    # Calculates base block height and base of block cavity zone
    rotatedPropStackedBlocks = \
        Obstacles.identifyBlockAndCavityBase(cursor, rotatedStackedBlocks,
                                             prefix = prefix)
    
    # Calculates obstacles properties
    obstaclePropertiesTable = \
        CalculatesIndicators.obstacleProperties(cursor = cursor,
                                                obstaclesTable = rotatedPropStackedBlocks,
                                                prefix = prefix)
    
    # Calculates obstacle zone properties
    zonePropertiesTable = \
        CalculatesIndicators.zoneProperties(cursor = cursor,
                                            obstaclePropertiesTable = obstaclePropertiesTable,
                                            prefix = prefix)
    
    # Init the upwind facades
    upwindInitedTable = \
        Obstacles.initUpwindFacades(cursor = cursor,
                                    obstaclesTable = zonePropertiesTable,
                                    prefix = prefix)
    # Update base height of upwind facades (if shared with the building below)
    upwindTable = \
        Obstacles.updateUpwindFacadeBase(cursor = cursor,
                                        upwindTable = upwindInitedTable,
                                        prefix = prefix)
        
    # Calculates downwind facades 
    downwindTable = \
        Obstacles.initDownwindFacades(cursor = cursor,
                                      obstaclesTable = zonePropertiesTable,
                                      prefix = prefix)
    
    # Calculates roughness properties of the study area
    z0, d, Hr, H_ob_max, lambda_f = \
        CalculatesIndicators.studyAreaProperties(cursor = cursor, 
                                                 upwindTable = upwindInitedTable, 
                                                 stackedBlockTable = rotatedStackedBlocks, 
                                                 vegetationTable = rotatedVegetation)

    # Print some of the roughness properties for the zone
    if feedback:
        feedback.setProgressText(f"""Roughness zone properties are:\n
                                     - z0: {z0}
                                     - d: {d}
                                     - Hr: {Hr}
                                     - H_ob_max: {H_ob_max}
                                     - lambda_f: {lambda_f}""")
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}


    # Save the rotated obstacles and facades as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                                  , tableName = rotatedPropStackedBlocks,
                           filedir = outputDataAbs["rotated_stacked_blocks"], delete = True)
        saveData.saveTable(cursor = cursor                         , tableName = rotatedVegetation,
                           filedir = outputDataAbs["rotated_vegetation"]    , delete = True)
        saveData.saveTable(cursor = cursor                      , tableName = upwindTable,
                           filedir = outputDataAbs["upwind_facades"]   , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
        saveData.saveTable(cursor = cursor                      , tableName = downwindTable,
                           filedir = outputDataAbs["downwind_facades"]   , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    # Used later to save results
    saveData.saveTable(cursor = cursor                          , tableName = rotatedPropStackedBlocks,
                       rotationCenterCoordinates = rotationCenterCoordinates,
                       rotateAngle = - windDirection,
                       filedir = outputDataAbs["stacked_blocks"], delete = True)
    
    
    # -----------------------------------------------------------------------------------
    # 4. CREATES THE 2D ROCKLE ZONES ----------------------------------------------------
    # -----------------------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Creates the 2D Röckle zones')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Creates the displacement zone (upwind)
    # displacementZonesTable, displacementVortexZonesTable = \
    #     Zones.displacementZones(cursor = cursor,
    #                             upwindTable = upwindTable,
    #                             zonePropertiesTable = zonePropertiesTable,
    #                             srid = srid,
    #                             prefix = prefix)
    displacementZonesTable, displacementVortexZonesTable = \
        Zones.displacementZones2(cursor = cursor,
                                  upwindWithPropTable = upwindTable,
                                  srid = srid,
                                  prefix = prefix)
    # # Creates the displacement zone (upwind)
    # displacementZonesTable = \
    #     Zones.displacementZones(cursor = cursor,
    #                             upwindTable = upwindTable,
    #                             zonePropertiesTable = zonePropertiesTable,
    #                             srid = srid,
    #                             prefix = prefix)[0]
    # displacementVortexZonesTable = \
    #     Zones.displacementZonesVortex(cursor = cursor,
    #                                   upwindWithPropTable = upwindTable,
    #                                   srid = srid,
    #                                   prefix = prefix)[0]
    
    
    # Save the resulting displacement zones as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                      , tableName = displacementZonesTable,
                  filedir = outputDataAbs["displacement"]       , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
        saveData.saveTable(cursor = cursor                          , tableName = displacementVortexZonesTable,
                  filedir = outputDataAbs["displacement_vortex"]    , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    # Creates the cavity and wake zones
    cavityZonesTable, wakeZonesTable = \
        Zones.cavityAndWakeZones(cursor = cursor, 
                                downwindWithPropTable = downwindTable,
                                srid = srid,
                                ellipseResolution = meshSize/3,
                                prefix = prefix).values()
    
    # Save the resulting displacement zones as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor             , tableName = cavityZonesTable,
                  filedir = outputDataAbs["cavity"]    , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
        saveData.saveTable(cursor = cursor           , tableName = wakeZonesTable,
                  filedir = outputDataAbs["wake"]    , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    
    # Creates the street canyon zones
    streetCanyonTable = \
        Zones.streetCanyonZones(cursor = cursor,
                                cavityZonesTable = cavityZonesTable,
                                zonePropertiesTable = zonePropertiesTable,
                                upwindTable = upwindTable,
                                downwindTable = downwindTable,
                                srid = srid,
                                prefix = prefix)
    
    # Save the resulting street canyon zones as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                    , tableName = streetCanyonTable,
                  filedir = outputDataAbs["street_canyon"]    , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    # Creates the rooftop zones
    rooftopPerpendicularZoneTable, rooftopCornerZoneTable = \
        Zones.rooftopZones(cursor = cursor,
                           upwindTable = upwindTable,
                           zonePropertiesTable = zonePropertiesTable,
                           prefix = prefix)
    # Save the resulting rooftop zones as fgb
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                              , tableName = rooftopPerpendicularZoneTable,
                  filedir = outputDataAbs["rooftop_perpendicular"]      , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
        saveData.saveTable(cursor = cursor                      , tableName = rooftopCornerZoneTable,
                  filedir = outputDataAbs["rooftop_corner"]     , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    # Creates the vegetation zones
    vegetationBuiltZoneTable, vegetationOpenZoneTable = \
        Zones.vegetationZones(cursor = cursor,
                              vegetationTable = rotatedVegetation,
                              wakeZonesTable = wakeZonesTable,
                              prefix = prefix)
    if debug or saveRockleZones:
        saveData.saveTable(cursor = cursor                              , tableName = vegetationBuiltZoneTable,
                           filedir = outputDataAbs["vegetation_built"]  , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
        saveData.saveTable(cursor = cursor                              , tableName = vegetationOpenZoneTable,
                           filedir = outputDataAbs["vegetation_open"]   , delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    # Define a dictionary of all building Rockle zones and same for veg
    dicOfBuildRockleZoneTable = {DISPLACEMENT_NAME       : displacementZonesTable,
                                DISPLACEMENT_VORTEX_NAME: displacementVortexZonesTable,
                                CAVITY_NAME             : cavityZonesTable,
                                WAKE_NAME               : wakeZonesTable,
                                STREET_CANYON_NAME      : streetCanyonTable,
                                ROOFTOP_PERP_NAME       : rooftopPerpendicularZoneTable,
                                ROOFTOP_CORN_NAME       : rooftopCornerZoneTable}
    dicOfVegRockleZoneTable = {VEGETATION_BUILT_NAME   : vegetationBuiltZoneTable,
                               VEGETATION_OPEN_NAME    : vegetationOpenZoneTable}    
    
    if outputRaster:
        # Creates a table with a polygon covering the raster zone envelope
        smallStudyZone = "SMALL_STUDY_ZONE"
        outputRasterExtent = outputRaster.extent()
        cursor.execute("""
           DROP TABLE IF EXISTS {0};
           CREATE TABLE {0}({5} GEOMETRY)
               AS SELECT ST_SETSRID(ST_ROTATE(ST_ENVELOPE('MULTIPOINT({1} {2},
                                               {3} {4})'),
                                              {6},
                                              {7},
                                              {8}), {9})
           """.format(smallStudyZone,
                       outputRasterExtent.xMinimum(),
                       outputRasterExtent.yMinimum(),
                       outputRasterExtent.xMaximum(),
                       outputRasterExtent.yMaximum(),
                       GEOM_FIELD,
                       DataUtil.degToRad(windDirection),
                       rotationCenterCoordinates[0],
                       rotationCenterCoordinates[1],
                       srid))
        # Identify the stacked blocks, blocks potentially impacting the
        # impacted zone and their corresponding Röckle zones 
        dicOfBuildRockleZoneTable, dicOfVegRockleZoneTable, rotatedPropStackedBlocks,\
        rotatedVegetation = \
            Zones.identifyImpactingStackedBlocks(cursor = cursor,
                                                 dicOfBuildRockleZoneTable = dicOfBuildRockleZoneTable,
                                                 dicOfVegRockleZoneTable = dicOfVegRockleZoneTable,
                                                 impactedZone = smallStudyZone,
                                                 stackedBlocksTable = rotatedPropStackedBlocks,
                                                 vegetationTable = rotatedVegetation,
                                                 crossWindExtend = crossWindZoneExtend,                                                 
                                                 prefix = prefix)
    # ----------------------------------------------------------------------
    # 5. SET THE 2D GRID IN THE ROCKLE ZONES -------------------------------
    # ----------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Creates the 2D grid')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
        
    # Creates the grid of points
    gridPoint = InitWindField.createGrid(cursor = cursor, 
                                         dicOfInputTables = dict(dicOfBuildRockleZoneTable,
                                                                 **dicOfVegRockleZoneTable),
                                         srid = srid,
                                         alongWindZoneExtend = alongWindZoneExtend, 
                                         crossWindZoneExtend = crossWindZoneExtend, 
                                         meshSize = meshSize,
                                         prefix = prefix)
    
    # Affects each 2D point to a build Rockle zone and calculates needed variables for 3D wind speed factors
    dicOfInitBuildZoneGridPoint, verticalLineTable = \
        InitWindField.affectsPointToBuildZone(  cursor = cursor, 
                                                gridTable = gridPoint,
                                                dicOfBuildRockleZoneTable = dicOfBuildRockleZoneTable,
                                                prefix = prefix)
        
    # Same for vegetation Röckle zones
    dicOfVegZoneGridPoint = \
        InitWindField.affectsPointToVegZone(cursor = cursor, 
                                            gridTable = gridPoint,
                                            dicOfVegRockleZoneTable = dicOfVegRockleZoneTable,
                                            prefix = prefix)
    
    # Remove some of the Röckle points where building Röckle zones overlap
    dicOfBuildZoneGridPoint = \
        InitWindField.removeBuildZonePoints(cursor = cursor, 
                                            dicOfInitBuildZoneGridPoint = dicOfInitBuildZoneGridPoint,
                                            prefix = prefix)
    
    # Manage backward cavity and wake zones in the leeward zone of tall buildings
    dicOfBuildZoneGridPoint, facadeWithinCavity =\
        InitWindField.manageBackwardZones(cursor = cursor, 
                                          dicOfBuildZoneGridPoint = dicOfBuildZoneGridPoint,
                                          cavity2dInitPoints = dicOfInitBuildZoneGridPoint[CAVITY_NAME],
                                          wake2dInitPoints = dicOfInitBuildZoneGridPoint[WAKE_NAME],
                                          streetCanyonTable = streetCanyonTable,
                                          gridTable = gridPoint,
                                          meshSize = meshSize,
                                          dz = dz,
                                          prefix = prefix)
    
    
    # -----------------------------------------------------------------------------------
    # 6. INITIALIZE THE 3D WIND FACTORS IN THE ROCKLE ZONES -------------------------------
    # -----------------------------------------------------------------------------------   
    if feedback:
        feedback.setProgressText('Initializes the 3D grid within Röckle zones')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Calculates the 3D wind speed factors for each building Röckle zone
    dicOfBuildZone3DWindFactor, maxBuildZoneHeight = \
        InitWindField.calculates3dBuildWindFactor(cursor = cursor,
                                                  dicOfBuildZoneGridPoint = dicOfBuildZoneGridPoint,
                                                  dz = dz,
                                                  prefix = prefix)
    if debug or saveRockleZones:
        for t in dicOfBuildZone3DWindFactor:
            cursor.execute("""
               DROP TABLE IF EXISTS point3D_Buildzone_{0};
               {5};
               {6};
               CREATE TABLE point3D_Buildzone_{0}
                   AS SELECT   a.{2}, b.*
                   FROM {3} AS a RIGHT JOIN {4} AS b
                       ON a.{1} = b.{1}
                   WHERE b.{1} IS NOT NULL
               """.format( t                            , ID_POINT,
                           GEOM_FIELD                   , gridPoint, 
                           dicOfBuildZone3DWindFactor[t], DataUtil.createIndex(tableName=gridPoint, 
                                                                               fieldName=ID_POINT,
                                                                               isSpatial=False),
                           DataUtil.createIndex(tableName=dicOfBuildZone3DWindFactor[t], 
                                                fieldName=ID_POINT,
                                                isSpatial=False)))
            saveData.saveTable(cursor = cursor,
                               tableName = "point3D_Buildzone_"+t,
                               filedir = outputDataAbs["point3D_BuildZone"]+t+".geojson",
                               delete = True,
                               rotationCenterCoordinates = rotationCenterCoordinates,
                               rotateAngle = - windDirection)
        
    # Calculates the 3D wind speed factors of the vegetation (considering all zone types)
    # after calculation of the top of the "sketch"
    maxHeight = H_ob_max
    if maxBuildZoneHeight: 
        if maxBuildZoneHeight > H_ob_max:
            maxHeight = maxBuildZoneHeight
    sketchHeight = maxHeight + verticalExtend
    vegetationWeightFactorTable = \
        InitWindField.calculates3dVegWindFactor(cursor = cursor,
                                                dicOfVegZoneGridPoint = dicOfVegZoneGridPoint,
                                                sketchHeight = sketchHeight,
                                                z0 = z0,
                                                d = d,
                                                dz = dz,
                                                prefix = prefix)
    if debug or saveRockleZones:
        cursor.execute("""
           DROP TABLE IF EXISTS point3D_AllVegZone;
           {4};
           {5};
           CREATE TABLE point3D_AllVegZone
               AS SELECT   a.{1}, b.*
               FROM {2} AS a RIGHT JOIN {3} AS b
                           ON a.{0} = b.{0}
               WHERE b.{0} IS NOT NULL
           """.format( ID_POINT                     , GEOM_FIELD, 
                       gridPoint                    , vegetationWeightFactorTable,
                       DataUtil.createIndex(tableName=gridPoint, 
                                            fieldName=ID_POINT,
                                            isSpatial=False),
                       DataUtil.createIndex(tableName=vegetationWeightFactorTable, 
                                            fieldName=ID_POINT,
                                            isSpatial=False)))
        saveData.saveTable(cursor = cursor,
                           tableName = "point3D_AllVegZone",
                           filedir = outputDataAbs["point3D_VegZone"]+".geojson",
                           delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)
    
    
    # ----------------------------------------------------------------
    # 7. DEALS WITH SUPERIMPOSED ZONES -------------------------------
    # ----------------------------------------------------------------
    # Calculates the final weighting factor for each point, dealing with duplicates (superimposition)
    dicAllWeightFactorsTables = dicOfBuildZone3DWindFactor.copy()
    dicAllWeightFactorsTables[ALL_VEGETATION_NAME] = vegetationWeightFactorTable
    allZonesPointFactor = \
        InitWindField.manageSuperimposition(cursor = cursor,
                                            dicAllWeightFactorsTables = dicAllWeightFactorsTables,
                                            facadeWithinCavity = facadeWithinCavity,
                                            upstreamPriorityTables = UPSTREAM_PRIORITY_TABLES,
                                            upstreamWeightingTables = UPSTREAM_WEIGHTING_TABLES,
                                            upstreamWeightingInterRules = UPSTREAM_WEIGHTING_INTER_RULES,
                                            upstreamWeightingIntraRules = UPSTREAM_WEIGHTING_INTRA_RULES,
                                            downstreamWeightingTable = DOWNSTREAM_WEIGTHING_TABLE,
                                            prefix = prefix,
                                            feedback = feedback)
    if debug or saveRockleZones:
        cursor.execute("""
            DROP TABLE IF EXISTS point3D_All;
            {4};
            {5};
            CREATE TABLE point3D_All
                AS SELECT   a.{1}, b.*
                FROM {2} AS a RIGHT JOIN {3} AS b
                            ON a.{0} = b.{0}
                WHERE b.{0} IS NOT NULL
            """.format( ID_POINT                    , GEOM_FIELD,
                        gridPoint                   , allZonesPointFactor,
                        DataUtil.createIndex(tableName=gridPoint, 
                                             fieldName=ID_POINT,
                                             isSpatial=False),
                        DataUtil.createIndex(tableName=allZonesPointFactor, 
                                             fieldName=ID_POINT,
                                             isSpatial=False)))
        saveData.saveTable(cursor = cursor,
                           tableName = "point3D_All",
                           filedir = outputDataAbs["point3D_All"]+".geojson",
                           delete = True,
                           rotationCenterCoordinates = rotationCenterCoordinates,
                           rotateAngle = - windDirection)    
    
    
    # -------------------------------------------------------------------
    # 8. 3D WIND SPEED INITIALIZATION -----------------------------------
    # -------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Initialize the 3D wind in the grid')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Identify 3D grid points intersected by buildings
    df_gridBuil = \
        InitWindField.identifyBuildPoints(cursor = cursor,
                                          gridPoint = gridPoint,
                                          stackedBlocksWithBaseHeight = rotatedPropStackedBlocks,
                                          dz = dz,
                                          tempoDirectory = tempoDirectory)
    
    # Set the initial 3D wind speed field
    df_wind0, nPoints, verticalWindProfile = \
        InitWindField.setInitialWindField(cursor = cursor, 
                                          initializedWindFactorTable = allZonesPointFactor,
                                          gridPoint = gridPoint,
                                          df_gridBuil = df_gridBuil,
                                          z0 = z0,
                                          sketchHeight = sketchHeight,
                                          profileType = profileType,
                                          meshSize = meshSize,
                                          dz = dz, 
                                          z_ref = z_ref,
                                          V_ref = v_ref, 
                                          tempoDirectory = tempoDirectory,
                                          d = d,
                                          H = Hr,
                                          lambda_f = lambda_f,
                                          verticalProfileFile = verticalProfileFile)
    
    # -------------------------------------------------------------------
    # 9. "RASTERIZE" THE DATA - PREPARE MATRICES FOR WIND CALCULATION ---
    # -------------------------------------------------------------------
    if feedback:
        feedback.setProgressText('Rasterize the data')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    # Set the ground as "building" (understand solid wall) - after getting grid size
    nx, ny, nz = nPoints.values()
    df_gridBuil = df_gridBuil.reindex(df_gridBuil.index.append(pd.MultiIndex.from_product([range(1,nx-1),
                                                                                          range(1,ny-1),
                                                                                          [0]])))

    # Set the buildGrid3D object to zero when a cell intersect a building 
    buildGrid3D = pd.Series(1, index = df_wind0.index, dtype = np.int32)
    buildGrid3D.loc[df_gridBuil.index] = 0
    
    # Convert building coordinates and wind speeds to numpy matrix...
    # (note that v axis direction is changed since we first use Röckle schemes
    # considering wind speed coming from North thus axis facing South)
    buildGrid3D = np.array([buildGrid3D.xs(i, level = 0).unstack().values for i in range(0,nx)])
    u0 = np.array([df_wind0[U].xs(i, level = 0).unstack().values for i in range(0,nx)])
    v0 = -np.array([df_wind0[V].xs(i, level = 0).unstack().values for i in range(0,nx)])
    w0 = np.array([df_wind0[W].xs(i, level = 0).unstack().values for i in range(0,nx)])
    
    # Identify all cells needing to be updated by the wind solver and store
    # their coordinates in a 1D array
    # (exclude buildings and sketch boundaries)
    cells4Solver = np.transpose(np.where(buildGrid3D == 1))
    cells4Solver = cells4Solver[cells4Solver[:, 0] > 0]
    cells4Solver = cells4Solver[cells4Solver[:, 1] > 0]
    cells4Solver = cells4Solver[cells4Solver[:, 2] > 0]
    cells4Solver = cells4Solver[cells4Solver[:, 0] < nx - 1]
    cells4Solver = cells4Solver[cells4Solver[:, 1] < ny - 1]
    cells4Solver = cells4Solver[cells4Solver[:, 2] < nz - 1]
    cells4Solver = cells4Solver.astype(np.int32)   
    
    # Identify building 3D coordinates
    buildingCoordinates = np.stack(np.where(buildGrid3D==0)).astype(np.int32)
    
    # Interpolation is made in order to have wind speed located on the face of
    # each grid cell
    u0[1:nx, :, :] =   (u0[0:nx-1, :, :] + u0[1:nx, :, :])/2
    v0[:, 1:ny, :] =   (v0[:, 0:ny-1, :] + v0[:,1:ny,:])/2
    w0[:, :, 1:nz] =   (w0[:, :, 0:nz-1] + w0[:, :, 1:nz])/2
    
    # Reset input and output wind speed to zero for building cells
    u0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    u0[buildingCoordinates[0]+1,buildingCoordinates[1],buildingCoordinates[2]]=0
    v0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    v0[buildingCoordinates[0],buildingCoordinates[1]+1,buildingCoordinates[2]]=0
    w0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    w0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]+1]=0
    
    # Create local grid space coordinates (x, y, z)
    Lz = (nz-1) * dz
    Lx = (nx-1) * meshSize
    Ly = (ny-1) * meshSize
    x = np.linspace(0, Lx, nx)  
    y = np.linspace(0, Ly, ny)
    z = np.linspace(0, Lz, nz)
    
    print("Time spent for wind speed initialization: {0} s".format(time.time()-timeStartCalculation))
    print("Shape: " + str(u0.shape) + " - " + "Nb cells: " + str(u0.shape[0] * u0.shape[1] * u0.shape[2]))
    # -------------------------------------------------------------------
    # 10. WIND SOLVER APPLICATION ----------------------------------------
    # ------------------------------------------------------------------- 
    if feedback:
        feedback.setProgressText('Apply the wind solver equations')
        if feedback.isCanceled():
            cursor.close()
            feedback.setProgressText("Calculation cancelled by user")
            return {}
    if not onlyInitialization:
        # Apply a mass-flow balance to have a more physical 3D wind speed field
        u, v, w = \
            WindSolver.solver(  x = x                       , y = y                 , z = z,
                                dx = meshSize               , dy = meshSize         , dz = dz,
                                u0 = u0                     , v0 = v0               , w0 = w0, cursor = cursor,
                                buildingCoordinates = buildingCoordinates   , cells4Solver = cells4Solver,
                                maxIterations = maxIterations, thresholdIterations = thresholdIterations,
                                feedback = feedback)
    else:
        u = u0
        v = v0
        w = w0
        
    # Wind speed values are recentered to the middle of the cells
    u[0:nx-1 ,0:ny-1 ,0:nz-1]=   (u[0:nx-1, 0:ny-1, 0:nz-1] + u[1:nx, 0:ny-1, 0:nz-1])/2
    v[0:nx-1 ,0:ny-1, 0:nz-1]=   (v[0:nx-1, 0:ny-1, 0:nz-1] + v[0:nx-1, 1:ny, 0:nz-1])/2
    w[0:nx-1, 0:ny-1, 0:nz-1]=   (w[0:nx-1, 0:ny-1, 0:nz-1] + w[0:nx-1, 0:ny-1, 1:nz])/2
    u0[0:nx-1 ,0:ny-1 ,0:nz-1]=   (u0[0:nx-1, 0:ny-1, 0:nz-1] + u0[1:nx, 0:ny-1, 0:nz-1])/2
    v0[0:nx-1 ,0:ny-1, 0:nz-1]=   (v0[0:nx-1, 0:ny-1, 0:nz-1] + v0[0:nx-1, 1:ny, 0:nz-1])/2
    w0[0:nx-1, 0:ny-1, 0:nz-1]=   (w0[0:nx-1, 0:ny-1, 0:nz-1] + w0[0:nx-1, 0:ny-1, 1:nz])/2
    
    # Reset input and output wind speed to zero for building cells
    u[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    v[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    w[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    u0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    v0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    w0[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    
    # -------------------------------------------------------------------
    # 11. ROTATE THE WIND FIELD TO THE INITIAL DISPOSITION --------------
    # ------------------------------------------------------------------- 
    # Get the relative position of the upper right corner of the grid from
    # the center of rotation used to rotate the grid
    cursor.execute(
        """{0};{1}
        """.format(DataUtil.createIndex(tableName=gridPoint, 
                                        fieldName=ID_POINT_X,
                                        isSpatial=False),
                    DataUtil.createIndex(tableName=gridPoint, 
                                         fieldName=ID_POINT_Y,
                                         isSpatial=False)))
    cursor.execute(
        """
        SELECT  {3}-ST_X(a.{0}) AS DIST_ROT_X,
                {4}-ST_Y(a.{0}) AS DIST_ROT_Y
        FROM {5} AS a
        WHERE   a.{1} = (SELECT MAX({1}) FROM {5})
                AND a.{2} = (SELECT MAX({2}) FROM {5})
        """.format(GEOM_FIELD                   , ID_POINT_X,
                   ID_POINT_Y                   , rotationCenterCoordinates[0],
                   rotationCenterCoordinates[1] , gridPoint))
    dist_rot_x, dist_rot_y = cursor.fetchall()[0]
    x += dist_rot_x
    y += dist_rot_y
    
    x_rot = np.zeros((nx, ny))
    y_rot = np.zeros((nx, ny))
    x_rot, y_rot, u_rot, v_rot = rotateData(theta = -windDirection*np.pi/180, nx = nx, 
                                            ny = ny                         , nz = nz, 
                                            x = x                           , y = y,
                                            x_rot = x_rot                   , y_rot = y_rot,
                                            u = u                           , v = v)

    x_rot, y_rot, u0_rot, v0_rot = rotateData(theta = -windDirection*np.pi/180  , nx = nx, 
                                              ny = ny                           , nz = nz, 
                                              x = x                             , y = y,
                                              x_rot = x_rot                     , y_rot = y_rot,
                                              u = u0                            , v = v0)
    # Set the real (x,y) grid coordinates
    x_rot += rotationCenterCoordinates[0]
    y_rot += rotationCenterCoordinates[1]
    
    # -------------------------------------------------------------------
    # 12. SAVE EACH OF THE UROCK OUTPUT ---------------------------------
    # ------------------------------------------------------------------- 
    # First rotate the coordinates of the grid of points
    rotated_grid = Obstacles.windRotation(cursor = cursor,
                                          dicOfInputTables = {gridPoint: gridPoint},
                                          rotateAngle = - windDirection,
                                          rotationCenterCoordinates = rotationCenterCoordinates)[0][gridPoint]
    
    dicVectorTables, netcdf_path =\
        saveData.saveBasicOutputs(cursor = cursor                , z_out = z_out,
                                  dz = dz                        , u = u_rot,
                                  v = v_rot                      , w = w, 
                                  gridName = rotated_grid        , verticalWindProfile = verticalWindProfile,
                                  outputFilePath = outputFilePath, outputFilename = outputFilename,
                                  meshSize = meshSize            , outputRaster = outputRaster,
                                  saveRaster = saveRaster        , saveVector = saveVector,
                                  saveNetcdf = saveNetcdf        , prefix_name = prefix,
                                  stacked_blocks = outputDataAbs["stacked_blocks"])
    
    # Save also the initialisation field if needed
    if debug:
        dicVectorTables_ini, netcdf_path_ini =\
            saveData.saveBasicOutputs(cursor = cursor                , z_out = z_out,
                                      dz = dz                        , u = u0_rot,
                                      v = v0_rot                     , w = w0, 
                                      gridName = rotated_grid        , verticalWindProfile = verticalWindProfile,
                                      outputFilePath = tempoDirectory, outputFilename = "wind_initiatlisation",
                                      meshSize = meshSize            , outputRaster = outputRaster,
                                      saveRaster = saveRaster        , saveVector = saveVector,
                                      saveNetcdf = saveNetcdf        , prefix_name = prefix,
                                      stacked_blocks = outputDataAbs["stacked_blocks"])  
    else:
        dicVectorTables_ini = None
        netcdf_path_ini = None

    # Last save the 2D grid for each Röckle zone
    saveData.saveRockleZones(cursor = cursor,
                             outputDataAbs = outputDataAbs,
                             dicOfBuildZoneGridPoint = dicOfBuildZoneGridPoint,
                             dicOfVegZoneGridPoint = dicOfVegZoneGridPoint,
                             gridPoint = gridPoint,
                             rotationCenterCoordinates = rotationCenterCoordinates, 
                             windDirection = windDirection)
    
    # Close the Database connection and remove the file
    if not debug:
        H2gisConnection.closeAndRemoveH2gisInstance(localH2InstanceDir = localH2InstanceDir, 
                                                    conn = conn,
                                                    cur = cursor)

    return  u_rot, v_rot, w, u0_rot, v0_rot, w0, x_rot, y_rot, z,\
            buildingCoordinates, cursor, rotated_grid, rotationCenterCoordinates,\
            verticalWindProfile, dicVectorTables, netcdf_path, netcdf_path_ini

@jit(nopython=True)
def rotateData(theta, nx, ny, nz, x, y, x_rot, y_rot, u, v):
    u_rot = np.zeros(u.shape)
    v_rot = np.zeros(v.shape)
    rot = np.array([[math.cos(theta), -math.sin(theta)],
                    [math.sin(theta), math.cos(theta)]])
    xmax = x.max()
    ymax = y.max()
    for i in range(nx):
        for j in range(ny):
            x_rot[i,j] = (xmax-x[i]) * rot[0,0] + (ymax-y[j]) * rot[0,1]
            y_rot[i,j] = (xmax-x[i]) * rot[1,0] + (ymax-y[j]) * rot[1,1]
            for k in range(nz):
                u_rot[i, j, k] = u[i,j,k] * rot[0,0] + v[i,j,k] * rot[0,1]
                v_rot[i, j, k] = u[i,j,k] * rot[1,0] + v[i,j,k] * rot[1,1]
    
    return x_rot, y_rot, u_rot, v_rot