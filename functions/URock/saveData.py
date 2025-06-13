#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 13:59:31 2021

@author: Jérémy Bernard, University of Gothenburg
"""
import pandas as pd
import geopandas as gpd
import shutil
import numpy as np
# from scipy.interpolate import griddata
# from rasterio.transform import from_origin
# import rasterio

from .DataUtil import radToDeg, windDirectionFromXY, createIndex, prefix
from .Obstacles import windRotation
from osgeo.osr import SpatialReference
from osgeo.gdal import Grid, GridOptions, FillNodata, Open, GA_Update, GetDriverByName
from .GlobalVariables import HORIZ_WIND_DIRECTION, HORIZ_WIND_SPEED, WIND_SPEED,\
    ID_POINT, TEMPO_DIRECTORY, TEMPO_HORIZ_WIND_FILE, VERT_WIND_SPEED, GEOM_FIELD,\
    OUTPUT_DIRECTORY, MESH_SIZE, OUTPUT_FILENAME, DELETE_OUTPUT_IF_EXISTS,\
    OUTPUT_RASTER_EXTENSION, OUTPUT_VECTOR_EXTENSION, OUTPUT_NETCDF_EXTENSION,\
    WIND_GROUP, WINDSPEED_PROFILE, RLON, RLAT, LON, LAT, LEVELS, WINDSPEED_X,\
    WINDSPEED_Y, WINDSPEED_Z, VERT_WIND, Z, OUTPUT_FILENAME, PREFIX_NAME,\
    BASE_HEIGHT_FIELD, HEIGHT_FIELD
from datetime import datetime
import netCDF4 as nc4
import os
import processing

def saveBasicOutputs(cursor, z_out, dz, u, v, w, gridName,
                     verticalWindProfile, outputFilePath, meshSize,
                     stacked_blocks, outputFilename = OUTPUT_FILENAME,
                     outputRaster = None, saveRaster = True,
                     saveVector = True, saveNetcdf = True,
                     prefix_name = PREFIX_NAME, tmp_dir = TEMPO_DIRECTORY):

    # Get the srid of the input geometry
    cursor.execute(""" SELECT ST_SRID({0}) AS srid FROM {1} LIMIT 1
                   """.format( GEOM_FIELD,
                               gridName))
    srid = cursor.fetchall()[0][0]
    
    # -------------------------------------------------------------------
    # SAVE NETCDF -------------------------------------------------------
    # ------------------------------------------------------------------- 
    final_netcdf_path = None
    if saveNetcdf:    
        # Get the coordinate in lat/lon of each point 
        # WARNING : for now keep the data in local coordinates)
        cursor.execute(""" 
           SELECT ST_X({0}) AS LON, ST_Y({0}) AS LAT FROM 
           (SELECT ST_TRANSFORM(ST_SETSRID({0},{2}), 4326) AS {0} FROM {1})
           """.format( GEOM_FIELD,
                       gridName,
                       srid))
        coord = np.array(cursor.fetchall())
        # Convert to a 2D (X, Y) array
        nx = u.shape[0]
        ny = u.shape[1]
        longitude = np.array([[coord[i * nx + j, 0] for i in range(ny)] for j in range(nx)])
        latitude = np.array([[coord[i * nx + j, 1] for i in range(ny)] for j in range(nx)])
        
    
        # Save the data into a NetCDF file
        # If delete = False, add a suffix to the file
        netcdf_base_dir_name = os.path.join(outputFilePath, 
                                            prefix(outputFilename, prefix_name))
        if os.path.isfile(netcdf_base_dir_name + OUTPUT_NETCDF_EXTENSION):
            if DELETE_OUTPUT_IF_EXISTS:
                os.remove(netcdf_base_dir_name + OUTPUT_NETCDF_EXTENSION)
            else:
                netcdf_base_dir_name = renameFileIfExists(filedir = netcdf_base_dir_name,
                                                          extension = OUTPUT_NETCDF_EXTENSION)  
        final_netcdf_path = saveToNetCDF(longitude = longitude,
                                         latitude = latitude,
                                         x = range(nx),
                                         y = range(ny),
                                         u = u,
                                         v = v,
                                         w = w,
                                         verticalWindProfile = verticalWindProfile,
                                         path = netcdf_base_dir_name,
                                         urock_srid = srid,
                                         horizontal_res = meshSize,
                                         vertical_res = dz)

    horizOutputUrock = {z_i : "HORIZ_OUTPUT_UROCK_{0}".format(str(z_i).replace(".","_")) for z_i in z_out}
    for z_i in z_out:
        # Keep only wind field for a single horizontal plan (and convert carthesian
        # wind speed into polar at least for horizontal)
        tempoTable = "TEMPO_HORIZ"
        if z_i % dz % (dz / 2) == 0:
            n_lev = int(z_i / dz) + 1
            ufin = u[:,:,n_lev]
            vfin = v[:,:,n_lev]
            wfin = w[:,:,n_lev]
        else:
            n_lev = int((z_i + dz /2) / dz)
            n_lev1 = n_lev + 1
            weight1 = (z_i - (n_lev - 0.5) * dz) / dz
            weight = 1 - weight1
            ufin = (weight * u[:,:,n_lev] + weight1 * u[:,:,n_lev1])
            vfin = (weight * v[:,:,n_lev] + weight1 * v[:,:,n_lev1])
            wfin = (weight * w[:,:,n_lev] + weight1 * w[:,:,n_lev1])
        df = pd.DataFrame({HORIZ_WIND_SPEED: ((ufin ** 2 + vfin ** 2) ** 0.5).flatten("F"),
                           WIND_SPEED: ((ufin ** 2 + vfin ** 2 + wfin ** 2) ** 0.5).flatten("F"), 
                           HORIZ_WIND_DIRECTION: radToDeg(windDirectionFromXY(ufin, vfin)).flatten("F"), 
                           VERT_WIND_SPEED: wfin.flatten("F")}).rename_axis(ID_POINT)
        
        # Load horizontal wind speed, wind direction and
        # vertical wind speed in a file containing a geometry field
        csv_tmp_file = os.path.join(tmp_dir, TEMPO_HORIZ_WIND_FILE)
        df.to_csv(csv_tmp_file)
        cursor.execute(
            """
            DROP TABLE IF EXISTS {9};
            CREATE TABLE {9}({3} INTEGER, {5} DOUBLE, {6} DOUBLE, {7} DOUBLE, {11} DOUBLE)
                AS SELECT {3}, {5}, {6}, {7}, {11} FROM CSVREAD('{10}');
            {0}{1}
            DROP TABLE IF EXISTS {2};
            CREATE TABLE {2}
                AS SELECT   a.{3}, {4}, b.{5}, 
                            b.{6}, b.{7}, b.{11}
                FROM {8} AS a
                LEFT JOIN {9} AS b
                ON a.{3} = b.{3}
            """.format(createIndex(tableName=gridName, 
                                            fieldName=ID_POINT,
                                            isSpatial=False),
                        createIndex(tableName=tempoTable, 
                                             fieldName=ID_POINT,
                                             isSpatial=False),
                        horizOutputUrock[z_i]       , ID_POINT,
                        GEOM_FIELD                  , HORIZ_WIND_SPEED,
                        HORIZ_WIND_DIRECTION        , VERT_WIND_SPEED,
                        gridName                    , tempoTable,
                        csv_tmp_file,
                        WIND_SPEED))
        
        # -------------------------------------------------------------------
        # SAVE VECTOR -------------------------------------------------------
        # ------------------------------------------------------------------- 
        if saveVector or saveRaster:
            outputDir_zi = os.path.join(outputFilePath, 
                                        "z" + str(z_i).replace(".","_"))
            if not os.path.exists(outputDir_zi):
                os.mkdir(outputDir_zi)
            outputVectorFile = saveTable(cursor = cursor,
                                         tableName = horizOutputUrock[z_i],
                                         filedir = os.path.join(outputDir_zi,
                                                                prefix(outputFilename, prefix_name)+\
                                                                OUTPUT_VECTOR_EXTENSION),
                                         delete = DELETE_OUTPUT_IF_EXISTS)
                
            # -------------------------------------------------------------------
            # SAVE RASTER -------------------------------------------------------
            # -------------------------------------------------------------------     
            if saveRaster:
                # Save the all direction, horizontal and vertical wind speeds into a a different raster
                for var in [WIND_SPEED, HORIZ_WIND_SPEED, VERT_WIND_SPEED]:
                    saveRasterFile(cursor = cursor, 
                                   outputVectorFile = outputVectorFile,
                                   outputFilePathAndNameBase = os.path.join(outputDir_zi,
                                                                            prefix(outputFilename, prefix_name)),
                                   horizOutputUrock = horizOutputUrock,
                                   outputRaster = outputRaster, 
                                   z_i = z_i, 
                                   meshSize = meshSize,
                                   var2save = var,
                                   stacked_blocks = stacked_blocks,
                                   srid = srid,
                                   tmp_dir = tmp_dir)   

    return horizOutputUrock, final_netcdf_path
    
def saveToNetCDF(longitude,
                 latitude,
                 x,
                 y,
                 u,
                 v,
                 w,
                 verticalWindProfile,
                 path,
                 urock_srid,
                 horizontal_res,
                 vertical_res):
    """
    Create a netCDF file and save wind speed, direction and initial 
    vertical wind profile in it (based on https://pyhogs.github.io/intro_netcdf4.html )
    
    Parameters
    _ _ _ _ _ _ _ _ _ _ 
        longitude: np.array (2D - X, Y)
            Longitude of each of the (X, Y) points
        latitude: np.array (2D - X, Y)
            Longitude of each of the (X, Y) points
        x: np.array (1D)
            X grid coordinates in local referential
        y: np.array (1D)
            Y grid coordinates in local referential
        u: np.array (3D)
            Wind speed along East axis
        v: 2D (X, Y) array
            Wind speed along North axis
        w: 2D (X, Y) array
            Wind speed along vertical axis
        verticalWindProfile: pd.DataFrame
            Initial wind speed profile for each each z from ground (2 columns)
        path: String
            Path and filename to save NetCDF file
        urock_srid: int
            EPSG code initially used for the URock calculations
    
    Returns
    -------
        String being the path, filename and extension of the netCdf file where
        are stored the results
    """    
    # Opens a netCDF file in writing mode ('w')
    f = nc4.Dataset(path + OUTPUT_NETCDF_EXTENSION,'w', format='NETCDF4')
    
    # 3D WIND SPEED DATA
    # Creates a group within this file for the 3D wind speed
    wind3dGrp = f.createGroup(WIND_GROUP)
    
    # Creates dimensions within this group
    wind3dGrp.createDimension('rlon', len(x))
    wind3dGrp.createDimension('rlat', len(y))
    wind3dGrp.createDimension('z', verticalWindProfile.index.size)
    
    # Build the variables
    rlon = wind3dGrp.createVariable(RLON, 'i4', 'rlon')
    rlat = wind3dGrp.createVariable(RLAT, 'i4', 'rlat')
    z = wind3dGrp.createVariable(Z, 'f4', 'z')
    lon = wind3dGrp.createVariable(LON, 'f8', ('rlon', 'rlat'))
    lat = wind3dGrp.createVariable(LAT, 'f8', ('rlon', 'rlat'))
    windSpeed_x = wind3dGrp.createVariable(WINDSPEED_X, 'f4', ('rlon', 'rlat', 'z'))
    windSpeed_y = wind3dGrp.createVariable(WINDSPEED_Y, 'f4', ('rlon', 'rlat', 'z'))  
    windSpeed_z = wind3dGrp.createVariable(WINDSPEED_Z, 'f4', ('rlon', 'rlat', 'z'))
    
    # Fill the variables
    rlon[:] = x
    rlat[:] = y
    z[:] = verticalWindProfile[Z].values
    lon[:,:] = longitude
    lat[:,:] = latitude
    windSpeed_x[:,:,:] = u
    windSpeed_y[:,:,:] = v
    windSpeed_z[:,:,:] = w
    
    # VERTICAL WIND PROFILE DATA
    # Creates a group within this file for the vertical wind profile
    vertWindProfGrp = f.createGroup(VERT_WIND)
    
    # Creates dimensions within this group
    vertWindProfGrp.createDimension('z', verticalWindProfile.index.size)
    
    # Build the variables 
    z_profile = vertWindProfGrp.createVariable(Z, 'f4', 'z')
    WindSpeed = vertWindProfGrp.createVariable(WINDSPEED_PROFILE, 'f4', 'z')
    
    # Fill the variables
    z_profile[:] = verticalWindProfile[Z].values
    WindSpeed[:] = verticalWindProfile[HORIZ_WIND_SPEED].values
    
    
    # ADD METADATA
    #Add local attributes to variable instances
    lon.units = 'degrees east'
    lat.units = 'degrees north'
    windSpeed_x.units = 'meter per second'
    windSpeed_y.units = 'meter per second'
    windSpeed_z.units = 'meter per second'
    z.units = 'meters'
    WindSpeed.units = 'meter per second'
    z_profile.units = 'meters'

    #Add global attributes
    f.description = "URock dataset containing one group of 3D wind field value and one group of input vertical wind speed profile"
    f.history = "Created " + datetime.today().strftime("%y-%m-%d")
    
    # Add the srid (epsg code) used for the URock processing calculation
    f.urock_srid = urock_srid
    
    # Add horizontal and vertical resolution into the metadata
    f.horizontal_res = horizontal_res
    f.vertical_res = vertical_res
    
    f.close()
    
    return path + OUTPUT_NETCDF_EXTENSION
    
def saveTable(cursor, tableName, filedir, delete = False, 
              rotationCenterCoordinates = None, rotateAngle = None):
    """ Save a table in .geojson or .shp (the table can be rotated before saving if needed).
    
    Parameters
	_ _ _ _ _ _ _ _ _ _ 
        cursor: conn.cursor
            A cursor object, used to perform spatial SQL queries
		tableName : String
			Name of the table to save
        filedir: String
            Directory (including filename and extension) of the file where to 
            store the table
        delete: Boolean, default False
            Whether or not the file is delete if exist
        rotationCenterCoordinates: tuple of float, default None
            x and y values of the point used as center of rotation
        rotateAngle: float, default None
            Counter clock-wise rotation angle (in degree)

    
    Returns
	_ _ _ _ _ _ _ _ _ _ 	
		output_filedir: String
            Directory (including filename and extension) of the saved file
            (could be different from input 'filedir' since the file may 
             have been renamed if exists)"""
    # Rotate the table if needed
    if rotationCenterCoordinates is not None and rotateAngle is not None:
        tableName = windRotation(cursor = cursor,
                                 dicOfInputTables = {tableName: tableName},
                                 rotateAngle = rotateAngle,
                                 rotationCenterCoordinates = rotationCenterCoordinates)[0][tableName]
    
    # Get extension
    extension = "." + filedir.split(".")[-1]
    filedirWithoutExt = ".".join(filedir.split(".")[0:-1])
    
    # Define the H2GIS function depending on extension
    if extension.upper() == ".GEOJSON":
        h2_function = "GEOJSONWRITE"
    elif extension.upper() == ".SHP":
        h2_function = "SHPWRITE"
    elif extension.upper() == ".FGB":
        h2_function = "FGBWRITE"
    else:
        print("The extension should be .geojson, .shp or .fgb")
    # Delete files if exists and delete = True
    if delete and os.path.isfile(filedir):
        output_filedir = filedir
        os.remove(filedir)
        if extension.upper() == ".SHP":
            os.remove(filedirWithoutExt+".dbf")
            os.remove(filedirWithoutExt+".shx")
            if os.path.isfile(filedirWithoutExt+".prj"):
                os.remove(filedirWithoutExt+".prj")
    # If delete = False, add a suffix to the file
    elif os.path.isfile(filedir):
        output_filedir = renameFileIfExists(filedir = filedirWithoutExt,
                                            extension = extension) + extension
    else:
        output_filedir = filedir
    # Write files
    cursor.execute("""CALL {0}('{1}','{2}')""".format(h2_function,
                                                      output_filedir,
                                                      tableName))
    return output_filedir

def renameFileIfExists(filedir, extension):
    """ Rename a file with a numbering prefix if exists.
    
    Parameters
	_ _ _ _ _ _ _ _ _ _ 
        filedir: String
            Directory (including filename but without extension) of the file
    
    Returns
	_ _ _ _ _ _ _ _ _ _ 	
		newFileDir: String
            Directory with renamed file"""
    i = 1
    newFileDir = filedir
    while(os.path.isfile(newFileDir + extension)):
        newFileDir = filedir + "({0})".format(i)
        i += 1
    return newFileDir


def saveRasterFile(cursor, outputVectorFile, outputFilePathAndNameBase, 
                   horizOutputUrock, outputRaster, z_i, meshSize, var2save,
                   stacked_blocks, srid, tmp_dir):
    """ Save results in a raster file.
    
    Parameters
	_ _ _ _ _ _ _ _ _ _ 
        outputFilePathAndNameBase: String
            Directory (including filename but without extension) of the file
    
    Returns
	_ _ _ _ _ _ _ _ _ _ 	
		None"""
    # Define output path name
    outputFilePathAndNameBaseRaster = outputFilePathAndNameBase + var2save
    # If delete = False, add a suffix to the filename
    if (os.path.isfile(outputFilePathAndNameBaseRaster + OUTPUT_RASTER_EXTENSION)) \
        and (not DELETE_OUTPUT_IF_EXISTS):
        outputFilePathAndNameBaseRaster = renameFileIfExists(filedir = outputFilePathAndNameBaseRaster,
                                                             extension = OUTPUT_RASTER_EXTENSION)
    
    # Whether or not a raster output is given as input, the rasterization process is slightly different
    if outputRaster:
        outputRasterExtent = outputRaster.extent()
        resX = (outputRasterExtent.xMaximum() - outputRasterExtent.xMinimum()) / outputRaster.width()
        resY = (outputRasterExtent.yMaximum() - outputRasterExtent.yMinimum()) / outputRaster.height()
        xmin = outputRasterExtent.xMinimum()
        ymax = outputRasterExtent.yMaximum()
        xmax = outputRasterExtent.xMaximum()
        ymin = outputRasterExtent.yMinimum()
        tmp_file = os.path.join(tmp_dir, f"interp_before_fillna_{var2save}.tif")
        # If a single output raster cell contains more than 4 points, average instead of interpolate
        if resX * resY > 4 * meshSize**2:
            Grid(destName = tmp_file,
                 srcDS = outputVectorFile,
                 options = GridOptions(format = OUTPUT_RASTER_EXTENSION.split(".")[-1],
                                       zfield = var2save, 
                                       width = outputRaster.width(), 
                                       height = outputRaster.height(),
                                       outputBounds = [xmin,
                                                       ymax,
                                                       xmax,
                                                       ymin],
                                       algorithm = "average:radius1={0}:radius2={0}".format(1.1*meshSize)))
        else:
            # Interpolate with building constraints
            interp_vec_to_rast(outputVectorFile = outputVectorFile, 
                               stacked_blocks = stacked_blocks,
                               outputFilePathAndNameBaseRaster = ".".join(tmp_file.split(".")[0:-1]), 
                               extent = f'{xmin},{xmax},{ymin},{ymax} [EPSG:{srid}]',
                               resX = resX, 
                               resY = resY,
                               z_i = z_i,
                               colname = var2save,
                               tmp_dir = tmp_dir)
            
        # Interpolate values to fill no data
        ds = Open(tmp_file) 
        driver = GetDriverByName("GTiff")
        output_ds = driver.CreateCopy(outputFilePathAndNameBaseRaster + OUTPUT_RASTER_EXTENSION,
                                      ds)
        band = output_ds.GetRasterBand(1)
        _ = FillNodata(targetBand = band, maskBand = None, 
                       maxSearchDist = 9999, smoothingIterations = 0)
        
        # Set the srid
        vec_srid = gpd.read_file(outputVectorFile, rows = slice(0,)).crs.to_epsg()
        srs = SpatialReference()
        srs.ImportFromEPSG(vec_srid)
        output_ds.SetProjection(srs.ExportToWkt())
        
        # Release the datasets.
        ds = ds.FlushCache()
        output_ds = output_ds.FlushCache()
        
        # ET = Open(outputFilePathAndNameBaseRaster + OUTPUT_RASTER_EXTENSION, GA_Update) 
        # ETband = ET.GetRasterBand(1)
 
        # FillNodata(targetBand = ETband, maskBand = None, 
        #            maxSearchDist = 999, smoothingIterations = 0)
        # ET = None
    else:
        cursor.execute(
            """
            SELECT  ST_XMIN({0}) AS XMIN, ST_XMAX({0}) AS XMAX,
                    ST_YMIN({0}) AS YMIN, ST_YMAX({0}) AS YMAX
            FROM    (SELECT ST_ACCUM({0}) AS {0} FROM {1})
            """.format(GEOM_FIELD            , horizOutputUrock[z_i]))
        vectorBounds = cursor.fetchall()[0]
        width = int((vectorBounds[1] - vectorBounds[0]) / meshSize) + 1
        height = int((vectorBounds[3] - vectorBounds[2]) / meshSize) + 1
        xmin = vectorBounds[0] - float(meshSize) / 2
        ymax = vectorBounds[3] + float(meshSize) / 2
        xmax = vectorBounds[0] + meshSize * (width - 0.5)
        ymin = vectorBounds[3] - meshSize * (height + 0.5)
        
        # Interpolate with building constraints
        interp_vec_to_rast(outputVectorFile = outputVectorFile, 
                           stacked_blocks = stacked_blocks,
                           outputFilePathAndNameBaseRaster = outputFilePathAndNameBaseRaster, 
                           extent = f'{xmin},{xmax},{ymin},{ymax} [EPSG:{srid}]',
                           resX = meshSize, 
                           resY = meshSize,
                           z_i = z_i,
                           colname = var2save,
                           tmp_dir = tmp_dir)

def interp_vec_to_rast(outputVectorFile, stacked_blocks, outputFilePathAndNameBaseRaster, 
                       extent, resX, resY, z_i, colname, tmp_dir):
    """ Interpolate and save wind data saved in a vector to a raster.
       
    Parameters
   	_ _ _ _ _ _ _ _ _ _ 
           outputVectorFile: String
               Directory (including filename and extension) of the file containing the 
               wind field to interpolate from vector to raster
           stacked_blocks: String
               Directory (including filename only - no extension) of the file containing stacked blocks 
               (building footprint + where it starts and ends vertically)
           outputFilePathAndNameBaseRaster: String
               Directory (including filename and extension) of the file that will 
               be used to save the results
           extent: String
               QGIS extent in the shape '{xmin},{xmax},{ymin},{ymax} [EPSG:epsgcode]'
           resX: float
               Resolution of the output raster along the X axis
           resY: float
               Resolution of the output raster along the Y axis
           z_i: float
               Height of the output wind field
           colname: string
               Name of the attribute column in teh vector table to convert to raster
           tmp_dir: String
               Path of the directory used for saving temporary results (should be unique)
           
       
   Returns
   	_ _ _ _ _ _ _ _ _ _ 	
           output_file_path: String
               Path and name of the file containing the resulting raster"""
    # gdf = gpd.read_file(outputVectorFile)
    
    # coordinates = np.array([point.coords[0] for point in gdf.geometry])  # [x, y]
    # values = gdf[colname].values  # The field 'V' to interpolate
    
    # # Define the raster grid (3m resolution)
    # resolution = min([resX, resY])
    # xmin, ymin, xmax, ymax = gdf.total_bounds
    # xmin -= resolution / 2
    # ymin -= resolution / 2
    # xmax -= resolution / 2
    # ymax += resolution / 2
    
    # # Create a grid of coordinates for interpolation
    # x_grid = np.arange(xmin, xmax, resolution)
    # y_grid = np.arange(ymin, ymax, resolution)
    # x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    
    # # Perform bilinear interpolation on the grid
    # grid_values = griddata(coordinates, values, (x_grid, y_grid), method='linear')
    # grid_values = grid_values[::-1, :]
    
    # # Define raster metadata
    # transform = from_origin(xmin, ymax, resolution, resolution)  # origin is top-left
    
    # # Save the interpolated grid to a raster file
    # interp_out = os.path.join(TEMPO_DIRECTORY,"interp_out.tif")
    # with rasterio.open(interp_out, 'w', driver='GTiff', 
    #                     height=grid_values.shape[0], width=grid_values.shape[1],
    #                     count=1, dtype=grid_values.dtype, crs=gdf.crs, transform=transform) as dst:
    #     dst.write(grid_values, 1)
    
    # Change the order of the points to make the TIN interpolation faster and working for all conditions
    order_changed = processing.run("native:orderbyexpression", 
                                   {'INPUT':outputVectorFile,
                                    'EXPRESSION':'randf(0,1)',
                                    'ASCENDING':False,
                                    'NULLS_FIRST':False,
                                    'OUTPUT':os.path.join(tmp_dir,
                                                          f"order_changed_{colname}")})["OUTPUT"]
    
    # Get the column number corresponding to the column name
    colnb = gpd.read_file(order_changed, rows = slice(0,)).columns.get_loc(colname) + 1
    
    # Interpolate the results without constraints
    interp_out = processing.run("qgis:tininterpolation",
                               {'INTERPOLATION_DATA':f'{order_changed}::~::0::~::{colnb}::~::0',
                                'METHOD':0,
                                'EXTENT':extent,
                                'PIXEL_SIZE':min(resX, resY),
                                'OUTPUT':os.path.join(tmp_dir,
                                                      f"interp_out_{colname}.tif")})["OUTPUT"]
    
    # If there are buildings in the study area, need to set wind speed = 0
    if os.stat(stacked_blocks).st_size > 0:
        # Rasterize the stacked blocks keeping the value of each stacked block base 
        block_base = processing.run("gdal:rasterize",
                                    {'INPUT':stacked_blocks,
                                     'FIELD':BASE_HEIGHT_FIELD,'BURN':0,'USE_Z':False,
                                     'UNITS':1,'WIDTH':resX,'HEIGHT':resY,
                                     'EXTENT':extent,
                                     'NODATA':None,'OPTIONS':'','DATA_TYPE':5,
                                     'INIT':-9999,'INVERT':False,'EXTRA':'',
                                     'OUTPUT':os.path.join(tmp_dir,
                                                           f"block_base_{colname}.tif")})["OUTPUT"]
        
        # Rasterize the stacked blocks keeping the value of each stacked block top 
        block_top = processing.run("gdal:rasterize",
                                   {'INPUT':stacked_blocks,
                                    'FIELD':HEIGHT_FIELD,'BURN':0,'USE_Z':False,
                                    'UNITS':1,'WIDTH':resX,'HEIGHT':resY,
                                    'EXTENT':extent,
                                    'NODATA':None,'OPTIONS':'','DATA_TYPE':5,
                                    'INIT':-9999,'INVERT':False,'EXTRA':'',
                                    'OUTPUT':os.path.join(tmp_dir,
                                                          f"block_top_{colname}.tif")})["OUTPUT"]
        
        # Keep the values only when there is no building at this position
        output_file_path = processing.run("gdal:rastercalculator",
                                          {'INPUT_A':block_base,'BAND_A':1,
                                           'INPUT_B':block_top,'BAND_B':1,
                                           'INPUT_C':interp_out,'BAND_C':1,
                                           'INPUT_D':None,'BAND_D':None,
                                           'INPUT_E':None,'BAND_E':None,
                                           'INPUT_F':None,'BAND_F':None,
                                           'FORMULA':f'((A == -9999) + (A < {z_i})) * ((B == -9999) + (B < {z_i})) * C',
                                           'NO_DATA':None,'EXTENT_OPT':0,'PROJWIN':None,
                                           'RTYPE':5,'OPTIONS':'','EXTRA':'',
                                           'OUTPUT':outputFilePathAndNameBaseRaster + OUTPUT_RASTER_EXTENSION})["OUTPUT"]   
    # Else directly save the result of the interpolation
    else:
        output_file_path = outputFilePathAndNameBaseRaster + OUTPUT_RASTER_EXTENSION
        shutil.copy2(src = interp_out, dst = output_file_path)
    
    return output_file_path
        
def saveRockleZones(cursor, outputDataAbs, dicOfBuildZoneGridPoint, dicOfVegZoneGridPoint,
                    gridPoint, rotationCenterCoordinates, windDirection):
    """ Save the 2D Röckle zones (building and vegetation) as points.
    
    Parameters
	_ _ _ _ _ _ _ _ _ _ 
        cursor: conn.cursor
            A cursor object, used to perform spatial SQL queries
		outputDataAbs : Dictionary
			Object containing the absolute path where should be saved the Röckle points
        dicOfBuildZoneGridPoint: Dictionary
            Dictionary containing all building table names to be saved
        dicOfVegZoneGridPoint: Dictionary
            Dictionary containing all vegetation table names to be saved
        gridPoint: String
            Grid table name
        rotationCenterCoordinates: tuple of float, default None
            x and y values of the point used as center of rotation
        windDirection: float, default None
            Wind direction used for calculation (° clock-wise from North)
        
    
    Returns
	_ _ _ _ _ _ _ _ _ _ 	
		None"""
    # Creates a folder if not exist
    if not os.path.exists(outputDataAbs["point_2DRockleZone"]):
        os.mkdir(outputDataAbs["point_2DRockleZone"])
    # Save Building Röckle zones
    for t in dicOfBuildZoneGridPoint:
        cursor.execute("""
           DROP TABLE IF EXISTS point_Buildzone_{0};
           {5};
           {6};
           CREATE TABLE point_Buildzone_{0}
               AS SELECT   a.{2}, b.*
               FROM {3} AS a RIGHT JOIN {4} AS b
                   ON a.{1} = b.{1}
               WHERE b.{1} IS NOT NULL
           """.format( t                            , ID_POINT, 
                       GEOM_FIELD                   , gridPoint, 
                       dicOfBuildZoneGridPoint[t]   , createIndex(tableName=gridPoint, 
                                                                  fieldName=ID_POINT,
                                                                  isSpatial=False),
                       createIndex(tableName=dicOfBuildZoneGridPoint[t], 
                                   fieldName=ID_POINT,
                                   isSpatial=False)))
        saveTable(cursor = cursor,
                  tableName = "point_Buildzone_"+t,
                  filedir = os.path.join(outputDataAbs["point_2DRockleZone"], t+".geojson"),
                  delete = True,
                  rotationCenterCoordinates = rotationCenterCoordinates,
                  rotateAngle = - windDirection)
    
    # Save vegetation Röckle zones
    for t in dicOfVegZoneGridPoint:
        saveTable(cursor = cursor,
                  tableName = dicOfVegZoneGridPoint[t],
                  filedir = os.path.join(outputDataAbs["point_2DRockleZone"], t+".geojson"),
                  delete = True,
                  rotationCenterCoordinates = rotationCenterCoordinates,
                  rotateAngle = - windDirection)