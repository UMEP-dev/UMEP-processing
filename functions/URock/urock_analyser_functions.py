#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 13:09:33 2022

@author: Jérémy Bernard, University of Gothenburg (Sweden)
"""

import xarray as xr
import os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
from pathlib import Path

from . import H2gisConnection
from .loadData import loadFile
from . import DataUtil
from .GlobalVariables import WIND_GROUP, TEMPO_DIRECTORY,\
    RLON, RLAT, GEOM_FIELD, LON , LAT, WINDSPEED_X, WINDSPEED_Y, WINDSPEED_Z,\
    Z, HORIZ_WIND_SPEED, WIND_SPEED

idx = pd.IndexSlice

# Properties for plotting
STREAM_DENSITY = 2
HEAD_WIDTH = 3
HEAD_LENGTH = 1.5
HEAD_AXIS_LENGTH = 1.5
WIDTH = 0.2

def plotSectionalViews(pluginDirectory, inputWindFile, lines_file='', srid_lines=None,
                       idLines='', isStream = False, savePlot = False,
                       polygons_file='', srid_polygons=None, idPolygons='', 
                       outputDirectory = None, simulationName = "",
                       fig = None, ax = None, scale = None, color = None,
                       feedback = None):

    if savePlot:
        plt.ioff()
    # Set the SRID for work (by default use the srid_lines if exists)
    if srid_polygons and not srid_lines:
        srid_all = srid_polygons
    else:
        srid_all = srid_lines
        
    # Create temporary table names (for tables that will be removed at the end of the process)
    allPointsTab = DataUtil.postfix("ALL_POINTS")
    linesTab = DataUtil.postfix("LINES")
    pointsIntersecTab = DataUtil.postfix("POINTS_INTERSECT")
    polygonsTab = DataUtil.postfix("POLYGONS")
    polygonsMeanTab = DataUtil.postfix("POLYGONS_MEAN")
    
    # Temporary files are declared
    pointsDir = os.path.join(TEMPO_DIRECTORY, "urock_allPoints.csv")
    outputPointsDir = os.path.join(TEMPO_DIRECTORY, "urock_selectedPoints.csv")
    outputPolygonsDir = os.path.join(TEMPO_DIRECTORY, "urock_selectedPolygons.csv")
    
    if feedback:
        feedback.setProgressText('Load NetCDF file in Python and save as csv file...')
        
    # Load the group of the NetCDF file containing the wind speed field
    ds = xr.open_dataset(inputWindFile, group = WIND_GROUP)
    
    # Get the SRID that has been used in the URock processing calculation
    urock_srid = xr.open_dataset(inputWindFile).urock_srid    
    
    # Send to a csv file
    ds.to_dataframe().to_csv(pointsDir, index_label = ['rlat', 'rlon', 'zlev'])
    
    if feedback:
        feedback.setProgressText('Load csv file into H2GIS Database...')    
    # Initialize an H2GIS database connection
    dBDir = os.path.join(Path(pluginDirectory).parent, 'functions','URock')
    cursor = H2gisConnection.startH2gisInstance(dbDirectory = dBDir,
                                                dbInstanceDir = TEMPO_DIRECTORY)
    
    # Load coordinates in a H2GIS table
    cursor.execute("""
       DROP TABLE IF EXISTS {0};
       CREATE TABLE {0}(ID_POINT SERIAL,
                        {3} INTEGER,
                        {4} INTEGER,
                        {8} DOUBLE,
                        {9} DOUBLE,
                        {10} DOUBLE,
                        {11} DOUBLE,
                        {5} GEOMETRY) AS
            SELECT  NULL, {3}, {4}, {8}, {9}, {10}, {11},
                    ST_TRANSFORM(ST_SETSRID(ST_MakePoint({6}, {7}), 4326), {1}) AS {5}
            FROM CSVREAD('{2}')
        """.format(allPointsTab             , urock_srid, 
                    pointsDir               , RLON,
                    RLAT                    , GEOM_FIELD,
                    LON                     , LAT,
                    Z                       , WINDSPEED_X,
                    WINDSPEED_Y             , WINDSPEED_Z))
    
    # DEAL WITH POLYGON (MEAN WIND PROFILES)
    fig_poly = None
    ax_poly = None
    if polygons_file and srid_polygons and idPolygons:
        if feedback:
            feedback.setProgressText('Calculates average wind profile (within polygons)...')    
        # Load polygons
        loadFile(cursor = cursor, 
                 filePath = polygons_file, 
                 tableName = polygonsTab, 
                 srid = srid_polygons,
                 srid_repro = urock_srid)    
        
        # Calculates horizontal mean wind speed within each polygon
        cursor.execute("""
           {0}{1}{2}{3}
           DROP TABLE IF EXISTS {4};
           CREATE TABLE {4}
               AS SELECT b.{5},
                         a.{6},
                         AVG(a.{7}) AS {7},
                         AVG(a.{8}) AS {8},
                         AVG(a.{9}) AS {9},
                         AVG(POWER(POWER(a.{7},2) + POWER(a.{8},2), 0.5)) AS {10},
                         AVG(POWER(POWER(a.{7},2) + POWER(a.{8},2) + POWER(a.{9},2), 0.5)) AS {11}
               FROM {12} AS a, {13} AS b
               WHERE    a.{14} && b.{14} AND ST_INTERSECTS(a.{14}, b.{14}) AND
                        a.{7} <> 0 AND a.{8} <> 0 AND a.{9} <> 0
               GROUP BY a.{6}, b.{5};
           CALL CSVWrite('{15}', 'SELECT * FROM {4}');
           """.format(  DataUtil.createIndex(tableName=allPointsTab, 
                                             fieldName=GEOM_FIELD,
                                             isSpatial=True),
                        DataUtil.createIndex(tableName=polygonsTab, 
                                             fieldName=GEOM_FIELD,
                                             isSpatial=True),
                        DataUtil.createIndex(tableName=allPointsTab, 
                                             fieldName=Z,
                                             isSpatial=False),
                        DataUtil.createIndex(tableName=polygonsTab, 
                                             fieldName=idPolygons,
                                             isSpatial=False),
                        polygonsMeanTab         , idPolygons,
                        Z                       , WINDSPEED_X,
                        WINDSPEED_Y             , WINDSPEED_Z,
                        HORIZ_WIND_SPEED        , WIND_SPEED,
                        allPointsTab            , polygonsTab,
                        GEOM_FIELD              , outputPolygonsDir))
    
        # Back to dataframe for the plotting the mean wind speed for each level
        df_selectedPolygons = pd.read_csv(outputPolygonsDir, 
                                          index_col = None, header = 0)
        windList = pd.Series(["Wind speed along x-axis (m/s)",
                              "Wind speed along y-axis (m/s)",
                              "Wind speed along z-axis (m/s)",
                              "Horizontal wind speed (m/s)",
                              "Wind speed (m/s)"],
                             index = [WINDSPEED_X, WINDSPEED_Y, WINDSPEED_Z, 
                                      HORIZ_WIND_SPEED, WIND_SPEED])
        polygonsList = df_selectedPolygons[idPolygons.upper()].unique()
        fig_poly = {}
        ax_poly = {}
        for w in windList.index:
            fig_poly[w], ax_poly[w] = plt.subplots()
            for p in polygonsList:
                data2plot = df_selectedPolygons[df_selectedPolygons[idPolygons.upper()] == p].sort_values(Z)
                ax_poly[w].plot(data2plot[w.upper()], data2plot[Z], label = "Polygon {0}".format(p))
                ax_poly[w].set_xlabel(windList[w]),
                ax_poly[w].set_ylabel("Height above ground (m)")
                plt.legend()
            if savePlot:
                fig_poly[w].savefig(os.path.join(outputDirectory,
                                                 simulationName + "_" + w + ".png"))
            
    
    # DEAL WITH LINES (ALONG LINE VERTICAL PROFILES)
    fig = None
    ax = None
    scale = None
    if lines_file and srid_lines and idLines:
        if feedback:
            feedback.setProgressText('Calculates vertical sectional plot (along lines)...')    
        # Get the resolution of the wind speed data
        cursor.execute("""
           SELECT ST_DISTANCE(a.{0}, b.{0}) AS dist 
           FROM {1} AS a, {1} AS b 
           WHERE a.{2} = 0 AND a.{3} = 0 AND b.{2} = 1 AND b.{3} = 0
           """.format(GEOM_FIELD            , allPointsTab,
                       RLON                 , RLAT))
        horiz_res = round(cursor.fetchall()[0][0])
        dist_max = horiz_res * (2 ** 0.5) / 2
    
        # Load lines
        loadFile(cursor = cursor, 
                 filePath = lines_file, 
                 tableName = linesTab, 
                 srid = srid_lines,
                 srid_repro = urock_srid)
    
        # Calculates a buffer of about half the horizontal resolution of the
        # wind data around the lines and project the points contained in this buffer
        # on the lines and calculate the distance to this point to the begining of the line
        # NOTE : ONLY LINES HAVING TWO POINTS ARE USED (SEGMENTS) 
        cursor.execute("""
           {7}{8}
           DROP TABLE IF EXISTS {0};
           CREATE TABLE {0}
               AS SELECT a.ID_POINT,
                         a.{9},
                         a.{10},
                         a.{11},
                         a.{12},
                         ST_DISTANCE(ST_PROJECTPOINT(a.{1}, b.{1}), ST_STARTPOINT(b.{1})) AS DIST,
                         ST_AZIMUTH(ST_STARTPOINT(b.{1}), ST_ENDPOINT(b.{1})) AS AZIMUTH,
                         b.{2}
               FROM {3} AS a, {4} AS b
               WHERE ST_NPOINTS(b.{1}) = 2 AND a.{1} && ST_EXPAND(b.{1}, {5}) AND ST_DWITHIN(a.{1}, b.{1}, {5});
           CALL CSVWrite('{6}', 'SELECT * FROM {0}');
           """.format(pointsIntersecTab         , GEOM_FIELD,
                       idLines                  , allPointsTab,
                       linesTab                 , dist_max,
                       outputPointsDir          , DataUtil.createIndex(tableName=allPointsTab, 
                                                                       fieldName=GEOM_FIELD,
                                                                       isSpatial=True),
                       DataUtil.createIndex(tableName=linesTab, 
                                            fieldName=GEOM_FIELD,
                                            isSpatial=True),
                       Z                        , WINDSPEED_X,
                       WINDSPEED_Y              , WINDSPEED_Z))
    
        # Back to dataframe for the plotting of the wind speed for each level
        df_selectedPoints = pd.read_csv(outputPointsDir, index_col = None, header = 0)
        df_selectedPoints["PROJECTED_HORIZ_WIND"] = \
            df_selectedPoints[WINDSPEED_X.upper()] * np.cos(df_selectedPoints["AZIMUTH"] - np.pi / 2) +\
                df_selectedPoints[WINDSPEED_Y.upper()] * np.cos(df_selectedPoints["AZIMUTH"])
    
        # Where wind speed equal to 0, we assume it is buildings
        buildIndexAll = df_selectedPoints[(df_selectedPoints[WINDSPEED_X.upper()]==0) &\
                                          (df_selectedPoints[WINDSPEED_Y.upper()]==0) &\
                                          (df_selectedPoints[WINDSPEED_Z.upper()]==0)].index
        df_selectedPoints.loc[buildIndexAll, [WINDSPEED_X.upper(), WINDSPEED_Y.upper(), WINDSPEED_Z.upper()]] = np.nan
            
        uniques_z = pd.unique(df_selectedPoints[Z.upper()])
        # Need to reindex regularly values for stream plot
        if isStream:
            uniques_z[uniques_z == 0] = 0 - float(horiz_res) / 2
            df_selectedPoints[Z.upper()] = df_selectedPoints[Z.upper()].replace(0, 0 - float(horiz_res) / 2)
            dic_all = {id_line : {zval : df_selectedPoints[(df_selectedPoints[idLines.upper()] == id_line) &
                                                            (df_selectedPoints[Z.upper()] == zval)].groupby("DIST").mean()
                                      for zval in uniques_z}
                           for id_line in pd.unique(df_selectedPoints[idLines.upper()])}
            dic_all = {id_line : {zval : conditional_interpolate(dic_all[id_line][zval].reindex(pd.Index(np.arange(0, 
                                                                                                                   dic_all[id_line][zval].index.max(), 
                                                                                                                   horiz_res))\
                                                                                                .union(dic_all[id_line][zval].index)).sort_index(),
                                                                 cols = [WINDSPEED_X.upper(),
                                                                         WINDSPEED_Y.upper(),
                                                                         WINDSPEED_Z.upper()],
                                                                 limit = 2).reindex(pd.Index(np.arange(0, 
                                                                                                       dic_all[id_line][zval].index.max(), 
                                                                                                       horiz_res)))
                                      for zval in dic_all[id_line]}
                           for id_line in dic_all}
                                               
        if not fig and not ax:
            fig = {}
            ax = {}
        if not scale:
            scale = {}
        for line in sorted(set(df_selectedPoints[idLines.upper()])):
            if not fig.get(line) and not ax.get(line):
                fig[line], ax[line] = plt.subplots(figsize = (15,7))
            df_plot = df_selectedPoints[df_selectedPoints[idLines.upper()] == line].groupby([Z.upper(), "DIST"]).mean()
            
            if isStream:
                uniques_dist = dic_all[line][uniques_z[0]].index.unique()
                D = np.array([[d for d in uniques_dist] for z in uniques_z])
                z = np.array([[z for d in uniques_dist] for z in uniques_z])
                wind_d = np.array([[dic_all[line][zval].loc[d, "PROJECTED_HORIZ_WIND"] for d in uniques_dist] for zval in uniques_z])
                wind_z = np.array([[dic_all[line][zval].loc[d, WINDSPEED_Z.upper()] for d in uniques_dist] for zval in uniques_z])
                ax[line].streamplot(D, z, wind_d, wind_z, density = STREAM_DENSITY,
                                    color = color)
            else:
                uniques_dist = df_plot.index.unique(level = "DIST")
                D = np.array([[d for d in uniques_dist] for z in uniques_z])
                z = np.array([[z for d in uniques_dist] for z in uniques_z])
                wind_d = np.array([[df_plot.loc[zval, d]["PROJECTED_HORIZ_WIND"] for d in uniques_dist] for zval in uniques_z])
                wind_z = np.array([[df_plot.loc[zval, d][WINDSPEED_Z.upper()] for d in uniques_dist] for zval in uniques_z])
                
                if not scale.get(line):
                    if np.max(np.abs(wind_d)) > 3 * np.median(np.abs(wind_d)):
                        scale[line] = np.max(np.abs(wind_d)) / (1.5 * horiz_res)
                    else:
                        scale[line] = np.median(np.abs(wind_d)) / (1.5 * horiz_res)
                Q = ax[line].quiver(D, z, wind_d, wind_z, 
                                    units = 'xy', scale = scale[line],
                                    headwidth = HEAD_WIDTH, headlength = HEAD_LENGTH,
                                    headaxislength = HEAD_AXIS_LENGTH,
                                    width = WIDTH, color = color,
                                    edgecolor = "k", linewidth = 0.2)
                ax[line].quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
                                   coordinates='figure', color = color)
    
                
            # Set buildings using a given color
            buildIndexes = df_plot[np.isnan(df_plot[WINDSPEED_X.upper()]) &\
                                   np.isnan(df_plot[WINDSPEED_Y.upper()]) &\
                                   np.isnan(df_plot[WINDSPEED_Z.upper()])].sort_index(axis = 0, 
                                                                                      level = 1).index
            # Sort distances from start of the line and z from ground
            sorted_dist = pd.Index(np.sort(df_plot.index.unique(level = "DIST")))
            sorted_z = pd.Index(np.sort(z[:,0]))
            rect = {}
            for i, bcell in enumerate(buildIndexes):
                # Get starting height and height of building rectangle
                loc_z = sorted_z.get_loc(bcell[0])
                if loc_z == 0:
                    rec_z0 = sorted_z[loc_z]
                    rec_height = 0.5 * (sorted_z[loc_z + 1] - sorted_z[loc_z])
                elif loc_z == sorted_z.size - 1:
                    rec_z0 = 0.5 * (sorted_z[loc_z - 1] + sorted_z[loc_z])
                    rec_height = sorted_z[loc_z] - sorted_z[loc_z - 1]            
                else:
                    rec_z0 = sorted_z[loc_z] - 0.5 * (sorted_z[loc_z] - sorted_z[loc_z - 1])
                    rec_height = 0.5 * (sorted_z[loc_z + 1] - sorted_z[loc_z - 1])                
                # Get starting distance and width of building rectangle
                loc_d = sorted_dist.get_loc(bcell[1])    
                if loc_d == 0:
                    rec_d0 = sorted_dist[loc_d] - 0.5 * (sorted_dist[loc_d + 1] - sorted_dist[loc_d])
                    rec_width = sorted_dist[loc_d + 1] - sorted_dist[loc_d]
                elif loc_d == sorted_dist.size - 1:
                    rec_d0 = 0.5 * (sorted_dist[loc_d-1] + sorted_dist[loc_d])
                    rec_width = sorted_dist[loc_d] - sorted_dist[loc_d - 1] 
                else:
                    rec_d0 = 0.5 * (sorted_dist[loc_d-1] + sorted_dist[loc_d])
                    rec_width = 0.5 * (sorted_dist[loc_d + 1] - sorted_dist[loc_d - 1])
                
                # Define and plot the building rectangle
                rect[i] = Rectangle((rec_d0, rec_z0), rec_width, rec_height,
                                    color='grey')
                ax[line].add_patch(rect[i])
            if savePlot:
                fig[line].savefig(os.path.join(outputDirectory, simulationName + "_line" + str(line) + ".png"))
            else:
                ax[line].set_title("Line {0}".format(line))

    return fig, ax, scale, fig_poly, ax_poly
            
def conditional_interpolate(df, cols, limit = 2):
    s = df[cols].isna().prod(axis = 1) == 0
    s = s.ne(s.shift()).cumsum()
    m = df[cols].groupby([s, df[cols].isna().prod(axis = 1) == 0])[cols[0]].transform('size').where(df[cols].isnull().prod(axis = 1).astype(bool))
    
    result = df.interpolate(limit_area='inside', method='slinear').mask(m >= limit)
    
    return result