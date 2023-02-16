
# This is the main file for the SOLWEIG model, python version

# Input variables:
# dsm = digital surface model
# scale = height to pixel size (2m pixel gives scale = 0.5)
# header = ESRI Ascii Grid header
# sizey,sizex = no. of pixels in x and y
# row,col = point of interest. Used if all data from one should be calculated
# svf,svfN,svfW,svfE,svfS = SVFs for building and ground
# svfveg,svfNveg,svfEveg,svfSveg,svfWveg = Veg SVFs blocking sky
# svfaveg,svfEaveg,svfSaveg,svfWaveg,svfNaveg = Veg SVFs blocking buildings
# vegdsm = Vegetation canopy DSM
# vegdsm2 = Vegetation trunk zone DSM
# albedo_b = buildings
# albedo_g = ground (if landcover==0)
# absK = human absorption coefficient for shortwave radiation
# absL = human absorption coefficient for longwave radiation
# ewall = Emissivity of building walls
# eground = Emissivity of ground (if landcover==0)
# Fside = The angular factors between a person and the surrounding surfaces
# Fup = The angular factors between a person and the surrounding surfaces
# PA = Posture of a human
# met = meteorological inputdata
# YYYY = Year
# altitude = Sun altitude (degree)
# azimuth = Sun azimuth (degree)
# zen = Sun zenith angle (radians)
# jday = day of year
# showimage = show image during execuation
# usevegdem = use vegetation scheme
# onlyglobal = calculate dir and diff from global
# buildings = Boolena grid to identify building pixels
# location = geographic location
# height = height of measurments point
# trans Trensmissivity of shortwave theough vegetation
# output = output settings
# fileformat = fileformat of output grids
# landcover = use landcover scheme !!!NEW IN 2015a!!!
# sensorheight = Sensorheight of wind sensor
# leafon = foliated vegetation or not
# lcgrid = grid with landcoverclasses
# lc_class = table with landcover properties
# dectime = decimal time
# altmax = maximum sun altitude
# dirwalls = aspect of walls
# walls = one pixel row outside building footprints
# cyl = consider man as cyliner instead of cube
# elvis = old thing from Jonsson et al.
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()

import numpy as np
import matplotlib.pylab as plt
# from osgeo.gdalconst import *
from osgeo import gdal, osr
import os

from Tgmaps_v1 import Tgmaps_v1
import Solweig_2021a_calc as so
from Tgmaps_v1 import Tgmaps_v1
import Solweig_v2015_metdata_noload as metload
from clearnessindex_2013b import clearnessindex_2013b
import wallalgorithms as wa

from misc import saveraster

# Settings
mainfolder = './data_for_fredrik/'
infolder = mainfolder
insvffolder = infolder + 'svfs/'
metfilepath = infolder + 'SOLWEIG_meteo_DCEP_osz_20180726_10minute_umep_UTC+1.txt'
outfolder = mainfolder + 'Out_test/'
dsmname = 'Grid_BUI_DSM_1200m_Res_5m.tif'
demname = 'Grid_Ground_DEM_1200m_Res_5m.tif'
cdsmname = ''
tdsmname = ''
wallheightname = 'WallHeight_no_veg.tif'
wallaspectname = 'WallAspect_no_veg.tif'
landcovername = ''
perezdata = ''
walltest = r'C:\Users\xlinfr\Downloads\testfiles\wall_height.tif'
showimage = 0
usevegdem = 0
onlyglobal = 0
height = 1.1
trans = 0.03
landcover = 0
sensorheight = 2.0
cyl = 1
elvis = 0
canopyToTrunkRatio = 0.3
PA = 'STAND'
UTC = 1
alt = 3.0
lat = 52.50250799716754
lon = 13.317437543639814
poisxy = np.zeros((1, 3)) - 999
poisxy[0, 1] = 30.
poisxy[0, 2] = 60.
poiname = 'TestTg2021a'
anisdiff = 0

albedo_b = 0.20
albedo_g = 0.15
absK = 0.7
absL = 0.98
ewall = 0.95
eground = 0.95

header = 'yyyy id   it imin dectime altitude azimuth kdir kdiff kglobal kdown   kup    keast ksouth ' \
            'kwest knorth ldown   lup    least lsouth lwest  lnorth   Ta      Tg     RH    Esky   Tmrt    ' \
            'I0     CI   Shadow  SVF_b  SVF_bv KsideI PET UTCI'

if not poisxy is None:
    poi_save = []
    data_out = outfolder + '/POI_' + str(poiname) + '.txt'
    np.savetxt(data_out, poi_save,  delimiter=' ', header=header, comments='')  # fmt=numformat,

if not os.path.exists(outfolder):
    os.makedirs(outfolder)

# load surface grids
# dsm = np.loadtxt(infolder + "DSM_KRbig.asc", skiprows=6) # to load an ERSIascii grid
dataSetDSM = gdal.Open(infolder + dsmname)
dsm = dataSetDSM.ReadAsArray().astype(np.float)
dataSet = gdal.Open(infolder + demname)
dem = dataSet.ReadAsArray().astype(np.float)

geotransform = dataSetDSM.GetGeoTransform()
scale = 1 / geotransform[1]
alt = np.median(dsm)

# not woring. lat lon given manually
# old_cs = osr.SpatialReference()
# # dsm_ref = dataSetDSM.crs().toWkt()
# dsm_ref = dataSetDSM.GetProjection().toWkt()
# old_cs.ImportFromWkt(dsm_ref)

# wgs84_wkt = """
#     GEOGCS["WGS 84",
#         DATUM["WGS_1984",
#             SPHEROID["WGS 84",6378137,298.257223563,
#                 AUTHORITY["EPSG","7030"]],
#             AUTHORITY["EPSG","6326"]],
#         PRIMEM["Greenwich",0,
#             AUTHORITY["EPSG","8901"]],
#         UNIT["degree",0.01745329251994328,
#             AUTHORITY["EPSG","9122"]],
#         AUTHORITY["EPSG","4326"]]"""

# new_cs = osr.SpatialReference()
# new_cs.ImportFromWkt(wgs84_wkt)

# transform = osr.CoordinateTransformation(old_cs, new_cs)

# width = dataSetDSM.RasterXSize
# height = dataSetDSM.RasterYSize
# minx = geotransform[0]
# miny = geotransform[3] + width * geotransform[4] + height * geotransform[5]
# lonlat = transform.TransformPoint(minx, miny)

# gdalver = float(gdal.__version__[0])
# if gdalver == 3.:
#     lon = lonlat[1] #changed to gdal 3
#     lat = lonlat[0] #changed to gdal 3
# else:
#     lon = lonlat[0] #changed to gdal 2
#     lat = lonlat[1] #changed to gdal 2


rows = dsm.shape[0]
cols = dsm.shape[1]
walllimit = 3

if usevegdem == 1:
    dataSet = gdal.Open(infolder + cdsmname)
    vegdsm = dataSet.ReadAsArray().astype(np.float)
    vegdsm2 = vegdsm * canopyToTrunkRatio
else:
    vegdsm = np.zeros([rows, cols])
    vegdsm2 = np.zeros([rows, cols])
    bush = np.zeros([rows, cols])

dataSet = gdal.Open(infolder + wallheightname)
wallheight = dataSet.ReadAsArray().astype(np.float)
dataSet = gdal.Open(infolder + wallaspectname)
wallaspect = dataSet.ReadAsArray().astype(np.float)
# walls = wa.findwalls(dsm, walllimit)
# dirwalls = wa.filter1Goodwin_as_aspect_v3(walls, scale, dsm)

dataSet = gdal.Open(insvffolder + "svf.tif")
svf = dataSet.ReadAsArray().astype(np.float)
dataSet = gdal.Open(insvffolder + "svfN.tif")
svfN = dataSet.ReadAsArray().astype(np.float)
dataSet = gdal.Open(insvffolder + "svfS.tif")
svfS = dataSet.ReadAsArray().astype(np.float)
dataSet = gdal.Open(insvffolder + "svfE.tif")
svfE = dataSet.ReadAsArray().astype(np.float)
dataSet = gdal.Open(insvffolder + "svfW.tif")
svfW = dataSet.ReadAsArray().astype(np.float)

if usevegdem == 1:
    dataSet = gdal.Open(insvffolder + "svfveg.tif")
    svfveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfNveg.tif")
    svfNveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfSveg.tif")
    svfSveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfEveg.tif")
    svfEveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfWveg.tif")
    svfWveg = dataSet.ReadAsArray().astype(np.float)

    dataSet = gdal.Open(insvffolder + "svfaveg.tif")
    svfaveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfNaveg.tif")
    svfNaveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfSaveg.tif")
    svfSaveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfEaveg.tif")
    svfEaveg = dataSet.ReadAsArray().astype(np.float)
    dataSet = gdal.Open(insvffolder + "svfWaveg.tif")
    svfWaveg = dataSet.ReadAsArray().astype(np.float)
else:
    svfveg = np.ones((rows, cols))
    svfNveg = np.ones((rows, cols))
    svfSveg = np.ones((rows, cols))
    svfEveg = np.ones((rows, cols))
    svfWveg = np.ones((rows, cols))
    svfaveg = np.ones((rows, cols))
    svfNaveg = np.ones((rows, cols))
    svfSaveg = np.ones((rows, cols))
    svfEaveg = np.ones((rows, cols))
    svfWaveg = np.ones((rows, cols))

header = None

if PA == 'STAND':
    Fside = 0.22
    Fup = 0.06
    height = 1.1
    Fcyl = 0.28
else:
    Fside = 0.166666
    Fup = 0.166667
    height = 0.75
    Fcyl = 0.2

metfile = 1
if metfile == 1:
    met = np.loadtxt(metfilepath, skiprows=1, delimiter=' ')
else:
    met = np.zeros((1, 24)) - 999.

    year = 2011
    month = 6
    day = 6
    hour = 12
    minu = 30

    if (year % 4) == 0:
        if (year % 100) == 0:
            if (year % 400) == 0:
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

    doy = np.sum(dayspermonth[0:month - 1]) + day

    Ta = 25.
    RH = 50
    radG = 880.
    radD = 150.
    radI = 950.

    met[0, 0] = year
    met[0, 1] = doy
    met[0, 2] = hour
    met[0, 3] = minu
    met[0, 11] = Ta
    met[0, 10] = RH
    met[0, 14] = radG
    met[0, 21] = radD
    met[0, 22] = radI

location = {'longitude': lon, 'latitude': lat, 'altitude': alt}
YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = metload.Solweig_2015a_metdata_noload(met, location, UTC)

if landcover == 1:
    dataSet = gdal.Open(infolder + landcovername)
    lcgrid = dataSet.ReadAsArray().astype(np.float)
    sitein = "landcoverclasses_2018a.txt"
    f = open(sitein)
    lin = f.readlines()
    lc_class = np.zeros((lin.__len__() - 1, 6))
    for i in range(1, lin.__len__()):
        lines = lin[i].split()
        for j in np.arange(1, 7):
            lc_class[i - 1, j - 1] = float(lines[j])

    buildings = np.copy(lcgrid)
    buildings[buildings == 7] = 1
    buildings[buildings == 6] = 1
    buildings[buildings == 5] = 1
    buildings[buildings == 4] = 1
    buildings[buildings == 3] = 1
    buildings[buildings == 2] = 0
else:
    buildings = dsm - dem
    buildings[buildings < 2.] = 1.
    buildings[buildings >= 2.] = 0.
    lcgrid = None

# outputfolder = infolder + 'Out_v2/'

# From core but outside core loop
# if not row:
#     #% This is settings for calculating PET at POI
#     pet = []
#     pet.append({'body' : 75,'age':35,'height':1.80,'activity':80,'sex':1.,'clo':0.9})

# fn=fopen([outputfolder 'Output_v2015a_POI_' PA '.txt'],'w');
# fprintf(fn, '%9s', 'year', 'DOY', 'hour', 'dectime', 'altitude', 'azimuth', 'Kdirect', 'Kdiffuse', 'Kglobal', 'Kdown', 'Kup', 'KsideI', 'Knorth', 'Keast', 'Ksouth', 'Kwest', 'Ldown', 'Lup', 'Lnorth', 'Least', 'Lsouth', 'Lwest', 'Ta', 'Tg', 'RH', 'Ea', 'Esky', 'Sstr', 'Tmrt', 'I0', 'CI', 'gvf', 'CITg', 'Shadow', 'SVF_b', 'SVF_b+v', 'PET', 'UTCI')
# fprintf(fn, '\r\n')

# Surface to air temperature difference at sunrise
# Tstart=0;%3.41; % dynamic as from 2015a
# Initialization of maps
Knight = np.zeros((rows, cols))
Tmrtday = np.zeros((rows, cols))
Lupday = np.zeros((rows, cols))
Ldownday = np.zeros((rows, cols))
Kupday = np.zeros((rows, cols))
Kdownday = np.zeros((rows, cols))
gvfday = np.zeros((rows, cols))
Tmrtdiurn = np.zeros((rows, cols))
Lupdiurn = np.zeros((rows, cols))
Ldowndiurn = np.zeros((rows, cols))
Tgmap1 = np.zeros((rows, cols))
Tgmap1E = np.zeros((rows, cols))
Tgmap1S = np.zeros((rows, cols))
Tgmap1W = np.zeros((rows, cols))
Tgmap1N = np.zeros((rows, cols))
Tgmap1Athens = np.zeros((rows, cols))

tmp = svf+svfveg-1.
tmp[tmp < 0.] = 0.
# matlab crazyness around 0
svfalfa = np.arcsin(np.exp((np.log((1.-tmp))/2.)))

# Creating vectors from meteorological input
DOY = met[:, 1]
hours = met[:, 2]
minu = met[:, 3]
Ta = met[:, 11]
RH = met[:, 10]
radG = met[:, 14]
radD = met[:, 21]
radI = met[:, 22]
P = met[:, 12]
Ws = met[:, 9]
Twater = []

#%Number of daytime hours
# Initialisation of time related variables
if Ta.__len__() == 1:
    timestepdec = 0.0
else:
    timestepdec = dectime[1] - dectime[0]
timeadd = 0.
timeaddE = 0.
timeaddS = 0.
timeaddW = 0.
timeaddN = 0.
firstdaytime = 1.
# bugfix so that model can start during daytime
# Parameterisarion for Lup
if not height:
    height = 1.1

# Radiative surface influence, Rule of thumb by Schmid et al. (1990).
first = np.round(height)

if first == 0.:
    first = 1.

second = np.round((height*20.))

if usevegdem == 1:
    # Vegetation transmittivity of shortwave radiation
    psi = leafon * trans
    psi[leafon == 0] = 0.5
    # amaxvalue
    vegmax = vegdsm.max()
    amaxvalue = dsm.max() - dsm.min()
    amaxvalue = np.maximum(amaxvalue, vegmax)

    # Elevation vegdsms if buildingDEM includes ground heights
    vegdsm = vegdsm + dsm
    vegdsm[vegdsm == dsm] = 0
    vegdsm2 = vegdsm2 + dsm
    vegdsm2[vegdsm2 == dsm] = 0

    # % Bush separation
    bush = np.logical_not((vegdsm2 * vegdsm)) * vegdsm
    svfbuveg = (svf - (1. - svfveg) * (1. - trans))  # major bug fixed 20141203
else:
    psi = leafon * 0. + 1.
    svfbuveg = svf
    amaxvalue = dsm.max() - dsm.min()

# Ts parameterisation maps
if landcover == 1.:
    [TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST, TmaxLST_wall] = Tgmaps_v1(lcgrid, lc_class)
else:
    TgK = Knight+0.37
    Tstart = Knight-3.41
    alb_grid = Knight+albedo_g
    emis_grid = Knight+eground
    TgK_wall = 0.37
    Tstart_wall = -3.41
    TmaxLST = 15.
    TmaxLST_wall = 15.
    lcgrid = None
    lc_class = None

# Import shadow matrices (Anisotropic sky)
if anisdiff == 1:
    ani = 1
    data = np.load(perezdata)
    shmat = data['shadowmat']
    vegshmat = data['vegshadowmat']
    if usevegdem == 1:
        diffsh = np.zeros((rows, cols, 145))
        for i in range(0, 145):
            diffsh[:, :, i] = shmat[:, :, i] - (1 - vegshmat[:, :, i]) * (1 - trans)
    else:
        diffsh = shmat
else:
    ani = 0
    diffsh = None

# If metfile starts at night
CI = 1.

poitmrt = np.zeros((Ta.__len__()))
i = 0

# New code from UMEP
tmrtplot = np.zeros((rows, cols))
TgOut1 = np.zeros((rows, cols))

numformat = '%3d %2d %3d %2d %6.5f ' + '%6.2f ' * 28

for i in np.arange(0, Ta.__len__()):
    print(i)
    # Daily water body temperature
    if landcover == 1:
        if (dectime[i] - np.floor(dectime[i])) == 0 or (i == 0):
            Twater = np.mean(Ta[jday[0] == np.floor(dectime[i])])

    # Nocturnal cloudfraction from Offerle et al. 2003
    if (dectime[i] - np.floor(dectime[i])) == 0:

        daylines = np.where(np.floor(dectime) == dectime[i])
        alt = altitude[0][daylines]
        alt2 = np.where(alt > 1)
        rise = alt2[0][0]
        [_, CI, _, _, _] = clearnessindex_2013b(zen[0, i + rise + 1], jday[0, i + rise + 1],
                                                Ta[i + rise + 1],
                                                RH[i + rise + 1] / 100., radG[i + rise + 1], location,
                                                P[i + rise + 1])  # i+rise+1 to match matlab code. correct?
        if (CI > 1) or (CI == np.inf):
            CI = 1

    Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, timeadd, \
        Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, Keast, Ksouth, Kwest, Knorth, Least, \
        Lsouth, Lwest, Lnorth, KsideI, TgOut1, TgOut, radIout, radDout = so.Solweig_2021a_calc(
            i, dsm, scale, rows, cols, svf, svfN, svfW, svfE, svfS, svfveg,
            svfNveg, svfEveg, svfSveg, svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg,
            vegdsm, vegdsm2, albedo_b, absK, absL, ewall, Fside, Fup, Fcyl, altitude[0][i],
            azimuth[0][i], zen[0][i], jday[0][i], usevegdem, onlyglobal, buildings, location,
            psi[0][i], landcover, lcgrid, dectime[i], altmax[0][i], wallaspect,
            wallheight, cyl, elvis, Ta[i], RH[i], radG[i], radD[i], radI[i], P[i], amaxvalue,
            bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
            TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timestepdec, 
            Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI, TgOut1, diffsh, ani)
    
    # Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, timeadd, \
    # Tgmap1, timeaddE, Tgmap1E, timeaddS, Tgmap1S, timeaddW, Tgmap1W, timeaddN, Tgmap1N, \
    # Keast, Ksouth, Kwest, Knorth, Least, Lsouth, Lwest, Lnorth, KsideI, Tgmap1Athens, TgAthens \
    #     = so.Solweig_2015a_calc(i, dsm, scale, rows, cols, svf, svfN, svfW, svfE, svfS, svfveg,
    #         svfNveg, svfEveg, svfSveg, svfWveg, svfaveg, svfEaveg, svfSaveg, svfWaveg, svfNaveg,
    #         vegdsm, vegdsm2, albedo_b, absK, absL, ewall, Fside, Fup, altitude[0][i],
    #         azimuth[0][i], zen[0][i], jday[0][i], usevegdem, onlyglobal, buildings, location,
    #         psi[0][i], landcover, lcgrid, dectime[i], altmax[0][i], wallaspect,
    #         wallheight, cyl, elvis, Ta[i], RH[i], radG[i], radD[i], radI[i], P[i], amaxvalue,
    #         bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
    #         TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timeaddE, timeaddS,
    #         timeaddW, timeaddN, timestepdec, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI, Tgmap1Athens)

    tmrtplot = tmrtplot + Tmrt

    if altitude[0][i] > 0:
        w = 'D'
    else:
        w = 'N'

    # Write to POIs
    if not poisxy is None:
        for k in range(0, poisxy.shape[0]):
            poi_save = np.zeros((1, 33))
            poi_save[0, 0] = YYYY[0][i]
            poi_save[0, 1] = jday[0][i]
            poi_save[0, 2] = hours[i]
            poi_save[0, 3] = minu[i]
            poi_save[0, 4] = dectime[i]
            poi_save[0, 5] = altitude[0][i]
            poi_save[0, 6] = azimuth[0][i]
            poi_save[0, 7] = radI[i]
            poi_save[0, 8] = radD[i]
            poi_save[0, 9] = radG[i]
            poi_save[0, 10] = Kdown[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 11] = Kup[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 12] = Keast[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 13] = Ksouth[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 14] = Kwest[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 15] = Knorth[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 16] = Ldown[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 17] = Lup[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 18] = Least[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 19] = Lsouth[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 20] = Lwest[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 21] = Lnorth[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 22] = Ta[i]
            poi_save[0, 23] = Tg[int(poisxy[k, 2]), int(poisxy[k, 1])] + Ta[i]
            poi_save[0, 24] = RH[i]
            poi_save[0, 25] = esky
            poi_save[0, 26] = Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 27] = I0
            poi_save[0, 28] = CI
            poi_save[0, 29] = shadow[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 30] = svf[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 31] = svfbuveg[int(poisxy[k, 2]), int(poisxy[k, 1])]
            poi_save[0, 32] = KsideI[int(poisxy[k, 2]), int(poisxy[k, 1])]

            data_out = outfolder + '/POI_' + str(poiname) + '.txt'
            f_handle = open(data_out, 'ab')
            np.savetxt(f_handle, poi_save, fmt=numformat)
            f_handle.close()

    if hours[i] < 10:
        XH = '0'
    else:
        XH = ''
    if minu[i] < 10:
        XM = '0'
    else:
        XM = ''
    
    # saveraster(dataSetDSM, outfolder + '/Tmrt_' + str(int(YYYY[0, i])) + '_' + str(int(DOY[i])) +
    #                '_' + XH + str(int(hours[i])) + XM + str(int(minu[i])) + '.tif', Tmrt)
    # if altitude[0][i] > 0:
    #     saveraster(dataSetDSM, outfolder + '/Tg_' + str(int(YYYY[0, i])) + '_' + str(int(DOY[i])) +
    #                '_' + XH + str(int(hours[i])) + XM + str(int(minu[i])) + '.tif', TgOut)
#
# # Loop through time series
# for i in np.arange(0, Ta.__len__()):
#     print i
#     # Daily water body temperature
#     if landcover == 1:
#         if ((dectime[i] - np.floor(dectime[i]))) == 0 or (i == 0):
#             Twater = np.mean(Ta[jday[0] == np.floor(dectime[i])])
#     else:
#         Twater = []
#
#     # Nocturnal cloudfraction from Offerle et al. 2003
#     if (dectime[i] - np.floor(dectime[i])) == 0:
#         # alt = altitude[i:altitude.__len__()]
#         daylines = np.where(np.floor(dectime) == dectime[i])
#         alt = altitude[0][daylines]
#         alt2 = np.where(alt > 1)
#         rise = alt2[0][0]
#         [_, CI, _, _, _] = clearnessindex_2013b(zen[0, i + rise + 1], jday[0, i + rise + 1], Ta[i + rise + 1],
#                                                 RH[i + rise + 1] / 100., radG[i + rise + 1], location, P[i + rise + 1])  # i+rise+1 to match matlab code. correct?
#         if (CI > 1) or (CI == np.inf):
#             CI = 1
#
#     # Main calcualtions
#     Tmrt, Kdown, Kup, Ldown, Lup, Tg, ea, esky, I0, CI, shadow, firstdaytime, timestepdec, timeadd, Tgmap1, timeaddE, \
#     Tgmap1E, timeaddS, Tgmap1S, timeaddW, Tgmap1W, timeaddN, Tgmap1N = so.Solweig_2015a_calc(i, dsm, scale, rows,
#                         cols, svf, svfN, svfW, svfE, svfS, svfveg, svfNveg, svfEveg, svfSveg, svfWveg, svfaveg,
#                         svfEaveg, svfSaveg, svfWaveg, svfNaveg, vegdem, vegdem2, albedo_b,
#                         absK, absL, ewall, Fside, Fup, altitude[0][i], azimuth[0][i], zen[0][i], jday[0][i],
#                         usevegdem, onlyglobal, buildings, location, psi[0][i], landcover, lcgrid, dectime[i],
#                         altmax[0][i], dirwalls, walls, cyl, elvis, Ta[i], RH[i], radG[i], radD[i], radI[i], P[i],
#                         amaxvalue, bush, Twater, TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST,
#                         TmaxLST_wall, first, second, svfalfa, svfbuveg, firstdaytime, timeadd, timeaddE, timeaddS,
#                         timeaddW, timeaddN, timestepdec, Tgmap1, Tgmap1E, Tgmap1S, Tgmap1W, Tgmap1N, CI)
#
#     filename = 'TmrtoutsideUMEP_' + str(int(YYYY[0, i])) + '_' + str(int(DOY[i])) + '_' + str(int(hours[i])) + str(int(minu[i])) + '.tif'
#     print(azimuth[0][i])
#     print(altitude[0][i])
#     print(hours[i])
#
#     # plt.matshow(Tmrt)
#     # plt.colorbar()
#     # plt.show()
#     # poitmrt[i] = Tmrt[51, 118]
#     k=4
#
#     # plt.plot(poitmrt)
#     # plt.show()
#
#     # def saveraster(self, gdal_data, filename, raster):
#     rows = dataSet.RasterYSize
#     cols = dataSet.RasterXSize
#
#     outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
#     outBand = outDs.GetRasterBand(1)
#
#     # write the data
#     outBand.WriteArray(Tmrt, 0, 0)
#     # flush data to disk, set the NoData value and calculate stats
#     outBand.FlushCache()
#     outBand.SetNoDataValue(-9999)
#
#     # georeference the image and set the projection
#     outDs.SetGeoTransform(dataSet.GetGeoTransform())
#     outDs.SetProjection(dataSet.GetProjection())
