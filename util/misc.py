__author__ = 'xlinfr'

import numpy as np
from osgeo import gdal, osr
from osgeo.gdalconst import GDT_Float32


# Slope and aspect used in SEBE and Wall aspect
def get_ders(dsm, scale):
    # dem,_,_=read_dem_grid(dem_file)
    dx = 1/scale
    # dx=0.5
    fy, fx = np.gradient(dsm, dx, dx)
    asp, grad = cart2pol(fy, fx, 'rad')
    slope = np.arctan(grad)
    asp = asp * -1
    asp = asp + (asp < 0) * (np.pi * 2)
    return slope, asp


def cart2pol(x, y, units='deg'):
    radius = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    if units in ['deg', 'degs']:
        theta = theta * 180 / np.pi
    return theta, radius


def saveraster(gdal_data, filename, raster):
    rows = gdal_data.RasterYSize
    cols = gdal_data.RasterXSize

    outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
    outBand = outDs.GetRasterBand(1)

    # write the data
    outBand.WriteArray(raster, 0, 0)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    outBand.SetNoDataValue(-9999)

    # georeference the image and set the projection
    outDs.SetGeoTransform(gdal_data.GetGeoTransform())
    outDs.SetProjection(gdal_data.GetProjection())

def saverasternd(gdal_data, filename, raster):
    rows = gdal_data.RasterYSize
    cols = gdal_data.RasterXSize

    outDs = gdal.GetDriverByName("GTiff").Create(filename, cols, rows, int(1), GDT_Float32)
    outBand = outDs.GetRasterBand(1)

    # write the data
    outBand.WriteArray(raster, 0, 0)
    # flush data to disk, set the NoData value and calculate stats
    outBand.FlushCache()
    # outBand.SetNoDataValue(-9999)

    # georeference the image and set the projection
    outDs.SetGeoTransform(gdal_data.GetGeoTransform())
    outDs.SetProjection(gdal_data.GetProjection())

def xy2latlon(crsWtkIn, x, y):
    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(crsWtkIn)

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

    lonlat = transform.TransformPoint(x, y)

    gdalver = float(gdal.__version__[0])
    if gdalver >= 3.:
        lon = lonlat[1] #changed to gdal 3
        lat = lonlat[0] #changed to gdal 3
    else:
        lon = lonlat[0] #changed to gdal 2
        lat = lonlat[1] #changed to gdal 2

    return lat, lon



def createTSlist():
    import pytz
    from datetime import datetime

    # Get the current time in naive UTC
    now = datetime.utcnow()

    # Dictionary to store timezones grouped by UTC offset
    timezones_by_offset = {}

    # Iterate through all timezones
    for tz in pytz.all_timezones:
        timezone = pytz.timezone(tz)
        localized_time = timezone.localize(now, is_dst=None)
        offset = localized_time.utcoffset()
        if offset is not None:
            offset_hours = offset.total_seconds() / 3600
            offset_str = f"UTC{offset_hours:+03.0f}:00"

            if offset_str not in timezones_by_offset:
                timezones_by_offset[offset_str] = {
                    "utc_offset": offset_hours,
                    "timezones": []
                }
            # Add up to 3 example timezones per offset
            if len(timezones_by_offset[offset_str]["timezones"]) < 3:
                timezones_by_offset[offset_str]["timezones"].append(tz)

    # Convert the dictionary to a list of dictionaries
    timezones_list = [
        {
            "utc_offset_str": offset,
            "utc_offset": data["utc_offset"],
            "timezones": data["timezones"]
        }
        for offset, data in timezones_by_offset.items()
    ]

    return timezones_list, timezones_by_offset
