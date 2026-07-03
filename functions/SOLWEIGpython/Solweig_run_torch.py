# This is the main function of the SOLWEIG model
# 2025-Mar-21
# Fredrik Lindberg, fredrikl@gvc.gu.se
# Goteborg Urban Climate Group
# Gothenburg University

# sommon imports

from __future__ import absolute_import

try:
    import torch
except ImportError:
    pass

import numpy as np
from ...util.umep_solweig_export_component import read_solweig_config
from ...util.SEBESOLWEIGCommonFiles.Solweig_v2015_metdata_noload_torch import (
    Solweig_2015a_metdata_noload,
)
from ...util.SEBESOLWEIGCommonFiles.clearnessindex_2013b_torch import (
    clearnessindex_2013b,
)

from . import Solweig_2026a_calc_forprocessing_torch as so

# from ...functions.SOLWEIGpython import WriteMetadataSOLWEIG # Not needed anymore?
from . import PET_calculations as p
from . import UTCI_calculations as utci
from .CirclePlotBar import PolarBarPlot
from .wall_surface_temperature import load_walls
from .wallOfInterest import pointOfInterest
from .patch_characteristics import hemispheric_image
from .wallsAsNetCDF import walls_as_netcdf
from .Tgmaps_v1_torch import Tgmaps_v1
from .. import wallalgorithms_torch as wa
from .ground_surface_torch import initiate_groundScheme
import json
import zipfile
import pandas as pd
import matplotlib.pylab as plt
from shutil import copyfile

# imports from osgeo/qgis dependency
try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    from ...util.misc import saveraster, xy2latlon_fromraster
    from qgis.core import QgsRasterLayer
except:
    pass

# imports for standalone
try:
    from umep import common
    from rasterio.transform import xy, rowcol
    import pyproj
    from tqdm import tqdm
    import geopandas as gpd
except:
    pass


def solweig_run(configPath, feedback):
    """
    Input:
    configPath : config file including geodata paths and settings.
    feedback : To communicate with qgis gui. Set to None if standalone
    """

    # Load config file
    configDict = read_solweig_config(configPath)

    # Load parameters settings for SOLWEIG
    with open(configDict["para_json_path"], "r") as jsn:
        param = json.load(jsn)

    with torch.no_grad():

        # reading variables from config and parameters that is not yet presented
        standAlone = int(configDict["standalone"])
        cyl = int(configDict["cyl"])
        albedo_b = param["Albedo"]["Effective"]["Value"]["Walls"]
        ewall = param["Emissivity"]["Value"]["Walls"]
        onlyglobal = int(configDict["onlyglobal"])
        elvis = 0.0
        absK = param["Tmrt_params"]["Value"]["absK"]
        absL = param["Tmrt_params"]["Value"]["absL"]

        # --- Load on CPU or GPU config
        device = torch.device("cpu")
        if (
            configDict["calculation_mode"] == "gpu"
            and torch.cuda.is_available()
        ):
            device = torch.device("cuda")
            if feedback is not None:
                feedback.setProgressText(
                    "PyTorch and NVIDIA/AMD GPU found. Initiating CUDA mode..."
                )
            else:
                print(
                    "PyTorch and NVIDIA/AMD GPU found. Initiating CUDA mode..."
                )

        elif (
            configDict["calculation_mode"] == "gpu"
            and hasattr(torch, "xpu")
            and torch.xpu.is_available()
        ):
            device = torch.device("xpu")
            if feedback is not None:

                feedback.setProgressText(
                    "PyTorch and Intel GPU found. Initiating XPU mode..."
                )
            else:
                print("PyTorch and Intel GPU found. Initiating XPU mode...")

        else:
            if feedback is not None:
                feedback.setProgressText(
                    "PyTorch found but compatible GPU not found. Initiating CPU mode..."
                )
            else:
                print(
                    "PyTorch found but compatible GPU not found. Initiating CPU mode..."
                )
        standAlone = int(configDict["standalone"])

        # Load DSM
        if standAlone == 1:
            dsm, dsm_transf, dsm_crs = common.load_raster(
                configDict["filepath_dsm"], bbox=None
            )
            scale = 1 / dsm_transf.a
            # dsm_height, dsm_width = dsm.shape  # y rows by x cols
            # y is flipped - so return max for lower row
            minx, miny = xy(dsm_transf, dsm.shape[0], 0)
            # Define the source and target CRS
            source_crs = pyproj.CRS(dsm_crs)
            target_crs = pyproj.CRS(4326)  # WGS 84
            # Create a transformer object
            transformer = pyproj.Transformer.from_crs(
                source_crs, target_crs, always_xy=True
            )
            # Perform the transformation
            lon, lat = transformer.transform(minx, miny)
            nd = -9999  # TODO: extract nodatavalue from rasterio
        else:
            # dsmlayer = QgsRasterLayer(configDict['filepath_dsm'])
            dsm_wkt = QgsRasterLayer(configDict["filepath_dsm"]).crs().toWkt()
            gdal_dsm = gdal.Open(configDict["filepath_dsm"])
            lat, lon, scale, minx, miny = xy2latlon_fromraster(
                dsm_wkt, gdal_dsm
            )
            dsm = torch.from_numpy(gdal_dsm.ReadAsArray().astype(float)).to(
                device
            )
            nd = gdal_dsm.GetRasterBand(1).GetNoDataValue()

        rows = dsm.shape[0]
        cols = dsm.shape[1]

        # response to issue #85
        dsm[dsm == nd] = 0.0
        if dsm.min() < 0:
            dsmraise = torch.abs(dsm.min())
            dsm = dsm + dsmraise
        else:
            dsmraise = 0

        alt = torch.median(dsm)
        if alt < 0:
            alt = 3

        # Vegetation
        transVeg = param["Tree_settings"]["Value"]["Transmissivity"]
        trunkratio = param["Tree_settings"]["Value"]["Trunk_ratio"]
        usevegdem = int(configDict["usevegdem"])
        if usevegdem == 1:
            if standAlone == 0:
                vegdsm = torch.from_numpy(
                    (
                        gdal.Open(configDict["filepath_cdsm"])
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
            else:
                vegdsm, _, _ = common.load_raster(
                    configDict["filepath_cdsm"], bbox=None
                )
            if configDict["filepath_tdsm"] != "":
                if standAlone == 0:
                    vegdsm2 = torch.from_numpy(
                        (
                            gdal.Open(configDict["filepath_tdsm"])
                            .ReadAsArray()
                            .astype(float)
                        )
                    ).to(device)
                else:
                    vegdsm2, _, _ = common.load_raster(
                        configDict["filepath_tdsm"], bbox=None
                    )
            else:
                vegdsm2 = vegdsm * trunkratio
        else:
            vegdsm = 0
            vegdsm2 = 0

        # Clean up intermediate tensors after vegetation loading
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        # Land cover
        landcover = int(configDict["landcover"])
        if landcover == 1:
            if standAlone == 0:
                lcgrid = torch.from_numpy(
                    (
                        gdal.Open(configDict["filepath_lc"])
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
            else:
                lcgrid, _, _ = common.load_raster(
                    configDict["filepath_lc"], bbox=None
                )
        else:
            lcgrid = torch.ones_like(dsm, device=device)

        # DEM for buildings #TODO: fix nodata in standalone
        demforbuild = int(configDict["demforbuild"])
        if demforbuild == 1:
            if standAlone == 0:
                gdal_dem = gdal.Open(
                    configDict["filepath_dem"]
                )  # .ReadAsArray().astype(float)
                dem = torch.from_numpy(
                    gdal_dem.ReadAsArray().astype(float)
                ).to(device)
                nd = gdal_dem.GetRasterBand(1).GetNoDataValue()
            else:
                dem, _, _ = common.load_raster(
                    configDict["filepath_dem"], bbox=None
                )
                nd = -9999  # TODO: standAlone nd exposure

            # response to issue and #230
            dem[dem == nd] = 0.0
            if dem.min() < 0:
                demraise = torch.abs(dem.min())
                dem = dem + demraise
            else:
                demraise = 0

        # Clean up DEM intermediate tensors
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        # SVF
        zip = zipfile.ZipFile(configDict["input_svf"], "r")
        zip.extractall(configDict["working_dir"])
        zip.close()

        if standAlone == 0:
            svf = torch.from_numpy(
                (
                    gdal.Open(configDict["working_dir"] + "/svf.tif")
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
            svfN = torch.from_numpy(
                (
                    gdal.Open(configDict["working_dir"] + "/svfN.tif")
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
            svfS = torch.from_numpy(
                (
                    gdal.Open(configDict["working_dir"] + "/svfS.tif")
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
            svfE = torch.from_numpy(
                (
                    gdal.Open(configDict["working_dir"] + "/svfE.tif")
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
            svfW = torch.from_numpy(
                (
                    gdal.Open(configDict["working_dir"] + "/svfW.tif")
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
        else:
            svf, _, _ = common.load_raster(
                configDict["working_dir"] + "/svf.tif", bbox=None
            )
            svfN, _, _ = common.load_raster(
                configDict["working_dir"] + "/svfN.tif", bbox=None
            )
            svfS, _, _ = common.load_raster(
                configDict["working_dir"] + "/svfS.tif", bbox=None
            )
            svfE, _, _ = common.load_raster(
                configDict["working_dir"] + "/svfE.tif", bbox=None
            )
            svfW, _, _ = common.load_raster(
                configDict["working_dir"] + "/svfW.tif", bbox=None
            )

        if usevegdem == 1:
            if standAlone == 0:
                svfveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfNveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfNveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfSveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfSveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfEveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfEveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfWveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfWveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)

                svfaveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfaveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfNaveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfNaveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfSaveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfSaveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfEaveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfEaveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
                svfWaveg = torch.from_numpy(
                    (
                        gdal.Open(configDict["working_dir"] + "/svfWaveg.tif")
                        .ReadAsArray()
                        .astype(float)
                    )
                ).to(device)
            else:
                svfveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfveg.tif", bbox=None
                )
                svfNveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfNveg.tif", bbox=None
                )
                svfSveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfSveg.tif", bbox=None
                )
                svfEveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfEveg.tif", bbox=None
                )
                svfWveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfWveg.tif", bbox=None
                )

                svfaveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfaveg.tif", bbox=None
                )
                svfNaveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfNaveg.tif", bbox=None
                )
                svfSaveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfSaveg.tif", bbox=None
                )
                svfEaveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfEaveg.tif", bbox=None
                )
                svfWaveg, _, _ = common.load_raster(
                    configDict["working_dir"] + "/svfWaveg.tif", bbox=None
                )
        else:
            svfveg = torch.ones((rows, cols)).to(device)
            svfNveg = torch.ones((rows, cols)).to(device)
            svfSveg = torch.ones((rows, cols)).to(device)
            svfEveg = torch.ones((rows, cols)).to(device)
            svfWveg = torch.ones((rows, cols)).to(device)
            svfaveg = torch.ones((rows, cols)).to(device)
            svfNaveg = torch.ones((rows, cols)).to(device)
            svfSaveg = torch.ones((rows, cols)).to(device)
            svfEaveg = torch.ones((rows, cols)).to(device)
            svfWaveg = torch.ones((rows, cols)).to(device)

        # Clean up after SVF loading
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        tmp = svf + svfveg - 1.0
        tmp[tmp < 0.0] = 0.0
        # %matlab crazyness around 0
        svfalfa = torch.arcsin(torch.exp((torch.log((1.0 - tmp)) / 2.0)))

        if standAlone == 0:
            wallheight = torch.from_numpy(
                (
                    gdal.Open(configDict["filepath_wh"])
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
            wallaspect = torch.from_numpy(
                (
                    gdal.Open(configDict["filepath_wa"])
                    .ReadAsArray()
                    .astype(float)
                )
            ).to(device)
        else:
            wallheight, _, _ = common.load_raster(
                configDict["filepath_wh"], bbox=None
            )
            wallheight = torch.tensor(wallheight, device=device)
            wallaspect, _, _ = common.load_raster(
                configDict["filepath_wa"], bbox=None
            )
            wallaspect = torch.tensor(wallaspect, device=device)

        # Metdata
        headernum = 1
        delim = " "
        Twater = []

        metdata = torch.from_numpy(
            np.loadtxt(
                configDict["input_met"], skiprows=headernum, delimiter=delim
            )
        ).to(device)

        location = {"longitude": lon, "latitude": lat, "altitude": alt}
        YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax = (
            Solweig_2015a_metdata_noload(
                metdata, location, int(configDict["utc"])
            )
        )

        DOY = metdata[:, 1]
        hours = metdata[:, 2]
        minu = metdata[:, 3]
        Ta = metdata[:, 11]
        RH = metdata[:, 10]
        radG = metdata[:, 14]
        radD = metdata[:, 21]
        radI = metdata[:, 22]
        P = metdata[:, 12]
        Ws = metdata[:, 9]

        # POIs check
        if configDict["poi_file"] != "":  # usePOI:
            header = (
                "yyyy id   it imin dectime altitude azimuth kdir kdiff kglobal kdown   kup    keast ksouth "
                "kwest knorth ldown   lup    least lsouth lwest  lnorth   Ta      Tg     RH    Esky   Tmrt    "
                "I0     CI   Shadow  SVF_b  SVF_bv  KsideI   PET UTCI  CI_TgG  KsideD  Lside   diffDown    Kside  "
            )
            # poiname = []
            poi_field = configDict[
                "poi_field"
            ]  # self.parameterAsString(parameters, self.POI_FIELD, context)
            if standAlone == 0:

                poi_field = configDict[
                    "woi_field"
                ]  # self.parameterAsStrings(parameters, self.WOI_FIELD, context)
                poisxy, poiname = pointOfInterest(
                    configDict["poi_file"], poi_field, scale, gdal_dsm
                )

            else:
                pois_gdf = gpd.read_file(configDict["poi_file"])
                numfeat = pois_gdf.shape[0]
                poisxy = torch.zeros((numfeat, 3)).to(device) - 999
                for idx, row in pois_gdf.iterrows():
                    y, x = rowcol(
                        dsm_transf,
                        row["geometry"].centroid.x,
                        row["geometry"].centroid.y,
                    )  # TODO: This produce different result since no standalone round coordinates
                    poiname.append(row[configDict["poi_field"]])
                    poisxy[idx, 0] = idx
                    poisxy[idx, 1] = x
                    poisxy[idx, 2] = y

            for k in range(0, poisxy.shape[0]):
                poi_save = []  # torch.zeros((1, 33))
                data_out = (
                    configDict["output_dir"]
                    + "/POI_"
                    + str(poiname[k])
                    + ".txt"
                )
                np.savetxt(
                    data_out,
                    (
                        poi_save.cpu().numpy()
                        if isinstance(poi_save, torch.Tensor)
                        else poi_save
                    ),
                    delimiter=" ",
                    header=header,
                    comments="",
                )

            # Num format for POI output
            numformat = "%d %d %d %d %.5f " + "%.2f " * 35

            # Other PET variables
            sensorheight = param["Wind_Height"]["Value"]["magl"]
            age = param["PET_settings"]["Value"]["Age"]
            mbody = param["PET_settings"]["Value"]["Weight"]
            ht = param["PET_settings"]["Value"]["Height"]
            clo = param["PET_settings"]["Value"]["clo"]
            activity = param["PET_settings"]["Value"]["Activity"]
            sex = param["PET_settings"]["Value"]["Sex"]
        else:
            poisxy = None

        # Posture settings
        if param["Tmrt_params"]["Value"]["posture"] == "Standing":
            Fside = param["Posture"]["Standing"]["Value"]["Fside"]
            Fup = param["Posture"]["Standing"]["Value"]["Fup"]
            height = torch.tensor(
                param["Posture"]["Standing"]["Value"]["height"]
            ).to(device)
            Fcyl = param["Posture"]["Standing"]["Value"]["Fcyl"]
            pos = 1
        else:
            Fside = param["Posture"]["Sitting"]["Value"]["Fside"]
            Fup = param["Posture"]["Sitting"]["Value"]["Fup"]
            height = torch.tensor(
                param["Posture"]["Sitting"]["Value"]["height"]
            ).to(device)
            Fcyl = param["Posture"]["Sitting"]["Value"]["Fcyl"]
            pos = 0

        # Radiative surface influence, Rule of thumb by Schmid et al. (1990).
        first = torch.round(height)
        if first == 0.0:
            first = 1.0
        second = torch.round((height * 20.0))

        if usevegdem == 1:
            # Conifer or deciduous
            if configDict["conifer_bool"]:
                leafon = torch.ones((1, DOY.shape[0])).to(device)
            else:
                leafon = torch.zeros((1, DOY.shape[0])).to(device)
                if (
                    param["Tree_settings"]["Value"]["First_day_leaf"]
                    > param["Tree_settings"]["Value"]["Last_day_leaf"]
                ):
                    leaf_bool = (
                        DOY > param["Tree_settings"]["Value"]["First_day_leaf"]
                    ) | (
                        DOY < param["Tree_settings"]["Value"]["Last_day_leaf"]
                    )
                else:
                    leaf_bool = (
                        DOY > param["Tree_settings"]["Value"]["First_day_leaf"]
                    ) & (
                        DOY < param["Tree_settings"]["Value"]["Last_day_leaf"]
                    )
                leafon[0, leaf_bool] = 1

            # % Vegetation transmittivity of shortwave radiation
            psi = leafon * transVeg
            psi[leafon == 0] = 0.5
            # amaxvalue
            vegmax = vegdsm.max()
            amaxvalue = dsm.max() - dsm.min()
            amaxvalue = torch.maximum(amaxvalue, vegmax)

            # Elevation vegdsms if buildingDEM includes ground heights
            vegdsm = vegdsm + dsm
            vegdsm[vegdsm == dsm] = 0
            vegdsm2 = vegdsm2 + dsm
            vegdsm2[vegdsm2 == dsm] = 0

            # % Bush separation
            bush = torch.logical_not((vegdsm2 * vegdsm)) * vegdsm

            svfbuveg = svf - (1.0 - svfveg) * (
                1.0 - transVeg
            )  # % major bug fixed 20141203
        else:
            psi = leafon * 0.0 + 1.0
            svfbuveg = svf
            bush = torch.zeros([rows, cols]).to(device)
            amaxvalue = 0

        # Clean up vegetation processing tensors
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        # Initialization of maps
        Knight = torch.zeros((rows, cols)).to(device)
        Tgmap1 = torch.zeros((rows, cols)).to(device)
        Tgmap1E = torch.zeros((rows, cols)).to(device)
        Tgmap1S = torch.zeros((rows, cols)).to(device)
        Tgmap1W = torch.zeros((rows, cols)).to(device)
        Tgmap1N = torch.zeros((rows, cols)).to(device)

        # Create building boolean raster from either land cover or height rasters
        if demforbuild == 0:
            buildings = lcgrid.clone()
            buildings[buildings == 7] = 1
            buildings[buildings == 6] = 1
            buildings[buildings == 5] = 1
            buildings[buildings == 4] = 1
            buildings[buildings == 3] = 1
            buildings[buildings == 2] = 0
        else:
            buildings = dsm - dem
            buildings[buildings < 2.0] = 1.0
            buildings[buildings >= 2.0] = 0.0

        if int(configDict["savebuild"]) == 1:
            if standAlone == 0:
                saveraster(
                    gdal_dsm,
                    configDict["output_dir"] + "/buildings.tif",
                    buildings.detach().cpu().numpy(),
                )
            else:
                common.save_raster(
                    configDict["output_dir"] + "/buildings.tif",
                    buildings.detach().cpu().numpy(),
                    dsm_transf,
                    dsm_crs,
                )

        # Import shadow matrices (Anisotropic sky)
        anisotropic_sky = int(configDict["aniso"])
        if anisotropic_sky == 1:  # UseAniso
            data = torch.load(
                configDict["input_aniso"],
                map_location=device,
                weights_only=True,
            )
            shmat = data["shadowmat"]
            vegshmat = data["vegshadowmat"]
            vbshvegshmat = data["vbshmat"]
            if usevegdem == 1:
                diffsh = torch.zeros((rows, cols, shmat.shape[2])).to(device)
                for i in range(0, shmat.shape[2]):
                    diffsh[:, :, i] = shmat[:, :, i] - (
                        1 - vegshmat[:, :, i]
                    ) * (
                        1 - transVeg
                    )  # changes in psi not implemented yet
            else:
                diffsh = shmat

            # Estimate number of patches based on shadow matrices
            if shmat.shape[2] == 145:
                patch_option = 1  # patch_option = 1 # 145 patches
            elif shmat.shape[2] == 153:
                patch_option = 2  # patch_option = 2 # 153 patches
            elif shmat.shape[2] == 306:
                patch_option = 3  # patch_option = 3 # 306 patches
            elif shmat.shape[2] == 612:
                patch_option = 4  # patch_option = 4 # 612 patches

            # asvf to calculate sunlit and shaded patches
            asvf = torch.arccos(torch.sqrt(svf))

            # Empty array for steradians
            steradians = torch.zeros((shmat.shape[2])).to(device)
        else:
            # anisotropic_sky = 0
            diffsh = None
            shmat = None
            vegshmat = None
            vbshvegshmat = None
            asvf = None
            patch_option = 0
            steradians = 0
        shadow = torch.zeros_like(dsm).to(device)

        # % Ts parameterisation maps
        if landcover == 1.0:
            # Get land cover properties for Tg wave (land cover scheme based on Bogren et al. 2000, explained in Lindberg et al., 2008 and Lindberg, Onomura & Grimmond, 2016)
            [
                TgK,
                Tstart,
                alb_grid,
                emis_grid,
                TgK_wall,
                Tstart_wall,
                TmaxLST,
                TmaxLST_wall,
            ] = Tgmaps_v1(lcgrid.clone(), param)
        else:
            TgK = Knight + param["Ts_deg"]["Value"]["Cobble_stone_2014a"]
            Tstart = Knight - param["Tstart"]["Value"]["Cobble_stone_2014a"]
            TmaxLST = param["TmaxLST"]["Value"]["Cobble_stone_2014a"]
            alb_grid = (
                Knight
                + param["Albedo"]["Effective"]["Value"]["Cobble_stone_2014a"]
            )
            emis_grid = (
                Knight + param["Emissivity"]["Value"]["Cobble_stone_2014a"]
            )
            TgK_wall = param["Ts_deg"]["Value"]["Walls"]
            Tstart_wall = param["Tstart"]["Value"]["Walls"]
            TmaxLST_wall = param["TmaxLST"]["Value"]["Walls"]

        # Parameterization for the 2026 ground scheme
        groundSurface = int(configDict["groundmodel"])

        if groundSurface == 1:
            # Initiate the maps if the surface temperature is available
            if configDict["input_surf"] != "":
                surfData = pd.read_csv(configDict["input_surf"])
                Tg = surfData["Tg"]
                Tm = torch.mean(surfData["Tg"])
                (
                    _,
                    _,
                    Rn,
                    Rn_past,
                    G,
                    cap_grid,
                    diff_grid,
                    a1_grid,
                    a2_grid,
                    a3_grid,
                ) = initiate_groundScheme(
                    lcgrid.copy(), param, DOY[0], Ta, location, device
                )
            else:
                (
                    Tg,
                    Tm,
                    Rn,
                    Rn_past,
                    G,
                    cap_grid,
                    diff_grid,
                    a1_grid,
                    a2_grid,
                    a3_grid,
                ) = initiate_groundScheme(
                    lcgrid.clone(), param, DOY[0], Ta, location, device
                )
        else:
            pass

        # Replace the ground view factors with integration of solid angles
        outgoingLW = int(configDict["outgoinglongwave"])

        # Import data for wall temperature parameterization TODO: fix for standalone
        wallScheme = int(configDict["wallscheme"])
        if wallScheme == 1:
            wallData = torch.load(
                configDict["input_wall"],
                map_location=device,
                weights_only=True,
            )
            voxelMaps = wallData["voxelId"]
            voxelTable = wallData["voxelTable"]
            # Get wall type from standalone
            if standAlone == 1:
                wall_type_standalone = {
                    "Brick_wall": "100",
                    "Concrete_wall": "101",
                    "Wood_wall": "102",
                }
                wall_type = wall_type_standalone[configDict["walltype"]]
            else:
                # Get wall type set in GUI
                wall_type = configDict[
                    "walltype"
                ]  # str(100 + int(self.parameterAsString(parameters, self.WALL_TYPE, context))) #TODO

            # Calculate wall height for wall scheme, i.e. include corners (thicker walls)
            walls_scheme = wa.findwalls_sp(
                dsm,
                2,
                device,
                torch.tensor([[1, 1, 1], [1, 0, 1], [1, 1, 1]]).to(device),
            )
            # Calculate wall aspect for wall scheme, i.e. include corners (thicker walls)
            dirwalls_scheme = wa.filter1Goodwin_as_aspect_v3(
                walls_scheme.copy(),
                scale,
                dsm,
                feedback,
                100.0 / 180.0,
                device,
            )

            # Used in wall temperature parameterization scheme
            first_timestep = (
                pd.to_datetime(int(YYYY[0][0].item()), format="%Y")
                + pd.to_timedelta(DOY[0].item() - 1, unit="d")
                + pd.to_timedelta(hours[0].item(), unit="h")
                + pd.to_timedelta(minu[0].item(), unit="m")
            )
            second_timestep = (
                pd.to_datetime(int(YYYY[0][1].item()), format="%Y")
                + pd.to_timedelta(DOY[1].item() - 1, unit="d")
                + pd.to_timedelta(hours[1].item(), unit="h")
                + pd.to_timedelta(minu[1].item(), unit="m")
            )

            timeStep = (second_timestep - first_timestep).seconds

            # Load voxelTable as Pandas DataFrame
            voxelTable, dirwalls_scheme = load_walls(
                voxelTable,
                param,
                wall_type,
                dirwalls_scheme,
                Ta[0],
                timeStep,
                albedo_b,
                ewall,
                alb_grid,
                landcover,
                lcgrid,
                dsm,
            )

            # Use wall of interest
            woi_file = configDict["woi_file"]
            if woi_file:
                # (dsm_minx, dsm_x_size, dsm_x_rotation, dsm_miny, dsm_y_rotation, dsm_y_size) = gdal_dsm.GetGeoTransform() #TODO: fix for standalone
                if standAlone == 0:
                    woi_field = configDict[
                        "woi_field"
                    ]  # self.parameterAsStrings(parameters, self.WOI_FIELD, context)
                    woisxy, woiname = pointOfInterest(
                        configDict["woi_file"], woi_field, scale, gdal_dsm
                    )
                else:
                    pois_gdf = gpd.read_file(configDict["poi_file"])
                    numfeat = pois_gdf.shape[0]
                    poisxy = torch.zeros((numfeat, 3)).to(device) - 999
                    for idx, row in pois_gdf.iterrows():
                        y, x = rowcol(
                            dsm_transf,
                            row["geometry"].centroid.x,
                            row["geometry"].centroid.y,
                        )  # TODO: This produce different result since no standalone round coordinates
                        poiname.append(row[configDict["poi_field"]])
                        poisxy[idx, 0] = idx
                        poisxy[idx, 1] = x
                        poisxy[idx, 2] = y

            # Create pandas datetime object to be used when createing an xarray DataSet where wall temperatures/radiation is stored and eventually saved as a NetCDf
            if configDict["wallnetcdf"] == 1:
                met_for_xarray = (
                    pd.to_datetime(YYYY[0][:], format="%Y")
                    + pd.to_timedelta(DOY - 1, unit="d")
                    + pd.to_timedelta(hours, unit="h")
                    + pd.to_timedelta(minu, unit="m")
                )
        else:
            wallScheme = 0
            voxelMaps = 0
            voxelTable = 0
            timeStep = 0
            # thermal_effusivity = 0
            walls_scheme = torch.ones((rows, cols)).to(device) * 10.0
            dirwalls_scheme = torch.ones((rows, cols)).to(device) * 10.0

        # Clean up wall scheme setup tensors
        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        # Initialisation of time related variables
        if Ta.__len__() == 1:
            timestepdec = 0
        else:
            timestepdec = dectime[1] - dectime[0]
        timeadd = 0.0
        # timeaddE = 0.
        # timeaddS = 0.
        # timeaddW = 0.
        # timeaddN = 0.
        firstdaytime = 1.0

        # Save hemispheric image
        if anisotropic_sky == 1:
            if not poisxy is None:
                patch_characteristics = hemispheric_image(
                    poisxy,
                    shmat,
                    vegshmat,
                    vbshvegshmat,
                    voxelMaps,
                    wallScheme,
                )

        # If metfile starts at night
        CI = 1.0

        # Main loop
        tmrtplot = torch.zeros((rows, cols)).to(device)

        # Initiate array for I0 values
        if torch.unique(DOY).shape[0] > 1:
            unique_days = torch.unique(DOY)
            first_unique_day = DOY[DOY == unique_days[0]]
            I0_array = torch.zeros((first_unique_day.shape[0])).to(device)
        else:
            first_unique_day = DOY.clone()
            I0_array = torch.zeros((DOY.shape[0])).to(device)

        if standAlone == 1:
            progress = tqdm(total=Ta.__len__())
        else:
            progress = None

        for i in torch.arange(0, Ta.__len__()):
            if feedback is not None:
                feedback.setProgress(
                    int(i * (100.0 / Ta.__len__()))
                )  # move progressbar forward
                if feedback.isCanceled():
                    feedback.setProgressText("Calculation cancelled")
                    break
            elif progress is not None:
                progress.update(1)

            # Daily water body temperature
            if landcover == 1:
                if ((dectime[i] - torch.floor(dectime[i]))) == 0 or (i == 0):
                    Twater = torch.mean(Ta[jday[0] == torch.floor(dectime[i])])
            # Nocturnal cloudfraction from Offerle et al. 2003
            if (dectime[i] - torch.floor(dectime[i])) == 0:
                daylines = torch.where(torch.floor(dectime) == dectime[i])
                if daylines.__len__() > 1:
                    alt = altitude[0][daylines]
                    alt2 = torch.where(alt > 1)
                    rise = alt2[0][0]
                    [_, CI, _, _, _] = clearnessindex_2013b(
                        zen[0, i + rise + 1],
                        jday[0, i + rise + 1],
                        Ta[i + rise + 1],
                        RH[i + rise + 1] / 100.0,
                        radG[i + rise + 1],
                        location,
                        P[i + rise + 1],
                    )
                    if (CI > 1.0) or (CI == torch.inf):
                        CI = 1.0
                else:
                    CI = 1.0

            # Timestep of the simulation used in the ground scheme calculation
            first_timestep = (
                pd.to_datetime(int(YYYY[0][0].item()), format="%Y")
                + pd.to_timedelta(DOY[0].item() - 1, unit="d")
                + pd.to_timedelta(hours[0].item(), unit="h")
                + pd.to_timedelta(minu[0].item(), unit="m")
            )
            second_timestep = (
                pd.to_datetime(int(YYYY[0][1].item()), format="%Y")
                + pd.to_timedelta(DOY[1].item() - 1, unit="d")
                + pd.to_timedelta(hours[1].item(), unit="h")
                + pd.to_timedelta(minu[1].item(), unit="m")
            )

            timeStep = (second_timestep - first_timestep).seconds

            (
                Tmrt,
                Kdown,
                Kup,
                Ldown,
                Lup,
                Tg,
                ea,
                esky,
                I0,
                CI,
                shadow,
                firstdaytime,
                timestepdec,
                timeadd,
                Tgmap1,
                Tgmap1E,
                Tgmap1S,
                Tgmap1W,
                Tgmap1N,
                Keast,
                Ksouth,
                Kwest,
                Knorth,
                Least,
                Lsouth,
                Lwest,
                Lnorth,
                KsideI,
                radIout,
                radDout,
                Lside,
                Lsky_patch_characteristics,
                CI_TgG,
                KsideD,
                dRad,
                Kside,
                steradians,
                voxelTable,
                Rn,
                Rn_past,
                Tm,
                G,
            ) = so.Solweig_2026a_calc(
                i,
                dsm,
                scale,
                rows,
                cols,
                svf,
                svfN,
                svfW,
                svfE,
                svfS,
                svfveg,
                svfNveg,
                svfEveg,
                svfSveg,
                svfWveg,
                svfaveg,
                svfEaveg,
                svfSaveg,
                svfWaveg,
                svfNaveg,
                vegdsm,
                vegdsm2,
                albedo_b,
                absK,
                absL,
                ewall,
                Fside,
                Fup,
                Fcyl,
                altitude[0][i],
                azimuth[0][i],
                zen[0][i],
                jday[0][i],
                usevegdem,
                onlyglobal,
                buildings,
                location,
                psi[0][i],
                landcover,
                lcgrid,
                dectime[i],
                altmax[0][i],
                wallaspect,
                wallheight,
                cyl,
                elvis,
                Ta[i],
                RH[i],
                radG[i],
                radD[i],
                radI[i],
                P[i],
                amaxvalue,
                bush,
                Twater,
                TgK,
                Tstart,
                alb_grid,
                emis_grid,
                TgK_wall,
                Tstart_wall,
                TmaxLST,
                TmaxLST_wall,
                first,
                second,
                svfalfa,
                svfbuveg,
                firstdaytime,
                timeadd,
                timestepdec,
                Tgmap1,
                Tgmap1E,
                Tgmap1S,
                Tgmap1W,
                Tgmap1N,
                CI,
                diffsh,
                shmat,
                vegshmat,
                vbshvegshmat,
                anisotropic_sky,
                asvf,
                patch_option,
                voxelMaps,
                voxelTable,
                Ws[i],
                wallScheme,
                timeStep,
                steradians,
                walls_scheme,
                dirwalls_scheme,
                groundSurface,
                outgoingLW,
                Tg,
                Rn,
                Rn_past,
                G,
                Tm,
                cap_grid,
                diff_grid,
                a1_grid,
                a2_grid,
                a3_grid,
                shadow,
            )

            # Save I0 for I0 vs. Kdown output plot to check if UTC is off
            if i < first_unique_day.shape[0]:
                I0_array[i] = I0
            elif i == first_unique_day.shape[0]:
                # Output I0 vs. Kglobal plot
                radG_for_plot = radG[DOY == first_unique_day[0]]
                # hours_for_plot = hours[DOY == first_unique_day[0]]
                dectime_for_plot = dectime[DOY == first_unique_day[0]]
                fig, ax = plt.subplots()
                ax.plot(dectime_for_plot, I0_array, label="I0")
                ax.plot(dectime_for_plot, radG_for_plot, label="Kglobal")
                ax.set_ylabel("Shortwave radiation [$Wm^{-2}$]")
                ax.set_xlabel("Decimal time")
                ax.set_title("UTC" + str(configDict["utc"]))
                ax.legend()
                fig.savefig(
                    configDict["output_dir"] + "/metCheck.png", dpi=150
                )

            tmrtplot = tmrtplot + Tmrt

            if altitude[0][i] > 0:
                w = "D"
            else:
                w = "N"

            # Write to POIs
            if not poisxy is None:
                for k in range(0, poisxy.shape[0]):
                    poi_save = torch.zeros((1, 40)).to(device)
                    poi_save[0, 0] = YYYY[0][i]
                    poi_save[0, 1] = jday[0][i]
                    poi_save[0, 2] = hours[i]
                    poi_save[0, 3] = minu[i]
                    poi_save[0, 4] = dectime[i]
                    poi_save[0, 5] = altitude[0][i]
                    poi_save[0, 6] = azimuth[0][i]
                    poi_save[0, 7] = radIout
                    poi_save[0, 8] = radDout
                    poi_save[0, 9] = radG[i]
                    poi_save[0, 10] = Kdown[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 11] = Kup[int(poisxy[k, 2]), int(poisxy[k, 1])]
                    poi_save[0, 12] = Keast[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 13] = Ksouth[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 14] = Kwest[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 15] = Knorth[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 16] = Ldown[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 17] = Lup[int(poisxy[k, 2]), int(poisxy[k, 1])]
                    poi_save[0, 18] = Least[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 19] = Lsouth[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 20] = Lwest[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 21] = Lnorth[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 22] = Ta[i]
                    poi_save[0, 23] = Tg[int(poisxy[k, 2]), int(poisxy[k, 1])]
                    poi_save[0, 24] = RH[i]
                    poi_save[0, 25] = esky
                    poi_save[0, 26] = Tmrt[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 27] = I0
                    poi_save[0, 28] = CI
                    poi_save[0, 29] = shadow[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 30] = svf[int(poisxy[k, 2]), int(poisxy[k, 1])]
                    poi_save[0, 31] = svfbuveg[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 32] = KsideI[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    # Recalculating wind speed based on powerlaw
                    WsPET = (1.1 / sensorheight) ** 0.2 * Ws[i]
                    WsUTCI = (10.0 / sensorheight) ** 0.2 * Ws[i]
                    resultPET = p._PET(
                        Ta[i],
                        RH[i],
                        Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])],
                        WsPET,
                        mbody,
                        age,
                        ht,
                        activity,
                        clo,
                        sex,
                    )
                    poi_save[0, 33] = resultPET
                    resultUTCI = utci.utci_calculator(
                        Ta[i],
                        RH[i],
                        Tmrt[int(poisxy[k, 2]), int(poisxy[k, 1])],
                        WsUTCI,
                    )
                    poi_save[0, 34] = resultUTCI
                    poi_save[0, 35] = CI_TgG
                    poi_save[0, 36] = KsideD[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 37] = Lside[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 38] = dRad[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    poi_save[0, 39] = Kside[
                        int(poisxy[k, 2]), int(poisxy[k, 1])
                    ]
                    data_out = (
                        configDict["output_dir"]
                        + "/POI_"
                        + str(poiname[k])
                        + ".txt"
                    )
                    # f_handle = file(data_out, 'a')
                    f_handle = open(data_out, "ab")
                    np.savetxt(
                        f_handle,
                        (
                            poi_save.cpu().numpy()
                            if isinstance(poi_save, torch.Tensor)
                            else poi_save
                        ),
                        fmt=numformat,
                    )
                    f_handle.close()

            # Clean up POI tensor after writing
            if device.type == "cuda":
                torch.cuda.empty_cache()
            elif device.type == "xpu":
                torch.xpu.empty_cache()
            # If wall temperature parameterization scheme is in use
            if (
                configDict["wallscheme"] == 1
            ):  # folderWallScheme: TODO: Fix for standalone
                # Store wall data for output
                if not woisxy is None:
                    for k in range(0, woisxy.shape[0]):
                        temp_wall = voxelTable.loc[
                            (
                                (voxelTable["ypos"] == woisxy[k, 2])
                                & (voxelTable["xpos"] == woisxy[k, 1])
                            ),
                            "wallTemperature",
                        ].to_numpy()
                        K_in = voxelTable.loc[
                            (
                                (voxelTable["ypos"] == woisxy[k, 2])
                                & (voxelTable["xpos"] == woisxy[k, 1])
                            ),
                            "K_in",
                        ].to_numpy()
                        L_in = voxelTable.loc[
                            (
                                (voxelTable["ypos"] == woisxy[k, 2])
                                & (voxelTable["xpos"] == woisxy[k, 1])
                            ),
                            "L_in",
                        ].to_numpy()
                        wallShade = voxelTable.loc[
                            (
                                (voxelTable["ypos"] == woisxy[k, 2])
                                & (voxelTable["xpos"] == woisxy[k, 1])
                            ),
                            "wallShade",
                        ].to_numpy()
                        temp_all = torch.concatenate(
                            [temp_wall, K_in, L_in, wallShade]
                        ).to(device)
                        # temp_all = torch.concatenate([temp_wall])
                        # wall_data = torch.zeros((1, 7 + temp_wall.shape[0]))
                        wall_data = torch.zeros((1, 7 + temp_all.shape[0])).to(
                            device
                        )
                        # Part of file name (wallid), i.e. WOI_wallid.txt
                        data_out = (
                            configDict["output_dir"]
                            + "/WOI_"
                            + str(woiname[k])
                            + ".txt"
                        )
                        if i == 0:
                            # Output file header
                            # header = 'yyyy id   it imin dectime Ta  SVF Ts'
                            header = (
                                "yyyy id   it imin dectime Ta  SVF"
                                + " Ts" * temp_wall.shape[0]
                                + " Kin" * K_in.shape[0]
                                + " Lin" * L_in.shape[0]
                                + " shade" * wallShade.shape[0]
                            )
                            # Part of file name (wallid), i.e. WOI_wallid.txt
                            # woiname = voxelTable.loc[((voxelTable['ypos'] == woisxy[k, 2]) & (voxelTable['xpos'] == woisxy[k, 1])), 'wallId'].to_numpy()[0]
                            woi_save = []  #
                            np.savetxt(
                                data_out,
                                woi_save,
                                delimiter=" ",
                                header=header,
                                comments="",
                            )
                        # Fill wall_data with variables
                        wall_data[0, 0] = YYYY[0][i]
                        wall_data[0, 1] = jday[0][i]
                        wall_data[0, 2] = hours[i]
                        wall_data[0, 3] = minu[i]
                        wall_data[0, 4] = dectime[i]
                        wall_data[0, 5] = Ta[i]
                        wall_data[0, 6] = svf[
                            int(woisxy[k, 2]), int(woisxy[k, 1])
                        ]
                        wall_data[0, 7:] = temp_all

                        # Num format for output file data
                        woi_numformat = (
                            "%d %d %d %d %.5f %.2f %.2f"
                            + " %.2f" * temp_all.shape[0]
                        )
                        # Open file, add data, save
                        f_handle = open(data_out, "ab")
                        np.savetxt(
                            f_handle,
                            (
                                wall_data.cpu().numpy()
                                if isinstance(wall_data, torch.Tensor)
                                else wall_data
                            ),
                            fmt=woi_numformat,
                        )
                        f_handle.close()
                        del wall_data, temp_all  # Clean up wall data tensors
                        if device.type == "cuda":
                            torch.cuda.empty_cache()
                        elif device.type == "xpu":
                            torch.xpu.empty_cache()
                # Save wall temperature/radiation as NetCDF TODO: fix for standAlone?
                if configDict["wallnetcdf"] == "1":  # wallNetCDF:
                    netcdf_output = configDict["outputDir"] + "/walls.nc"
                    walls_as_netcdf(
                        voxelTable,
                        rows,
                        cols,
                        met_for_xarray,
                        i,
                        dsm,
                        configDict["filepath_dsm"],
                        netcdf_output,
                    )

            if hours[i] < 10:
                XH = "0"
            else:
                XH = ""
            if minu[i] < 10:
                XM = "0"
            else:
                XM = ""

            time_code = (
                str(int(YYYY[0, i]))
                + "_"
                + str(int(DOY[i]))
                + "_"
                + XH
                + str(int(hours[i]))
                + XM
                + str(int(minu[i]))
                + w
            )

            if configDict["outputtmrt"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Tmrt_"
                        + time_code
                        + ".tif",
                        Tmrt.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Tmrt_"
                        + time_code
                        + ".tif",
                        Tmrt.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )
            if configDict["outputkup"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Kup_"
                        + time_code
                        + ".tif",
                        Kup.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Kup_"
                        + time_code
                        + ".tif",
                        Kup.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )
            if configDict["outputkdown"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Kdown_"
                        + time_code
                        + ".tif",
                        Kdown.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Kdown_"
                        + time_code
                        + ".tif",
                        Kdown.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )
            if configDict["outputlup"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Lup_"
                        + time_code
                        + ".tif",
                        Lup.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Lup_"
                        + time_code
                        + ".tif",
                        Lup.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )
            if configDict["outputldown"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Ldown_"
                        + time_code
                        + ".tif",
                        Ldown.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Ldown_"
                        + time_code
                        + ".tif",
                        Ldown.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )
            if configDict["outputsh"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Shadow_"
                        + time_code
                        + ".tif",
                        shadow.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Shadow_"
                        + time_code
                        + ".tif",
                        shadow.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )

            if configDict["outputkdiff"] == "1":
                if standAlone == 0:
                    saveraster(
                        gdal_dsm,
                        configDict["output_dir"]
                        + "/Kdiff_"
                        + time_code
                        + ".tif",
                        dRad.detach().cpu().numpy(),
                    )
                else:
                    common.save_raster(
                        configDict["output_dir"]
                        + "/Kdiff_"
                        + time_code
                        + ".tif",
                        dRad.detach().cpu().numpy(),
                        dsm_transf,
                        dsm_crs,
                    )

            # Clean up iteration tensors after output
            if device.type == "cuda" or device.type == "xpu":
                del (
                    Tmrt,
                    Kdown,
                    Kup,
                    Ldown,
                    Lup,
                    Keast,
                    Ksouth,
                    Kwest,
                    Knorth,
                    Least,
                    Lsouth,
                    Lwest,
                    Lnorth,
                    KsideI,
                    radIout,
                    radDout,
                    Lside,
                    KsideD,
                    dRad,
                    Kside,
                )
                if device.type == "cuda":
                    torch.cuda.empty_cache()
                elif device.type == "xpu":
                    torch.xpu.empty_cache()

            # Sky view image of patches
            if (anisotropic_sky == 1) & (i == 0) & (not poisxy is None):
                for k in range(poisxy.shape[0]):
                    Lsky_patch_characteristics[:, 2] = patch_characteristics[
                        :, k
                    ]
                    skyviewimage_out = (
                        configDict["output_dir"]
                        + "/POI_"
                        + str(poiname[k])
                        + ".png"
                    )
                    PolarBarPlot(
                        Lsky_patch_characteristics,
                        altitude[0][i],
                        azimuth[0][i],
                        "Hemisphere partitioning",
                        skyviewimage_out,
                        0,
                        5,
                        0,
                    )

        # Save files for Tree Planter
        if configDict["outputtreeplanter"] == "1":  # outputTreeplanter:
            if feedback is not None:
                feedback.setProgressText("Saving files for Tree Planter tool")
            # Save DSM
            copyfile(
                configDict["filepath_dsm"],
                configDict["output_dir"] + "/DSM.tif",
            )

            # Save CDSM
            if usevegdem == 1:
                copyfile(
                    configDict["filepath_cdsm"],
                    configDict["output_dir"] + "/CDSM.tif",
                )

            albedo_g = param["Albedo"]["Effective"]["Value"][
                "Cobble_stone_2014a"
            ]
            eground = param["Emissivity"]["Value"]["Cobble_stone_2014a"]

            # Saving settings from SOLWEIG for SOLWEIG1D in TreePlanter
            settingsHeader = "UTC, posture, onlyglobal, landcover, anisotropic, cylinder, albedo_walls, albedo_ground, emissivity_walls, emissivity_ground, absK, absL, elevation, patch_option"
            settingsFmt = (
                "%i",
                "%i",
                "%i",
                "%i",
                "%i",
                "%i",
                "%1.2f",
                "%1.2f",
                "%1.2f",
                "%1.2f",
                "%1.2f",
                "%1.2f",
                "%1.2f",
                "%i",
            )
            settingsData = torch.tensor(
                [
                    [
                        int(configDict["utc"]),
                        pos,
                        onlyglobal,
                        landcover,
                        anisotropic_sky,
                        cyl,
                        albedo_b,
                        albedo_g,
                        ewall,
                        eground,
                        absK,
                        absL,
                        alt,
                        patch_option,
                    ]
                ]
            ).to(device)

            np.savetxt(
                configDict["output_dir"] + "/treeplantersettings.txt",
                (
                    settingsData.cpu().numpy()
                    if isinstance(settingsData, torch.Tensor)
                    else settingsData
                ),
                fmt=settingsFmt,
                header=settingsHeader,
                delimiter=" ",
            )

        # Copying met file for SpatialTC
        copyfile(
            configDict["input_met"],
            configDict["output_dir"] + "/metforcing.txt",
        )

        tmrtplot = (
            tmrtplot / Ta.__len__()
        )  # fix average Tmrt instead of sum, 20191022
        if standAlone == 0:
            saveraster(
                gdal_dsm,
                configDict["output_dir"] + "/Tmrt_average.tif",
                tmrtplot.detach().cpu().numpy(),
            )
        else:
            common.save_raster(
                configDict["output_dir"] + "/Tmrt_average.tif",
                tmrtplot,
                dsm_transf,
                dsm_crs,
            )

        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
