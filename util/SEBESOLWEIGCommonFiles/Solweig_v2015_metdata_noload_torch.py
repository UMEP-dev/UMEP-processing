from __future__ import absolute_import

# from importdata import importdata
from . import sun_position_torch as sp

# import sun_position as sp
import datetime
import calendar

try:
    import torch
except:
    pass


def Solweig_2015a_metdata_noload(inputdata, location, UTC):
    """
    This function is used to process the input meteorological file.
    It also calculates Sun position based on the time specified in the met-file

    :param inputdata:
    :param location:
    :param UTC:
    :return:
    """

    met = inputdata
    device = (
        met.device if isinstance(met, torch.Tensor) else torch.device("cpu")
    )
    data_len = len(met[:, 0])
    dectime = torch.tensor(
        met[:, 1] + met[:, 2] / 24 + met[:, 3] / (60 * 24.0)
    )
    dectimemin = met[:, 3] / (60 * 24.0)
    if data_len == 1:
        halftimestepdec = 0
    else:
        halftimestepdec = (dectime[1] - dectime[0]) / 2.0
    time = dict()
    time["sec"] = 0
    time["UTC"] = UTC
    sunmaximum = 0.0
    leafon1 = 97  # TODO this should change
    leafoff1 = 300  # TODO this should change

    # initialize matrices
    altitude = torch.empty(size=(1, data_len), device=device)
    azimuth = torch.empty(size=(1, data_len), device=device)
    zen = torch.empty(size=(1, data_len), device=device)
    jday = torch.empty(size=(1, data_len), device=device)
    YYYY = torch.empty(size=(1, data_len), device=device)
    leafon = torch.empty(size=(1, data_len), device=device)
    altmax = torch.empty(size=(1, data_len), device=device)

    sunmax = dict()

    for i, row in enumerate(met[:, 0]):
        if met[i, 1] == 221:
            test = 4
        YMD = datetime.datetime(int(met[i, 0]), 1, 1) + datetime.timedelta(
            int(met[i, 1]) - 1
        )
        # Finding maximum altitude in 15 min intervals (20141027)
        if (i == 0) or (
            torch.remainder(dectime[i], torch.floor(dectime[i])) == 0
        ):
            fifteen = 0.0
            sunmaximum = -90.0
            sunmax["zenith"] = 90.0
            while sunmaximum <= 90.0 - sunmax["zenith"]:
                sunmaximum = 90.0 - sunmax["zenith"]
                fifteen = fifteen + 15.0 / 1440.0
                HM = datetime.timedelta(days=(60 * 10) / 1440.0 + fifteen)
                YMDHM = YMD + HM
                time["year"] = YMDHM.year
                time["month"] = YMDHM.month
                time["day"] = YMDHM.day
                time["hour"] = YMDHM.hour
                time["min"] = YMDHM.minute
                sunmax = sp.sun_position(time, location)
        altmax[0, i] = sunmaximum

        half = datetime.timedelta(days=float(halftimestepdec))
        H = datetime.timedelta(hours=int(met[i, 2]))
        M = datetime.timedelta(minutes=int(met[i, 3]))
        YMDHM = YMD + H + M - half
        time["year"] = YMDHM.year
        time["month"] = YMDHM.month
        time["day"] = YMDHM.day
        time["hour"] = YMDHM.hour
        time["min"] = YMDHM.minute
        sun = sp.sun_position(time, location)
        sun_zenith = sun["zenith"]
        sun_azimuth = sun["azimuth"]

        if (sun_zenith > 89.0) & (sun_zenith <= 90.0):
            sun_zenith = 89.0

        altitude[0, i] = 90.0 - sun_zenith
        zen[0, i] = sun_zenith * (torch.pi / 180.0)
        azimuth[0, i] = sun_azimuth

        YYYY[0, i] = met[i, 0]
        doy = YMD.timetuple().tm_yday
        jday[0, i] = doy
        if (doy > leafon1) | (doy < leafoff1):
            leafon[0, i] = 1
        else:
            leafon[0, i] = 0

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()

    return YYYY, altitude, azimuth, zen, jday, leafon, dectime, altmax
