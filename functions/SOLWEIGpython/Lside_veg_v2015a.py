from __future__ import absolute_import
import numpy as np
from .Lvikt_veg import Lvikt_veg
import torch


def Lside_veg_v2015a(
    svfS,
    svfW,
    svfN,
    svfE,
    svfEveg,
    svfSveg,
    svfWveg,
    svfNveg,
    svfEaveg,
    svfSaveg,
    svfWaveg,
    svfNaveg,
    azimuth,
    altitude,
    Ta,
    Tw,
    SBC,
    ewall,
    Ldown,
    esky,
    t,
    F_sh,
    CI,
    LupE,
    LupS,
    LupW,
    LupN,
):
    # Load device
    device = (
        Ldown.device
        if isinstance(Ldown, torch.Tensor)
        else (
            Ta.device
            if isinstance(Ta, torch.Tensor)
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )
    )
    # This m-file is the current one that estimates L from the four cardinal points 20100414

    # Building height angle from svf
    svfalfaE = torch.arcsin(torch.exp((torch.log(1 - svfE)) / 2))
    svfalfaS = torch.arcsin(torch.exp((torch.log(1 - svfS)) / 2))
    svfalfaW = torch.arcsin(torch.exp((torch.log(1 - svfW)) / 2))
    svfalfaN = torch.arcsin(torch.exp((torch.log(1 - svfN)) / 2))

    vikttot = 4.4897
    aziW = azimuth + t
    aziN = azimuth - 90 + t
    aziE = azimuth - 180 + t
    aziS = azimuth - 270 + t

    F_sh = 2 * F_sh - 1  # (cylindric_wedge scaled 0-1)

    c = 1 - CI
    Lsky_allsky = esky * SBC * ((Ta + 273.15) ** 4) * (1 - c) + c * SBC * (
        (Ta + 273.15) ** 4
    )

    ## Least
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(svfE, svfEveg, svfEaveg, vikttot)

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaE)
        betaB = torch.arctan(torch.tan((svfalfaE) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180 - t)) and (azimuth <= (360 - t)):
            Lwallsun = (
                SBC
                * ewall
                * ((Ta + 273.15 + Tw * torch.sin(aziE * (torch.pi / 180))) ** 4)
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    Lsky = ((svfE + svfEveg - 1) * Lsky_allsky) * viktsky * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
    Lground = LupE * 0.5
    Lrefl = (Ldown + LupE) * (viktrefl) * (1 - ewall) * 0.5
    Least = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lsouth
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(svfS, svfSveg, svfSaveg, vikttot)

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaS)
        betaB = torch.arctan(torch.tan((svfalfaS) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaS)*(1+F_sh))
        if (azimuth <= (90 - t)) or (azimuth > (270 - t)):
            Lwallsun = (
                SBC
                * ewall
                * ((Ta + 273.15 + Tw * torch.sin(aziS * (torch.pi / 180))) ** 4)
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    Lsky = ((svfS + svfSveg - 1) * Lsky_allsky) * viktsky * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
    Lground = LupS * 0.5
    Lrefl = (Ldown + LupS) * (viktrefl) * (1 - ewall) * 0.5
    Lsouth = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lwest
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(svfW, svfWveg, svfWaveg, vikttot)

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaW)
        betaB = torch.arctan(torch.tan((svfalfaW) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaW)*(1+F_sh))
        if (azimuth > (360 - t)) or (azimuth <= (180 - t)):
            Lwallsun = (
                SBC
                * ewall
                * ((Ta + 273.15 + Tw * torch.sin(aziW * (torch.pi / 180))) ** 4)
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    Lsky = ((svfW + svfWveg - 1) * Lsky_allsky) * viktsky * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
    Lground = LupW * 0.5
    Lrefl = (Ldown + LupW) * (viktrefl) * (1 - ewall) * 0.5
    Lwest = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lnorth
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(svfN, svfNveg, svfNaveg, vikttot)

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaN)
        betaB = torch.arctan(torch.tan((svfalfaN) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaN)*(1+F_sh))
        if (azimuth > (90 - t)) and (azimuth <= (270 - t)):
            Lwallsun = (
                SBC
                * ewall
                * ((Ta + 273.15 + Tw * torch.sin(aziN * (torch.pi / 180))) ** 4)
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    Lsky = ((svfN + svfNveg - 1) * Lsky_allsky) * viktsky * 0.5
    Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
    Lground = LupN * 0.5
    Lrefl = (Ldown + LupN) * (viktrefl) * (1 - ewall) * 0.5
    Lnorth = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    return Least, Lsouth, Lwest, Lnorth
