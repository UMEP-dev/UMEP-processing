from __future__ import absolute_import
from .Lvikt_veg_torch import Lvikt_veg

try:
    import torch

except:
    pass


def Lside_veg_v2022a(
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
    anisotropic_longwave,
):

    # This m-file is the current one that estimates L from the four cardinal points 20100414
    device = (
        Ldown.device
        if isinstance(Ldown, torch.Tensor)
        else (
            Ta.device
            if isinstance(Ta, torch.Tensor)
            else torch.device(
                "cuda"
                if torch.cuda.is_available()
                else (
                    "xpu"
                    if (hasattr(torch, "xpu") and torch.xpu.is_available())
                    else "cpu"
                )
            )
        )
    )

    # Convert all SVF inputs to the same device
    svfS = torch.as_tensor(svfS, device=device)
    svfW = torch.as_tensor(svfW, device=device)
    svfN = torch.as_tensor(svfN, device=device)
    svfE = torch.as_tensor(svfE, device=device)
    svfEveg = torch.as_tensor(svfEveg, device=device)
    svfSveg = torch.as_tensor(svfSveg, device=device)
    svfWveg = torch.as_tensor(svfWveg, device=device)
    svfNveg = torch.as_tensor(svfNveg, device=device)
    svfEaveg = torch.as_tensor(svfEaveg, device=device)
    svfSaveg = torch.as_tensor(svfSaveg, device=device)
    svfWaveg = torch.as_tensor(svfWaveg, device=device)
    svfNaveg = torch.as_tensor(svfNaveg, device=device)

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
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
        svfE, svfEveg, svfEaveg, vikttot
    )

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaE)
        betaB = torch.arctan(torch.tan((svfalfaE) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
        if (azimuth > (180 - t)) and (azimuth <= (360 - t)):
            Lwallsun = (
                SBC
                * ewall
                * (
                    (Ta + 273.15 + Tw * torch.sin(aziE * (torch.pi / 180)))
                    ** 4
                )
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = (
                SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
            )
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    # Longwave from ground (see Lcyl_v2022a for remaining fluxes)
    if anisotropic_longwave == 1:
        Lground = LupE * 0.5
        Least = Lground
    else:
        Lsky = ((svfE + svfEveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lground = LupE * 0.5
        Lrefl = (Ldown + LupE) * (viktrefl) * (1 - ewall) * 0.5
        Least = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lsouth
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
        svfS, svfSveg, svfSaveg, vikttot
    )

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaS)
        betaB = torch.arctan(torch.tan((svfalfaS) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaS)*(1+F_sh))
        if (azimuth <= (90 - t)) or (azimuth > (270 - t)):
            Lwallsun = (
                SBC
                * ewall
                * (
                    (Ta + 273.15 + Tw * torch.sin(aziS * (torch.pi / 180)))
                    ** 4
                )
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = (
                SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
            )
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    # Longwave from ground (see Lcyl_v2022a for remaining fluxes)
    if anisotropic_longwave == 1:
        Lground = LupS * 0.5
        Lsouth = Lground
    else:
        Lsky = ((svfS + svfSveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lground = LupS * 0.5
        Lrefl = (Ldown + LupS) * (viktrefl) * (1 - ewall) * 0.5
        Lsouth = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lwest
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
        svfW, svfWveg, svfWaveg, vikttot
    )

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaW)
        betaB = torch.arctan(torch.tan((svfalfaW) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaW)*(1+F_sh))
        if (azimuth > (360 - t)) or (azimuth <= (180 - t)):
            Lwallsun = (
                SBC
                * ewall
                * (
                    (Ta + 273.15 + Tw * torch.sin(aziW * (torch.pi / 180)))
                    ** 4
                )
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = (
                SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
            )
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    # Longwave from ground (see Lcyl_v2022a for remaining fluxes)
    if anisotropic_longwave == 1:
        Lground = LupW * 0.5
        Lwest = Lground
    else:
        Lsky = ((svfW + svfWveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lground = LupW * 0.5
        Lrefl = (Ldown + LupW) * (viktrefl) * (1 - ewall) * 0.5
        Lwest = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

    ## Lnorth
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
        svfN, svfNveg, svfNaveg, vikttot
    )

    if altitude > 0:  # daytime
        alfaB = torch.arctan(svfalfaN)
        betaB = torch.arctan(torch.tan((svfalfaN) * F_sh))
        betasun = ((alfaB - betaB) / 2) + betaB
        # betasun = torch.arctan(0.5*torch.tan(svfalfaN)*(1+F_sh))
        if (azimuth > (90 - t)) and (azimuth <= (270 - t)):
            Lwallsun = (
                SBC
                * ewall
                * (
                    (Ta + 273.15 + Tw * torch.sin(aziN * (torch.pi / 180)))
                    ** 4
                )
                * viktwall
                * (1 - F_sh)
                * torch.cos(betasun)
                * 0.5
            )
            Lwallsh = (
                SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
            )
        else:
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
    else:  # nighttime
        Lwallsun = 0
        Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

    # Longwave from ground (see Lcyl_v2022a for remaining fluxes)
    if anisotropic_longwave == 1:
        Lground = LupN * 0.5
        Lnorth = Lground
    else:
        Lsky = ((svfN + svfNveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lground = LupN * 0.5
        Lrefl = (Ldown + LupN) * (viktrefl) * (1 - ewall) * 0.5
        Lnorth = Lsky + Lwallsun + Lwallsh + Lveg + Lground + Lrefl

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return Least, Lsouth, Lwest, Lnorth


def Lside_veg_v2026(
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
    anisotropic_longwave,
):

    # This m-file is the current one that estimates L from the four cardinal points 20100414
    device = (
        Ldown.device
        if isinstance(Ldown, torch.Tensor)
        else (
            Ta.device
            if isinstance(Ta, torch.Tensor)
            else torch.device(
                "cuda"
                if torch.cuda.is_available()
                else (
                    "xpu"
                    if (hasattr(torch, "xpu") and torch.xpu.is_available())
                    else "cpu"
                )
            )
        )
    )

    # Convert all SVF inputs to the same device
    svfS = torch.as_tensor(svfS, device=device)
    svfW = torch.as_tensor(svfW, device=device)
    svfN = torch.as_tensor(svfN, device=device)
    svfE = torch.as_tensor(svfE, device=device)
    svfEveg = torch.as_tensor(svfEveg, device=device)
    svfSveg = torch.as_tensor(svfSveg, device=device)
    svfWveg = torch.as_tensor(svfWveg, device=device)
    svfNveg = torch.as_tensor(svfNveg, device=device)
    svfEaveg = torch.as_tensor(svfEaveg, device=device)
    svfSaveg = torch.as_tensor(svfSaveg, device=device)
    svfWaveg = torch.as_tensor(svfWaveg, device=device)
    svfNaveg = torch.as_tensor(svfNaveg, device=device)

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
    [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
        svfE, svfEveg, svfEaveg, vikttot
    )

    if anisotropic_longwave == 1:
        Least = torch.zeros_like(Ldown, device=device)
        Lnorth = torch.zeros_like(Ldown, device=device)
        Lwest = torch.zeros_like(Ldown, device=device)
        Lsouth = torch.zeros_like(Ldown, device=device)

        return Least, Lsouth, Lwest, Lnorth
    else:
        if altitude > 0:  # daytime
            alfaB = torch.arctan(svfalfaE)
            betaB = torch.arctan(torch.tan((svfalfaE) * F_sh))
            betasun = ((alfaB - betaB) / 2) + betaB
            # betasun = torch.arctan(0.5*torch.tan(svfalfaE)*(1+F_sh)) #TODO This should be considered in future versions
            if (azimuth > (180 - t)) and (azimuth <= (360 - t)):
                Lwallsun = (
                    SBC
                    * ewall
                    * (
                        (Ta + 273.15 + Tw * torch.sin(aziE * (torch.pi / 180)))
                        ** 4
                    )
                    * viktwall
                    * (1 - F_sh)
                    * torch.cos(betasun)
                    * 0.5
                )
                Lwallsh = (
                    SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
                )
            else:
                Lwallsun = 0
                Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
        else:  # nighttime
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

        Lsky = ((svfE + svfEveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lrefl = Ldown * (viktrefl) * (1 - ewall) * 0.5
        Least = Lsky + Lwallsun + Lwallsh + Lveg + Lrefl

        # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

        ## Lsouth
        [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
            svfS, svfSveg, svfSaveg, vikttot
        )

        if altitude > 0:  # daytime
            alfaB = torch.arctan(svfalfaS)
            betaB = torch.arctan(torch.tan((svfalfaS) * F_sh))
            betasun = ((alfaB - betaB) / 2) + betaB
            # betasun = torch.arctan(0.5*torch.tan(svfalfaS)*(1+F_sh))
            if (azimuth <= (90 - t)) or (azimuth > (270 - t)):
                Lwallsun = (
                    SBC
                    * ewall
                    * (
                        (Ta + 273.15 + Tw * torch.sin(aziS * (torch.pi / 180)))
                        ** 4
                    )
                    * viktwall
                    * (1 - F_sh)
                    * torch.cos(betasun)
                    * 0.5
                )
                Lwallsh = (
                    SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
                )
            else:
                Lwallsun = 0
                Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
        else:  # nighttime
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

        Lsky = ((svfS + svfSveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lrefl = Ldown * (viktrefl) * (1 - ewall) * 0.5
        Lsouth = Lsky + Lwallsun + Lwallsh + Lveg + Lrefl

        # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

        ## Lwest
        [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
            svfW, svfWveg, svfWaveg, vikttot
        )

        if altitude > 0:  # daytime
            alfaB = torch.arctan(svfalfaW)
            betaB = torch.arctan(torch.tan((svfalfaW) * F_sh))
            betasun = ((alfaB - betaB) / 2) + betaB
            # betasun = torch.arctan(0.5*torch.tan(svfalfaW)*(1+F_sh))
            if (azimuth > (360 - t)) or (azimuth <= (180 - t)):
                Lwallsun = (
                    SBC
                    * ewall
                    * (
                        (Ta + 273.15 + Tw * torch.sin(aziW * (torch.pi / 180)))
                        ** 4
                    )
                    * viktwall
                    * (1 - F_sh)
                    * torch.cos(betasun)
                    * 0.5
                )
                Lwallsh = (
                    SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
                )
            else:
                Lwallsun = 0
                Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
        else:  # nighttime
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

        Lsky = ((svfW + svfWveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lrefl = Ldown * (viktrefl) * (1 - ewall) * 0.5
        Lwest = Lsky + Lwallsun + Lwallsh + Lveg + Lrefl

        # clear alfaB betaB betasun Lsky Lwallsh Lwallsun Lveg Lground Lrefl viktveg viktwall viktsky

        ## Lnorth
        [viktveg, viktwall, viktsky, viktrefl] = Lvikt_veg(
            svfN, svfNveg, svfNaveg, vikttot
        )

        if altitude > 0:  # daytime
            alfaB = torch.arctan(svfalfaN)
            betaB = torch.arctan(torch.tan((svfalfaN) * F_sh))
            betasun = ((alfaB - betaB) / 2) + betaB
            # betasun = torch.arctan(0.5*torch.tan(svfalfaN)*(1+F_sh))
            if (azimuth > (90 - t)) and (azimuth <= (270 - t)):
                Lwallsun = (
                    SBC
                    * ewall
                    * (
                        (Ta + 273.15 + Tw * torch.sin(aziN * (torch.pi / 180)))
                        ** 4
                    )
                    * viktwall
                    * (1 - F_sh)
                    * torch.cos(betasun)
                    * 0.5
                )
                Lwallsh = (
                    SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * F_sh * 0.5
                )
            else:
                Lwallsun = 0
                Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5
        else:  # nighttime
            Lwallsun = 0
            Lwallsh = SBC * ewall * ((Ta + 273.15) ** 4) * viktwall * 0.5

        Lsky = ((svfN + svfNveg - 1) * Lsky_allsky) * viktsky * 0.5
        Lveg = SBC * ewall * ((Ta + 273.15) ** 4) * viktveg * 0.5
        Lrefl = Ldown * (viktrefl) * (1 - ewall) * 0.5
        Lnorth = Lsky + Lwallsun + Lwallsh + Lveg + Lrefl

        if device.type == "cuda":
            torch.cuda.empty_cache()
        elif device.type == "xpu":
            torch.xpu.empty_cache()
        return Least, Lsouth, Lwest, Lnorth
