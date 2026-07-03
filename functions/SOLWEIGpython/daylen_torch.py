try:
    import torch

except:
    pass


def daylen(DOY, XLAT):
    # Calculation of declination of sun (Eqn. 16). Amplitude= +/-23.45
    # deg. Minimum = DOY 355 (DEC 21), maximum = DOY 172.5 (JUN 21/22).
    # Sun angles.  SOC limited for latitudes above polar circles.
    # Calculate daylength, sunrise and sunset (Eqn. 17)

    device = (
        DOY.device
        if isinstance(DOY, torch.Tensor)
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
    if not isinstance(XLAT, torch.Tensor):
        XLAT = torch.tensor(XLAT, device=device)

    RAD = torch.tensor(torch.pi / 180.0, device=device)

    DEC = -23.45 * torch.cos(2.0 * torch.pi * (DOY + 10.0) / 365.0)

    SOC = torch.tan(RAD * DEC) * torch.tan(RAD * XLAT)
    SOC = torch.clamp(SOC, -1.0, 1.0)
    # SOC=alt

    DAYL = 12.0 + 24.0 * torch.arcsin(SOC) / torch.pi
    SNUP = 12.0 - DAYL / 2.0
    SNDN = 12.0 + DAYL / 2.0
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return DAYL, DEC, SNDN, SNUP
