from __future__ import division

try:
    import torch
except:
    pass


def diffusefraction(radG, altitude, Kt, Ta, RH):
    """
    This function estimates diffuse and directbeam radiation according to
    Reindl et al (1990), Solar Energy 45:1

    :param radG:
    :param altitude:
    :param Kt: # radiation at the top of the atmosphere
    :param Ta:
    :param RH:
    :return:
    """

    alfa = altitude * (torch.pi / 180)

    if Ta <= -999.00 or RH <= -999.00 or torch.isnan(Ta) or torch.isnan(RH):
        if Kt <= 0.3:
            radD = radG * (1.020 - 0.248 * Kt)
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG * (1.45 - 1.67 * Kt)
        else:
            radD = radG * 0.147
    else:
        RH = RH / 100
        if Kt <= 0.3:
            radD = radG * (
                1
                - 0.232 * Kt
                + 0.0239 * torch.sin(alfa)
                - 0.000682 * Ta
                + 0.0195 * RH
            )
        elif Kt > 0.3 and Kt < 0.78:
            radD = radG * (
                1.329
                - 1.716 * Kt
                + 0.267 * torch.sin(alfa)
                - 0.00357 * Ta
                + 0.106 * RH
            )
        else:
            radD = radG * (
                0.426 * Kt
                - 0.256 * torch.sin(alfa)
                + 0.00349 * Ta
                + 0.0734 * RH
            )

    radI = (radG - radD) / (torch.sin(alfa))

    # Corrections for low sun altitudes (20130307)
    if radI < 0:
        radI = 0

    if altitude < 1 and radI > radG:
        radI = radG

    if radD > radG:
        radD = radG

    return radI, radD
