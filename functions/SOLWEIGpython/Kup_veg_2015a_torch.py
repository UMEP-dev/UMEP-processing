try:
    import torch
except:
    pass


def Kup_veg_2015a(
    radI,
    radD,
    radG,
    altitude,
    svfbuveg,
    albedo_b,
    F_sh,
    gvfalb,
    gvfalbE,
    gvfalbS,
    gvfalbW,
    gvfalbN,
    gvfalbnosh,
    gvfalbnoshE,
    gvfalbnoshS,
    gvfalbnoshW,
    gvfalbnoshN,
):

    Kup = (gvfalb * radI * torch.sin(altitude * (torch.pi / 180.0))) + (
        radD * svfbuveg
        + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
    ) * gvfalbnosh

    KupE = (gvfalbE * radI * torch.sin(altitude * (torch.pi / 180.0))) + (
        radD * svfbuveg
        + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
    ) * gvfalbnoshE

    KupS = (gvfalbS * radI * torch.sin(altitude * (torch.pi / 180.0))) + (
        radD * svfbuveg
        + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
    ) * gvfalbnoshS

    KupW = (gvfalbW * radI * torch.sin(altitude * (torch.pi / 180.0))) + (
        radD * svfbuveg
        + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
    ) * gvfalbnoshW

    KupN = (gvfalbN * radI * torch.sin(altitude * (torch.pi / 180.0))) + (
        radD * svfbuveg
        + albedo_b * (1 - svfbuveg) * (radG * (1 - F_sh) + radD * F_sh)
    ) * gvfalbnoshN

    return Kup, KupE, KupS, KupW, KupN
