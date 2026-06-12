try:
    import torch

except:
    pass


def Lvikt_veg(svf, svfveg, svfaveg, vikttot):
    device = None
    if isinstance(svf, torch.Tensor):
        device = svf.device
    elif isinstance(svfveg, torch.Tensor):
        device = svfveg.device
    elif isinstance(svfaveg, torch.Tensor):
        device = svfaveg.device
    else:
        torch.device(
            "cuda"
            if torch.cuda.is_available()
            else (
                "xpu"
                if (hasattr(torch, "xpu") and torch.xpu.is_available())
                else "cpu"
            )
        )

    svf = torch.as_tensor(svf, device=device)
    svfveg = torch.as_tensor(svfveg, device=device)
    svfaveg = torch.as_tensor(svfaveg, device=device)

    # Least
    viktonlywall = (
        vikttot
        - (
            63.227 * svf**6
            - 161.51 * svf**5
            + 156.91 * svf**4
            - 70.424 * svf**3
            + 16.773 * svf**2
            - 0.4863 * svf
        )
    ) / vikttot

    viktaveg = (
        vikttot
        - (
            63.227 * svfaveg**6
            - 161.51 * svfaveg**5
            + 156.91 * svfaveg**4
            - 70.424 * svfaveg**3
            + 16.773 * svfaveg**2
            - 0.4863 * svfaveg
        )
    ) / vikttot

    viktwall = viktonlywall - viktaveg

    svfvegbu = svfveg + svf - 1  # Vegetation plus buildings
    viktsky = (
        63.227 * svfvegbu**6
        - 161.51 * svfvegbu**5
        + 156.91 * svfvegbu**4
        - 70.424 * svfvegbu**3
        + 16.773 * svfvegbu**2
        - 0.4863 * svfvegbu
    ) / vikttot
    viktrefl = (
        vikttot
        - (
            63.227 * svfvegbu**6
            - 161.51 * svfvegbu**5
            + 156.91 * svfvegbu**4
            - 70.424 * svfvegbu**3
            + 16.773 * svfvegbu**2
            - 0.4863 * svfvegbu
        )
    ) / vikttot
    viktveg = (
        vikttot
        - (
            63.227 * svfvegbu**6
            - 161.51 * svfvegbu**5
            + 156.91 * svfvegbu**4
            - 70.424 * svfvegbu**3
            + 16.773 * svfvegbu**2
            - 0.4863 * svfvegbu
        )
    ) / vikttot
    viktveg = viktveg - viktwall

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return viktveg, viktwall, viktsky, viktrefl
