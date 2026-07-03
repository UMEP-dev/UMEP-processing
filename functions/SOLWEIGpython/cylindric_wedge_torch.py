try:
    import torch
except:
    pass


def cylindric_wedge(zen, svfalfa, rows, cols):

    # Fraction of sunlit walls based on sun altitude and svf wieghted building angles
    # input:
    # sun zenith angle "beta"
    # svf related angle "alfa"
    try:
        import torch
    except:
        Exception("Error, pytorch must be imported.")

    device = (
        svfalfa.device
        if isinstance(svfalfa, torch.Tensor)
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
    beta = zen
    # alfa=svfalfa
    alfa = torch.zeros((rows, cols), device=device) + svfalfa
    # measure the size of the image
    # sizex=size(svfalfa,2)
    # sizey=size(svfalfa,1)

    xa = 1 - 2.0 / (torch.tan(alfa) * torch.tan(beta))
    ha = 2.0 / (torch.tan(alfa) * torch.tan(beta))
    ba = 1.0 / torch.tan(alfa)
    hkil = 2.0 * ba * ha

    qa = torch.zeros((rows, cols), device=device)
    # qa(length(svfalfa),length(svfalfa))=0;
    qa[xa < 0] = torch.tan(beta) / 2

    Za = torch.zeros((rows, cols), device=device)
    # Za(length(svfalfa),length(svfalfa))=0;
    Za[xa < 0] = (((ba[xa < 0] ** 2) - ((qa[xa < 0] ** 2) / 4)) ** 0.5).float()

    phi = torch.zeros((rows, cols), device=device)
    # phi(length(svfalfa),length(svfalfa))=0;
    phi[xa < 0] = torch.arctan(Za[xa < 0] / qa[xa < 0])

    A = torch.zeros((rows, cols), device=device)
    # A(length(svfalfa),length(svfalfa))=0;
    A[xa < 0] = (
        torch.sin(phi[xa < 0]) - phi[xa < 0] * torch.cos(phi[xa < 0])
    ) / (1 - torch.cos(phi[xa < 0]))

    ukil = torch.zeros((rows, cols), device=device)
    # ukil(length(svfalfa),length(svfalfa))=0
    ukil[xa < 0] = (2 * ba[xa < 0] * xa[xa < 0] * A[xa < 0]).float()

    Ssurf = hkil + ukil

    F_sh = (2 * torch.pi * ba - Ssurf) / (2 * torch.pi * ba)  # Xa
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return F_sh


def cylindric_wedge_voxel(zen, svfalfa):

    torch.seterr(divide="ignore", invalid="ignore")

    # Fraction of sunlit walls based on sun altitude and svf wieghted building angles
    # input:
    # sun zenith angle "beta"
    # svf related angle "alfa"

    beta = zen

    xa = 1 - 2.0 / (torch.tan(svfalfa) * torch.tan(beta))
    ha = 2.0 / (torch.tan(svfalfa) * torch.tan(beta))
    ba = 1.0 / torch.tan(svfalfa)
    hkil = 2.0 * ba * ha

    device = (
        svfalfa.device
        if isinstance(svfalfa, torch.Tensor)
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

    qa = torch.zeros((svfalfa.shape[0]), device=device)
    # qa(length(svfalfa),length(svfalfa))=0;
    qa[xa < 0] = torch.tan(beta) / 2

    Za = torch.zeros((svfalfa.shape[0]), device=device)
    # Za(length(svfalfa),length(svfalfa))=0;
    Za[xa < 0] = ((ba[xa < 0] ** 2) - ((qa[xa < 0] ** 2) / 4)) ** 0.5

    phi = torch.zeros((svfalfa.shape[0]), device=device)
    # phi(length(svfalfa),length(svfalfa))=0;
    phi[xa < 0] = torch.arctan(Za[xa < 0] / qa[xa < 0])

    A = torch.zeros((svfalfa.shape[0]), device=device)
    # A(length(svfalfa),length(svfalfa))=0;
    A[xa < 0] = (
        torch.sin(phi[xa < 0]) - phi[xa < 0] * torch.cos(phi[xa < 0])
    ) / (1 - torch.cos(phi[xa < 0]))

    ukil = torch.zeros((svfalfa.shape[0]), device=device)
    # ukil(length(svfalfa),length(svfalfa))=0
    ukil[xa < 0] = 2 * ba[xa < 0] * xa[xa < 0] * A[xa < 0]

    Ssurf = hkil + ukil

    F_sh = (2 * torch.pi * ba - Ssurf) / (2 * torch.pi * ba)  # Xa
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return F_sh
