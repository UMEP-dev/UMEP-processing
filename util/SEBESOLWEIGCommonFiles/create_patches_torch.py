try:
    import torch
except:
    pass


def create_patches(patch_option, device):

    deg2rad = torch.pi / 180

    # patch_option = 1 = 145 patches (Robinson & Stone, 2004)
    # patch_option = 2 = 153 patches (Wallenberg et al., 2022)
    # patch_option = 3 = 306 patches -> test
    # patch_option = 4 = 612 patches -> test

    skyvaultalt = torch.atleast_2d([])
    skyvaultazi = torch.atleast_2d([])

    # Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
    # Patch option 1, 145 patches, Original Robinson & Stone (2004) after Tregenza (1987)/Tregenza & Sharples (1993)
    if patch_option == 1:
        annulino = torch.tensor(
            [0, 12, 24, 36, 48, 60, 72, 84, 90], device=device
        )
        skyvaultaltint = torch.tensor(
            [6, 18, 30, 42, 54, 66, 78, 90], device=device
        )  # Robinson & Stone (2004)
        azistart = torch.tensor(
            [0, 4, 2, 5, 8, 0, 10, 0], device=device
        )  # Fredrik/Nils
        patches_in_band = torch.tensor(
            [30, 30, 24, 24, 18, 12, 6, 1], device=device
        )  # Robinson & Stone (2004)
    # Patch option 2, 153 patches, Wallenberg et al. (2022)
    elif patch_option == 2:
        annulino = torch.tensor(
            [0, 12, 24, 36, 48, 60, 72, 84, 90], device=device
        )
        skyvaultaltint = torch.tensor(
            [6, 18, 30, 42, 54, 66, 78, 90], device=device
        )  # Robinson & Stone (2004)
        azistart = torch.tensor(
            [0, 4, 2, 5, 8, 0, 10, 0], device=device
        )  # Fredrik/Nils
        patches_in_band = torch.tensor(
            [31, 30, 28, 24, 19, 13, 7, 1], device=device
        )  # Nils
    # Patch option 3, 306 patches, test
    elif patch_option == 3:
        annulino = torch.tensor(
            [0, 12, 24, 36, 48, 60, 72, 84, 90], device=device
        )
        skyvaultaltint = torch.tensor(
            [6, 18, 30, 42, 54, 66, 78, 90], device=device
        )  # Robinson & Stone (2004)
        azistart = torch.tensor(
            [0, 4, 2, 5, 8, 0, 10, 0], device=device
        )  # Fredrik/Nils
        patches_in_band = torch.tensor(
            [31 * 2, 30 * 2, 28 * 2, 24 * 2, 19 * 2, 13 * 2, 7 * 2, 1],
            device=device,
        )  # Nils
    # Patch option 4, 612 patches, test
    elif patch_option == 4:
        annulino = torch.tensor(
            [0, 4.5, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 90],
            device=device,
        )  # Nils
        skyvaultaltint = torch.tensor(
            [3, 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 90],
            device=device,
        )  # Nils
        patches_in_band = torch.tensor(
            [
                31 * 2,
                31 * 2,
                30 * 2,
                30 * 2,
                28 * 2,
                28 * 2,
                24 * 2,
                24 * 2,
                19 * 2,
                19 * 2,
                13 * 2,
                13 * 2,
                7 * 2,
                7 * 2,
                1,
            ],
            device=device,
        )  # Nils
        azistart = torch.tensor(
            [0, 0, 4, 4, 2, 2, 5, 5, 8, 8, 0, 0, 10, 10, 0], device=device
        )  # Nils

    skyvaultaziint = torch.tensor(
        [360 / patches for patches in patches_in_band]
    )

    for j in range(0, skyvaultaltint.shape[0]):
        for k in range(0, patches_in_band[j]):
            skyvaultalt = skyvaultalt + (skyvaultaltint[j],)
            skyvaultazi = skyvaultazi + (k * skyvaultaziint[j] + azistart[j],)

    del deg2rad
    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()

    return (
        skyvaultalt,
        skyvaultazi,
        annulino,
        skyvaultaltint,
        patches_in_band,
        skyvaultaziint,
        azistart,
    )
