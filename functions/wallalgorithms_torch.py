from builtins import range

# -*- coding: utf-8 -*-
__author__ = "xlinfr"

import math

try:
    import torch
    import torch.nn.functional as F
except:
    pass

# import scipy.misc as sc
import scipy.ndimage as sc


def findwalls_sp(arr_dsm, walllimit, device, footprint=None):
    """
    Identifie les murs de manière ultra-optimisée en mémoire (sans F.unfold).
    """
    # 1. S'assurer que l'entrée est un tenseur PyTorch
    if isinstance(arr_dsm, torch.Tensor):
        dsm_tensor = arr_dsm
    else:
        dsm_tensor = torch.tensor(arr_dsm, device=device)

    # 2. Définir le footprint par défaut (forme de diamant / points cardinaux)
    if footprint is None or footprint is False:
        footprint = torch.tensor(
            [[0, 1, 0], [1, 1, 1], [0, 1, 0]], device=device
        )

    fh, fw = footprint.shape
    pad_h, pad_w = fh // 2, fw // 2

    # Padding adaptatif basé sur la taille du filtre
    padded_a = dsm_tensor.unsqueeze(0).unsqueeze(0)
    padded_a = F.pad(
        padded_a, pad=(pad_w, pad_w, pad_h, pad_h), mode="replicate"
    )
    padded_a = padded_a.squeeze(0).squeeze(0)

    # Initialisation de la matrice des maximums avec une valeur minimale (-infini)
    max_neighbors = torch.full_like(dsm_tensor, float("-inf"))

    # Trouver les coordonnées où le filtre est actif (égal à 1)
    y_indices, x_indices = torch.where(footprint == 1)

    # Utilisation du glissement par vue (0 copie mémoire)
    H, W = dsm_tensor.shape
    for dy, dx in zip(y_indices, x_indices):
        # Cette ligne crée une "vue" virtuelle sans allouer de RAM
        shifted_view = padded_a[dy : dy + H, dx : dx + W]
        # Comparaison élément par élément optimisée
        max_neighbors = torch.maximum(max_neighbors, shifted_view)

    # 3. Identification des pixels de murs
    walls = max_neighbors - dsm_tensor

    # Appliquer la limite de hauteur des murs
    walls[walls < walllimit] = 0

    # 4. Remise à zéro des bordures extérieures
    walls[0, :] = 0
    walls[-1, :] = 0
    walls[:, 0] = 0
    walls[:, -1] = 0

    return walls


def findwalls(a, walllimit, feedback, total):
    # This function identifies walls based on a DSM and a wall-height limit
    # Walls are represented by outer pixels within building footprints
    #
    # Fredrik Lindberg, Goteborg Urban Climate Group
    # fredrikl@gvc.gu.se
    # 20150625

    col, row = a.shape
    walls = torch.zeros((col, row))
    domain = torch.tensor([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    index = 0
    for i in torch.arange(1, row - 1):
        if feedback.isCanceled():
            feedback.setProgressText("Calculation cancelled")
            break
        for j in torch.arange(1, col - 1):
            dom = a[j - 1 : j + 2, i - 1 : i + 2]
            walls[j, i] = torch.max(
                dom[torch.where(domain == 1)]
            )  # new 20171006
            index = index + 1
            feedback.setProgress(int(index * total))

    walls = torch.clone(walls - a)  # new 20171006
    walls[(walls < walllimit)] = 0

    walls[0 : walls.shape[0], 0] = 0
    walls[0 : walls.shape[0], walls.shape[1] - 1] = 0
    walls[0, 0 : walls.shape[0]] = 0
    walls[walls.shape[0] - 1, 0 : walls.shape[1]] = 0

    return walls


def filter1Goodwin_as_aspect_v3(
    walls_for_dir, scale, a, feedback, total, device, tile_size=2048
):
    """Applies the Goodwin et al.

    (2010) filter processing to calculate wall aspect based on a wall pixel
    grid, a DSM (a), and a scale factor.

    Optimized for GPU using 2D Convolutions to eliminate nested pixel loops
    and CPU-GPU transfer bottlenecks.

    :param walls_for_dir: Tensor (H, W), Binary grid of walls
    :param scale: Float, Scale factor determining filter size
    :param a: Tensor (H, W), Digital Surface Model (DSM)
    :param feedback: Optional, QGIS/Processing feedback object for progress
    :param total: Float, Progress step multiplier
    :return: Tensor (H, W), Wall aspect directions
    """
    row, col = a.shape

    # 1. Compute kernel footprint based on scale factor
    filtersize = torch.floor((scale + 0.0000000001) * 9)
    if filtersize <= 2:
        filtersize = 3
    elif filtersize != 9 and filtersize % 2 == 0:
        filtersize = filtersize + 1

    filtersize = int(filtersize)
    filthalveceil = int(torch.ceil(torch.tensor(filtersize) / 2.0))
    filthalvefloor = int(torch.floor(torch.tensor(filtersize) / 2.0))
    n = filtersize - 1

    # Initialize base structural elements on CPU for rotation
    filtmatrix = torch.zeros((filtersize, filtersize))
    buildfilt = torch.zeros((filtersize, filtersize))

    filtmatrix[:, filthalveceil - 1] = 1
    buildfilt[filthalveceil - 1, 0:filthalvefloor] = 1
    buildfilt[filthalveceil - 1, filthalveceil:filtersize] = 2

    filtmatrix_list = []
    buildfilt1_list = []
    buildfilt2_list = []

    # 2. Pre-calculate all 180 directional filters on CPU
    with torch.no_grad():
        for h in range(180):
            filtmatrix1temp = sc.rotate(
                filtmatrix.numpy(), h, order=1, reshape=False, mode="nearest"
            )
            filtmatrix1 = torch.round(torch.from_numpy(filtmatrix1temp))

            filtmatrixbuildtemp = sc.rotate(
                buildfilt.numpy(), h, order=0, reshape=False, mode="nearest"
            )
            filtmatrixbuild = torch.round(
                torch.from_numpy(filtmatrixbuildtemp)
            )

            index = 270 - h
            if h in (150, 30):
                filtmatrixbuild[:, n] = 0
            if index == 225:
                filtmatrix1[0, 0] = 1
                filtmatrix1[n, n] = 1
            if index == 135:
                filtmatrix1[0, n] = 1
                filtmatrix1[n, 0] = 1

            filtmatrix_list.append(filtmatrix1.unsqueeze(0).unsqueeze(0))
            buildfilt1_list.append(
                (filtmatrixbuild == 1).float().unsqueeze(0).unsqueeze(0)
            )
            buildfilt2_list.append(
                (filtmatrixbuild == 2).float().unsqueeze(0).unsqueeze(0)
            )

    # Stacking kernels on CPU first - Now perfectly defined!
    all_kernels_walls = torch.cat(filtmatrix_list, dim=0)
    all_kernels_dsm1 = torch.cat(buildfilt1_list, dim=0)
    all_kernels_dsm2 = torch.cat(buildfilt2_list, dim=0)

    # 3. Setup global output allocation arrays
    final_y = torch.zeros((row, col), device=device)
    final_x = torch.zeros((row, col), device=device)

    walls_binary = (walls_for_dir > 0).float().to(device)
    a_device = a.float().to(device)

    # Calculate total tiles for feedback tracking
    num_tiles_y = math.ceil(row / tile_size)
    num_tiles_x = math.ceil(col / tile_size)
    total_tiles = num_tiles_y * num_tiles_x
    tile_count = 0

    # 4. Loop through spatial tiles
    for r_start in range(0, row, tile_size):
        r_end = min(r_start + tile_size, row)

        for c_start in range(0, col, tile_size):
            c_end = min(c_start + tile_size, col)

            if feedback is not None and feedback.isCanceled():
                return final_y

            # Determine padded coordinates (the halo zone)
            pad_top = min(r_start, filthalvefloor)
            pad_bottom = min(row - r_end, filthalvefloor)
            pad_left = min(c_start, filthalvefloor)
            pad_right = min(col - c_end, filthalvefloor)

            # Slice tile out with padding included
            tile_a = (
                a_device[
                    (r_start - pad_top) : (r_end + pad_bottom),
                    (c_start - pad_left) : (c_end + pad_right),
                ]
                .unsqueeze(0)
                .unsqueeze(0)
            )

            tile_walls = (
                walls_binary[
                    (r_start - pad_top) : (r_end + pad_bottom),
                    (c_start - pad_left) : (c_end + pad_right),
                ]
                .unsqueeze(0)
                .unsqueeze(0)
            )

            # Local allocation sizes for this specific tile (including pads)
            tile_rows, tile_cols = tile_a.shape[2], tile_a.shape[3]

            # Running Maximum Setup for this tile
            z_max = torch.full((tile_rows, tile_cols), -1.0, device=device)
            h_best = torch.zeros(
                (tile_rows, tile_cols), dtype=torch.long, device=device
            )
            dsm_best1 = torch.zeros((tile_rows, tile_cols), device=device)
            dsm_best2 = torch.zeros((tile_rows, tile_cols), device=device)

            # 5. Process angles in small chunks inside this single spatial tile
            chunk_size = 10
            for idx in range(0, 180, chunk_size):
                end_idx = min(idx + chunk_size, 180)

                k_walls = all_kernels_walls[idx:end_idx].to(device)
                k_dsm1 = all_kernels_dsm1[idx:end_idx].to(device)
                k_dsm2 = all_kernels_dsm2[idx:end_idx].to(device)

                # Notice: padding=0 because we manually padded our spatial tiles!
                walls_conv = F.conv2d(tile_walls, k_walls, padding=0).squeeze(
                    0
                )
                dsm_conv1 = F.conv2d(tile_a, k_dsm1, padding=0).squeeze(0)
                dsm_conv2 = F.conv2d(tile_a, k_dsm2, padding=0).squeeze(0)

                # Account for spatial shrinkage since F.conv2d with padding=0 trims edges
                c_rows, c_cols = walls_conv.shape[1], walls_conv.shape[2]

                # Align running arrays dynamically to conv output window
                z_max_crop = z_max[
                    filthalvefloor : filthalvefloor + c_rows,
                    filthalvefloor : filthalvefloor + c_cols,
                ]

                chunk_max, chunk_h_local = torch.max(walls_conv, dim=0)
                is_new_max = chunk_max >= z_max_crop

                if is_new_max.any():
                    # Update local trackers
                    z_max[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ] = torch.where(is_new_max, chunk_max, z_max_crop)

                    chunk_h_absolute = chunk_h_local + idx
                    h_best_crop = h_best[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ]
                    h_best[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ] = torch.where(is_new_max, chunk_h_absolute, h_best_crop)

                    h_local_unsqueeze = chunk_h_local.unsqueeze(0)
                    chunk_dsm1 = torch.gather(
                        dsm_conv1, dim=0, index=h_local_unsqueeze
                    ).squeeze(0)
                    chunk_dsm2 = torch.gather(
                        dsm_conv2, dim=0, index=h_local_unsqueeze
                    ).squeeze(0)

                    dsm_best1_crop = dsm_best1[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ]
                    dsm_best2_crop = dsm_best2[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ]

                    dsm_best1[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ] = torch.where(is_new_max, chunk_dsm1, dsm_best1_crop)
                    dsm_best2[
                        filthalvefloor : filthalvefloor + c_rows,
                        filthalvefloor : filthalvefloor + c_cols,
                    ] = torch.where(is_new_max, chunk_dsm2, dsm_best2_crop)

            # Un-pad results to extract the pure valid window of this tile
            tile_y = 270.0 - h_best.float()
            tile_x = torch.where(dsm_best1 > dsm_best2, 1, 2)

            valid_tile_y = tile_y[
                pad_top : tile_rows - pad_bottom,
                pad_left : tile_cols - pad_right,
            ]
            valid_tile_x = tile_x[
                pad_top : tile_rows - pad_bottom,
                pad_left : tile_cols - pad_right,
            ]

            # Write back cleanly into the global array without overlap seams
            final_y[r_start:r_end, c_start:c_end] = valid_tile_y
            final_x[r_start:r_end, c_start:c_end] = valid_tile_x

            # Progress handling
            tile_count += 1
            if feedback is not None:
                feedback.setProgress(
                    int((tile_count / total_tiles) * total * 0.9)
                )

    # 6. Global Post-processing calculations
    border_mask = torch.zeros((row, col), dtype=torch.bool, device=device)
    start = filthalveceil - 1
    end_row = row - filthalveceil - 1
    end_col = col - filthalveceil - 1
    border_mask[start:end_row, start:end_col] = True

    valid_mask = (walls_binary == 1) & border_mask
    final_y = torch.where(valid_mask, final_y, torch.zeros_like(final_y))
    final_x = torch.where(valid_mask, final_x, torch.zeros_like(final_x))

    final_y[final_x == 1] = final_y[final_x == 1] - 180
    final_y[final_y < 0] = final_y[final_y < 0] + 360

    # Incorporate derivative fallback values for flat results
    grad, asp = get_ders(a, scale)
    asp_device = (
        torch.from_numpy(asp).to(device)
        if not isinstance(asp, torch.Tensor)
        else asp.to(device)
    )

    final_y = final_y + ((walls_binary == 1) * 1) * ((final_y == 0) * 1) * (
        asp_device / (math.pi / 180.0)
    )

    if feedback is not None:
        feedback.setProgress(int(total))

    return final_y


def cart2pol(x, y, units="deg"):
    radius = torch.sqrt(x**2 + y**2)
    theta = torch.arctan2(y, x)
    if units in ["deg", "degs"]:
        theta = theta * 180 / torch.pi
    return theta, radius


def get_ders(dsm, scale):
    # dem,_,_=read_dem_grid(dem_file)
    dx = 1 / scale
    # dx=0.5
    fy, fx = torch.gradient(dsm, spacing=dx)
    asp, grad = cart2pol(fy, fx, "rad")
    grad = torch.arctan(grad)
    asp = asp * -1
    asp = asp + (asp < 0) * (torch.pi * 2)
    return grad, asp
