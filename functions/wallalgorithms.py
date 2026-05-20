from builtins import range

# -*- coding: utf-8 -*-
__author__ = "xlinfr"

import math
import numpy as np
import torch
import torch.nn.functional as F


# import scipy.misc as sc
import scipy.ndimage as sc
from scipy.ndimage import maximum_filter


def findwalls_sp(arr_dsm, walllimit, footprint=None):
    """
    This function identifies walls based on a DSM and a wall height limit.
    
    Parameters:
    arr_dsm (Tensor or array-like): Digital Surface Model
    walllimit (float): Wall height limit
    footprint (Tensor, optional): Structuring element for the maximum filter.
                                  Defaults to a diamond shape (cardinal neighbors).
    """
    # Ensure input is a PyTorch tensor and preserve its device (CPU or GPU)
    if isinstance(arr_dsm, torch.Tensor):
        dsm_tensor = arr_dsm
    else:
        dsm_tensor = torch.tensor(arr_dsm)

    # 1. Padding: 'edge' mode in NumPy becomes 'replicate' in PyTorch
    # PyTorch requires at least a 3D tensor for 'replicate', so we unsqueeze(0)
    padded_a = dsm_tensor.unsqueeze(0)
    padded_a = F.pad(padded_a, pad=(1, 1, 1, 1), mode="replicate")
    padded_a = padded_a.squeeze(0)

    # 2. Maximum Filter using Max Pooling
    # Default footprint is a 3x3 diamond shape (cardinal points)
    if footprint is None or footprint is False:
        # For max_pool2d, we simulate the footprint by looking at a 3x3 window.
        # To strictly replicate a diamond footprint in pure PyTorch, we can use 
        # a 2D morphological dilation with a custom kernel:
        kernel = torch.tensor([[0., 1., 0.],
                               [1., 1., 1.], # Including center for max filter
                               [0., 1., 0.]], device=dsm_tensor.device)
        
        # Reshape padded_a to (Batch=1, Channel=1, H, W) for 2D convolutions
        x = padded_a.unsqueeze(0).unsqueeze(0)
        
        # We use unfolds or a specialized max_pool, but since your footprint has 0s,
        # the cleanest PyTorch way is to mask the 3x3 windows:
        windows = F.unfold(x, kernel_size=3, padding=0) # Shape: (1, 9, L)
        # Multiply by kernel flattened, replacing 0s with -inf so they are ignored in max
        kernel_flat = kernel.flatten().unsqueeze(1)
        windows = windows + torch.where(kernel_flat == 1, 0.0, float('-inf'))
        
        # Take the maximum over the 9 neighbors
        max_neighbors = windows.max(dim=1)[0] # Shape: (1, L)
        
        # Reshape back to original DSM size (col, row)
        max_neighbors = max_neighbors.view(dsm_tensor.shape)
    else:
        # If a custom footprint is provided, you can apply a similar unfold method
        x = padded_a.unsqueeze(0).unsqueeze(0)
        windows = F.unfold(x, kernel_size=footprint.shape, padding=0)
        footprint_flat = footprint.to(dsm_tensor.device).flatten().unsqueeze(1)
        windows = windows + torch.where(footprint_flat == 1, 0.0, float('-inf'))
        max_neighbors = windows.max(dim=1)[0].view(dsm_tensor.shape)

    # 3. Identify wall pixels
    walls = max_neighbors - dsm_tensor

    # Apply wall height limit
    walls[walls < walllimit] = 0

    # 4. Set the outermost boundaries to zero
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

    col = a.shape[0]
    row = a.shape[1]
    walls = torch.zeros((col, row))
    domain = torch.tensor([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    index = 0
    for i in torch.arange(1, row - 1):
        if feedback.isCanceled():
            feedback.setProgressText("Calculation cancelled")
            break
        for j in torch.arange(1, col - 1):
            dom = a[j - 1 : j + 2, i - 1 : i + 2]
            walls[j, i] = torch.max(dom[torch.where(domain == 1)])  # new 20171006
            index = index + 1
            feedback.setProgress(int(index * total))

    walls = torch.clone(walls - a)  # new 20171006
    walls[(walls < walllimit)] = 0

    walls[0 : walls.shape[0], 0] = 0
    walls[0 : walls.shape[0], walls.shape[1] - 1] = 0
    walls[0, 0 : walls.shape[0]] = 0
    walls[walls.shape[0] - 1, 0 : walls.shape[1]] = 0

    return walls


def filter1Goodwin_as_aspect_v3(walls_for_dir, scale, a, feedback, total):
    """
    tThis function applies the filter processing presented in Goodwin et al (2010) but instead for removing
    linear fetures it calculates wall aspect based on a wall pixels grid, a dsm (a) and a scale factor

    Fredrik Lindberg, 2012-02-14
    fredrikl@gvc.gu.se

    Translated: 2015-09-15

    :param walls:
    :param scale:
    :param a:
    :return: dirwalls
    """

    walls = walls_for_dir.clone().to(a.device)

    row = a.shape[0]
    col = a.shape[1]

    filtersize = torch.floor((scale + 0.0000000001) * 9)
    if filtersize <= 2:
        filtersize = 3
    else:
        if filtersize != 9:
            if filtersize % 2 == 0:
                filtersize = filtersize + 1

    filthalveceil = int(torch.ceil(filtersize / 2.0))
    filthalvefloor = int(torch.floor(filtersize / 2.0))

    filtmatrix = torch.zeros((int(filtersize), int(filtersize)))
    buildfilt = torch.zeros((int(filtersize), int(filtersize)))

    filtmatrix[:, filthalveceil - 1] = 1
    n = filtmatrix.shape[0] - 1
    buildfilt[filthalveceil - 1, 0:filthalvefloor] = 1
    buildfilt[filthalveceil - 1, filthalveceil : int(filtersize)] = 2

    y = torch.zeros((row, col), device=a.device)  # final direction
    z = torch.zeros((row, col), device=a.device)  # temporary direction
    x = torch.zeros((row, col), device=a.device)  # building side
    walls[walls > 0] = 1

    for h in range(
        0, 180
    ):  # =0:1:180 #%increased resolution to 1 deg 20140911
        if feedback is not None:
            feedback.setProgress(int(h * total))
            if feedback.isCanceled():
                feedback.setProgressText("Calculation cancelled")
                break
        filtmatrix1temp = sc.rotate(
            filtmatrix.cpu().numpy()
            if filtmatrix.device.type == "cuda"
            else filtmatrix.numpy(),
            h,
            order=1,
            reshape=False,
            mode="nearest",
        )  # bilinear
        filtmatrix1 = torch.round(
            torch.from_numpy(filtmatrix1temp).to(walls.device)
        )
        # filtmatrix1temp = sc.imrotate(filtmatrix, h, 'bilinear')
        # filtmatrix1 = torch.round(filtmatrix1temp / 255.)
        # filtmatrixbuildtemp = sc.imrotate(buildfilt, h, 'nearest')
        filtmatrixbuildtemp = sc.rotate(
            buildfilt.cpu().numpy()
            if buildfilt.device.type == "cuda"
            else buildfilt.numpy(),
            h,
            order=0,
            reshape=False,
            mode="nearest",
        )  # Nearest neighbor
        # filtmatrixbuild = torch.round(filtmatrixbuildtemp / 127.)
        filtmatrixbuild = torch.round(
            torch.from_numpy(filtmatrixbuildtemp).to(walls.device)
        )
        index = 270 - h
        if h == 150:
            filtmatrixbuild[:, n] = 0
        if h == 30:
            filtmatrixbuild[:, n] = 0
        if index == 225:
            # n = filtmatrix.shape[0] - 1  # length(filtmatrix);
            filtmatrix1[0, 0] = 1
            filtmatrix1[n, n] = 1
        if index == 135:
            # n = filtmatrix.shape[0] - 1  # length(filtmatrix);
            filtmatrix1[0, n] = 1
            filtmatrix1[n, 0] = 1

        for i in range(
            int(filthalveceil) - 1, row - int(filthalveceil) - 1
        ):  # i=filthalveceil:sizey-filthalveceil
            for j in range(
                int(filthalveceil) - 1, col - int(filthalveceil) - 1
            ):  # (j=filthalveceil:sizex-filthalveceil
                if walls[i, j] == 1:
                    wallscut = (
                        walls[
                            i - filthalvefloor : i + filthalvefloor + 1,
                            j - filthalvefloor : j + filthalvefloor + 1,
                        ]
                        * filtmatrix1
                    )
                    dsmcut = a[
                        i - filthalvefloor : i + filthalvefloor + 1,
                        j - filthalvefloor : j + filthalvefloor + 1,
                    ]
                    if z[i, j] < wallscut.sum():  # sum(sum(wallscut))
                        z[i, j] = wallscut.sum()  # sum(sum(wallscut));
                        if torch.sum(dsmcut[filtmatrixbuild == 1]) > torch.sum(
                            dsmcut[filtmatrixbuild == 2]
                        ):
                            x[i, j] = 1
                        else:
                            x[i, j] = 2

                        y[i, j] = index

    y[(x == 1)] = y[(x == 1)] - 180
    y[(y < 0)] = y[(y < 0)] + 360

    grad, asp = get_ders(a, scale)

    y = y + ((walls == 1) * 1) * ((y == 0) * 1) * (asp / (math.pi / 180.0))

    dirwalls = y

    return dirwalls


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
