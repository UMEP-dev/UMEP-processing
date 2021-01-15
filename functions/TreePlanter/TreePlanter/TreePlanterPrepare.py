import numpy as np
import time
from ..TreeGenerator import makevegdems
from ..TreePlanter.TreePlanterClasses import Treerasters
from ..TreePlanter.TreePlanterClasses import Position

def treeplanter(treeinput,treedata,treerasters,tmrt_1d):

    treeinput.tmrt_s = treeinput.tmrt_s * treeinput.buildings     # Remove all Tmrt values that are in shade or on top of buildings

    bld_copy = treeinput.buildings.copy()

    # Creating boolean for where it is possible to plant a tree
    bd_b = np.int_( np.ceil( (treedata.dia / 2) / treeinput.gt[1] ) )

    # Buffer on building raster so that trees can't be planted next to walls. Can be planted one radius from walls.
    for i1 in range(bd_b):
        walls = np.zeros((treeinput.rows, treeinput.cols))
        domain = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        for i in range(1, treeinput.cols - 1):
            for j in range(1, treeinput.rows - 1):
                dom = bld_copy[j - 1:j + 2, i - 1:i + 2]
                walls[j, i] = np.min(dom[np.where(domain == 1)])

        walls = bld_copy - walls
        bld_copy = 1 - bld_copy
        bld_copy = bld_copy + walls
        bld_copy = 1 - bld_copy

    # Remove all possible positions outside selected area
    bld_copy = bld_copy * treeinput.selected_area

    # Calculating sum of Tmrt in shade for each possible position in the Tmrt matrix
    sum_tmrt_tsh = np.zeros((treeinput.rows, treeinput.cols))   # Empty matrix for sum of Tmrt in tree shadow
    sum_tmrt = np.zeros((treeinput.rows, treeinput.cols))  # Empty matrix for sum of Tmrt in sun under tree shadow

    res_y, res_x = np.where(bld_copy == 1)  # Coordinates for where it is possible to plant a tree (buildings, area of interest excluded)

    pos_ls = np.zeros((res_y.__len__(), 6))      # Length of vectors with y and x positions. Will have x and y positions, tmrt in shade and in sun and an id for each position

    index = 0

    # This loop pastes the tree shadow into all possible positions and calculates the sum of Tmrt under the shadow and
    # the sum Tmrt for the same area but sunlit
    for i in range(res_y.__len__()):
        y1 = np.int_(res_y[i] - treerasters.buffer_y[0])
        y2 = np.int_(res_y[i] + treerasters.buffer_y[1])
        x1 = np.int_(res_x[i] - treerasters.buffer_x[0])
        x2 = np.int_(res_x[i] + treerasters.buffer_x[1])

        ts_temp1 = np.zeros((treeinput.rows, treeinput.cols))
        ts_temp1[y1:y2,x1:x2] = treerasters.treeshade

        # Klistra in skugga!! Fixa för TreePlanterTreeshade.py
        # Gör samma som i TreePlanterOptimizer (regional groups, etc)
        for j in range(tmrt_1d.__len__()):
            ts_temp2 = np.zeros((treeinput.rows, treeinput.cols))
            ts_temp2[y1:y2, x1:x2] = treerasters.treeshade_bool[:, :, j]
            sum_tmrt[res_y[i],res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * treeinput.shadow[:,:,j] * treeinput.tmrt_ts[:,:,j])
            sum_tmrt_tsh[res_y[i], res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * treeinput.shadow[:,:,j] * tmrt_1d[j,0])

        pos_ls[index, 1] = res_x[i]                         # X position of tree
        pos_ls[index, 2] = res_y[i]                         # Y position of tree
        pos_ls[index, 3] = sum_tmrt_tsh[res_y[i], res_x[i]] # Sum of Tmrt in tree shade - vector
        pos_ls[index, 4] = sum_tmrt[res_y[i], res_x[i]]     # Sum of Tmrt in same area as tree shade but sunlit - vector
        pos_ls[index, 5] = 1

        index += 1

    pos_bool = pos_ls[:,3] != 0
    pos_ls = pos_ls[pos_bool,:]
    # Gives a unique value ranging from 1 to length of pos_ls.shape[0]+1, for each position where it is possible to plant a tree
    pos_ls[:, 0] = np.arange(1,pos_ls.shape[0]+1)

    # Adding sum_tmrt and sum_tmrt_tsh to the Treerasters class as well as calculating the difference between sunlit and shaded
    treerasters.tmrt(sum_tmrt, sum_tmrt_tsh)

    # Adding pos_ls and adding all unique values from pos_ls to their respective positions in a 2D matrix and returns a Positions class
    positions = Position(pos_ls,treeinput.rows,treeinput.cols)

    return treerasters, positions
