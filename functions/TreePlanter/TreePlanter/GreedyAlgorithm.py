import numpy as np
import time
from ..TreeGenerator import makevegdems
from ..TreePlanter.TreePlanterClasses import Treerasters
from ..TreePlanter.TreePlanterClasses import Position

def greedyplanter(treeinput,treedata,treerasters,tmrt_1d,trees,feedback):

    treeinput.tmrt_s = treeinput.tmrt_s * treeinput.buildings     # Remove all Tmrt values that are in shade or on top of buildings

    bld_copy = treeinput.buildings.copy()

    # Creating boolean for where it is possible to plant a tree
    bd_b = np.int_(np.ceil(treedata.dia / 2))

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

    # This loop pastes the tree shadow into all possible positions and calculates the sum of Tmrt under the shadow and
    # the sum Tmrt for the same area but sunlit
    best_y = np.zeros((trees))
    best_x = np.zeros((trees))
    # index = 0

    tmrt_copy = treeinput.tmrt_ts.copy()
    shadows_copy = treeinput.shadow.copy()

    # Calculating sum of Tmrt in shade for each possible position in the Tmrt matrix
    sum_tmrt_tsh = np.zeros((treeinput.rows, treeinput.cols))   # Empty matrix for sum of Tmrt in tree shadow
    sum_tmrt = np.zeros((treeinput.rows, treeinput.cols))  # Empty matrix for sum of Tmrt in sun under tree shadow

    res_y, res_x = np.where(bld_copy == 1)  # Coordinates for where it is possible to plant a tree (buildings, area of interest excluded)

    for tree in range(trees):
        for i in range(res_y.__len__()):
            y1 = np.int_(res_y[i] - treerasters.buffer_y[0])
            y2 = np.int_(res_y[i] + treerasters.buffer_y[1])
            x1 = np.int_(res_x[i] - treerasters.buffer_x[0])
            x2 = np.int_(res_x[i] + treerasters.buffer_x[1])

            ts_temp1 = np.zeros((treeinput.rows, treeinput.cols))
            ts_temp1[y1:y2,x1:x2] = treerasters.treeshade

            # Estimating Tmrt in tree shade and in sun
            for j in range(tmrt_1d.__len__()):
                ts_temp2 = np.zeros((treeinput.rows, treeinput.cols))
                ts_temp2[y1:y2, x1:x2] = treerasters.treeshade_bool[:, :, j]
                sum_tmrt[res_y[i],res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * shadows_copy[:,:,j] * tmrt_copy[:,:,j])
                sum_tmrt_tsh[res_y[i], res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * shadows_copy[:,:,j] * tmrt_1d[j,0])

        # Adding sum_tmrt and sum_tmrt_tsh to the Treerasters class as well as calculating the difference between sunlit and shaded
        treerasters.tmrt(sum_tmrt, sum_tmrt_tsh)

        if tree == 0:
            possible_locations = np.sum(treerasters.d_tmrt > 0)
            feedback.setProgressText(str(possible_locations) + " possible locations for trees...")

        temp_y, temp_x = np.where(treerasters.d_tmrt == np.max(treerasters.d_tmrt))
        
        best_y[tree] = temp_y[0]
        best_x[tree] = temp_x[0]
        
        y1 = np.int_(temp_y[0] - treerasters.buffer_y[0])
        y2 = np.int_(temp_y[0] + treerasters.buffer_y[1])
        x1 = np.int_(temp_x[0] - treerasters.buffer_x[0])
        x2 = np.int_(temp_x[0] + treerasters.buffer_x[1])

        for j in range(tmrt_1d.__len__()):
            temp_shadow = np.zeros((treeinput.rows,treeinput.cols))
            temp_shadow[y1:y2,x1:x2] = treerasters.treeshade_bool[:,:,j]
            temp_shadow = 1 - temp_shadow
            shadows_copy[:,:,j] = shadows_copy[:,:,j] * temp_shadow
            tmrt_copy[:,:,j] = tmrt_copy[:,:,j] * temp_shadow

        # Test
        y1 = np.int_(temp_y[0] - treerasters.buffer_y[0] - treerasters.buffer_y[1])
        if y1 < 0:
            y1 = 0
        y2 = np.int_(temp_y[0] + treerasters.buffer_y[1] + treerasters.buffer_y[0])
        if y2 < 0:
            y2 = 0
        x1 = np.int_(temp_x[0] - treerasters.buffer_x[0] - treerasters.buffer_x[1])
        if x1 < 0:
            x1 = 0
        x2 = np.int_(temp_x[0] + treerasters.buffer_x[1] + treerasters.buffer_x[0])
        if x2 < 0:
            x2 = 0

        recalc_positions = np.zeros((bld_copy.shape[0], bld_copy.shape[1])) # Alla noll
        recalc_positions[y1:y2,x1:x2] = 1                                   # Runt bästa position = 1
        bld_copy[temp_y[0]-bd_b:temp_y[0]+bd_b+1,temp_x[0]-bd_b:temp_x[0]+bd_b+1] = 0 # Där trädet står noll
        treerasters.d_tmrt[y1:y2,x1:x2] = 0
        sum_tmrt[y1:y2,x1:x2] = 0
        sum_tmrt_tsh[y1:y2,x1:x2] = 0
        recalc_positions = recalc_positions * bld_copy
        res_y, res_x = np.where(recalc_positions == 1)

        # Remove one radian of canopy diameter surrounding the position of the tree
        # bld_copy[temp_y[0]-bd_b:temp_y[0]+bd_b+1,temp_x[0]-bd_b:temp_x[0]+bd_b+1] = 0

        feedback.setProgress(int(tree * (100 / trees)))

        # index += 1


    return best_y, best_x
