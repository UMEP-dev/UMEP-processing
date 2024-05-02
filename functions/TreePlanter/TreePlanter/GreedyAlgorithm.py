import numpy as np
import time
from ...TreeGenerator import makevegdems
# from ..TreeGenerator import makevegdems
from ..TreePlanter.TreePlanterClasses import Treerasters
from ..TreePlanter.TreePlanterClasses import Position
from ..TreePlanter.TreePlanterTreeshade import tree_slice

def greedyplanter(treeinput,treedata,treerasters,tmrt_1d,trees,feedback):

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

    # This loop pastes the tree shadow into all possible positions and calculates the sum of Tmrt under the shadow and
    # the sum Tmrt for the same area but sunlit
    best_y = np.zeros((trees))
    best_x = np.zeros((trees))
    # index = 0

    tmrt_copy = treeinput.tmrt_ts.copy()
    shadows_copy = treeinput.shadow.copy()

    tmrt_max = 0

    # Calculating sum of Tmrt in shade for each possible position in the Tmrt matrix
    sum_tmrt_tsh = np.zeros((treeinput.rows, treeinput.cols))   # Empty matrix for sum of Tmrt in tree shadow
    sum_tmrt = np.zeros((treeinput.rows, treeinput.cols))  # Empty matrix for sum of Tmrt in sun under tree shadow

    res_y, res_x = np.where(bld_copy == 1)  # Coordinates for where it is possible to plant a tree (buildings, area of interest excluded)

    for tree in range(trees):
        if feedback.isCanceled():
            break
        for i in range(res_y.__len__()):

            if feedback.isCanceled():
                break

            y1 = np.int_(res_y[i] - treerasters.buffer_y[0])
            y2 = np.int_(res_y[i] + treerasters.buffer_y[1])
            x1 = np.int_(res_x[i] - treerasters.buffer_x[0])
            x2 = np.int_(res_x[i] + treerasters.buffer_x[1])

            yslice1, xslice1, yslice2, xslice2 = tree_slice(y1,y2,x1,x2,treeinput,treerasters)

            ts_temp1 = np.zeros((treeinput.rows, treeinput.cols))
            ts_temp1[yslice2, xslice2] = treerasters.treeshade[yslice1, xslice1]

            # Estimating Tmrt in tree shade and in sun
            for j in range(tmrt_1d.__len__()):
                ts_temp2 = np.zeros((treeinput.rows, treeinput.cols))
                ts_temp2[yslice2, xslice2] = treerasters.treeshade_bool[yslice1, xslice1, j]
                sum_tmrt[res_y[i],res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * shadows_copy[:,:,j] * tmrt_copy[:,:,j])
                sum_tmrt_tsh[res_y[i], res_x[i]] += np.sum(ts_temp2 * treeinput.buildings * shadows_copy[:,:,j] * tmrt_1d[j,0])

        # Adding sum_tmrt and sum_tmrt_tsh to the Treerasters class as well as calculating the difference between sunlit and shaded
        treerasters.tmrt(sum_tmrt, sum_tmrt_tsh)

        if tree == 0:
            possible_locations = np.sum(treerasters.d_tmrt > 0)
            feedback.setProgressText(str(possible_locations) + " possible locations for trees...")

        temp_y, temp_x = np.where(treerasters.d_tmrt == np.max(treerasters.d_tmrt))
        
        tmrt_max += np.max(treerasters.d_tmrt)

        best_y[tree] = temp_y[0]
        best_x[tree] = temp_x[0]
        
        # Add new tree shade and remove Tmrt from Tmrt.SOLWEIG where there is shade from new tree
        y1 = np.int_(temp_y[0] - treerasters.buffer_y[0])
        y2 = np.int_(temp_y[0] + treerasters.buffer_y[1])
        x1 = np.int_(temp_x[0] - treerasters.buffer_x[0])
        x2 = np.int_(temp_x[0] + treerasters.buffer_x[1])

        yslice1, xslice1, yslice2, xslice2 = tree_slice(y1,y2,x1,x2,treeinput,treerasters)

        for j in range(tmrt_1d.__len__()):
            temp_shadow = np.zeros((treeinput.rows,treeinput.cols))
            temp_shadow[yslice2, xslice2] = treerasters.treeshade_bool[yslice1 , xslice1, j]
            temp_shadow = 1 - temp_shadow
            shadows_copy[:,:,j] = shadows_copy[:,:,j] * temp_shadow
            tmrt_copy[:,:,j] = tmrt_copy[:,:,j] * temp_shadow

        # Determine where to recalcaulate d_tmrt
        y1 = np.int_(temp_y[0] - treerasters.buffer_y[0] - treerasters.buffer_y[1])
        y2 = np.int_(temp_y[0] + treerasters.buffer_y[1] + treerasters.buffer_y[0])
        x1 = np.int_(temp_x[0] - treerasters.buffer_x[0] - treerasters.buffer_x[1])
        x2 = np.int_(temp_x[0] + treerasters.buffer_x[1] + treerasters.buffer_x[0])

        _, __, yslice2, xslice2 = tree_slice(y1,y2,x1,x2,treeinput,treerasters)
   
        recalc_positions = np.zeros((bld_copy.shape[0], bld_copy.shape[1])) 
        recalc_positions[yslice2, xslice2] = 1                                   
        treerasters.d_tmrt[yslice2, xslice2] = 0
        sum_tmrt[yslice2, xslice2] = 0
        sum_tmrt_tsh[yslice2, xslice2] = 0

        # Remove position and one radian of tree canopy of added tree
        yt1 = np.int_(temp_y[0] - treerasters.buffer_y[0])
        yt2 = np.int_(temp_y[0] + treerasters.buffer_y[1])
        xt1 = np.int_(temp_x[0] - treerasters.buffer_x[0])
        xt2 = np.int_(temp_x[0] + treerasters.buffer_x[1])

        yslice1, xslice1, yslice2, xslice2 = tree_slice(yt1,yt2,xt1,xt2,treeinput,treerasters)

        added_tree = np.zeros((bld_copy.shape[0], bld_copy.shape[1]))
        added_tree[yslice2, xslice2] = treerasters.cdsm[yslice1, xslice1]
        added_tree = added_tree > 0
        added_tree = 1 - added_tree

        for i1 in range(bd_b):
            walls = np.ones((treeinput.rows, treeinput.cols))
            domain = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
            for i in range(1, treeinput.cols - 1):
                for j in range(1, treeinput.rows - 1):
                    dom = added_tree[j - 1:j + 2, i - 1:i + 2]
                    walls[j, i] = np.min(dom[np.where(domain == 1)])
            
            added_tree = added_tree * walls

        bld_copy = bld_copy * added_tree

        recalc_positions = recalc_positions * bld_copy
        res_y, res_x = np.where(recalc_positions == 1)
            
        if (np.max(recalc_positions) == 0) & (tree != trees-1):
            best_bool = (best_y > 0) & (best_x > 0)
            best_y = best_y[best_bool]
            best_x = best_x[best_bool]
            feedback.setProgressText('Not enough space for all trees! Found locations for ' + str(int(best_y.shape[0])) + ' out of ' + str(int(trees)) + ' trees.')
            break

        feedback.setProgress(int(tree * (100 / trees)))

    return best_y, best_x, tmrt_max
