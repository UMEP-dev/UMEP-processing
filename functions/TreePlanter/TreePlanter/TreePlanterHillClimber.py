import numpy as np
import itertools
from scipy.ndimage import label
import time
import datetime
from ..TreePlanter import HillClimberAlgorithm
from ..TreePlanter.TreePlanterTreeshade import tsh_gen
from ..TreePlanter.TreePlanterTreeshade import tsh_gen_ts
from ..TreePlanter.adjustments import treenudge
from ..TreePlanter import StartingPositions

def combine(tup, t):
    return tuple(itertools.combinations(tup, t))

def treeoptinit(treerasters, treeinput, positions, treedata, shadow_rg, tmrt_1d, trees, r_iters, sa, feedback):

    #time_sum = 0

    dia = treedata.dia  # Diameter of tree canopy

    i_tmrt = np.zeros((r_iters)) # Empty vector to be filled with Tmrt values for each tree
    i_y = np.zeros((r_iters, trees))    # Empty vector to be filled with corresponding y position of the above
    i_x = np.zeros((r_iters, trees))    # Empty vector to be filled with corresponding x position of the above

    tree_pos_all = np.zeros((r_iters,trees))

    # Pad for kernel for neighboring trees
    pos_m_pad_t = np.pad(positions.pos_m, pad_width=((1, 1), (1, 1)), mode='constant',
                       constant_values=0)

    break_loop = 0

    #sa = 0  # Hill-climbing with random restarts
    #sa = 1  # Genetic x or y random
    #sa = 2  # Genetic parents

    percentage_progress = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    iters_progress = np.int_(percentage_progress * r_iters)
    progress_counter = itertools.cycle(range(percentage_progress.shape[0]))
    progress = progress_counter.__next__()

    if (sa == 1):
        tree_pos_c = np.zeros((trees))
        tree_pos_y = 0
        tree_pos_x = 0
        d_tmrt_p = np.zeros((trees))
        i_tmrt_max = np.max(i_tmrt)

    # Iterate for r_iters number of iterations. Will find optimal positions for trees and return the positions and decrease in tmrt
    for counter in range(r_iters):
        # Check if plugin is cancelled
        if feedback.isCanceled():
            break

        # Printing progress
        if counter == iters_progress[progress]:
            feedback.setProgressText(str(percentage_progress[progress] * 100) + " percent of iterations finished...")
            progress = progress_counter.__next__()

        r_count = 0
        while r_count < 1:
            # Check if plugin is cancelled
            if feedback.isCanceled():
                break

            # Creating starting positions.
            # If sa = 0, random restart, i
            # If sa = 1 evolutionary restart, i.e. y or x is random, the other is kept from previous run (previous local optimum)
            if sa == 0:
                tree_pos, tp_c, break_loop = StartingPositions.random_start(positions.pos[:,0], trees, tree_pos_all, r_iters)
            elif sa == 1:
                tree_pos_y, tree_pos_x, tree_pos, tree_pos_c, tp_c, break_loop = \
                    StartingPositions.genetic_start(tree_pos_x, tree_pos_y, tree_pos_c, positions, trees, tree_pos_all,
                                                   r_iters, counter, dia)

            if (tp_c == 100):
                feedback.setProgressText('Possibly too many trees to fit in planting area. Try a lower number.')
                break

            if ((counter == 0) | (sa == 0)):
                tree_pos_y = np.zeros((trees), dtype=int)  # Random y-positions for trees
                tree_pos_x = np.zeros((trees), dtype=int)  # Random x-positions for trees

                # Y and X positions of starting positions
                for i2 in range(tree_pos.__len__()):
                    tree_pos_x[i2] = positions.pos[positions.pos[:, 0] == tree_pos[i2], 1]
                    tree_pos_y[i2] = positions.pos[positions.pos[:, 0] == tree_pos[i2], 2]

            # Euclidean distance between random positions so that trees are not too close to each other
            it_comb = combine(tree_pos, 2)
            eucl_dist = np.zeros((it_comb.__len__(), 1))

            for i3 in range(it_comb.__len__()):
                ax = positions.pos[positions.pos[:, 0] == it_comb[i3][0], 1]
                ay = positions.pos[positions.pos[:, 0] == it_comb[i3][0], 2]
                a = np.array(([ay[0], ax[0]]))
                bx = positions.pos[positions.pos[:, 0] == it_comb[i3][1], 1]
                by = positions.pos[positions.pos[:, 0] == it_comb[i3][1], 2]
                b = np.array(([by[0], bx[0]]))
                eucl_dist[i3, 0] = np.linalg.norm(a - b)

            if (np.min(eucl_dist[:, 0]) >= dia):
                r_count = 1
                # Lägg till alla positioner för eucl >= dia

            #r_count = 1
            tree_pos_all[counter] = np.sort(tree_pos)

        if (break_loop == (r_iters + 1)):
            break

        tp_nc = np.zeros((trees, 1))    # 1 if tree is in the same position as in previous iteration, otherwise 0
        tp_nc_a = np.zeros((trees))     # 1 if tree is stuck but have been checked for better position and none was found, otherwise 0

        ti = itertools.cycle(range(trees)) # Iterator to move between trees moving around in the study area

        # Create matrices for tree paths and add starting positions
        tree_paths_temp = np.zeros((treerasters.rows, treerasters.cols, trees))
        for t_i in range(trees):
            tree_paths_temp[tree_pos_y[t_i],tree_pos_x[t_i], t_i] = 1

        t1 = np.zeros((1,5))

        # Moving trees, i.e. optimization
        while np.sum(tp_nc[:,0]) < trees:

            i = ti.__next__()

            # Running optimizer
            # t1 = best shading position
            t1, nc, y_out, x_out = HillClimberAlgorithm.topt(tree_pos_y, tree_pos_x, treerasters, treeinput, dia, shadow_rg, tmrt_1d, positions, i, pos_m_pad_t, t1)
            tp_nc[i, 0] = nc

            if (tp_nc[i,0] == 0):
                tree_pos_y = y_out
                tree_pos_x = x_out

                tp_nc_a[i] = 0

            # Possibly moving trees that are stuck, where tree shadows intersect
            elif ((tp_nc_a[i] == 0) & (tp_nc[i,0] == 1)):
                y_out, x_out, tp_nc, tp_nc_a, t1 = treenudge(y_out, x_out, tp_nc, tp_nc_a, i, t1, counter, i_tmrt, treerasters, treeinput, positions,
                          tmrt_1d)
                if ((tp_nc_a[i] == 1) & (tp_nc[i, 0] == 0)):
                    tree_pos_y = y_out
                    tree_pos_x = x_out

            if (tp_nc[i,0] == 0): # Tree paths of current random run
                tree_paths_temp[tree_pos_y[i], tree_pos_x[i], i] = 1

            # Changing position of tree
            if (t1[0, 2] > i_tmrt[counter]):
                i_tmrt[counter] = t1[0, 2]
                i_x[counter, :] = tree_pos_x[:]
                i_y[counter, :] = tree_pos_y[:]

        # Check whether there is a new equal or worse position for each individual tree.
        if (sa == 1):
            if i_tmrt[counter] > i_tmrt_max:
                i_tmrt_max = np.max(i_tmrt)
                for idx in range(trees):
                    d_tmrt_p[idx] = treerasters.d_tmrt[tree_pos_y[idx], tree_pos_x[idx]]
            else:
                d_tmrt_temp = np.zeros((trees))
                for idx in range(trees):
                    d_tmrt_temp[idx] = treerasters.d_tmrt[tree_pos_y[idx], tree_pos_x[idx]]
                low_p = d_tmrt_temp <= d_tmrt_p
                tree_pos_c[low_p] += 1
                high_p = d_tmrt_temp > d_tmrt_p
                tree_pos_c[high_p] = 0

        if (i_tmrt[counter] == np.max(i_tmrt)):   # Path of trees with best positions from starting to ending
            tree_paths = tree_paths_temp.copy()

        if (tp_c == 100):
            break

        # Progress bar
        feedback.setProgress(int(counter * (100 / r_iters)))

    # Save locations for occurrence map
    i_y_all = i_y.copy()
    i_x_all = i_x.copy()

    # Finding best position from all r_iters iteration, i.e. if r_iters = 1000 then best position out of 1000 runs
    t_max = np.max(i_tmrt)
    y = np.where(i_tmrt == t_max)

    # Calculate heat map for best 33 % positions

    # Returns all the unique positions found by the algorithm and the potential decrease in tmrt
    #from misc import max_tmrt
    #unique_tmrt, unique_tmrt_max = max_tmrt(i_y, i_x, i_tmrt, trees, positions.pos)

    # Optimal positions of trees
    i_y = i_y[y[0][0], :]
    i_x = i_x[y[0][0], :]

    return i_y, i_x, t_max, i_y_all, i_x_all
    #return i_y, i_x, unique_tmrt, unique_tmrt_max, tree_paths
