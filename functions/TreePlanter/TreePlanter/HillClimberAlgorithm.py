import numpy as np
import time
import timeit
from scipy.ndimage import label
from ..TreePlanter.TreePlanterTreeshade import tsh_gen
from ..TreePlanter.TreePlanterTreeshade import tsh_gen_ts
from ..TreePlanter.TreePlanterTreeshade import tsh_gen_mt1, tsh_gen_mt2
import itertools

# This function will look for better shading position for a tree by looking one pixel east, west, north, south of it's
# current position. It will compare it's position to the other trees in the study area, i.e. other moving trees.
def topt(y, x, treerasters, treeinput, dia, shadow_rg, tmrt_1d, positions, ti, pos_m_pad, p_tmrt):

    # x = x positions of trees
    # y = y positions of trees
    # treerasters = rasters with tmrt, shadows, etc
    # dia = diameter of tree canopy
    # buildings = raster with buildings
    # shadow_rg = vector with which timesteps shade which pixels
    # tmrt from SOLWEIG
    # i = which tree to move

    t_y = y[ti]    # y-position of tree to move
    t_x = x[ti]    # x-position of tree to move

    p_tmrt[0, 3] = t_y
    p_tmrt[0, 4] = t_x

    # Array with movement one step south, one step north, one step west, one step east
    k_d = 1

    yt = t_y + k_d
    xt = t_x + k_d
    domain = np.array([[1,1,1],
                       [1,0,1],
                       [1,1,1]])
    # domain = np.array([[0,1,0],
    #                    [1,0,1],
    #                    [0,1,0]])
    t_pos = pos_m_pad[yt - k_d:yt + k_d+1, xt - k_d:xt + k_d+1]
    t_pos = t_pos * domain
    t_pos = t_pos.flatten()
    t_pos = t_pos[t_pos != 0]
    t_yx = np.zeros((t_pos.shape[0],2))
    for i in range(t_pos.shape[0]):
        t_yx[i,0] = positions.pos[positions.pos[:, 0] == t_pos[i], 2]
        t_yx[i,1] = positions.pos[positions.pos[:, 0] == t_pos[i], 1]

    t_yx = np.int_(t_yx)

    # Check if it is possible to move to all positions
    pos_bool = np.zeros((t_yx.shape[0]))
    for i in range(t_yx.shape[0]):
        if np.any(positions.pos[((positions.pos[:, 1] == t_yx[i,1]) & (positions.pos[:, 2] == t_yx[i,0])), 0]):
            pos_bool[i] = 1

    pos_bool = np.bool_(pos_bool)
    t_yx = t_yx[pos_bool,:]

    # Where not to move, i.e. positions of all other trees
    not_t = ((x[:] != t_x) | (y[:] != t_y))

    y_n = y[not_t]    # y position of trees that are not moving
    x_n = x[not_t]    # x position of trees that are not moving

    # Check euclidean distance between  moving tree and none-moving trees, i.e. where possible to move.
    # Also checking for canopy diameter / 2 to see if tree fits
    eucl = np.zeros((y_n.shape[0],t_yx.shape[0]))
    e_bool = np.ones((t_yx.shape[0]))
    e_bool_sh = np.zeros((t_yx.shape[0]))
    for i in range(y_n.shape[0]):
        yx_temp = np.array((y_n[i], x_n[i]))
        for j in range(t_yx.shape[0]):
            eucl[i,j] = np.linalg.norm(t_yx[j,:] - yx_temp)
            if (eucl[i,j] < dia):
                e_bool[j] = 0
            if (eucl[i,j] < treerasters.euclidean_d):
                e_bool_sh[j] = 1

    e_bool = np.bool_(e_bool)   # Boolean of where it's possible and not to move
    t_yx = t_yx[e_bool,:]  # Positions where it's possible to move to
    e_bool_sh = np.bool_(e_bool_sh)

    # Checking euclidean distance between none-moving trees
    yx_nmt = tuple(itertools.combinations(np.arange(y_n.shape[0]), 2))
    e_bool_nmt = np.zeros((yx_nmt.__len__()))
    e_nmt = np.zeros((yx_nmt.__len__()))
    for i in range(e_bool_nmt.shape[0]):
        ax = x_n[yx_nmt[i][0]]
        ay = y_n[yx_nmt[i][0]]
        a = np.array(([ay, ax]))
        bx = x_n[yx_nmt[i][1]]
        by = y_n[yx_nmt[i][1]]
        b = np.array(([by, bx]))
        e_nmt[i] = np.linalg.norm(a - b)
        if (e_nmt[i] < treerasters.euclidean_d):
            e_bool_nmt[i] = 1

    e_bool_nmt = np.bool_(e_bool_nmt)

    tree_tmrt = np.zeros((t_yx.shape[0] + 1, 5))  # y, x, tmrt shade, tmrt sun, tmrt diff, sum tmrt all trees
    tree_tmrt[-1, :] = p_tmrt

    compare_mt = np.zeros((t_yx.shape[0]))
    compare = 0
    if ((np.any(e_bool_nmt == 1)) | (np.any(e_bool_sh == 1))):

    # Check if tree shadows overlap
        treesh_ts_bool_pad, treesh_ts_bool_pad_large, compare = tsh_gen_ts(y_n, x_n, treerasters, treeinput)  # Boolean tree shadows for trees that are not moving

        if (compare == 0):    # If none of the none-moving trees overlap, go in here
            treesh_bool_mt = tsh_gen_mt1(t_yx[:,0], t_yx[:,1], treerasters, treeinput)
            if np.any(((treesh_bool_mt == 1) & (treesh_ts_bool_pad_large == 1))):
                for j in range(t_yx.shape[0]):
                    if (e_bool_sh[j] == 1):   # If boolean distance between coming position of the moving tree possibly have overlapping shadows with the other trees, continue
                        y_t = np.array([t_yx[j, 0]])    # Y-position of the current position of the currently moving tree
                        x_t = np.array([t_yx[j, 1]])    # X-position of the current position of the currently moving tree
                        treesh_bool_mt_pad, treesh_bool_mt_pad_large, _ = tsh_gen_ts(y_t, x_t, treerasters, treeinput)    # Boolean tree shadows for moving tree
                        if np.any(((treesh_ts_bool_pad_large == 1) & (treesh_bool_mt_pad_large == 1))):
                            for ij in range(treesh_bool_mt_pad.shape[2]):   # Check if timesteps overlap
                                temp_bool = ((treesh_ts_bool_pad[:,:,ij] == 1) | (treesh_bool_mt_pad[:,:,ij] == 1)) * (treeinput.shadow[:,:,ij] == 1) * (treeinput.buildings == 1)
                                tree_tmrt[j,0] += np.sum(temp_bool * tmrt_1d[ij, 0])
                                tree_tmrt[j,1] += np.sum(temp_bool * treeinput.tmrt_ts[:,:,ij])
                                compare_mt[j] = 1

                    tree_tmrt[j,2] = tree_tmrt[j,1] - tree_tmrt[j,0]

        else:
            for j in range(t_yx.shape[0]):
                if (e_bool_sh[j] == 1):
                    y_t = np.array([t_yx[j, 0]])
                    x_t = np.array([t_yx[j, 1]])
                    treesh_bool_mt_pad = tsh_gen_mt2(y_t, x_t, treerasters, treeinput)    # Boolean tree shadows for moving tree
                    for ij in range(treesh_bool_mt_pad.shape[2]):
                        temp_bool = ((treesh_ts_bool_pad[:, :, ij] == 1) | (treesh_bool_mt_pad[:, :, ij] == 1)) * (
                                    treeinput.shadow[:, :, ij] == 1) * (treeinput.buildings == 1)
                        tree_tmrt[j, 0] += np.sum(temp_bool * tmrt_1d[ij, 0])
                        tree_tmrt[j, 1] += np.sum(temp_bool * treeinput.tmrt_ts[:, :, ij])
                        compare_mt[j] = 1
                tree_tmrt[j, 2] = tree_tmrt[j, 1] - tree_tmrt[j, 0]

    # Calculation of shadows for the currently moving tree
    for i in range(t_yx.shape[0]):
        y_t = np.array([t_yx[i,0]])
        x_t = np.array([t_yx[i,1]])

        tree_tmrt[i, 3] = y_t               # y position of currently moving tree
        tree_tmrt[i, 4] = x_t               # x position of  currently moving tree

        # Comparison between the currently moving tree and other trees if their shadows overlap
        # If moving tree is overlapping with any of the other trees

        # If any of the none-moving trees are overlapping with each other but the moving tree is not overlapping with them
        if ((compare == 1) & (compare_mt[i] == 0)):
            start_time = time.time()
            for j in range(treesh_ts_bool_pad.shape[2]):
                tree_bool_temp = treesh_ts_bool_pad[:,:,j] * treeinput.shadow[:,:,j] * treeinput.buildings
                tree_tmrt[i, 0] += np.sum(tree_bool_temp * tmrt_1d[j, 0])
                tree_tmrt[i, 1] += np.sum(tree_bool_temp * treeinput.tmrt_ts[:,:,j])
            tree_tmrt[i, 0] += treerasters.tmrt_shade[y_t, x_t]
            tree_tmrt[i, 1] += treerasters.tmrt_sun[y_t, x_t]
            tree_tmrt[i, 2] = tree_tmrt[i, 1] - tree_tmrt[i, 0]

        # If no trees are overlapping
        elif ((compare == 0) & (compare_mt[i] == 0)):
            start_time = time.time()
            for j in range(y_n.shape[0]):
                tree_tmrt[i, 0] += treerasters.tmrt_shade[y_n[j], x_n[j]]   # Tree shade
                tree_tmrt[i, 1] += treerasters.tmrt_sun[y_n[j], x_n[j]]     # Sunlit
                tree_tmrt[i, 2] += treerasters.d_tmrt[y_n[j], x_n[j]]       # Sunlit - Tree shade
            tree_tmrt[i,0] += treerasters.tmrt_shade[y_t,x_t]    # Potential decrease in Tmrt from shade due to other trees
            tree_tmrt[i,1] += treerasters.tmrt_sun[y_t,x_t]      # Tmrt in sun for the tree shadow
            tree_tmrt[i,2] += treerasters.d_tmrt[y_t,x_t]        # Difference in Tmrt between sunlit and shaded

    nc = 0 # no change

    if (np.sum(e_bool) > 0):
        t_max = np.where(tree_tmrt[:, 2] == np.max(tree_tmrt[:, 2]))          # Which position has highest tmrt
        if (t_max[0].__len__() > 1):
            t_rand = np.zeros((t_max[0].__len__()))
            for iy in range(t_max[0].__len__()):
                t_out = tree_tmrt[[t_max[0][iy]],:]
                if ((t_out[0, 3] == t_y) & (t_out[0, 4] == t_x)):
                    t_rand[iy] = 1
                    y[ti] = t_out[0, 3]
                    x[ti] = t_out[0, 4]
                    nc = 1
                    break
            if (np.sum(t_rand) == 0):
                t_max_rand = np.random.choice(t_max[0], 1)
                t_out = tree_tmrt[t_max_rand, :]
                y[ti] = t_out[0, 3]
                x[ti] = t_out[0, 4]
        else:
            t_out = tree_tmrt[t_max[0],:]                                       # Position of where tmrt is highest
            if ((np.int_(t_out[0,3]) == t_y) & (np.int_(t_out[0,4]) == t_x)):   # If starting position, i.e. input position is best, don't move
                nc = 1
            else:
                y[ti] = t_out[0,3]
                x[ti] = t_out[0,4]
    else:
        t_out = tree_tmrt                             # If tree can't move because of other trees, no move
        nc = 1

    return t_out, nc, y, x
