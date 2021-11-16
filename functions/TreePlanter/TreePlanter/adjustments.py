import numpy as np
from scipy.ndimage import label

from ..TreePlanter.TreePlanterTreeshade import tsh_gen_ts

def treenudge(y_out, x_out, tp_nc, tp_nc_a, i, t1, counter, i_tmrt, treerasters, treeinput, positions, tmrt_1d):

    tree_adjusted = 0

    if (np.sum(tp_nc[:, 0]) > 1):  # If more than one tree is standing still, check if shadows are next to each other
        nc_r = np.where(tp_nc[:, 0] == 1)  # Rows of trees standing still
        nc_y = np.squeeze(y_out[nc_r])  # Y-positions of trees
        nc_x = np.squeeze(x_out[nc_r])  # X-positions of trees

        tsh_, tsh_bool, compare, = tsh_gen_ts(nc_y, nc_x, treerasters, treeinput)  # Shadows of trees

        if (compare == 1):  # If there are less regional groups than trees standing still, some must be side by side
            label_nc, num_feat_nc = label(tsh_bool)  # Regional groups

            nc_y_c = np.squeeze(y_out.copy())
            nc_x_c = np.squeeze(x_out.copy())

            # Find out label, i.e. regional group of shadows
            nc_label = np.zeros((y_out.shape[0]))
            for iy in range(nc_y.shape[0]):
                nc_y_t = np.int_(y_out[iy])
                nc_x_t = np.int_(x_out[iy])
                nc_label[iy] = label_nc[nc_y_t, nc_x_t]

            nc_i_y = np.int_(np.squeeze(y_out[i]))
            nc_i_x = np.int_(np.squeeze(x_out[i]))
            nc_i = label_nc[nc_i_y, nc_i_x]
            # Only look at the regional group of the currently moving tree
            nc_i_r = np.where(nc_label == nc_i)

            nc_sum = np.sum(tp_nc[nc_i_r, 0])  # Calculate the sum stuck trees in rg nc_i (should be all if adjustment should be made)

            if ((nc_i_r[0].__len__() > 1) & (nc_sum == nc_i_r[0].__len__())):
                y_dir = np.array([1, -1, 0, 0])  # New temp y-position
                x_dir = np.array([0, 0, 1, -1])  # New temp x-position

                dir_bool = np.ones((y_dir.shape[0]))
                # Check if it is possible to move in all directions, i.e. are these positions possible for trees
                for iy in range(y_dir.shape[0]):
                    y_dir_temp = np.squeeze(y_out[nc_i_r]) + y_dir[iy]
                    x_dir_temp = np.squeeze(x_out[nc_i_r]) + x_dir[iy]
                    for iy1 in range(y_dir_temp.shape[0]):
                        if not (np.any(positions.pos[(positions.pos[:, 1] == x_dir_temp[iy1]) & (
                                positions.pos[:, 2] == y_dir_temp[iy1]), 0])):
                            dir_bool[iy] = 0

                dir_bool = np.bool_(dir_bool)
                y_dir = y_dir[dir_bool]  # Remove positions where it is not possible to put a tree
                x_dir = x_dir[dir_bool]  # Remove positions where it is not possible to put a tree

                tmrt_nc = np.zeros((y_dir.shape[0], 3))  # New tmrt
                if (np.any(y_dir.shape[0])):  # Only go in to if if there are any possible locations for trees
                    for iy in range(y_dir.shape[0]):  # Loop for shadows for new positions
                        nc_y_c[nc_i_r] = np.squeeze(y_out[nc_i_r]) + y_dir[iy]  # New y-position
                        nc_x_c[nc_i_r] = np.squeeze(x_out[nc_i_r]) + x_dir[iy]  # New x-position

                        tsh_bool_nc_t, tsh_bool_nc_t_large, comp_ = tsh_gen_ts(nc_y_c, nc_x_c, treerasters,
                                                                               treeinput)  # Shadows and rg for new position
                        # Calculate Tmrt for intersecting shadows
                        for j in range(tmrt_1d.shape[0]):
                            tsh_bool_nc_t[:, :, j] = tsh_bool_nc_t[:, :, j] * treeinput.shadow[:, :,
                                                                              j] * treeinput.buildings
                            tmrt_nc[iy, 0] += np.sum(tsh_bool_nc_t[:, :, j] * tmrt_1d[j, 0])  # Tree shade
                            tmrt_nc[iy, 1] += np.sum(tsh_bool_nc_t[:, :, j] * treeinput.tmrt_ts[:, :, j])  # Sunlit
                        tmrt_nc[iy, 2] = tmrt_nc[iy, 1] - tmrt_nc[iy, 0]  # New Tmrt (sunlit - tree shade)

                    if (np.around(np.max(tmrt_nc[:, 2]), decimals=1) > np.around(i_tmrt[counter], decimals=1)):  # If any new Tmrt decrease is higher, continue
                        nc_max_r = np.where(tmrt_nc[:, 2] == np.max(tmrt_nc[:, 2]))  # Where is new Tmrt decrease highest

                        if (nc_max_r[0].__len__() > 1):  # If more than one positions has max tmrt, choose a random of the positions
                            nc_max_r = np.random.choice(nc_max_r[0], 1)

                        nc_y_out = np.squeeze(y_out[nc_i_r]) + y_dir[nc_max_r]  # New y-positions
                        nc_x_out = np.squeeze(x_out[nc_i_r]) + x_dir[nc_max_r]  # New x-positions
                        nc_tmrt_out = np.max(tmrt_nc[nc_max_r, 2])  # New Tmrt

                        y_out[nc_i_r] = nc_y_out  # New y-positions
                        x_out[nc_i_r] = nc_x_out  # New x-positions

                        t1[0, 2] = nc_tmrt_out  # New Tmrt
                        t1[0, 3] = y_out[i]
                        t1[0, 4] = x_out[i]

                        tp_nc[nc_i_r, 0] = 0  # Reset tp_nc, i.e. trees have moved
                        tp_nc_a[nc_i_r] = 1  # Reset tp_nc_a, i.e. trees have been adjusted

    return y_out, x_out, tp_nc, tp_nc_a, t1

## Adjusting tree with lowest Tmrt
def tree_adjust(i_y, i_x, i_tmrt, counter, trees, treerasters, treeinput):

    a_counter = 0   # Adjustment counter

    while a_counter < trees:
        tmrt_temp = np.zeros((trees))
        for ix in range(trees):
            tmrt_temp[ix] = treerasters.d_tmrt[np.int(i_y[counter,ix]),np.int(i_x[counter,ix])]

        y_min = tmrt_temp[:] == np.min(tmrt_temp[:])    # Bool of position of tree with least decrease in Tmrt
        if y_min.shape[0] > 1:
            y_min = y_min.cumsum(axis=0).cumsum(axis=0) == 1    # If more than one index with np.min, return only first
        y_max = tmrt_temp[:] != np.min(tmrt_temp[:])    # Bool of positions of the other trees (highest)
        y_te = i_y[counter,:]   # Current y-positions of all trees for current run
        x_te = i_x[counter,:]   # Current x-positions of all trees for current run
        y_high = np.int_(y_te[y_max])   # y-positions of trees with highest decrease in Tmrt
        x_high = np.int_(x_te[y_max])   # x-positions of trees with highest decrease in Tmrt
        y_low = np.int_(y_te[y_min])   # y-position of tree with lowest decrease in Tmrt
        x_low = np.int_(x_te[y_min])   # x-position of tree with lowest decrease in Tmrt

        tsh_rg_temp, tsh_bool_temp, comp_ = tsh_gen_ts(y_high, x_high, treerasters, treeinput)  # Shadows for high trees
        tsh_rg_min, tsh_bool_min, comp__ = tsh_gen_ts(y_low, x_low, treerasters, treeinput) # Shadows for low tree

        d_tmrt_pad = np.pad(treerasters.d_tmrt, pad_width=((treerasters.tpy[0], treerasters.tpy[0]), (treerasters.tpx[0], treerasters.tpx[0])), mode='constant',
                                  constant_values=0)    # Padded sun vs shade tmrt rasters

        tsh_bool_all = tsh_bool_temp == 0   # All pixels that are sunlit

        d_tmrt_pad = treeinput.d_tmrt * tsh_bool_all  # Remaining pixels, i.e. pixels that are not shaded for sun vs. shade

        tsh_bool_all = np.bool_(1 - tsh_bool_all)   # Reverse tsh_bool_all

        # See if old trees overlap, i.e. the tree that is to be adjusted overlaps with the other trees
        if np.any((tsh_bool_all == 1) & (tsh_bool_min == 1)):  # If not, proceed. If they overlap, adjustment for overlap needs to be made, etc.

            d_tmrt_vec = d_tmrt_pad.flatten()   # Flatten d_tmrt_pad and create vector

            d_tmrt_vec_s = -np.sort(-d_tmrt_vec)    # Sort vector to find pixels with highest difference in Tmrt (sun vs. shade)

            d_bool = d_tmrt_vec_s > np.min(tmrt_temp[:])    # All pixels that have higher difference in Tmrt compared to position of tree with lowest difference (decrease)
            d_tmrt_vec_s = d_tmrt_vec_s[d_bool] # Remove all other pixels

            d_tmrt_vec_s = np.unique(d_tmrt_vec_s)  # Only save unique values

            a_nc = 0    # Adjustment change parameter
            for ix in range(d_tmrt_vec_s.shape[0]):
                tmrt_adjust = d_tmrt_vec_s[ix]  # Trying values from d_tmrt_vec_s
                y_adjust, x_adjust = np.where(d_tmrt_pad == tmrt_adjust)    # Find coordinates for tmrt_adjust
                for iy in range(y_adjust.shape[0]):     # If there are more than one position with tmrt_adjust
                    y_temp = np.array([y_adjust[iy] - treerasters.tpy[0]])  # y-Position of tmrt_adjust
                    x_temp = np.array([x_adjust[iy] - treerasters.tpx[0]])  # x-Position of tmrt_adjust
                    tsh_rg_adjust, tsh_bool_adjust, comp___ = tsh_gen_ts(y_temp, x_temp, treerasters, treeinput)  # Create boolean shadow from y_temp, x_temp
                    if not np.any((tsh_bool_all == 1) & (tsh_bool_adjust == 1)):  # If the shadow of the new position does not overlap with the the other trees, proceed
                        y_te[y_min] = y_temp
                        x_te[y_min] = x_temp
                        i_y[counter, :] = y_te  # New y-position
                        i_x[counter, :] = x_te  # New x-position
                        i_tmrt[counter] += (tmrt_adjust - np.min(tmrt_temp))    # New adjusted Tmrt
                        a_nc = 1
                        print('Adjusted')
                        break
                    else:
                        break
                if a_nc == 1:
                    break
        else:
            a_counter = trees
            a_nc = 1
            break
        a_counter += 1

    return i_y, i_x, i_tmrt