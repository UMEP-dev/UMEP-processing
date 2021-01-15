import numpy as np
import time

# This function returns a raster with a boolean shadow and the regional group shadow for the tree in position y,x
def tsh_gen(y,x,treerasters,treeinput):
    tsh_pos_pad = np.zeros((treeinput.rows, treeinput.cols, y.__len__()))
    tsh_pos_bool_pad = np.zeros((treeinput.rows, treeinput.cols, y.__len__()))
    tsh_pos_bool_pad_all = np.zeros((treeinput.rows, treeinput.cols))
    tsh_pos_bool_pad_all_ = np.zeros((treeinput.rows, treeinput.cols))

    for i in range(y.__len__()):
        y1 = np.int_(y[i] - treerasters.buffer_y[0])
        y2 = np.int_(y[i] + treerasters.buffer_y[1])
        x1 = np.int_(x[i] - treerasters.buffer_x[0])
        x2 = np.int_(x[i] + treerasters.buffer_x[1])

        tsh_pos_pad[y1:y2,x1:x2,i] = treerasters.treeshade_rg
        tsh_pos_pad[:, :, i] = tsh_pos_pad[:,:,i] * treeinput.buildings_pad
        tsh_pos_bool_pad[:,:,i] = tsh_pos_pad[:,:,i] > 0
        tsh_pos_bool_pad_all += tsh_pos_bool_pad[:,:,i]

    compare = 0
    if y.__len__() > 1:
        if np.any(tsh_pos_bool_pad_all > 1):
            tsh_pos_bool_pad_all_ = tsh_pos_bool_pad_all > 1
            compare = 1

    tsh_pos_bool_pad_all = tsh_pos_bool_pad_all > 0

    return tsh_pos_pad, tsh_pos_bool_pad, tsh_pos_bool_pad_all, tsh_pos_bool_pad_all_, compare

def tsh_gen_ts(y,x,treerasters,treeinput):
    tsh_bool_pad = np.zeros((treeinput.rows, treeinput.cols, treerasters.treeshade_bool.shape[2]))
    tsh_bool_pad_large = np.zeros((treeinput.rows, treeinput.cols))
    compare = 0

    for i in range(y.__len__()):
        y1 = np.int_(y[i] - treerasters.buffer_y[0])
        y2 = np.int_(y[i] + treerasters.buffer_y[1])
        x1 = np.int_(x[i] - treerasters.buffer_x[0])
        x2 = np.int_(x[i] + treerasters.buffer_x[1])

        tsh_bool_pad[y1:y2, x1:x2, :] += treerasters.treeshade_bool[:, :, :]
        tsh_bool_pad_large[y1:y2, x1:x2] += treerasters.treeshade_rg > 0

    if y.__len__() > 1:
        if np.any(tsh_bool_pad > 1):
            compare = 1

    tsh_bool_pad = tsh_bool_pad > 0
    tsh_bool_pad_large = tsh_bool_pad_large > 0

    return tsh_bool_pad, tsh_bool_pad_large, compare

def tsh_gen_mt1(y,x,treerasters,treeinput):
    tsh_bool_pad_large = np.zeros((treeinput.rows, treeinput.cols))

    for i in range(y.__len__()):
        y1 = np.int_(y[i] - treerasters.buffer_y[0])
        y2 = np.int_(y[i] + treerasters.buffer_y[1])
        x1 = np.int_(x[i] - treerasters.buffer_x[0])
        x2 = np.int_(x[i] + treerasters.buffer_x[1])

        tsh_bool_pad_large[y1:y2, x1:x2] += treerasters.treeshade_rg > 0

    tsh_bool_pad_large = tsh_bool_pad_large > 0

    return tsh_bool_pad_large

def tsh_gen_mt2(y,x,treerasters,treeinput):
    tsh_bool_pad = np.zeros((treeinput.rows, treeinput.cols, treerasters.treeshade_bool.shape[2]))

    for i in range(y.__len__()):
        y1 = np.int_(y[i] - treerasters.buffer_y[0])
        y2 = np.int_(y[i] + treerasters.buffer_y[1])
        x1 = np.int_(x[i] - treerasters.buffer_x[0])
        x2 = np.int_(x[i] + treerasters.buffer_x[1])

        tsh_bool_pad[y1:y2, x1:x2, :] += treerasters.treeshade_bool[:, :, :]

    tsh_bool_pad = tsh_bool_pad > 0

    return tsh_bool_pad