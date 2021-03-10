import numpy as np
import itertools

# Creating trees at random positions
def random_start(pos, trees, tree_pos_all, r_iters):

    tree_pos = np.random.choice(pos, trees)  # Random positions for trees
    tp_c = 0
    break_loop = 0

    while np.any(np.where((tree_pos_all == np.sort(tree_pos)).all(axis=1))):
        tree_pos = np.random.choice(pos, trees)
        tp_c += 1
        # print('re-random')
        if (tp_c == 100):
            break_loop = r_iters + 1
            break

    return tree_pos, tp_c, break_loop

# Creating trees evolutionary from previous starting position. First population is random.
def genetic_start(tree_pos_x, tree_pos_y, tree_pos_c, positions, trees, tree_pos_all, r_iters, counter, dia):
    pos = positions.pos[:, 0]
    pos_y = positions.pos[:, 2]
    pos_x = positions.pos[:, 1]

    tree_pos = 0
    tp_c = 0
    break_loop = 0

    tree_pos = np.zeros((trees))

    # Random for first run
    if counter == 0:
        tree_pos, tp_c, break_loop = random_start(pos, trees, tree_pos_all, r_iters)

    # Take either x or y position from previous local optimum for each tree, make the other coordinate random
    else:
        distance = 0
        tries = np.zeros((trees))
        while distance < 1:
            exists = np.zeros((trees))
            ti = itertools.cycle(range(trees))  # Iterator to move between trees moving around in the study area
            new_y = np.zeros((trees))
            new_x = np.zeros((trees))
            tree_pos_c_temp = tree_pos_c.copy()
            while np.sum(exists) < trees:
                i = ti.__next__()
                # Random y or x coordinate i.e. mutation
                if ((tree_pos_c[i] > 3) or (tries[i] > 50)):
                    # Decide which coordinate will be random; 0 = y, 1 = x
                    xy_random = np.random.choice([0, 1], 1)
                    # Random y-position
                    if xy_random == 0:
                        new_x[i] = tree_pos_x[i]
                        pos_y_temp = pos_y[pos_x == new_x[i]]
                        new_y[i] = np.random.choice(pos_y_temp, 1)
                    # Random x-position
                    else:
                        new_y[i] = tree_pos_y[i]
                        pos_x_temp = pos_x[pos_y == new_y[i]]
                        new_x[i] = np.random.choice(pos_x_temp, 1)
                    tree_pos[i] = positions.pos[
                        ((positions.pos[:, 1] == new_x[i]) & (positions.pos[:, 2] == new_y[i])), 0]
                    exists[i] = 1
                    tree_pos_c_temp[i] = 0
                    tries[i] = 0
                # Parents
                elif exists[i] == 0:
                    y_temp = np.random.choice(tree_pos_y, 1)
                    x_temp = np.random.choice(tree_pos_x, 1)
                    if np.any(positions.pos[((positions.pos[:, 1] == x_temp) & (positions.pos[:, 2] == y_temp)), 0]):
                        tree_pos[i] = positions.pos[
                            ((positions.pos[:, 1] == x_temp) & (positions.pos[:, 2] == y_temp)), 0]
                        new_y[i] = y_temp
                        new_x[i] = x_temp
                        exists[i] = 1

            # Euclidean distance between random positions so that trees are not too close to each other
            yxp = tuple(itertools.combinations(np.arange(tree_pos.shape[0]), 2))
            eucl_dist = np.zeros((yxp.__len__()))

            for j in range(yxp.__len__()):
                ax = new_x[yxp[j][0]]
                ay = new_y[yxp[j][0]]
                a = np.array(([ay, ax]))
                bx = new_x[yxp[j][1]]
                by = new_y[yxp[j][1]]
                b = np.array(([by, bx]))
                eucl_dist[j] = np.linalg.norm(a - b)

            if (np.min(eucl_dist) >= dia):
                distance = 1

            tries[i] += 1

        tree_pos_c = tree_pos_c_temp
        tree_pos_y = new_y
        tree_pos_x = new_x

    tree_pos_y = np.int_(tree_pos_y)
    tree_pos_x = np.int_(tree_pos_x)

    return tree_pos_y, tree_pos_x, tree_pos, tree_pos_c, tp_c, break_loop

### Creating trees evolutionary from previous starting position. First population is random.
##def geneticstart(tree_pos_x, tree_pos_y, tree_pos_c, positions, trees, tree_pos_all, r_iters, counter):
##
##    pos = positions.pos[:, 0]
##    pos_y = positions.pos[:, 2]
##    pos_x = positions.pos[:, 1]
##
##    tree_pos = 0
##    tp_c = 0
##    break_loop = 0
##
##    tree_pos = np.zeros((trees))
##
##    # Random for first run
##    if counter == 0:
##        tree_pos, tp_c, break_loop = randomstart(pos, trees, tree_pos_all, r_iters)
##
##    # Take either x or y position from previous local optimum for each tree, make the other coordinate random
##    else:
##        for i in range(trees):
##            if tree_pos_c[i] < 4:
##                # Decide which coordinate will be random; 0 = y, 1 = x
##                xy_random = np.random.choice([0,1], 1)
##                # Random y-position
##                if xy_random == 0:
##                    x_temp = tree_pos_x[i]
##                    pos_y_temp = pos_y[pos_x == x_temp]
##                    tree_pos_y[i] = np.random.choice(pos_y_temp, 1)
##                # Random x-position
##                else:
##                    y_temp = tree_pos_y[i]
##                    pos_x_temp = pos_x[pos_y == y_temp]
##                    tree_pos_x[i] = np.random.choice(pos_x_temp, 1)
##            # If the tree have ended up in local optimums with lower Tmrt x times, make both coordinates random
##            else:
##                tree_pos_temp = np.random.choice(pos, 1)
##                tree_pos_y[i] = pos_y[pos == tree_pos_temp]  # Random y-position
##                tree_pos_x[i] = pos_x[pos == tree_pos_temp]  # Random x-position
##                tree_pos_c[i] = 0
##
##            tree_pos[i] = positions.pos[((positions.pos[:,1] == tree_pos_x[i]) & (positions.pos[:,2] == tree_pos_y[i])), 0]
##
##        tree_pos_y = np.int_(tree_pos_y)
##        tree_pos_x = np.int_(tree_pos_x)
##
##    return tree_pos_y, tree_pos_x, tree_pos, tree_pos_c, tp_c, break_loop
