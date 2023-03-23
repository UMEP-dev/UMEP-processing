#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 14:57:25 2021

# Sandro Oswald, sandro.oswald@boku.ac.at
# Vienna Urban Climate Group
# University of natural sciences (BOKU)

# Jérémy Bernard, jeremy.bernard@zaclys.net
# University of Gothenburg
# Department of Earth Sciences
"""
import numpy as np
import time
from .GlobalVariables import MAX_ITERATIONS, THRESHOLD_ITERATIONS, DESCENDING_Y
from numba import jit
import pandas as pd

def solver(x, y, z, dx, dy, dz, u0, v0, w0, buildingCoordinates, cells4Solver, cursor,
           maxIterations = MAX_ITERATIONS, thresholdIterations = THRESHOLD_ITERATIONS,
           feedback = None):
    """ Use the mass-balance solver minimizing the modification of the initial
    wind speed field. The method used is based on Pardyjak and Brown (2003).
    
    References:
        Pardyjak, Eric R, et Michael Brown. « QUIC-URB v. 1.1: Theory and 
        User’s Guide ». Los Alamos National Laboratory, Los Alamos, NM, 2003.

    
    		Parameters
    		_ _ _ _ _ _ _ _ _ _ 
    
            x: 1D array
                X-axis cell coordinates in a local reference system (starting from 0)
            y: 1D array
                Y-axis cell coordinates in a local reference system (starting from 0)   
            z: 1D array
                Z-axis cell coordinates in a local reference system (starting from 0)
            dx: int
                Grid spacing along X-axis
            dy: int
                Grid spacing along Y-axis  
            dz: int
                Grid spacing along Z-axis
            u0: 3D array
                Initialized 3D wind speed value in X direction
            v0: 3D array
                Initialized 3D wind speed value in Y direction 
            w0: 3D array
                Initialized 3D wind speed value in Z direction
            buildingCoordinates: 3D array
                Building 3D coordinates
            cells4Solver: 1D array
                Array of 3D cell coordinates for which the wind solver is applied
            maxIterations: int, default MAX_ITERATIONS
                Maximum number of wind solver iterations (solver stops if reached)
            thresholdIterations: float, default THRESHOLD_ITERATIONS
                Threshold for stopping wind solver: when the relative 
                variation of lambda between 2 iterations goes under this 
                threshold, the wind solver stops
            feedback: Qgis.core class QgsProcessingFeedback
                Base class for providing feedback to QGIS from a processing algorithm (if not in standalone mode).
        
    		Returns
    		_ _ _ _ _ _ _ _ _ _ 
    
            u: 3D array
                Updated 3D wind speed value in X direction
            v: 3D array
                Updated 3D wind speed value in Y direction 
            w: 3D array
                Updated 3D wind speed value in Z direction"""    

    print("Start to apply the wind solver")
    timeStartCalculation = time.time()

    # Get number of cells in each 3 dimensions
    nx = x.size
    ny = y.size
    nz = z.size
    
    # Create empty matrix for the 3D wind speed calculation
    u = np.zeros((nx, ny, nz))
    v = np.zeros((nx, ny, nz))
    w = np.zeros((nx, ny, nz))

    # Preallocating lambda and lambda + 1 and set values to 0 on sketch boundaries
    lambdaN = np.ones([nx, ny, nz])
    lambdaN1 = np.ones([nx, ny, nz])
    lambdaN[0, :, :] = 0.
    lambdaN[:, 0, :] = 0.
    lambdaN[:, :, 0] = 0.
    lambdaN[-1, :, :] = 0.
    lambdaN[:, -1, :] = 0.
    lambdaN[:, :, -1] = 0.
    lambdaN1[0, :, :] = 0.
    lambdaN1[:, 0, :] = 0.
    lambdaN1[:, :, 0] = 0.
    lambdaN1[-1, :, :] = 0.
    lambdaN1[:, -1, :] = 0.
    lambdaN1[:, :, -1] = 0.

    Xi = ((np.cos(np.pi / nx) + (dx / dy) ** 2 * np.cos(np.pi / ny)) / (1 + (dx / dy) ** 2)) ** 2

    omega = 2. * ((1 - np.sqrt(1 - Xi)) / Xi)
    if (omega < 1) or (omega > 2):
        omega = 1.78

    # Coefficients such as defined in Pardyjak et Brown (2003) 
    alpha1 = 1.
    alpha2 = 1.
    eta = alpha1 / alpha2
    A = dx ** 2 / dy ** 2
    B = eta ** 2 * dx ** 2 / dz ** 2

    # Set coefficients according to table 1 (Pardyjak et Brown, 2003) 
    # to modify the Equation near obstacles
    e = np.ones([nx, ny, nz])
    f = np.ones([nx, ny, nz])
    g = np.ones([nx, ny, nz])
    h = np.ones([nx, ny, nz])
    m = np.ones([nx, ny, nz])
    n = np.ones([nx, ny, nz])
    o = np.ones([nx, ny, nz])
    p = np.ones([nx, ny, nz])
    q = np.ones([nx, ny, nz])
    
    # Identify index having wall below AND front, left, right or behind
    indBelow = pd.MultiIndex.from_tuples(list(zip(*[buildingCoordinates[0], 
                                                    buildingCoordinates[1], 
                                                    buildingCoordinates[2] + 1])))
    indBelowFront = indBelow.intersection(pd.MultiIndex.from_tuples(list(zip(*[buildingCoordinates[0], 
                                                                               buildingCoordinates[1] - 1, 
                                                                               buildingCoordinates[2]]))))
    indBelowBehind = indBelow.intersection(pd.MultiIndex.from_tuples(list(zip(*[buildingCoordinates[0], 
                                                                               buildingCoordinates[1] + 1, 
                                                                               buildingCoordinates[2]]))))
    indBelowLeft = indBelow.intersection(pd.MultiIndex.from_tuples(list(zip(*[buildingCoordinates[0] + 1, 
                                                                               buildingCoordinates[1], 
                                                                               buildingCoordinates[2]]))))
    indBelowRight = indBelow.intersection(pd.MultiIndex.from_tuples(list(zip(*[buildingCoordinates[0] - 1, 
                                                                               buildingCoordinates[1], 
                                                                               buildingCoordinates[2]]))))
    indBelowAnyAround = indBelowFront.union(indBelowBehind).union(indBelowLeft).union(indBelowRight)
    indBelowLeftRight = indBelowLeft.union(indBelowRight)
    indBelowFrontBehind = indBelowFront.union(indBelowBehind)  
    
    # Go descending order along y
    if DESCENDING_Y:
        e[buildingCoordinates[0] + 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.
        e[indBelowLeft.get_level_values(0), indBelowLeft.get_level_values(1), indBelowLeft.get_level_values(2)] = 0.
        f[buildingCoordinates[0] - 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.
        f[indBelowRight.get_level_values(0), indBelowRight.get_level_values(1), indBelowRight.get_level_values(2)] = 0. 
        g[buildingCoordinates[0], buildingCoordinates[1] + 1, buildingCoordinates[2]] = 0.
        g[indBelowBehind.get_level_values(0), indBelowBehind.get_level_values(1), indBelowBehind.get_level_values(2)] = 0.
        h[buildingCoordinates[0], buildingCoordinates[1] - 1, buildingCoordinates[2]] = 0.  
        h[indBelowFront.get_level_values(0), indBelowFront.get_level_values(1), indBelowFront.get_level_values(2)] = 0.
        m[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] + 1] = 0.
        n[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] - 1] = 0.
    else:    
        e[buildingCoordinates[0] - 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.
        e[indBelowRight.get_level_values(0), indBelowRight.get_level_values(1), indBelowRight.get_level_values(2)] = 0.        
        f[buildingCoordinates[0] + 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.
        f[indBelowLeft.get_level_values(0), indBelowLeft.get_level_values(1), indBelowLeft.get_level_values(2)] = 0.    
        g[buildingCoordinates[0], buildingCoordinates[1] - 1, buildingCoordinates[2]] = 0.
        g[indBelowFront.get_level_values(0), indBelowFront.get_level_values(1), indBelowFront.get_level_values(2)] = 0.
        h[buildingCoordinates[0], buildingCoordinates[1] + 1, buildingCoordinates[2]] = 0.
        h[indBelowBehind.get_level_values(0), indBelowBehind.get_level_values(1), indBelowBehind.get_level_values(2)] = 0.
        m[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] - 1] = 0.
        n[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] + 1] = 0.
    
    o[buildingCoordinates[0] - 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.5
    o[buildingCoordinates[0] + 1, buildingCoordinates[1], buildingCoordinates[2]] = 0.5
    p[buildingCoordinates[0], buildingCoordinates[1] - 1, buildingCoordinates[2]] = 0.5
    p[buildingCoordinates[0], buildingCoordinates[1] + 1, buildingCoordinates[2]] = 0.5
    q[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] + 1] = 0.5
    q[buildingCoordinates[0], buildingCoordinates[1], buildingCoordinates[2] - 1] = 0.5
    n[indBelowAnyAround.get_level_values(0), indBelowAnyAround.get_level_values(1), indBelowAnyAround.get_level_values(2)] = 0.
    o[indBelowLeftRight.get_level_values(0), indBelowLeftRight.get_level_values(1), indBelowLeftRight.get_level_values(2)] = 0.5
    p[indBelowFrontBehind.get_level_values(0), indBelowFrontBehind.get_level_values(1), indBelowFrontBehind.get_level_values(2)] = 0.5
    q[indBelowAnyAround.get_level_values(0), indBelowAnyAround.get_level_values(1), indBelowAnyAround.get_level_values(2)] = 0.5
       
    for N in range(maxIterations):
        print("Iteration {0} (max {1})".format( N + 1, 
                                                maxIterations))
        lambdaN = np.copy(lambdaN1)

        lambdaN1 = calcLambda(cells4Solver, lambdaN, lambdaN1, omega, alpha1,
                              u0, v0, w0, dx, dy, dz, e, f, g, h, m, n, o, p, q,
                              DESCENDING_Y, A, B)
        
        # Calculate how much lambda evolves between 2 consecutive iterations                                      
        eps = np.sum(np.abs(lambdaN1 - lambdaN)) / np.sum(np.abs(lambdaN1))
        
        # Check if the condition for ending process is reached
        if eps < thresholdIterations:
            break
        else:
            print("   eps = {0} >= {1}".format(np.round(eps,6),
                                               thresholdIterations))
            # Feedback to QGIS every 50 iterations
            if (N % 50 == 0) & (feedback is not None):
                textToSend = """Iteration {0} (max {1}) - eps = {2} >= {3}
                            """.format( N + 1, 
                                        maxIterations,
                                        np.round(eps,6),
                                        thresholdIterations)
                feedback.setProgressText(textToSend)
                if feedback.isCanceled():
                    feedback.setProgressText("Calculation cancelled by user")
                    break
            
    
    # Calculates the final wind speed

    # go descending order along y
    if DESCENDING_Y:
        u[1:nx, :, :] = u0[1:nx, :, :] + 0.5 * (
                1. / (alpha1 ** 2)) * (lambdaN1[0:nx-1, :, :] - lambdaN1[1:nx, :, :]) / dx
        v[:, 1:ny, :] = v0[:, 1:ny, :] + 0.5 * (
                1. / (alpha1 ** 2)) * (lambdaN1[:, 0:ny-1, :] - lambdaN1[:, 1:ny, :]) / dy
        w[:, :, 1:nz] = w0[:, :, 1:nz] + 0.5 * (
                1. / (alpha2 ** 2)) * (lambdaN1[:, :, 0:nz - 1] - lambdaN1[:, :, 1:nz]) / dz
    else:
        u[1:nx, :, :] = u0[1:nx, :, :] + 0.5 * (
                1. / (alpha1 ** 2)) * (lambdaN1[1:nx, :, :] - lambdaN1[0:nx - 1, :, :]) / dx
        v[:, 1:ny, :] = v0[:, 1:ny, :] + 0.5 * (
                1. / (alpha1 ** 2)) * (lambdaN1[:, 1:ny, :] - lambdaN1[:, 0:ny - 1, :]) / dy
        w[:, :, 1:nz] = w0[:, :, 1:nz] + 0.5 * (
                1. / (alpha2 ** 2)) * (lambdaN1[:, :, 1:nz] - lambdaN1[:, :, 0:nz - 1]) / dz

    # Reset input and output wind speed to zero for building cells
    u[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    u[buildingCoordinates[0]+1,buildingCoordinates[1],buildingCoordinates[2]]=0
    v[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    v[buildingCoordinates[0],buildingCoordinates[1]+1,buildingCoordinates[2]]=0
    w[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]] = 0
    w[buildingCoordinates[0],buildingCoordinates[1],buildingCoordinates[2]+1]=0

    print("Time spent by the wind speed solver: {0} s".format(time.time()-timeStartCalculation))
    
    return u, v, w

@jit(nopython=True)
def calcLambda(cells4Solver, lambdaN, lambdaN1, omega, alpha1, u0, v0, w0, dx, dy, dz, e, f, g, h, m, n, o, p, q, DESCENDING_Y, A, B):
    # Go descending order along y
    if DESCENDING_Y:
        for k, j, i in np.flip(cells4Solver):
            lambdaN1[i, j, k] = omega * (
                ((-1.) * (dx ** 2 * (-2. * alpha1 ** 2) * (((u0[i, j, k] - u0[i + 1, j, k]) / (dx) + (
                        v0[i, j, k] - v0[i, j + 1, k]) / (dy) +
                                                            (w0[i, j, k] - w0[i, j, k + 1]) / (dz)))) + (
                          e[i, j, k] * lambdaN[i - 1, j, k] + f[i, j, k] * lambdaN1[i + 1, j, k] + A * (
                          g[i, j, k] * lambdaN[i, j - 1, k] + h[i, j, k] * lambdaN1[i, j + 1, k]) + B * (
                                  m[i, j, k] * lambdaN[i, j, k - 1] + n[i, j, k] * lambdaN1[i, j, k + 1]))) / (
                        2. * (o[i, j, k] + A * p[i, j, k] + B * q[i, j, k]))) + (1 - omega) * lambdaN1[i, j, k]  
                                      
    else:
        for i, j, k in cells4Solver:
            lambdaN1[i, j, k] = omega * (
                ((-1.) * (dx ** 2 * (-2. * alpha1 ** 2) * (((u0[i + 1, j, k] - u0[i, j, k]) / (dx) + (
                        v0[i, j + 1, k] - v0[i, j, k]) / (dy) +
                                                            (w0[i, j, k + 1] - w0[i, j, k]) / (dz)))) + (
                          e[i, j, k] * lambdaN[i + 1, j, k] + f[i, j, k] * lambdaN1[i - 1, j, k] + A * (
                          g[i, j, k] * lambdaN[i, j + 1, k] + h[i, j, k] * lambdaN1[i, j - 1, k]) + B * (
                                  m[i, j, k] * lambdaN[i, j, k + 1] + n[i, j, k] * lambdaN1[i, j, k - 1]))) / (
                        2. * (o[i, j, k] + A * p[i, j, k] + B * q[i, j, k]))) + (1 - omega) * lambdaN1[i, j, k]  
                                      
    return lambdaN1