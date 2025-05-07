import numpy as np

def Tgmaps_v1(lc_grid, solweig_parameters):

    #Tgmaps_v1 Populates grids with cooeficients for Tg wave
    #   Detailed explanation goes here
    lc_grid[lc_grid >= 100] = 2
    id = np.unique(lc_grid)
    id = lc_grid[lc_grid <= 7].astype(int)
    TgK = np.copy(lc_grid)
    Tstart = np.copy(lc_grid)
    alb_grid = np.copy(lc_grid)
    emis_grid = np.copy(lc_grid)
    TmaxLST = np.copy(lc_grid)

    for i in id:
        # row = (lc_class[:, 0] == id[i])
        Tstart[Tstart == i] = solweig_parameters['Tstart']['Value'][solweig_parameters['Names']['Value'][str(i)]]
        alb_grid[alb_grid == i] = solweig_parameters['Albedo']['Effective']['Value'][solweig_parameters['Names']['Value'][str(i)]]
        emis_grid[emis_grid == i] = solweig_parameters['Emissivity']['Value'][solweig_parameters['Names']['Value'][str(i)]]
        TmaxLST[TmaxLST == i] = solweig_parameters['TmaxLST']['Value'][solweig_parameters['Names']['Value'][str(i)]]
        TgK[TgK == i] = solweig_parameters['Ts_deg']['Value'][solweig_parameters['Names']['Value'][str(i)]]

    TgK_wall = solweig_parameters['Ts_deg']['Value']['Walls']
    Tstart_wall = solweig_parameters['Tstart']['Value']['Walls']
    TmaxLST_wall = solweig_parameters['TmaxLST']['Value']['Walls']

    return TgK, Tstart, alb_grid, emis_grid, TgK_wall, Tstart_wall, TmaxLST, TmaxLST_wall
