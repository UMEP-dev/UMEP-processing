try:
    import torch
except:
    pass


def Tgmaps_v1(lc_grid, solweig_parameters):

    # Tgmaps_v1 Populates grids with cooeficients for Tg wave
    #   Detailed explanation goes here
    lc_grid[lc_grid >= 100] = 2
    id = torch.unique(lc_grid)
    id = lc_grid[lc_grid <= 7].to(int)
    TgK = torch.clone(lc_grid)
    Tstart = torch.clone(lc_grid)
    alb_grid = torch.clone(lc_grid)
    emis_grid = torch.clone(lc_grid)
    TmaxLST = torch.clone(lc_grid)

    for i in id:
        # row = (lc_class[:, 0] == id[i])
        Tstart[Tstart == i] = solweig_parameters["Tstart"]["Value"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ]
        alb_grid[alb_grid == i] = solweig_parameters["Albedo"]["Effective"][
            "Value"
        ][solweig_parameters["Names"]["Value"][str((int(i.item())))]]
        emis_grid[emis_grid == i] = solweig_parameters["Emissivity"]["Value"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ]
        TmaxLST[TmaxLST == i] = solweig_parameters["TmaxLST"]["Value"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ]
        TgK[TgK == i] = solweig_parameters["Ts_deg"]["Value"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ]

    TgK_wall = solweig_parameters["Ts_deg"]["Value"]["Walls"]
    Tstart_wall = solweig_parameters["Tstart"]["Value"]["Walls"]
    TmaxLST_wall = solweig_parameters["TmaxLST"]["Value"]["Walls"]

    return (
        TgK,
        Tstart,
        alb_grid,
        emis_grid,
        TgK_wall,
        Tstart_wall,
        TmaxLST,
        TmaxLST_wall,
    )
