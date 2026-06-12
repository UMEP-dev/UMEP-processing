"""
@author: Eliott Bridoux, University of Gothenburg
"""

import math

try:
    import torch
except:
    pass
# Stefan-Boltzmann s constant
SBC = 5.67e-8


def saturated_vp(T):
    """
    Used in the calculation of the surface temperature for the water bodies.

    :Parameters:
    :T: the temperature just above the water

    :Return:
    :qs: the saturated vapor pressure
    :slope: the slope of the tangent to the saturated vapor pressure-temperature curve
    """

    L = 2.260e6  # Latent heat of vaporisation
    R = 8.314  # Gas constant

    # August-Roche-Magnus approx. in Pa
    qs = 6109.4 * torch.exp(17.625 * T / (T + 243.04))

    # Clausius-Clapeyron equation
    slope = L * qs / R / (T + 273.15) ** 2

    return slope, qs


def initiate_groundScheme(
    lc_grid, solweig_parameters, day, Ta, location, device
):
    """
    Setup the maps used in the ground scheme calculations depending on the landcover

    :Parameters:
    :lc_grid: Array of landcover values
    :solweig_param: Dict of physical parameters
    :day: day of the year (int)
    :Ta: air temperature (float)
    :location: Dict containing latitude and longitude

    :Return: Array for the surface temperature, radiation flux and physical parameters
    """

    # Get the landcover data from lc_grid array
    lc_grid[lc_grid >= 100] = 2
    id = torch.unique(lc_grid)
    id = id.int()

    # Physical parameters grids
    cap_grid = torch.clone(lc_grid)  # Heat capacity
    diff_grid = torch.clone(lc_grid)  # Thermal diffusivity
    a1_grid = torch.clone(lc_grid)
    a2_grid = torch.clone(lc_grid)
    a3_grid = torch.clone(lc_grid)  # The 3 OHM coefficients

    # Initial fluxes and temperature grids
    Rn = torch.zeros_like(lc_grid)  # Net radiation
    Rn_past = torch.zeros_like(lc_grid)  # Stored net radiation
    G = torch.zeros_like(lc_grid)  # Ground heat flux
    Tg = torch.clone(lc_grid)  # Surface temperature
    Tm = torch.clone(lc_grid)  # Mean daily surface temperature

    for i in id:
        cap_grid[cap_grid == i] = solweig_parameters["Heat capacity"]["Value"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ]
        diff_grid[diff_grid == i] = solweig_parameters["Thermal_diffusivity"][
            "Value"
        ][solweig_parameters["Names"]["Value"][str((int(i.item())))]]

        # Coefficients of the OHM per land cover
        mean_a1 = solweig_parameters["OHM_coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][0]
        phi_a1 = solweig_parameters["OHM_coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][1]
        a1_grid[a1_grid == i] = mean_a1 * (
            1
            + 0.33
            * torch.sin(2 * torch.pi / 365.25 * day + phi_a1)
            * torch.sign(torch.tensor(location["latitude"], device=device))
        )
        a2_grid[a2_grid == i] = solweig_parameters["OHM_coefficients"][
            "Values"
        ][solweig_parameters["Names"]["Value"][str((int(i.item())))]][2]
        a3_grid[a3_grid == i] = solweig_parameters["OHM_coefficients"][
            "Values"
        ][solweig_parameters["Names"]["Value"][str((int(i.item())))]][3]

        # Initial ground surface temperature parameters
        offset_Tg = solweig_parameters["Tg_ini coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][0]

        slope_Tg = solweig_parameters["Tg_ini coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][0]

        ratio_Tg = float(
            solweig_parameters["Tg_ini coefficients"]["Values"][
                solweig_parameters["Names"]["Value"][str((int(i.item())))]
            ][1]
        )
        phi_Tg = 1.6

        # Correct the offset value given the latitude
        offset_Tg += slope_Tg * location["latitude"]

        # Mean daily soil temperature parameters
        ampl_Tm = solweig_parameters["Tm_ini coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][0]
        slope_Tm = solweig_parameters["Tm_ini coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][1]
        phi_Tm = 1.7
        offset_Tm = solweig_parameters["Tm_ini coefficients"]["Values"][
            solweig_parameters["Names"]["Value"][str((int(i.item())))]
        ][2]

        # Correct the offset value given the latitude
        offset_Tm += slope_Tm * location["latitude"]

        if i == 0 or i == 1:
            # For paved and asphalt landcover
            Tg[Tg == i] = (
                Ta[0]
                + offset_Tg
                * (
                    1
                    + ratio_Tg
                    * torch.sin(2 * torch.pi / 365.25 * day + phi_Tg)
                    * torch.sign(
                        torch.tensor(location["latitude"], device=device)
                    )
                )
                + 4
            )
            Tm[Tm == i] = (
                torch.mean(Ta)
                + ampl_Tm
                * torch.sin(2 * torch.pi / 365.25 * day + phi_Tm)
                * torch.sign(torch.tensor(location["latitude"], device=device))
                + offset_Tm
                + 4
            )

        elif i == 2:
            # For roofs
            Tg[Tg == i] = (
                Ta[0]
                + offset_Tg
                * (
                    1
                    + ratio_Tg
                    * torch.sin(2 * torch.pi / 365.25 * day + phi_Tg)
                    * torch.sign(
                        torch.tensor(location["latitude"], device=device)
                    )
                )
                + 4
            )
            Tm[Tm == i] = torch.mean(Ta) + offset_Tm

        elif i == 5:
            # For grass surfaces
            Tg[Tg == i] = Ta[0] + offset_Tg * (
                1
                + ratio_Tg
                * torch.sin(2 * torch.pi / 365.25 * day + phi_Tg)
                * torch.sign(torch.tensor(location["latitude"], device=device))
            )
            Tm[Tm == i] = (
                torch.mean(Ta)
                + ampl_Tm
                * torch.sin(2 * torch.pi / 365.25 * day + phi_Tm)
                * torch.sign(torch.tensor(location["latitude"], device=device))
                + offset_Tm
            )

        elif i == 6:
            # For bare soil landcover
            Tg[Tg == i] = (
                Ta[0]
                + offset_Tg
                * (
                    1
                    + ratio_Tg
                    * torch.sin(2 * torch.pi / 365.25 * day + phi_Tg)
                    * torch.sign(
                        torch.tensor(location["latitude"], device=device)
                    )
                )
                + 2
            )
            Tm[Tm == i] = (
                torch.mean(Ta)
                + ampl_Tm
                * torch.sin(2 * torch.pi / 365.25 * day + phi_Tm)
                * torch.sign(torch.tensor(location["latitude"], device=device))
                + offset_Tm
                + 2
            )

        elif i == 7:
            # For water bodies
            Tg[Tg == i] = Ta[0]

    if device.type == "cuda":
        torch.cuda.empty_cache()
    elif device.type == "xpu":
        torch.xpu.empty_cache()
    return (
        Tg,
        Tm,
        Rn,
        Rn_past,
        G,
        cap_grid,
        diff_grid,
        a1_grid,
        a2_grid,
        a3_grid,
    )


def surfaceTemperature_calc(
    Kdown,
    Ldown,
    Rn,
    Rn_past,
    G,
    Tg,
    Tm,
    alb,
    emis,
    cap,
    diff,
    lc_grid,
    a1,
    a2,
    a3,
    timestep,
    RH,
    shadow,
    shadow_past,
):
    """
    Calculation of the ground surface temperature

    Based on the force restore method with the ground heat flux modeled according to the heat storage flux formulation depicted in the OHM (Grimmond 1991).
    A simple model is implemented to assess the surface temperature for water landcover (lc==7) since the OHM fails to model the corresponding ground heat
    flux. The temporal integration is done using the Runge-Kutta 2nd order scheme

    :Parameters:
    :Kdown: the downwelling shortwave radiation
    :Ldown: the downwelling longwave radiation
    :Rn_past: net radiation stored to calculate the radiation rate
    :G: ground heat flux stored from the previous timestep
    :Tg: ground surface temperature grid
    :Tm: Temperature of the deep soil
    :alb: albedo grid of the ground surface
    :emis: emissivity grid of the surface
    :cap: heat capacity grid of the material
    :diff: thermal diffusivity grid of the material
    :lc_grid: landcover grid to identify the water bodies
    :a: coefficient grids related to the OHM
    :timestep: the timestep between 2 iteration of the simulation in min

    :Return:
    :Tg: Surface temperature
    :Rn: Net radiation for the current timestep
    :Rn_past: Stored net radiation for upcoming radiation rate calc
    :G: Ground heat flux for the current timestep
    """

    # Store the past ground surface temperature
    Tg_stored = Tg

    # Damping depths of the daily surface temperature wave
    D = torch.sqrt((2 * diff) / (2 * math.pi / 86400))

    ### Runge Kutta method for the surface temperature calc
    # First estimate of the surface temperature and of the deep soil temperature given past ground heat flux
    k1 = 2 * G / cap / D - 2 * math.pi / 86400 * (Tg - Tm)
    Tg_temp = Tg + k1 * timestep

    ### Estimate k2 the surface temperature step based on updated heat fluxes
    # The fluxes involved are calculated using the estimated Ts
    Lup_temp = SBC * emis * (Tg_temp + 273.15) ** 4 + Ldown * (
        1 - emis
    )  # Temporary outgoing longwave rad (W.m-2)
    Rn_temp = Kdown * (1 - alb) + Ldown - Lup_temp  # Temporary net rad (W.m-2)
    RnStar_temp = (Rn_temp - Rn) / 1  # Temporary radiation rate (W.m-2.h-1)
    G_temp = (
        a1 * Rn_temp + a2 * RnStar_temp + a3
    )  # Temporary ground heat flux (W.m-2)

    # Damping of the ground heat flux if it increases (or drops) too quickly
    deltaG = abs(G_temp - G)
    radCriterion = abs(
        a1 * (Rn_temp - Rn_past)
    )  # Criterion regarding the radiation step
    mask = torch.logical_and(
        deltaG > radCriterion, abs(shadow - shadow_past) > 0.5
    )  # Grid of the pixels where the ground heat flux spikes
    G_temp[mask] = G[mask] + torch.sign(G_temp - G)[mask] * radCriterion[mask]

    # Correction of the temperature estimates
    k2 = 2 * G_temp / cap / D - 2 * math.pi / 86400 * (Tg_temp - Tm)
    Tg += (k1 + k2) / 2 * timestep

    ### Finally calculation of the updated heat fluxes
    Rn_past = Rn
    G_past = G
    Lup = SBC * emis * (Tg + 273.15) ** 4 + (1 - emis) * Ldown
    Rn = (1 - alb) * Kdown + Ldown - Lup
    Rn_star = (Rn - Rn_past) / 1
    G = a1 * Rn + a2 * Rn_star + a3

    # Damping of the ground heat flux if it increases (or decreases) too quickly
    deltaG = abs(G - G_past)
    radCriterion = abs(
        a1 * (Rn - Rn_past)
    )  # Criterion regarding the radiation step
    mask = torch.logical_and(
        deltaG > radCriterion, abs(shadow - shadow_past) > 0.5
    )  # Grid of the pixels where the ground heat flux spikes
    G[mask] = G_past[mask] + torch.sign(G - G_past)[mask] * radCriterion[mask]

    ### Water bodies surface temperature estimate
    beta = 0.45  # Amount of shortwave rad absorbed by the water surface layer
    thickness = 1  # Depth of the water layer
    lamb = 2.260e6  # Latent heat of vaporisation
    rho = 1000  # Density of water (kg.m-3)
    Rn_water = (
        Kdown
        * (1 - alb)
        * (
            beta
            + (1 - beta) * (1 - torch.exp(torch.tensor(-1, device=Tg.device)))
        )
        + Ldown
        - Lup
    )  # Net radiation for the top water layer beta described the transmitted rad
    _, es = saturated_vp(Tg)
    E = 0.0858 * (es / 1000) * (1 - RH / 100) / 3600 / 1000 * rho * lamb
    deltaTg = torch.clone(lc_grid)
    deltaTg = (
        timestep
        / cap
        / thickness
        * (Rn_water - E - diff * cap / thickness * (Tg - Tm))
    )
    Tg[lc_grid == 7] = Tg_stored[lc_grid == 7] + deltaTg[lc_grid == 7]

    return Tg, Rn, Rn_past, G


def outgoingLongwave_calc(
    Tg,
    Tgwall,
    Ta,
    Ldown,
    emis,
    alb,
    buildings,
    shadow,
    sunwall,
    walls,
    rows,
    cols,
    sizepx,
):
    """
    Calculation of the outgoing longwave radiation from the ground,

    :Parameters:
    :Tg: ground surface temperature grid
    :Tgwall: wall surface temperature grid
    :Ta: air temperature grid
    :emis: emissivity grid of the surface
    :alb: albedo grid of the surface
    :emis_wall: emissivity of the wall (for now float = 0.9)
    :buildings: boolean grid, 0 if the landcover is roof and 1 if there is no building
    :shadow: boolean grid, 0 when the pixel is shadowed and 1 when sunlit
    :sunwall: Grid where non zero values indicate a sunlit wall and its height
    :walls: grid containing the heights of the walls
    :aspect: grid containing the angles of the normal dir to walls
    :rows: number of rows in the grids
    :cols: number of columns in the grids
    :sizepx: size of a pixel in m (1/scale)

    :Return:
    """

    # Assessment of the distance from a pixel at which most of the radiation are received (cf view factor Lambert)
    device = Tg.device
    factor = torch.tensor(
        0.99, device=device
    )  # Percentage of radiation accounted for
    zs = 1.1  # in m
    r_max = zs * torch.sqrt(
        factor / (1 - factor)
    )  # in m, maximum radius for the radiation calc

    # Emissivity of the wall
    emis_wall = 0.9
    alb_wall = 0.2

    # Copy of the sunlit wall grid and replacement of the wall height with 1 if sunlit
    sunlitwall = sunwall
    sunlitwall[sunlitwall > 0] = 1

    # Boolean array 1 if the pixel is a wall, 0 if not
    wallbol = walls > 0

    # The alb grids only take into account the sunlit surfaces in the alb calculation albnosh calculate it for all the surfaces
    albsunlit = alb * shadow

    # Boolean array 1 if the pixel is a wall, 0 if not
    wallbol = (walls > 0) * 1

    # step in meters between every iteration
    step = 1

    # Grid of the outgoing longwave radiation coming from the ground
    Lup = (SBC * emis * (Tg + 273.15) ** 4 + Ldown * (1 - emis)) * buildings

    ### Initialize the ground view factor grids as torch.zeros()
    # Upwelling longwave radiation
    gvfLup = torch.zeros((rows, cols), device=device)
    gvfLupE = torch.zeros((rows, cols), device=device)
    gvfLupS = torch.zeros((rows, cols), device=device)
    gvfLupW = torch.zeros((rows, cols), device=device)
    gvfLupN = torch.zeros((rows, cols), device=device)

    # Albedo of the sunlit surfaces
    gvfalbsun = torch.zeros((rows, cols), device=device)
    gvfalbsunE = torch.zeros((rows, cols), device=device)
    gvfalbsunS = torch.zeros((rows, cols), device=device)
    gvfalbsunW = torch.zeros((rows, cols), device=device)
    gvfalbsunN = torch.zeros((rows, cols), device=device)

    # Albedo complete
    gvfalbtot = torch.zeros((rows, cols), device=device)
    gvfalbtotE = torch.zeros((rows, cols), device=device)
    gvfalbtotS = torch.zeros((rows, cols), device=device)
    gvfalbtotW = torch.zeros((rows, cols), device=device)
    gvfalbtotN = torch.zeros((rows, cols), device=device)

    # Longwave radiation coming from the side
    gvfLsideE = torch.zeros((rows, cols), device=device)
    gvfLsideS = torch.zeros((rows, cols), device=device)
    gvfLsideW = torch.zeros((rows, cols), device=device)
    gvfLsideN = torch.zeros((rows, cols), device=device)

    # Add the radiation from the pixel directly below, only for the total gvf
    # Do not take the roofs into account for now
    view_factor = (sizepx / 2) ** 2 / ((sizepx / 2) ** 2 + zs**2)
    gvfLup = gvfLup + Lup * view_factor
    gvfalbsun = gvfalbsun + albsunlit * view_factor * buildings
    gvfalbtot = gvfalbtot + alb * view_factor * buildings

    # Division of the 360° field of view in 20 and convert the array in radian
    azimuths = torch.linspace(18, 360, steps=20, device=device)
    azimuths = azimuths * (torch.pi / 180)

    ### Loop for the number of azimuth values
    for azimuth in azimuths:
        # Copy of the building grid
        building_copy = buildings.clone()

        # Boolean array 1 if the pixel is (or was) a wall, 0 if not
        pastwalls = wallbol.clone()

        # Grid of the longwave radiation emitted by the walls
        Lwall = SBC * emis_wall * (Tgwall + Ta + 273.15) ** 4 * wallbol

        # Initialisation of the tables
        # First the ones containing the translated rasters (temporary)
        building_temp = buildings.clone()
        Lup_temp = Lup.clone()
        Lwall_temp = Lwall.clone()
        albsun_temp = albsunlit.clone()
        albtot_temp = alb.clone()
        walls_temp = torch.zeros((rows, cols), device=device)
        sunlitwall_temp = torch.zeros((rows, cols), device=device)
        onlywall_temp = torch.zeros((rows, cols), device=device)

        # Then the tables containing the sum of the radiations (or albedo) for this azimuth
        Lup_sum = torch.zeros((rows, cols), device=device)
        LsideE_sum = torch.zeros((rows, cols), device=device)
        LsideN_sum = torch.zeros((rows, cols), device=device)
        LsideW_sum = torch.zeros((rows, cols), device=device)
        LsideS_sum = torch.zeros((rows, cols), device=device)
        albsun_sum = torch.zeros((rows, cols), device=device)
        albtot_sum = torch.zeros((rows, cols), device=device)

        ### Shadow casting algorithm
        # Translation ranges from 1/2 a pixel to the max radius r_max
        for r in torch.arange(sizepx / 2, r_max, step=step):

            # Step of the raster translation
            dx = -torch.cos(azimuth)
            dy = -torch.sin(azimuth)

            # Scale so that the grid is at least translated from 1px
            if abs(dx) > abs(dy):
                dx = -r * torch.sign(torch.cos(azimuth))
                dy = (
                    -r
                    * abs(torch.tan(azimuth))
                    * torch.sign(torch.sin(azimuth))
                )
            else:
                dx = (
                    -r
                    / abs(torch.tan(azimuth))
                    * torch.sign(torch.cos(azimuth))
                )
                dy = -r * torch.sign(torch.sin(azimuth))

            # Select the interested part of the initial raster and the translated one from their four corners and
            # translating toward the direction azimuth = 0° for dx > 0
            if dx > 0:
                x_select_start = dx
                x_select_end = rows
                x_transl_start = 0
                x_transl_end = rows - dx
            else:
                x_select_start = 0
                x_select_end = rows + dx
                x_transl_start = -dx
                x_transl_end = rows

            # translating toward the direction azimuth = 90° for dy > 0
            if dy > 0:
                y_select_start = dy
                y_select_end = cols
                y_transl_start = 0
                y_transl_end = cols - dy
            else:
                y_select_start = 0
                y_select_end = cols + dy
                y_transl_start = -dy
                y_transl_end = cols

            # Copy the initial rasters and input inside translated raster temporary and changing every iteration
            # Building grid
            building_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = buildings[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Ground longwave radiation grid
            Lup_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = Lup[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Wall longwave radiation grid
            Lwall_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = Lwall[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Albedo grid for the sunlit area
            albsun_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = albsunlit[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Albedo grid for all the area
            albtot_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = alb[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Sunlit wall grid
            sunlitwall_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = sunlitwall[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # All walls grid
            walls_temp[
                int(x_transl_start) : math.ceil(x_transl_end),
                int(y_transl_start) : math.ceil(y_transl_end),
            ] = wallbol[
                int(x_select_start) : math.ceil(x_select_end),
                int(y_select_start) : math.ceil(y_select_end),
            ]

            # Change the boolean building grid, if the px was already a building it remains one (px value = 0)
            building_copy = torch.min(building_copy, building_temp)

            # For each pixel add the translated Lup to the received rad if there where there is no building
            view_factor = ((r + step) ** 2 / (zs**2 + (r + step) ** 2)) - (
                r**2 / (zs**2 + r**2)
            )
            Lup_sum += Lup_temp * view_factor * building_copy / 20
            albsun_sum += albsun_temp * view_factor * building_copy / 20
            albtot_sum += albtot_temp * view_factor * building_copy / 20

            # Create a boolean grid to assert that the sunlit walls are not inside a building
            onlywall_temp = torch.logical_and(
                walls_temp, torch.logical_not(pastwalls)
            )
            onlywall_temp = onlywall_temp * building_copy
            onlysunwall_temp = sunlitwall_temp * building_copy

            pastwalls = torch.logical_or(pastwalls, walls_temp)

            # Compute the view factor of the wall surface for a given translatio distance
            viewfactor_wall = (
                1
                / 2 ** (1 / 2)
                / 3
                * torch.sqrt(
                    1 + (r + step) / torch.sqrt((r + step) ** 2 + zs**2)
                )
                * (
                    2
                    + (r + sizepx / 2)
                    / torch.sqrt((r + sizepx / 2) ** 2 + zs**2)
                )
                * (
                    1
                    - (r + sizepx / 2)
                    / torch.sqrt((r + sizepx / 2) ** 2 + zs**2)
                )
                / zs
                * torch.sqrt((r + sizepx / 2) ** 2 + zs**2)
            )

            # Then add the radiation incoming from those walls
            Lup_sum += (
                onlywall_temp
                * Lwall_temp
                * viewfactor_wall
                * building_copy
                / 20
            )
            albsun_sum += (
                onlysunwall_temp
                * alb_wall
                * viewfactor_wall
                * building_copy
                / 20
            )
            albtot_sum += (
                onlywall_temp * alb_wall * viewfactor_wall * building_copy / 20
            )

            # Finally add the radiation in Lside
            dphi = torch.arctan((r + step) / zs) - torch.arctan(r / zs)
            dtrigo = zs / torch.sqrt(r**2 + zs**2) * r / torch.sqrt(
                r**2 + zs**2
            ) - zs / torch.sqrt((r + step) ** 2 + zs**2) * (
                r + step
            ) / torch.sqrt(
                (r + step) ** 2 + zs**2
            )

            # Calculation of the solid angle for each of the cardinal points
            # plus add the radiation from a potential wall
            steradiansW, steradiansS, steradiansE, steradiansN = 0, 0, 0, 0
            if (azimuth >= 0) and (azimuth < torch.pi):
                dthetaW = 2 * torch.pi / 20
                steradiansW += dthetaW * (dphi + dtrigo) / 2
                LsideW_sum += onlywall_temp * Lwall_temp * viewfactor_wall / 10

            if (azimuth >= torch.pi / 2) and (azimuth < 3 * torch.pi / 2):
                dthetaS = 2 * torch.pi / 20
                steradiansS += dthetaS * (dphi + dtrigo) / 2
                LsideS_sum += onlywall_temp * Lwall_temp * viewfactor_wall / 10

            if (azimuth >= torch.pi) and (azimuth < 2 * torch.pi):
                dthetaE = 2 * torch.pi / 20
                steradiansE += dthetaE * (dphi + dtrigo) / 2
                LsideE_sum += onlywall_temp * Lwall_temp * viewfactor_wall / 10

            if (azimuth >= 3 * torch.pi / 2) or (azimuth < torch.pi / 2):
                dthetaN = 2 * torch.pi / 20
                steradiansN += dthetaN * (dphi + dtrigo) / 2
                LsideN_sum += onlywall_temp * Lwall_temp * viewfactor_wall / 10

            LsideW_sum += Lup_temp / torch.pi * steradiansW * building_copy
            LsideS_sum += Lup_temp / torch.pi * steradiansS * building_copy
            LsideE_sum += Lup_temp / torch.pi * steradiansE * building_copy
            LsideN_sum += Lup_temp / torch.pi * steradiansN * building_copy

        ### Add the value for the computed part of the field of view
        gvfLup += Lup_sum
        gvfalbsun += albsun_sum
        gvfalbtot += albtot_sum

        # Add the value if the azimuth correspond to the side of the compass
        if (azimuth >= 0) and (azimuth < torch.pi):
            gvfLupW += Lup_sum
            gvfalbsunW += albsun_sum
            gvfalbtotW += albtot_sum
            gvfLsideW += LsideW_sum

        if (azimuth >= torch.pi / 2) and (azimuth < 3 * torch.pi / 2):
            gvfLupS += Lup_sum
            gvfalbsunS += albsun_sum
            gvfalbtotS += albtot_sum
            gvfLsideS += LsideS_sum

        if (azimuth >= torch.pi) and (azimuth < 2 * torch.pi):
            gvfLupE += Lup_sum
            gvfalbsunE += albsun_sum
            gvfalbtotE += albtot_sum
            gvfLsideE += LsideE_sum

        if (azimuth >= 3 * torch.pi / 2) or (azimuth < torch.pi / 2):
            gvfLupN += Lup_sum
            gvfalbsunN += albsun_sum
            gvfalbtotN += albtot_sum
            gvfLsideN += LsideN_sum

    # If the px is associated with a roof landcover, for now Lup = 0
    # Here their Lup value is allocated to those px
    gvfLup += (SBC * emis * (Tg + 273.15) ** 4) * (buildings * -1 + 1)
    gvfLsideE += (SBC * emis * (Tg + 273.15) ** 4) * 0.5 * (buildings * -1 + 1)
    gvfLsideN += (SBC * emis * (Tg + 273.15) ** 4) * 0.5 * (buildings * -1 + 1)
    gvfLsideW += (SBC * emis * (Tg + 273.15) ** 4) * 0.5 * (buildings * -1 + 1)
    gvfLsideS += (SBC * emis * (Tg + 273.15) ** 4) * 0.5 * (buildings * -1 + 1)
    gvfalbsun += albsunlit * (buildings * -1 + 1)
    gvfalbtot += alb * (buildings * -1 + 1)
    gvfalbsunE += albsunlit * 0.5 * (buildings * -1 + 1)
    gvfalbsunN += albsunlit * 0.5 * (buildings * -1 + 1)
    gvfalbsunW += albsunlit * 0.5 * (buildings * -1 + 1)
    gvfalbsunS += albsunlit * 0.5 * (buildings * -1 + 1)
    gvfalbtotE += alb * 0.5 * (buildings * -1 + 1)
    gvfalbtotN += alb * 0.5 * (buildings * -1 + 1)
    gvfalbtotW += alb * 0.5 * (buildings * -1 + 1)
    gvfalbtotS += alb * 0.5 * (buildings * -1 + 1)

    return (
        gvfLup,
        gvfalbsun,
        gvfalbtot,
        gvfLupE,
        gvfalbsunE,
        gvfalbtotE,
        gvfLupS,
        gvfalbsunS,
        gvfalbtotS,
        gvfLupW,
        gvfalbsunW,
        gvfalbtotW,
        gvfLupN,
        gvfalbsunN,
        gvfalbtotN,
        gvfLsideW,
        gvfLsideS,
        gvfLsideE,
        gvfLsideN,
    )
