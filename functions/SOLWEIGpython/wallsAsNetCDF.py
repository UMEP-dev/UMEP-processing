import xarray as xr
import numpy as np

import gc

def walls_as_netcdf(voxelTable, rows, cols, timeSlots, iteration, raster_path, output_path):
    '''This function creates a 4D NetCDF with wall temperatures and corresponding emitted longwave radiation'''
    # rows = number of rows (latitudinal position)
    # cols = number of columns (longitudinal position)
    # level = number of voxel levels (elevation position)
    # raster_path is used to load an existing .tif layer and create arrays with latitudes and longitudes to be used in xarray/NetCDF

    # Highest number of voxels used to determine z/height level of NetCDF
    levels = voxelTable.loc[voxelTable['voxelHeight'] == voxelTable['voxelHeight'].max(), 'voxelHeight'].to_numpy()[0].astype(int)

    # Range of height levels
    height_levels = np.arange(1, levels+1)

    # Create two empty numpy arrays to fill with wall temperatures and corresponding longwave radiation from current time step
    wallTemperature = np.ones((cols, rows, levels)) * np.nan
    longwaveRadiation = np.ones((cols, rows, levels)) * np.nan

    # Add current time step wall temperature and longwave radiation to numpy array, which will be used to update the NetCDF.
    for y, x, z, wallTemp, wallRad in zip(voxelTable['ypos'].astype(int), voxelTable['xpos'].astype(int), voxelTable['voxelHeight'].astype(int), voxelTable['wallTemperature'], voxelTable['LongwaveRadiation']):
        wallTemperature[x, y, z-1] = wallTemp
        longwaveRadiation[x, y, z-1] = wallRad

    # NetCDF compression
    # comp = dict(zlib=True, complevel=5)
    comp = dict(zlib=True)

    # If first time step, create an empty NetCDF to fill with values
    if iteration == 0:
        raster_file = xr.open_rasterio(raster_path)
        lat = raster_file.y.to_numpy()
        lon = raster_file.x.to_numpy()
        temp_data = np.zeros((cols, rows, levels, timeSlots.shape[0]))
        data_xr = xr.Dataset(
            data_vars=dict(
            wall_temperature=(["lon", "lat", "height", "time"], temp_data),
            #longwave_radiation=(["y", "x", "height", "time"], temp_data),
                ),
            coords=dict(
            lon=lon,
            lat=lat,
            height=height_levels,
            time=timeSlots,
            ),
            attrs=raster_file.attrs
        )
        
        # Update wall temperature and longwave radiation for current timestep (iteration)
        data_xr.wall_temperature[:, :, :, iteration] = wallTemperature
        # data_xr.longwave_radiation[:, :, :, iteration] = longwaveRadiation
        
        encodings = {var: comp for var in data_xr.data_vars}

        # Save as NetCDF
        data_xr.to_netcdf(output_path, encoding=encodings)

        data_xr.close()
    # If not first time step, load existing NetCDF as an xarray dataset
    else:
        data_xr = xr.load_dataset(output_path)
        
        # Update wall temperature and longwave radiation for current timestep (iteration)
        data_xr.wall_temperature[:, :, :, iteration] = wallTemperature
        # data_xr.longwave_radiation[:, :, :, iteration] = longwaveRadiation
        
        encodings = {var: comp for var in data_xr.data_vars}
        
        # Save as NetCDF
        data_xr.to_netcdf(output_path, encoding=encodings)

# with xr.open_dataset(output_path, engine="scipy") as data_xr:
# data_xr.to_netcdf(output_path, mode="w", engine="scipy")