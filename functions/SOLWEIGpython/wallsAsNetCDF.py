import xarray as xr
import rioxarray
import numpy as np

def walls_as_netcdf(voxelTable, rows, cols, timeSlots, iteration, dsm, raster_path, output_path):
    '''This function creates a 4D NetCDF with wall temperatures and corresponding emitted longwave radiation'''
    # rows = number of rows (latitudinal position)
    # cols = number of columns (longitudinal position)
    # level = number of voxel levels (elevation position)
    # raster_path is used to load an existing .tif layer and create arrays with latitudes and longitudes to be used in xarray/NetCDF

    # Highest number of voxels used to determine z/height level of NetCDF
    # levels = voxelTable.loc[voxelTable['voxelHeight'] == voxelTable['voxelHeight'].max(), 'voxelHeight'].to_numpy()[0].astype(int)

    levels = voxelTable.loc[voxelTable['voxelHeightMasl'] == voxelTable['voxelHeightMasl'].max(), 'voxelHeightMasl'].to_numpy()[0].astype(int)

    # Range of height levels
    height_levels = np.arange(1, levels+1)

    # Create empty numpy array to fill with wall temperatures from current time step
    wallTemperature = np.full((cols, rows, levels), np.nan, dtype=np.float32)

    # Add current time step wall temperature and longwave radiation to numpy array, which will be used to update the NetCDF.
    #for y, x, z, wallTemp in zip(voxelTable['ypos'].astype(int), voxelTable['xpos'].astype(int), voxelTable['voxelHeight'].astype(int), voxelTable['wallTemperature'].astype(np.float32)):
    for y, x, z, wallTemp in zip(voxelTable['ypos'].astype(int), voxelTable['xpos'].astype(int), voxelTable['voxelHeightMasl'].astype(int), voxelTable['wallTemperature'].astype(np.float32)):
        wallTemperature[x, y, z-1] = wallTemp

    # NetCDF compression
    comp = dict(zlib=True, 
                complevel=5,
                chunksizes=(100, 100, 10, 1))

    # If first time step, create an empty NetCDF to fill with values
    if iteration == 0:
        # raster_file = xr.open_rasterio(raster_path)
        # raster_file = xr.open_dataset(raster_path, decode_coords='all', engine='rasterio')
        raster_file = rioxarray.open_rasterio(raster_path)
        lat = raster_file.y.to_numpy()
        lon = raster_file.x.to_numpy()
        # temp_data = np.zeros((cols, rows, levels, timeSlots.shape[0]))
        temp_data = np.full((cols, rows, levels, timeSlots.shape[0]), np.nan, dtype=np.float32)
        data_xr = xr.Dataset(
            data_vars=dict(
            wall_temperature=(["lon", "lat", "height", "time"], temp_data),
                ),
            coords=dict(
            lon=lon,
            lat=lat,
            height=height_levels,
            time=timeSlots
            ),
            attrs={"crs": raster_file.rio.crs.to_string()}
        )
        
        # Update wall temperature and longwave radiation for current timestep (iteration)
        data_xr.wall_temperature[:, :, :, iteration] = wallTemperature
        
        encodings = {var: comp for var in data_xr.data_vars}

        # Save as NetCDF
        data_xr.to_netcdf(output_path, engine="netcdf4", encoding=encodings)

        data_xr.close()
    # If not first time step, load existing NetCDF as an xarray dataset
    else:
        with xr.open_dataset(output_path, mode='r+') as data_xr:
            # Update wall temperature for current timestep (iteration)
            data_xr.wall_temperature[:, :, :, iteration] = wallTemperature
            # Save changes
            data_xr.to_netcdf(output_path, engine="netcdf4", mode="a")    