from osgeo import gdal, osr
import numpy as np
from scipy.ndimage import label
import os

class Inputdata():
    '''Class containing input data for Tree planter'''
    __slots__ = ('dataSet', 'buildings', 'selected_area', 'dsm', 'cdsm', 'cdsm_b', 'shadow', 'tmrt_ts', 'tmrt_s', 'tmrt_avg', 'rows', 'cols', 'scale', 'lat', 'lon', 'gt')
    def __init__(self,r_range, sh_fl, tmrt_fl, infolder, inputPolygonlayer, feedback):

        self.dataSet = gdal.Open(infolder + '/buildings.tif')            # GIS data
        self.buildings = self.dataSet.ReadAsArray().astype(np.float)    # Building raster
        self.buildings = self.buildings == 1.0
        self.rows = self.buildings.shape[0]                             # Rows of input rasters from SOLWEIG
        self.cols = self.buildings.shape[1]                             # Cols of input rasters from SOLWEIG
        self.cdsm = np.zeros((self.rows,self.cols))                               # Canopy digital surface model
        self.cdsm_b = np.zeros((self.rows,self.cols))  # Canopy digital surface model
        self.shadow = np.zeros((self.rows,self.cols, r_range.__len__()))          # Shadow rasters
        self.tmrt_ts = np.zeros((self.rows,self.cols, r_range.__len__()))        # Tmrt for each timestep
        self.tmrt_s = np.zeros((self.rows,self.cols))                            # Sum of tmrt for all timesteps

        # Loading DEm, DSM (and CDSM) rasters
        dataSet = gdal.Open(infolder + '/DSM.tif')
        self.dsm = dataSet.ReadAsArray().astype(np.float)
        # dataSet = gdal.Open(infolder + '/DEM.tif')
        # self.dem = dataSet.ReadAsArray().astype(np.float)

        # Check if CDSM exists
        if os.path.exists(infolder + '/CDSM.tif'):
            dataSet = gdal.Open(infolder + '/CDSM.tif')
            self.cdsm = dataSet.ReadAsArray().astype(np.float)
            self.cdsm_b = self.cdsm == 0
            self.buildings = (self.buildings == True) & (self.cdsm_b == True)

        c = 0
        for iy in r_range:
            dataSet1 = gdal.Open(sh_fl[iy])
            feedback.setProgressText('Loading ' + sh_fl[iy] + '..')
            self.shadow[:, :, c] = dataSet1.ReadAsArray().astype(np.float)
            dataSet2 = gdal.Open(tmrt_fl[iy])
            feedback.setProgressText('Loading ' + tmrt_fl[iy] + '..')
            self.tmrt_ts[:, :, c] = np.around(dataSet2.ReadAsArray().astype(np.float), decimals=1) * self.shadow[:, :, c]
            self.tmrt_s = self.tmrt_s + self.tmrt_ts[:, :, c]

            c += 1

        self.tmrt_avg = (self.tmrt_s / c)

        # Find latlon etc for input data.
        old_cs = osr.SpatialReference()
        # dsm_ref = dsmlayer.crs().toWkt()
        dsm_ref = self.dataSet.GetProjection()
        old_cs.ImportFromWkt(dsm_ref)

        wgs84_wkt = """
            GEOGCS["WGS 84",
                DATUM["WGS_1984",
                    SPHEROID["WGS 84",6378137,298.257223563,
                        AUTHORITY["EPSG","7030"]],
                    AUTHORITY["EPSG","6326"]],
                PRIMEM["Greenwich",0,
                    AUTHORITY["EPSG","8901"]],
                UNIT["degree",0.01745329251994328,
                    AUTHORITY["EPSG","9122"]],
                AUTHORITY["EPSG","4326"]]"""

        new_cs = osr.SpatialReference()
        new_cs.ImportFromWkt(wgs84_wkt)

        transform = osr.CoordinateTransformation(old_cs, new_cs)

        width1 = self.dataSet.RasterXSize
        height1 = self.dataSet.RasterYSize
        self.gt = self.dataSet.GetGeoTransform()
        minx = self.gt[0]
        maxy = self.gt[3]
        miny = self.gt[3] + width1 * self.gt[4] + height1 * self.gt[5]
        maxx = minx + self.gt[1] * width1
        
        lonlat = transform.TransformPoint(minx, miny)
        geotransform = self.dataSet.GetGeoTransform()
        self.scale = 1 / self.gt[1]
        self.lon = lonlat[0]
        self.lat = lonlat[1]

        # Import Planting area
        rasterize_options = gdal.RasterizeOptions(options=[
            '-burn', '1',
            '-te', str(minx), str(miny), str(maxx), str(maxy),
            '-tr', str(self.gt[1]), str(self.gt[5])
        ])

        gdal.Rasterize(infolder + '/selected_area.tif', inputPolygonlayer, options=rasterize_options)
        dataSetSel = gdal.Open(infolder + '/selected_area.tif')
        self.selected_area = dataSetSel.ReadAsArray().astype(np.float)

        try:
            del dataSetSel
            os.remove(infolder + '/selected_area.tif')
            print('Successfully removed selected_area.tif')
        except:
            print('Could not remove selected_area.tif')

        # Buffer zone to remove potential edge effects
        buffer_percentage = 0.05
        buffer_y = np.int_(self.rows * buffer_percentage)
        buffer_x = np.int_(self.cols * buffer_percentage)
        buffer_zone = np.zeros((self.rows, self.cols))
        buffer_zone[buffer_y:-buffer_y, buffer_x:-buffer_x] = 1

        self.selected_area = self.selected_area * buffer_zone

class Treerasters():
    '''Class containing calculated shadows, regional grouping of shadows \
    if many timesteps, tmrt in shade, tmrt sunlit, difference between \
    shade and sunlit'''

    __slots__ = ('treeshade', 'treeshade_rg', 'treeshade_bool', 'cdsm', 'buffer_y', 'buffer_x', 'tpy', 'tpx', 'rows', 'cols', 'rows_s', 'cols_s', 'euclidean', 'euclidean_d', 'tmrt_sun', 'tmrt_shade', 'd_tmrt')
    def __init__(self, treeshade, treeshade_rg, treeshade_bool, cdsm, treedata):
        # Find min and max rows and cols where there are shadows
        shy, shx = np.where((treeshade > 0) | (cdsm > 0))
        shy_min = np.min(shy); shy_max = np.max(shy) + 1
        shx_min = np.min(shx); shx_max = np.max(shx) + 1

        y = np.int_(np.abs(treedata.treey - shy_min))
        x = np.int_(np.abs(treedata.treex - shx_min))

        # Cropping to only where there is a shadow
        self.treeshade = treeshade[shy_min:shy_max, shx_min:shx_max]
        self.treeshade_rg = treeshade_rg[shy_min:shy_max, shx_min:shx_max]
        self.treeshade_bool = 1-treeshade_bool[shy_min:shy_max, shx_min:shx_max,:]
        self.cdsm = cdsm[shy_min:shy_max, shx_min:shx_max]
        # y, x = np.where(cdsm_clip == treedata.height)  # Position of tree in clipped shadow image
        self.buffer_y = np.zeros((2))
        self.buffer_x = np.zeros((2))
        self.buffer_y[0] = np.int_(y)
        self.buffer_y[1] = np.int_(self.cdsm.shape[0] - y)
        self.buffer_x[0] = np.int_(x)
        self.buffer_x[1] = np.int(self.cdsm.shape[1] - x)
        self.tpy = y
        self.tpx = x
        self.rows = treeshade.shape[0]
        self.cols = treeshade.shape[1]
        self.rows_s = self.treeshade_rg.shape[0]
        self.cols_s = self.treeshade_rg.shape[1]

        a = np.array((self.tpy, self.tpx))
        b = np.zeros((4,2))
        b[0,:] = np.array((0, 0))   # Upper left corner
        b[1,:] = np.array((0, self.treeshade.shape[0] - 1)) # Lower left corner
        b[2,:] = np.array((self.treeshade.shape[1] - 1, 0)) # Upper right corner
        b[3,:] = np.array((self.treeshade.shape[0] - 1, self.treeshade.shape[1] - 1))   # Lower right corner
        eucl = np.zeros((b.shape[0],1))
        for i in range(b.shape[0]):
            eucl[i,0] = np.linalg.norm(a - b[i,:])
        eucl_d = np.array([eucl[0] + eucl[3], eucl[1] + eucl[2]])
        self.euclidean = np.max(eucl[:,0])
        self.euclidean_d = np.max(eucl_d)

        self.tmrt_sun = 0
        self.tmrt_shade = 0
        self.d_tmrt = 0

    def tmrt(self, tmrt_sun, tmrt_shade):
        '''Calculate difference in Tmrt between sun and shade'''
        nr_dec = 1      
        self.tmrt_sun = np.around(tmrt_sun, decimals=nr_dec)
        self.tmrt_shade = np.around(tmrt_shade, decimals=nr_dec)
        self.d_tmrt = (self.tmrt_sun - self.tmrt_shade)
        #print(self.tmrt_sun - self.tmrt_shade)
        # if filter == 1:
        #     import scipy as sp
        #     # Finding courtyards and small separate areas
        #     sun_vs_tsh_filtered = sp.ndimage.label(self.d_tmrt)
        #     sun_sh_d = sun_vs_tsh_filtered[0]
        #     for i in range(1, sun_sh_d.max() + 1):
        #         if np.sum(sun_sh_d[sun_sh_d == i]) < 2000:
        #             self.d_tmrt[sun_sh_d == i] = 0
        # self.tmrt_sun = np.around(self.tmrt_sun, decimals=nr_dec)
        # self.d_tmrt = np.around(self.d_tmrt, decimals=nr_dec)

class Position():
# Class containing y and x positions of trees and their corresponding sum of Tmrt in shade and sum of Tmrt in same area as shade but sunlit
# as well as a unique number for each position. Also a matrix with the unique number in each y,x position in the matrix.
    __slots__ = ('pos', 'pos_m')
    def __init__(self,vector,rows,cols):
        self.pos = vector

        self.pos_m = np.zeros((rows, cols))
        for idx in range(vector.shape[0]):
            y = np.int_(vector[idx, 2])
            x = np.int_(vector[idx, 1])
            self.pos_m[y, x] = vector[idx, 0]

class Treedata():
# Class containing data for the tree that is used in Tree planter, i.e. the tree that is being "planted" and studied
    __slots__ = ('ttype', 'height', 'trunk', 'dia', 'treey', 'treex')
    def __init__(self, ttype, height, trunk, dia, treey, treex):
        self.ttype = ttype
        self.height = height
        self.trunk = trunk
        self.dia = dia
        self.treey = treey
        self.treex = treex


class Regional_groups():
    # Class for creation of regional groups for shadows, i.e. which timesteps shade which pixels
    # Returns a matrix with regional groups and a vector with the corresponding timesteps for each regional group
    # range_ = between which timesteps to calculate regional groups
    # shadow_ = matrix with sum of shadows for all timesteps
    # shadow_ts = shadows for each timestep
    __slots__ = ('shadow', 'timesteps', 'tmrt_rg')
    def __init__(self, range_, shadow_, shadow_ts, tmrt):

        t_r = range(range_.__len__())
        t_l = t_r.__len__()

        shade_u = np.unique(shadow_)                            # Unique values in summation matrix for tree shadows
        shade_max = np.max(shade_u)                             # Maximum value of unique values

        for i in range(1, shade_u.shape[0]):                    # Loop over all unique values
            shade_b = (shadow_ == shade_u[i])                   # Boolean shadow for each timestep i
            shade_r = label(shade_b)                            # Create regional groups
            shade_r_u = np.unique(shade_r[0])                   # Find out how many regional groups, i.e. unique values
            if np.sum(shade_r_u) > 1:                           # If more than there groups, i.e. 0, 1, 2, ... , continue
                for j in range(2, shade_r_u.shape[0]):          # Loop over the unique values and give all but 1 new values
                    shade_b2 = (shade_r[0] == shade_r_u[j])     # Boolean of shadow for each unique value
                    shade_max += 1                              # Add +1 to the maximum value of unique values, continues (creates new unique values)
                    shadow_[shade_b2] = shade_max               # Add these to the building summation matrix

        shade_u_u = np.unique(shadow_)                          # New unique values of regional groups
        sh_vec_t = np.zeros((shade_u_u.shape[0], t_l + 1))      # Empty array for storing which timesteps are found in each regional group
        sh_vec_t[:, 0] = shade_u_u                              # Adding the unique regional groups to the first column
        tmrt_t = np.zeros((shade_u_u.shape[0], 2))
        tmrt_t[:,0] = shade_u_u

        for i in range(1, shade_u_u.shape[0]):                  # Loop over the unique values
            shade_b = (shadow_ == shade_u_u[i])                 # Boolean of each regional group
            for j in t_r:                                       # Loop over each timestep
                shade_b2 = (shadow_ts[:, :, j].copy() == 1)     # Boolean of shadow for each timestep
                shade_b3 = (shade_b) & (shade_b2)               # Find out where they overlap, i.e. which timesteps are found in each regional group
                if np.sum(shade_b3) > 0:                        # If they overlap, continue
                    sh_vec_t[i, 1 + j] = 1                      # Add 1 to timestep column
                    tmrt_t[i,1] += tmrt[j,0]

        sh_vec_unique = np.unique(sh_vec_t[:,1:], axis=0)
        for i in range(sh_vec_unique.shape[0]):
            sh_rg = sh_vec_t[:,1:] == sh_vec_unique[i,:]
            sh_b = np.all(sh_rg, axis=1)
            if np.sum(sh_b) > 1:
                sh_temp = sh_vec_t[sh_b,0]
                sh_vec_t[sh_b,0] = sh_temp[0]
                sh_temp2 = sh_temp[1:]
                for j in sh_temp2:
                    shadow_[shadow_ == j] = sh_temp[0]

        sh_vec_t, sh_idx, sh_inv = np.unique(sh_vec_t, return_index=True, return_inverse=True, axis=0)
        tmrt_t = tmrt_t[sh_idx,:]

        self.shadow = shadow_
        self.timesteps = sh_vec_t
        self.tmrt_rg = tmrt_t

class ClippedInputdata():
    '''This class clips the input rasters based on a buffer zone around the selected area. 
    This buffer zone is based on how far the tree shadow can reach outside the study area if a tree is at the edge '''
    #__slots__ = ('buildings', 'selected_area', 'dem', 'dsm', 'cdsm', 'cdsm_b', 'shadow', 'tmrt_ts', 'tmrt_s')
    __slots__ = ('dataSet', 'buildings', 'selected_area', 'dsm', 'cdsm', 'cdsm_b', 'shadow', 'tmrt_ts', 'tmrt_s', 'rows', 'cols', 'scale', 'lat', 'lon', 'gt', 'shadows_pad', 'tmrt_ts_pad', 'buildings_pad', 'rows_pad', 'cols_pad',
    'clip_rows', 'clip_cols')
    def __init__(self, treeinput, treerasters):
        # Estimate extent of selected area
        sa_rows, sa_cols = np.where(treeinput.selected_area == 1)

        # Empty arrays to store coordinates in
        self.clip_rows = np.zeros((2))
        self.clip_cols = np.zeros((2))

        # Calculate the buffer to clip from
        self.clip_rows[0] = sa_rows.min() - treerasters.buffer_y[0]
        if self.clip_rows[0] < 0:
            self.clip_rows[0] = 0
        self.clip_rows[1] = sa_rows.max() + treerasters.buffer_y[1]
        if self.clip_rows[1] > treeinput.selected_area.shape[0]:
            self.clip_rows[1] = treeinput.selected_area.shape[0]
        self.clip_cols[0] = sa_cols.min() - treerasters.buffer_x[0]
        if self.clip_cols[0] < 0:
            self.clip_cols[0] = 0
        self.clip_cols[1] = sa_cols.max() + treerasters.buffer_x[1]
        if self.clip_cols[1] > treeinput.selected_area.shape[1]:
            self.clip_cols[1] = treeinput.selected_area.shape[1]

        # self.clip_rows[0] = sa_rows.min(); self.clip_rows[1] = sa_rows.max()
        # self.clip_cols[0] = sa_cols.min(); self.clip_cols[1] = sa_cols.max()

        # Turn into integers to be able to use as indices
        self.clip_rows = np.int_(self.clip_rows); self.clip_cols = np.int_(self.clip_cols)

        # Clip input rasters
        self.buildings = treeinput.buildings[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.selected_area = treeinput.selected_area[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        # self.dem = treeinput.dem[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.dsm = treeinput.dsm[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.cdsm = treeinput.cdsm[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.cdsm_b = treeinput.cdsm_b[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.shadow = treeinput.shadow[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.tmrt_ts = treeinput.tmrt_ts[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]
        self.tmrt_s = treeinput.tmrt_s[self.clip_rows[0]:self.clip_rows[1], self.clip_cols[0]:self.clip_cols[1]]

        # Save other stuff from input rasters
        self.dataSet = treeinput.dataSet
        #self.aoi = treeinput.aoi
        self.scale = treeinput.scale
        self.lat = treeinput.lat
        self.lon = treeinput.lon
        self.gt = treeinput.gt

        # Rows and cols of new extent
        self.rows = self.buildings.shape[0]
        self.cols = self.buildings.shape[1]