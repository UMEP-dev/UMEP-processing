import re
from PIL import Image
import rasterio
import numpy as np

# Open the landcover file
lc_filepath = r'/home/lemap/Documents/suede/datasets/gotenburg/cut2/lc_cut.tif'
lc_filepathout = r'/home/lemap/Documents/suede/datasets/gotenburg/cut2/lc_cut_out.tif'
landcover = rasterio.open(lc_filepath)
crs = landcover.profile
lc_data = landcover.read()

lc_data[lc_data == 4] = 2
lc_data[lc_data == 3] = 2

new_lc = rasterio.open(lc_filepathout, 'w', **crs)
new_lc.write(lc_data)