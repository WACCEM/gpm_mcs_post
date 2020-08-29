import numpy as np
import os, glob
import time, datetime, calendar
from netCDF4 import Dataset, num2date, chartostring
# from scipy.ndimage import label, binary_dilation, generate_binary_structure
# from skimage.measure import regionprops
# from math import pi
# from scipy.stats import skew
import xarray as xr
# import pandas as pd
from multiprocessing import Pool
np.set_printoptions(threshold=np.inf)

#########################################################
# Load MCS track stats
print('Loading IR data')
print((time.ctime()))

# region = 'nam'
region = 'spac'
startdate = '20140101'
enddate = '20141231'
rr_min = 2.0  # [mm/h] rain rate threshold to define a PF
pixel_radius = 10.0  # [km] pixel size
pf_min_area_thresh = 400  # [km^2] minimum PF area threshold
# pf_link_area_thresh = 

# Convert minimum PF area to number of pixels
pf_min_npix = np.ceil(pf_min_area_thresh / (pixel_radius**2))

stats_path = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/'
pixelfile_path = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/mcstracking_ccs4_4h/{startdate}_{enddate}/'

mcsstats_filebase = 'robust_mcs_tracks_extc_'
pixel_filebase = 'mcstrack_'

# Robust MCS statistics file
robustmcs_file = f'{stats_path}{mcsstats_filebase}{startdate}_{enddate}.nc'
# Find all MCS pixel-level files
pixelfilelist = sorted(glob.glob(f'{pixelfile_path}{pixel_filebase}*.nc'))
nfiles = len(pixelfilelist)

print(robustmcs_file)


# Read Robust MCS statistics file
dsrobust = xr.open_dataset(robustmcs_file, decode_times=False)
ntracks = dsrobust.dims['tracks']
ntimes = dsrobust.dims['times']
nmaxpf = dsrobust.dims['nmaxpf']
rmcs_basetime = dsrobust.base_time.values

# Create variables for PF
pf_maxrate = np.full((ntracks, ntimes, nmaxpf), np.nan, dtype=float)

##############################################################
# Find precipitation feature in each mcs
print(f'Total Number of Tracks: {ntracks}')

from calc_pfstats_singlefile import calc_pfstats_singlefile

for filename in pixelfilelist:
    results = calc_pfstats_singlefile(filename, pixel_filebase, rmcs_basetime, rr_min, pf_min_npix, pixel_radius, nmaxpf)
    # import pdb; pdb.set_trace()

import pdb; pdb.set_trace()
