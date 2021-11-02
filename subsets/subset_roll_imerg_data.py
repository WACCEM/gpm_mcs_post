"""
Subset and roll the GPM Tb+IMERG combined global data from -180~180 to 0~360.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
04/16/2021 - Written.
"""

import numpy as np
import glob, os, sys
import xarray as xr
import pandas as pd
# import time, datetime, calendar, pytz
import dask
from dask.distributed import Client, LocalCluster


#-----------------------------------------------------------------------
def subset_data(infiles, out_dir, minlat, maxlat):
    # print('Reading input data ... ')
    ds = xr.open_dataset(infiles, decode_times=False)
    # Subset latitude
    ds = ds.sel(lat=slice(minlat, maxlat))
    lon_attrs = ds['lon'].attrs

    # Roll lon=1800 to make the data start from 0~360 instead of -180~180
    ds = ds.roll(lon=1800, roll_coords=True)
    # Convert longitude coordinates from -180~180 to 0~360
    lon360 = ds['lon'].values % 360
    ds = ds.assign_coords(lon=lon360)
    ds['lon'].attrs = lon_attrs
    ds['lon'].attrs['valid_min'] = 0.
    ds['lon'].attrs['valid_max'] = 360.

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in ds.data_vars}
    # Update time variable dtype as 'double' for better precision
    bt_dict = {'time': {'zlib':True, 'dtype':'float64'}}
    encoding.update(bt_dict)

    # Write to netcdf file
    filename_out = out_dir + os.path.basename(infiles)
    ds.to_netcdf(path=filename_out, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)
    print(f'Saved as: {filename_out}')

    # import pdb; pdb.set_trace()
    status = 1
    return status


if __name__=='__main__':

    start_date = sys.argv[1]
    end_date = sys.argv[2]
    run_parallel = int(sys.argv[3])
    # start_date = '2020-01-20'
    # end_date = '2020-02-28'
    minlat, maxlat = -60., 60.

    in_year = start_date[0:4]
    in_dir = f'/global/cscratch1/sd/feng045/waccem/IR_IMERG_Combined/{in_year}/'
    out_dir = f'/global/cscratch1/sd/feng045/waccem/IR_IMERG_Combined_60NS/{in_year}/'

    basename = 'merg_'
    out_basename = 'merg_'

    os.makedirs(out_dir, exist_ok=True)

    # Number of workers for Dask
    n_workers = 32
    # Threads per worker
    threads_per_worker = 1

    # Create a range of dates
    dates = pd.date_range(start=start_date, end=end_date, freq='D')
    dates = dates.strftime('%Y%m%d')
    # Find all files from the dates    
    infiles = []
    for ii in range(0, len(dates)):
        infiles.extend(sorted(glob.glob(f'{in_dir}{basename}{dates[ii]}*.nc')))
    nf1 = len(infiles)

    if run_parallel==0:
        # serial version
        for ifile in range(nf1):
            print(infiles[ifile])
            status = subset_data(infiles[ifile], out_dir, minlat, maxlat)

    elif run_parallel==1:
        # parallel version
            
        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

        results = []
        for ifile in range(nf1):
            print(infiles[ifile])

            # Call subset function
            status = dask.delayed(subset_data)(infiles[ifile], out_dir, minlat, maxlat)
            results.append(status)

        # Collect results from Dask
        # print("Precompute step")
        results = dask.compute(*results)