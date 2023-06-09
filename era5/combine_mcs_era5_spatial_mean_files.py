import os, sys
import glob
import time
import xarray as xr
import pandas as pd
import numpy as np

if __name__ == '__main__':

    root_dir = '/pscratch/sd/f/feng045/waccem/mcs_global/'
    dir3d = f'{root_dir}era5_3d/mean_300km/'
    dir2d = f'{root_dir}era5_2d/avg_space/'
    dir2d_derived = f'{root_dir}era5_2d_derived/mean_300km/'
    outdir = '/pscratch/sd/f/feng045/waccem/mcs_global/era5_envs/'
    outbasename = 'mcs_era5_mean_envs_'

    start_year = '2010-01-01'
    end_year = '2021-01-01'

    basename = 'mcs_era5_'
    
    # Generate time marks within the start/end datetime
    file_startdates = pd.date_range(start=start_year, end=end_year, freq='YS', inclusive='left')
    file_enddates = file_startdates + pd.offsets.DateOffset(years=1)
    # Convert to strings
    file_startdates = file_startdates.strftime('%Y%m%d.%H%M')
    file_enddates = file_enddates.strftime('%Y%m%d.%H%M')
    # For 2000
    # file_startdates = ['20000601.0000']
    # file_enddates = ['20010101.0000']

    # Loop over each date
    for ii in range(0, len(file_enddates)):
        # Track period
        track_period = f"{file_startdates[ii]}_{file_enddates[ii]}"
        print(track_period)

        # Find input files
        files3d = sorted(glob.glob(f'{dir3d}{basename}*{track_period}.nc'))
        nfiles3d = len(files3d)
        files2d = sorted(glob.glob(f'{dir2d}{basename}*{track_period}.nc'))
        nfiles2d = len(files2d)
        files2d_derived = sorted(glob.glob(f'{dir2d_derived}{basename}*{track_period}.nc'))
        nfiles2d_derived = len(files2d_derived)

        if (nfiles3d > 0) & (nfiles2d > 0) & (nfiles2d_derived > 0):
            # Read 3D variables
            ds3d = xr.open_mfdataset(files3d, combine='by_coords')

            # Read 2D variables
            ds2d = xr.open_mfdataset(files2d, combine='by_coords')
            # Rename Z to Z_sfc
            ds2d = ds2d.rename_vars({'Z': 'Z_sfc'})

            # Read 2D derived variables
            ds2dd = xr.open_mfdataset(files2d_derived, combine='by_coords')
            
            # Merged all DataSets
            dsout = xr.merge([ds3d, ds2d, ds2dd])

            # Output file
            outfilename = f'{outdir}{outbasename}{track_period}.nc'
            # Update global attributes
            dsout.attrs['Title'] = 'ERA5 domain mean environment data for MCS'
            dsout.attrs['Institution'] = 'Pacific Northwest National Laboratory'
            dsout.attrs['Contact'] = 'Zhe Feng, zhe.feng@pnnl.gov'
            dsout.attrs['Created_on'] = time.ctime(time.time())
            # Set encoding/compression for all variables
            comp = dict(zlib=True)
            encoding = {var: comp for var in dsout.data_vars}
            # Write to netcdf file
            dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                            unlimited_dims='tracks', encoding=encoding)
            print(f"Output saved as: {outfilename}")
        
        else:
            print(f'Missing files!')
            print(f'3D files: {nfiles3d}')
            print(f'2D files: {nfiles2d}')
            print(f'2D derived files: {nfiles2d_derived}')



        # import pdb; pdb.set_trace()