"""
Combine individual MCS ERA5 environmental files for the same year into a single file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "01-Sep-2021"

import xarray as xr
import sys, os, glob
import pandas as pd
import yaml

if __name__ == "__main__":

    config_file = sys.argv[1]
    start_year = sys.argv[2]
    end_year = sys.argv[3]

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    mcs_dir = config['stats_dir']
    data_dir = config['output_dir']
    envs_dir = config['envs_dir']
    mcsfile_basename = config['mcsfile_basename']
    in_basename = config['in_env_basename']
    out_basename = config['out_env_basename']

    # Make output directory
    os.makedirs(envs_dir, exist_ok=True)

    # Generate time marks within the start/end datetime
    # file_startdates = pd.date_range(start=start_year, end=end_year, freq='YS', closed='left').strftime('%Y%m%d')
    # file_enddates = pd.date_range(start=start_year, end=end_year, freq='Y', closed='left').strftime('%Y%m%d')
    file_startdates = [f'{start_year}0120.0000']
    file_enddates = [f'{end_year}0301.0000']

    # Loop over each time
    for ii in range(0, len(file_enddates)):
        
        track_period = f"{file_startdates[ii]}_{file_enddates[ii]}"
        print(track_period)

        # mcs_file = f"{mcs_dir}{mcsfile_basename}{track_period}.nc"
        outfilename = f"{envs_dir}{out_basename}{track_period}.nc"

        # Find all variable files        
        era5files = sorted(glob.glob(f"{data_dir}{in_basename}{track_period}.nc"))
        print(f'Reading ERA5 files: {len(era5files)}')

        # Read all variable files
        dsout = xr.open_mfdataset(era5files, combine='by_coords', compat='override')
        print(f'Done reading.')
        
        # Set encoding/compression for all variables
        comp = dict(zlib=True, dtype='float32')
        encoding = {var: comp for var in dsout.data_vars}

        # Write to netcdf file
        dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                        unlimited_dims='tracks', encoding=encoding)
        print(f"Output saved as: {outfilename}")
