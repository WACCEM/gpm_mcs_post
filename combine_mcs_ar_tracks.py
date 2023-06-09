"""
Combines the MCS track numbers in Atmospheric Rivers from two separate regions to the same file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "17-Mar-2022"

import numpy as np
import glob, os, sys
import xarray as xr
import yaml

if __name__ == '__main__':
    config_file = sys.argv[1]

    start_year = 2001
    end_year = 2020

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    stats_dir = config['stats_dir']
    output_dir = stats_dir

    basename1 = 'mcs_ar1_tracknumbers_'
    nchar = len(basename1)
    outbasename = 'mcs_ar_tracknumbers_'

    files1 = []
    files2 = []
    for iyear in range(start_year, end_year):
        files1.extend(glob.glob(f'{stats_dir}mcs_ar1_tracknumbers_{iyear}*.nc'))
        files2.extend(glob.glob(f'{stats_dir}mcs_ar2_tracknumbers_{iyear}*.nc'))
    nfiles1 = len(files1)
    nfiles2 = len(files2)

    if (nfiles1 != nfiles2):
        print(f'Number of files not matching!')
    else:
        print(f'Number of files: {nfiles1}')

    # Loop over each year
    for ii in range(0, nfiles1):
        # Open and combine two files
        ds = xr.open_mfdataset([files1[ii], files2[ii]], concat_dim='tracks', combine='nested')
        
        # Make output filename
        fname = os.path.basename(files1[ii])
        ftime = fname[nchar:nchar+27]
        outfilename = f'{output_dir}{outbasename}{ftime}.nc'

        # Write to netCDF file
        ds.to_netcdf(path=outfilename, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
        print(f'Output saved: {outfilename}')

        # import pdb; pdb.set_trace()