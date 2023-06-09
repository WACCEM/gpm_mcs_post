"""
Combine individual MCS center 2D files for the same year into a single file.
This code should be run on an interactive node.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "19-May-2023"

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
    data_dir = config['output_dir']
    combine_dir = config['combine_dir']
    in_basename = 'mcs_2d_mask_'
    out_basename = 'mcs_2d_mask_'

    # Make output directory
    os.makedirs(combine_dir, exist_ok=True)

    # Generate time marks within the start/end datetime
    file_startdates = pd.date_range(start=start_year, end=end_year, freq='YS', inclusive='left')
    file_enddates = file_startdates + pd.offsets.DateOffset(years=1)
    # Convert to strings
    file_startdates = file_startdates.strftime('%Y%m%d.%H%M')
    file_enddates = file_enddates.strftime('%Y%m%d.%H%M')

    # Loop over each time
    for ii in range(0, len(file_enddates)):
        
        track_period = f"{file_startdates[ii]}_{file_enddates[ii]}"
        print(track_period)

        # Output filename
        outfilename = f"{combine_dir}{out_basename}{track_period}.nc"

        # Find all variable files        
        infiles = sorted(glob.glob(f"{data_dir}{in_basename}{track_period}_t*.nc"))
        print(f'Reading input files: {len(infiles)}')

        # Read all variable files
        dsout = xr.open_mfdataset(infiles, combine='by_coords')
        # dsout = xr.open_mfdataset(infiles, combine='by_coords', compat='override')
        print(f'Done reading.')
        
        # Set encoding/compression for all variables
        comp = dict(zlib=True, dtype='float32')
        encoding = {var: comp for var in dsout.data_vars}

        # Write to netcdf file
        dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                        unlimited_dims='tracks', encoding=encoding)
        print(f"Output saved as: {outfilename}")
