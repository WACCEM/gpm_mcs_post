"""
This script subsets robust MCS tracks initiated within a region and time period and saves to a new netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "21-Jun-2021"

import sys
import os
import numpy as np
import xarray as xr
import pandas as pd

def subset_file(infile, outdir, lon_box, lat_box, start_time, end_time):

    status = 0

    # Read input MCS track file
    ds = xr.open_dataset(infile).load()

    # Get track initiation lat/lon location and time
    lon0 = ds['meanlon'].isel(times=0)
    lat0 = ds['meanlat'].isel(times=0)
    time0 = ds['base_time'].isel(times=0)

    # Get track index for MCS initiation location within a region and time within a period
    trackid = np.where(
        (lon0 >= np.min(lon_box)) & (lon0 <= np.max(lon_box)) & \
        (lat0 >= np.min(lat_box)) & (lat0 <= np.max(lat_box)) & \
        (time0 >= start_time) & (time0 <= end_time)
    )[0]

    # Subset the tracks from the dataset
    dsout = ds.sel(tracks=trackid, drop=True)
    # Add subset lat/lon box to global attribute
    dsout.attrs['lon_box'] = lon_box
    dsout.attrs['lat_box'] = lat_box
    dsout.attrs['subset_start_time'] = start_time.strftime('%Y-%m-%dT%H:%M:%S')
    dsout.attrs['subset_end_time'] = end_time.strftime('%Y-%m-%dT%H:%M:%S')

    # Strip the path to get input file name
    infilename = os.path.basename(infile)
    # Make output file name 
    outfile = f'{outdir}/{infilename}'

    # Write output file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')

    status = 1
    return status


def main():
    # Read from input
    lon_min = float(sys.argv[1])
    lon_max = float(sys.argv[2])
    lat_min = float(sys.argv[3])
    lat_max = float(sys.argv[4])
    sdate = sys.argv[5]
    edate = sys.argv[6]
    infile = sys.argv[7]
    outdir = sys.argv[8]

    # Subset region boundary
    lon_box = [lon_min, lon_max]
    lat_box = [lat_min, lat_max]

    # Convert start/end date to Pandas Timestamp
    # start_time = pd.to_datetime(start_time)
    # end_time = pd.to_datetime(end_time)

    # Get start date and time for the first day and hour of the month
    sdatetime = pd.date_range(sdate, periods=1, freq='MS')
    # Get end date and time for the last day and hour of the month
    edatetime = pd.date_range(edate, periods=1, freq='M') + pd.DateOffset(hours=23) + pd.DateOffset(minutes=59)
    # Get the time stamp values
    start_time = sdatetime[0]
    end_time = edatetime[0]

    # Define output directory
    os.makedirs(outdir, exist_ok=True)

    # Call function
    status = subset_file(infile, outdir, lon_box, lat_box, start_time, end_time)


if __name__ == "__main__":
    main()