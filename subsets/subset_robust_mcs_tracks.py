"""
This script subsets robust MCS tracks initiated within a region and saves to a new netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "08-Dec-2020"

import sys
import os
import numpy as np
import xarray as xr

def subset_file(infile, outdir, lon_box, lat_box):

    status = 0

    # Read input MCS track file
    ds = xr.open_dataset(infile)

    # Get track initiation lat/lon location
    lon0 = ds.meanlon.isel(times=0)
    lat0 = ds.meanlat.isel(times=0)

    # Get track index for MCS initiation location within a region
    trackid = np.where((lon0 >= np.min(lon_box)) & (lon0 <= np.max(lon_box)) & \
                        (lat0 >= np.min(lat_box)) & (lat0 <= np.max(lat_box)))[0]


    # Subset the tracks from the dataset
    dsout = ds.sel(tracks=trackid, drop=True)
    # Add subset lat/lon box to global attribute
    dsout.attrs['lon_box'] = lon_box
    dsout.attrs['lat_box'] = lat_box

    # Strip the path to get input file name
    infilename = os.path.basename(infile)
    # Make output file name 
    outfile = f'{outdir}/{infilename}'

    # Write output file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')

    status = 1
    return status


def main():
    # Read from input
    lon_min = float(sys.argv[1])
    lon_max = float(sys.argv[2])
    lat_min = float(sys.argv[3])
    lat_max = float(sys.argv[4])
    infile = sys.argv[5]
    outdir = sys.argv[6]

    # Subset region boundary
    lon_box = [lon_min, lon_max]
    lat_box = [lat_min, lat_max]
    # lon_box = [-15, 45]
    # lat_box = [-20, 30]
    # lon_box = [-15, 40]
    # lat_box = [-15, 25]

    # Define output directory
    # outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/spac/subset_africa_stats/'
    os.makedirs(outdir, exist_ok=True)

    # Call function
    status = subset_file(infile, outdir, lon_box, lat_box)


if __name__ == "__main__":
    main()