"""
This script filters both AR and TC tracks from robust MCS statistics files and 
saves to a new statistics netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "05-Mar-2024"

import numpy as np
import glob, os, sys
import xarray as xr
import yaml

if __name__ == '__main__':

    indates = sys.argv[1]
    config_file = sys.argv[2]

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    stats_dir = config['stats_dir']

    # Original robust MCS stats file
    statsfile = f'{stats_dir}mcs_tracks_final_{indates}.nc'

    # MCS tracknumbers in AR file
    # arfile = f'{stats_dir}mcs_ar_tracknumbers_{indates}.nc'
    tcfile = f'{stats_dir}mcs_tc_tracknumbers_{indates}.nc'

    # Output directory and file
    outdir = stats_dir
    outfile = f'{outdir}mcs_tracks_final_extc_{indates}.nc'

    # Read MCS stats file
    ds = xr.open_dataset(statsfile)

    # Read AR file
    # dsar = xr.open_dataset(arfile)
    dstc = xr.open_dataset(tcfile)

    # Subtrack tracknumber by 1 to get track indices
    # ar_trackid = dsar.mcs_tracknumber - 1
    # ar_nhours = dsar.mcs_nhours
    tc_trackid = dstc.mcs_tracknumber - 1
    tc_nhours = dstc.mcs_nhours

    # Find unique tracknumbers from AR and TC
    # artc_trackid = np.unique(np.concatenate((ar_trackid.data,tc_trackid.data)))
    artc_trackid = np.unique(tc_trackid.data)
    # import pdb; pdb.set_trace()

    # Get all MCS track indices
    alltracks = ds.tracks

    # Find track indices that are not in AR&TC
    # nonartc_trackid = alltracks[~np.isin(alltracks, artc_trackid)]
    nonartc_trackid = (alltracks[~np.isin(alltracks, artc_trackid)]).values

    # Select tracks not in AR
    # dsout = ds.sel(tracks=nonartc_trackid)
    dsout = ds.isel(tracks=nonartc_trackid, drop=True)
    # Update tracks coordinate
    nmcs = len(nonartc_trackid)
    tracks = np.arange(0, nmcs, 1)
    dsout['tracks'] = tracks
    # import pdb; pdb.set_trace()

    # Write to netCDF file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')
