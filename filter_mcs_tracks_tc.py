"""
This script filters TC tracks from robust MCS statistics files and saves to a new statistics netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "29-Mar-2020"

import numpy as np
import glob, os, sys
import xarray as xr

region = sys.argv[1]
date = sys.argv[2]
# region = 'npac'
# date = '20150101_20151231'

datadir = os.path.expandvars('$SCRATCH') + f'/waccem/mcs_region/{region}/stats_ccs4_4h/'
statsfile = f'{datadir}robust_mcs_tracks_{date}.nc'
print(f'Statsfile: {statsfile}')

# MCS tracknumbers in TC file
tcfile = f'{datadir}mcs_tc_{date}.nc'

# Output directory and file
outdir = datadir
outfile = f'{outdir}robust_mcs_tracks_extc_{date}.nc'

# Read MCS stats file
ds = xr.open_dataset(statsfile)

# Read TC file
dstc = xr.open_dataset(tcfile)

# Subtrack tracknumber by 1 to get track indices
tc_trackid = dstc.mcs_tracknumber - 1
tc_nhours = dstc.mcs_nhours

# Get all MCS track indices
alltracks = ds.tracks

# Find track indices that are not in TC
nontc_trackid = (alltracks[~np.isin(alltracks, tc_trackid)]).values

# Select tracks not in TC
dsout = ds.isel(tracks=nontc_trackid, drop=True)
# Update tracks coordinate
nmcs = len(nontc_trackid)
tracks = np.arange(0, nmcs, 1)
dsout['tracks'] = tracks

# import pdb; pdb.set_trace()

# Write to netCDF file
dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
print(f'Output saved: {outfile}')
