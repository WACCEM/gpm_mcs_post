import numpy as np
import glob, os, sys
import xarray as xr

region = sys.argv[1]
date = sys.argv[2]
# region = 'npac'
# date = '20150101_20151231'

datadir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/'
# Original robust MCS stats file
statsfile = f'{datadir}robust_mcs_tracks_{date}.nc'

# MCS tracknumbers in AR file
arfile = f'{datadir}mcs_ar_tracknumbers_{date}.nc'

# Output directory and file
outdir = datadir
outfile = f'{outdir}robust_mcs_tracks_exar_{date}.nc'

# Read MCS stats file
ds = xr.open_dataset(statsfile)

# Read AR file
dsar = xr.open_dataset(arfile)

# Subtrack tracknumber by 1 to get track indices
ar_trackid = dsar.mcs_tracknumber - 1
ar_nhours = dsar.mcs_nhours

# Get all MCS track indices
alltracks = ds.tracks

# Find track indices that are not in AR
nonar_trackid = alltracks[~np.isin(alltracks, ar_trackid)]

# Select tracks not in AR
dsout = ds.sel(tracks=nonar_trackid)
# Update tracks coordinate
nmcs = len(nonar_trackid)
tracks = np.arange(0, nmcs, 1)
dsout['tracks'] = tracks

# Write to netCDF file
dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
print(f'Output saved: {outfile}')