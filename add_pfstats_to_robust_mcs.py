"""
Calculate MCS PF land fraction using track statistics and pixel-level files 
and updates the land fraction in the MCS track statistics file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "02-Sep-2021"

import numpy as np
import os, glob, sys
import time, datetime, calendar
from pytz import utc
import xarray as xr
import yaml
import dask
from dask.distributed import Client, LocalCluster
from calc_pfstats_singlefile import calc_pfstats_singlefile


if __name__ == "__main__":

    config_file = sys.argv[1]
    year = int(sys.argv[2])
    region = sys.argv[3]
    # config_file = config_gpm_region.yml
    # year = 2002
    # region = 'asia'

    startdate = f'{year}0101'
    enddate = f'{year}1231'
    
    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)    
    n_workers = config['n_workers']
    threads_per_worker = config['threads_per_worker']
    rr_min = config['rr_min']
    pixel_radius = config['pixel_radius']
    pf_min_area_thresh = config['pf_min_area_thresh']
    landfrac_thresh = config['landfrac_thresh']
    stats_path = config['stats_path'].replace('REGION', region)
    pixelfile_path = config['pixelfile_path'].replace('REGION', region) + f'/{startdate}_{enddate}/'
    landmask_file = config['landmask_file'].replace('REGION', region)
    mcsstats_filebase = config['mcsstats_filebase']
    pixel_filebase = config['pixel_filebase']
    out_path = config['out_path'].replace('REGION', region)
    out_file = f'{out_path}{mcsstats_filebase}{startdate}_{enddate}.nc'
  
    # Robust MCS statistics file
    robustmcs_file = f'{stats_path}{mcsstats_filebase}{startdate}_{enddate}.nc'

    # Convert minimum PF area to number of pixels
    pf_min_npix = np.ceil(pf_min_area_thresh / (pixel_radius**2))

    # Find all MCS pixel-level files
    pixelfilelist = sorted(glob.glob(f'{pixelfile_path}{pixel_filebase}*.nc'))
    nfiles = len(pixelfilelist)

    # Compute base time from pixel file names
    nleadingchar = np.array(len(pixel_filebase)).astype(int)
    pixelfile_basetime = np.full(nfiles, np.nan, dtype=float)
    for ifile in range(nfiles):
        fname = os.path.basename(pixelfilelist[ifile])
        pixelfile_basetime[ifile] = calendar.timegm(
                                    datetime.datetime(
                                    int(fname[nleadingchar:nleadingchar+4]),
                                    int(fname[nleadingchar+4:nleadingchar+6]),
                                    int(fname[nleadingchar+6:nleadingchar+8]),
                                    int(fname[nleadingchar+9:nleadingchar+11]),
                                    int(fname[nleadingchar+11:nleadingchar+13]),
                                    0, tzinfo=utc).timetuple())

    # Read Robust MCS statistics file
    print(robustmcs_file)
    dsrobust = xr.open_dataset(robustmcs_file, decode_times=False)
    ntracks = dsrobust.dims['tracks']
    ntimes = dsrobust.dims['times']
    nmaxpf = dsrobust.dims['nmaxpf']
    rmcs_basetime = dsrobust.base_time.values

    print(f'Total Number of Tracks: {ntracks}')

    # Initialize dask
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
    client = Client(cluster)

    # Create a list to store matchindices for each pixel file
    trackindices_all = []
    timeindices_all = []
    results = []
    # for ifile in range(10):
    for ifile in range(nfiles):
        filename = pixelfilelist[ifile]
        print(filename)

        # Find all matching time indices from robust MCS stats file to the current pixel file
        matchindices = np.array(np.where(np.abs(rmcs_basetime - pixelfile_basetime[ifile]) < 1))
        # The returned match indices are for [tracks, times] dimensions respectively
        idx_track = matchindices[0]
        idx_time = matchindices[1]
        # Save matchindices for the current pixel file to the overall list
        trackindices_all.append(idx_track)
        timeindices_all.append(idx_time)

        # Call function to calculate PF stats
        result = dask.delayed(calc_pfstats_singlefile)(
                    filename, 
                    pixel_filebase, 
                    idx_track,
                    rr_min, 
                    pf_min_npix, 
                    pixel_radius, 
                    nmaxpf, 
                    landmask_file, 
                    landfrac_thresh
                    )
        results.append(result)

    # Trigger dask computation
    final_result = dask.compute(*results)

    # Create variables for PF
    missing_val = -999
    pf_npf = np.full((ntracks, ntimes), missing_val, dtype=int)
    pf_landfrac = np.full((ntracks, ntimes), np.nan, dtype=float)

    # Put the PF results to output track stats variables
    # Loop over each pixel file (parallel return results)
    for ifile in range(nfiles):
    # for ifile in range(10):
        # Get the return results for this pixel file
        iVAR = final_result[ifile]
        if iVAR is not None:
            trackindices = trackindices_all[ifile]
            timeindices = timeindices_all[ifile]
            # Put variable to the MCS track array
            pf_npf[trackindices, timeindices] = iVAR['pf_npf']
            pf_landfrac[trackindices, timeindices] = iVAR['pf_landfrac']

    # Replace landfrac in the robust statistics file
    pf_landfrac = xr.DataArray(pf_landfrac, coords={'tracks':dsrobust.tracks, 'times':dsrobust.times}, dims=('tracks', 'times'))
    dsrobust['pf_landfrac'] = pf_landfrac
    dsrobust.attrs['updated_on'] = time.ctime(time.time())

    # Write output
    dsrobust.to_netcdf(path=out_file, mode="w", format="NETCDF4_CLASSIC")
    print(f'Output saved as: {out_file}')