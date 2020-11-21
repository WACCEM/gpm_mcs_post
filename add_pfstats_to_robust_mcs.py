import numpy as np
import os, glob, sys
import time, datetime, calendar
from pytz import utc
import xarray as xr
import dask
from dask.distributed import Client, LocalCluster
from calc_pfstats_singlefile import calc_pfstats_singlefile


if __name__ == "__main__":
    year = int(sys.argv[1])
    region = sys.argv[2]
    # year = 2002
    # region = 'asia'
    
    # Set up dask workers and threads
    n_workers = 32
    threads_per_worker = 2

    startdate = f'{year}0101'
    enddate = f'{year}1231'
    rr_min = 2.0  # [mm/h] rain rate threshold to define a PF
    pixel_radius = 10.0  # [km] pixel size
    pf_min_area_thresh = 0  # [km^2] minimum PF area threshold
    # Threshold of pixel-level land fraction to define a land pixel
    landfrac_thresh = 90    # [%]

    # Convert minimum PF area to number of pixels
    pf_min_npix = np.ceil(pf_min_area_thresh / (pixel_radius**2))

    stats_path = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/robust/filtered/'
    pixelfile_path = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/mcstracking_ccs4_4h/{startdate}_{enddate}/'
    mcsstats_filebase = 'robust_mcs_tracks_extc_'
    pixel_filebase = 'mcstrack_'
    landmask_file = f'/global/cscratch1/sd/feng045/waccem/mcs_region/map-data/IMERG_landmask_{region}.nc'

    out_path = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/'
    out_file = f'{out_path}{mcsstats_filebase}{startdate}_{enddate}.nc'

    
    # Robust MCS statistics file
    robustmcs_file = f'{stats_path}{mcsstats_filebase}{startdate}_{enddate}.nc'
    # Find all MCS pixel-level files
    pixelfilelist = sorted(glob.glob(f'{pixelfile_path}{pixel_filebase}*.nc'))
    nfiles = len(pixelfilelist)

    # Compute base time from pixel file names
    nleadingchar = np.array(len(pixel_filebase)).astype(int)
    pixelfile_basetime = np.full(nfiles, np.nan, dtype=np.float)
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

    # Create variables for PF
    missing_val = -999
    pf_npf = np.full((ntracks, ntimes), missing_val, dtype=int)
    pf_landfrac = np.full((ntracks, ntimes), np.nan, dtype=float)
    # pf_pflon = np.full((ntracks, ntimes, nmaxpf), np.nan, dtype=float)
    # pf_pflat = np.full((ntracks, ntimes, nmaxpf), np.nan, dtype=float)
    # pf_pfnpix = np.full((ntracks, ntimes, nmaxpf), missing_val, dtype=int)
    # pf_maxrate = np.full((ntracks, ntimes, nmaxpf), np.nan, dtype=float)

    # Initialize dask
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
    client = Client(cluster)

    # Create a list to store matchindices for each pixel file
    nmatchcloud_all = []
    matchindices_all = []
    results = []
    # for ifile in range(10):
    for ifile in range(nfiles):
        filename = pixelfilelist[ifile]
        # Find all matching time indices from robust MCS stats file to the current pixel file
        matchindices = np.array(np.where(np.abs(rmcs_basetime - pixelfile_basetime[ifile]) < 1))
        # The returned match indices are for [tracks, times] dimensions respectively
        idx_track = matchindices[0]
        idx_time = matchindices[1]
        # Save matchindices for the current pixel file to the overall list
        nmatchcloud_all.append(len(idx_track))
        matchindices_all.append(matchindices)
        # import pdb; pdb.set_trace()

        # Call function to calculate PF stats
        result = dask.delayed(calc_pfstats_singlefile)(
                    filename, 
                    pixel_filebase, 
                    # rmcs_basetime, 
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

    # Put the PF results to output track stats variables
    # Loop over each pixel file (parallel return results)
    for ifile in range(nfiles):
    # for ifile in range(10):
        # Get the return results for this pixel file
        tmp = final_result[ifile]
        if tmp is not None:
            # Get the match track indices (matchindicestmp contains: [track_index, time_index]) for this pixel file
            # nmatchpftmp = tmp[0]
            # matchindicestmp = tmp[1]
            nmatchpftmp = nmatchcloud_all[ifile]
            matchindicestmp = matchindices_all[ifile]
            # Loop over each match cloud in this pixel file
            for imatch in range(nmatchpftmp):
                pf_npf[matchindicestmp[0, imatch], matchindicestmp[1, imatch]] = tmp[0][imatch]
                pf_landfrac[matchindicestmp[0, imatch], matchindicestmp[1, imatch]] = tmp[1][imatch]
                # pf_pflon[matchindicestmp[0, imatch], matchindicestmp[1, imatch], :] = tmp[3][imatch, :]
                # pf_pflat[matchindicestmp[0, imatch], matchindicestmp[1, imatch], :] = tmp[4][imatch, :]
                # pf_pfnpix[matchindicestmp[0, imatch], matchindicestmp[1, imatch], :] = tmp[5][imatch, :]
                # pf_maxrate[matchindicestmp[0, imatch], matchindicestmp[1, imatch], :] = tmp[6][imatch, :]
            # import pdb; pdb.set_trace()

    # import pdb; pdb.set_trace()

    # Replace landfrac in the robust statistics file
    pf_landfrac = xr.DataArray(pf_landfrac, coords={'tracks':dsrobust.tracks, 'times':dsrobust.times}, dims=('tracks', 'times'))
    dsrobust['pf_landfrac'] = pf_landfrac

    # Write output
    dsrobust.to_netcdf(path=out_file, mode="w", format="NETCDF4_CLASSIC")
    print(f'Output saved as: {out_file}')
    
    # import pdb; pdb.set_trace()