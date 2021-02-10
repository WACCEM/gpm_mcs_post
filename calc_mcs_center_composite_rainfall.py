"""
Composites MCS PF-center rainfall within a specified region and period and saves to a netCDF file.
All individual frames of MCS precipitation are saved in the output.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
12/03/2020 - Written.
"""

import numpy as np
import os, glob, sys
import time, datetime, calendar
import pytz
import xarray as xr
import pandas as pd
import dask
from dask.distributed import Client, LocalCluster
from composite_mcs_center import composite_mcs_center

if __name__ == "__main__":
    start_time = sys.argv[1]
    end_time = sys.argv[2]
    # start_time = '2010-06-01T00'
    # end_time = '2010-09-01T00'
    region = 'asia'

    # Domain boundary
    lon_box = [70, 95]
    lat_box = [0, 30]
    # Composite window size in one direction [km]
    # Output window size is 2 x window_size_km
    window_size_km = 700
    
    run_parallel = 1
    # Set up dask workers and threads
    n_workers = 30
    threads_per_worker = 1

    # Pixel size [km]
    pixel_radius = 10

    # Tracking start/end date
    year = start_time[0:4]
    startdate = f'{year}0101'
    enddate = f'{year}1231'

    # Format output date (yyyymmdd)
    d0 = pd.to_datetime(start_time).strftime("%Y%m%d")
    d1 = pd.to_datetime(end_time).strftime("%Y%m%d")

    # Input and output data directory
    stats_path = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/robust/filtered/'
    pixelfile_path = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/mcstracking_ccs4_4h/{startdate}_{enddate}/'
    # out_path = f'/global/cscratch1/sd/feng045/E3SM/GPM_IMERG/{region}/'
    out_path = f'/global/cscratch1/sd/feng045/E3SM/GPM_IMERG/{region}/mean_direction/'
    out_file = f'{out_path}mcs_center_composite_{d0}_{d1}.nc'
    
    mcsstats_filebase = 'robust_mcs_tracks_extc_'
    pixel_filebase = 'mcstrack_'

    os.makedirs(out_path, exist_ok=True)
    
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
                                    0, tzinfo=pytz.UTC).timetuple())

    # Read Robust MCS statistics file
    print(robustmcs_file)
    dsrobust = xr.open_dataset(robustmcs_file, decode_times=False)
    ntracks = dsrobust.dims['tracks']
    ntimes = dsrobust.dims['times']
    nmaxpf = dsrobust.dims['nmaxpf']

    # Calculate lifetime mean largest PF lat/lon for each track as proxy for MCS location
    pflon_avg = dsrobust['pf_lon'].isel(nmaxpf=0).mean(dim='times')
    pflat_avg = dsrobust['pf_lat'].isel(nmaxpf=0).mean(dim='times')

    # Get MCS initiation time
    init_time = dsrobust['base_time'].isel(times=0)

    # Calculate base time for start/end time window (divide datetime value by 1e9 to get seconds)
    bt0 = pd.to_datetime(start_time, utc=True).value//1e9
    bt1 = pd.to_datetime(end_time, utc=True).value//1e9

    # Find track indices within time and space window
    trackid = np.where((init_time >= bt0) & (init_time <= bt1) &
                       (pflon_avg >= np.min(lon_box)) & (pflon_avg <= np.max(lon_box)) &
                       (pflat_avg >= np.min(lat_box)) & (pflat_avg <= np.max(lat_box)))[0]

    # Subset the tracks from the dataset
    dsrobust = dsrobust.isel(tracks=trackid, drop=True)
    # Keep the original tracknumbers (this is important to match pixel-level files!)
    tracknumbers = dsrobust['tracks'].values
    rmcs_basetime = dsrobust['base_time'].values
    mcs_status = dsrobust['mcs_status'].values
    # Get the largest PF center lat/lon
    pf_lon = dsrobust['pf_lon'].isel(nmaxpf=0).values
    pf_lat = dsrobust['pf_lat'].isel(nmaxpf=0).values
    uspeed = dsrobust['uspeed'].values
    vspeed = dsrobust['vspeed'].values
    uspeed_avg = dsrobust['uspeed'].mean(dim='times').values
    vspeed_avg = dsrobust['vspeed'].mean(dim='times').values
    # print(f'Total Number of Tracks: {ntracks}')

    print(f'Total number of selected tracks: {len(tracknumbers)}')

    # Get track min/max base time
    rmcs_bt_min = np.nanmin(rmcs_basetime)
    rmcs_bt_max = np.nanmax(rmcs_basetime)

    # Find pixel files within the tracks time window
    subset_idx = np.where((pixelfile_basetime >= rmcs_bt_min) & (pixelfile_basetime <= rmcs_bt_max))[0]
    # Convert file list to numpy array then subset
    pixelfilelist = np.array(pixelfilelist)[subset_idx]
    pixelfile_basetime = pixelfile_basetime[subset_idx]
    nfiles = len(pixelfilelist)

    # Create variables to store composite output
    # Count total number of valid base time (i.e., total number of MCS PF)
    npf_all = np.count_nonzero(~np.isnan(rmcs_basetime) & (mcs_status == 1))
    # Get size of the MCS-center cut window
    nx = np.ceil(window_size_km / pixel_radius).astype(int)
    ny = nx
    rain_save = np.full((npf_all, ny*2+1, nx*2+1), 0, dtype=np.float)
    rain_save_ne = np.full((npf_all, ny*2+1, nx*2+1), 0, dtype=np.float)
    rain_save_se = np.full((npf_all, ny*2+1, nx*2+1), 0, dtype=np.float)
    rain_save_sw = np.full((npf_all, ny*2+1, nx*2+1), 0, dtype=np.float)
    rain_save_nw = np.full((npf_all, ny*2+1, nx*2+1), 0, dtype=np.float)
    basetime_save = np.full(npf_all, np.NaN, dtype=np.float32)
    basetime_save_ne = np.full(npf_all, np.NaN, dtype=np.float32)
    basetime_save_se = np.full(npf_all, np.NaN, dtype=np.float32)
    basetime_save_sw = np.full(npf_all, np.NaN, dtype=np.float32)
    basetime_save_nw = np.full(npf_all, np.NaN, dtype=np.float32)
    

    # Loop over subset files
    npf_save = 0
    npf_save_ne = 0
    npf_save_se = 0
    npf_save_sw = 0
    npf_save_nw = 0
    results = []
    final_result = []

    if run_parallel==0:

        for ifile in range(nfiles):
        # for ifile in range(60):
            filename = pixelfilelist[ifile]
            # Find all matching time indices from robust MCS stats file to the current pixel file
            # and require the MCS status == 1 (MCS size threshold is reached)
            # and must have a PF
            matchindices = np.array(np.where( (np.abs(rmcs_basetime - pixelfile_basetime[ifile]) < 1) & 
                                            (mcs_status == 1) & (~np.isnan(pf_lon)) ) )
            # The returned match indices are for [tracks, times] dimensions respectively
            idx_track = matchindices[0]
            idx_time = matchindices[1]
            nmatch = len(idx_track)
            if (nmatch > 0):            
                # Get the track numbers
                itracknum = tracknumbers[idx_track]
                # print(itracknum)
                # Get the PF center lat/lon
                ipflon = pf_lon[idx_track, idx_time]
                ipflat = pf_lat[idx_track, idx_time]
                # iuspeed = uspeed[idx_track, idx_time]
                # ivspeed = vspeed[idx_track, idx_time]
                iuspeed = uspeed_avg[idx_track]
                ivspeed = vspeed_avg[idx_track]

                # Call composite function
                result = composite_mcs_center(
                            filename,
                            itracknum,
                            ipflon,
                            ipflat,
                            iuspeed,
                            ivspeed,
                            nx,
                            ny,
                            )
                final_result.append(result)

    elif run_parallel==1:
        print(f'Parallel version by dask')

        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

        for ifile in range(nfiles):
        # for ifile in range(60):
            filename = pixelfilelist[ifile]
            # print(pixelfile_basetime[ifile])
            # Find all matching time indices from robust MCS stats file to the current pixel file
            # and require the MCS status == 1 (MCS size threshold is reached)
            # and must have a PF
            matchindices = np.array(np.where( (np.abs(rmcs_basetime - pixelfile_basetime[ifile]) < 1) & 
                                            (mcs_status == 1) & (~np.isnan(pf_lon)) ) )
            # The returned match indices are for [tracks, times] dimensions respectively
            idx_track = matchindices[0]
            idx_time = matchindices[1]
            nmatch = len(idx_track)
            if (nmatch > 0):            
                # Get the track numbers
                itracknum = tracknumbers[idx_track]
                # print(itracknum)
                # Get the PF center lat/lon
                ipflon = pf_lon[idx_track, idx_time]
                ipflat = pf_lat[idx_track, idx_time]
                # iuspeed = uspeed[idx_track, idx_time]
                # ivspeed = vspeed[idx_track, idx_time]
                iuspeed = uspeed_avg[idx_track]
                ivspeed = vspeed_avg[idx_track]

                # Call composite function
                result = dask.delayed(composite_mcs_center)(
                            filename,
                            itracknum,
                            ipflon,
                            ipflat,
                            iuspeed,
                            ivspeed,
                            nx,
                            ny,
                            )
                results.append(result)
        
        # Trigger dask comptation
        final_result = dask.compute(*results)

    # Get the number of returned items
    n_result = len(final_result)

    # Loop over the list of returned items to save to output arrays
    for ifile in range(n_result):
        # Get the return results for this pixel file
        tmp = final_result[ifile]
        if tmp is not None:
            rain_cut = tmp[0]
            npf = tmp[1]
            bt = tmp[2]

            rain_ne = tmp[3]
            rain_se = tmp[4]
            rain_sw = tmp[5]
            rain_nw = tmp[6]
            npf_ne = tmp[7]
            npf_se = tmp[8]
            npf_sw = tmp[9]
            npf_nw = tmp[10]
            bt_ne = tmp[11]
            bt_se = tmp[12]
            bt_sw = tmp[13]
            bt_nw = tmp[14]

            # Save the rain cut field
            if npf > 0:
                rain_save[npf_save:(npf_save+npf), :, :] = rain_cut
                # Add number of PFs
                npf_save += npf
                # Save base time for all returned PFs for this file
                basetime_save[npf_save:(npf_save+npf)] = bt

            if npf_ne > 0:
                rain_save_ne[npf_save_ne:(npf_save_ne+npf_ne), :, :] = rain_ne
                # Add number of PFs
                npf_save_ne += npf_ne
                # Save base time for all returned PFs for this file
                basetime_save_ne[npf_save_ne:(npf_save_ne+npf_ne)] = bt_ne
            if npf_se > 0:
                rain_save_se[npf_save_se:(npf_save_se+npf_se), :, :] = rain_se
                # Add number of PFs
                npf_save_se += npf_se
                # Save base time for all returned PFs for this file
                basetime_save_se[npf_save_se:(npf_save_se+npf_se)] = bt_se
            if npf_sw > 0:
                rain_save_sw[npf_save_sw:(npf_save_sw+npf_sw), :, :] = rain_sw
                # Add number of PFs
                npf_save_sw += npf_sw
                # Save base time for all returned PFs for this file
                basetime_save_sw[npf_save_sw:(npf_save_sw+npf_sw)] = bt_sw
            if npf_nw > 0:
                rain_save_nw[npf_save_nw:(npf_save_nw+npf_nw), :, :] = rain_nw
                # Add number of PFs
                npf_save_nw += npf_nw
                # Save base time for all returned PFs for this file
                basetime_save_nw[npf_save_nw:(npf_save_nw+npf_nw)] = bt_nw

    # Find indices where base time is valid
    idx_valid = np.nonzero(~np.isnan(basetime_save))[0]
    # idx_valid_ne = np.nonzero(~np.isnan(basetime_save_ne))[0]
    # idx_valid_se = np.nonzero(~np.isnan(basetime_save_se))[0]
    # idx_valid_sw = np.nonzero(~np.isnan(basetime_save_sw))[0]
    # idx_valid_nw = np.nonzero(~np.isnan(basetime_save_nw))[0]
    # ntimes = len(idx_valid)
    # Remove invalid times
    basetime_save = basetime_save[idx_valid]
    basetime_save_ne = basetime_save_ne[idx_valid]
    basetime_save_se = basetime_save_se[idx_valid]
    basetime_save_sw = basetime_save_sw[idx_valid]
    basetime_save_nw = basetime_save_nw[idx_valid]
    rain_save = rain_save[idx_valid, :, :]
    rain_save_ne = rain_save_ne[idx_valid, :, :]
    rain_save_se = rain_save_se[idx_valid, :, :]
    rain_save_sw = rain_save_sw[idx_valid, :, :]
    rain_save_nw = rain_save_nw[idx_valid, :, :]

    # Create x, y coordinates
    xcoord = (np.arange(nx*2+1)-nx) * pixel_radius
    ycoord = (np.arange(ny*2+1)-ny) * pixel_radius

    # Define xarray dataset for output
    dsout = xr.Dataset({'precipitation': (['time', 'y', 'x'], rain_save), \
                        'precipitation_ne': (['time', 'y', 'x'], rain_save_ne), \
                        'precipitation_se': (['time', 'y', 'x'], rain_save_se), \
                        'precipitation_sw': (['time', 'y', 'x'], rain_save_sw), \
                        'precipitation_nw': (['time', 'y', 'x'], rain_save_nw), \
                        }, \
                        coords={'time': (['time'], basetime_save), \
                                'time_ne': (['time'], basetime_save_ne), \
                                'time_se': (['time'], basetime_save_se), \
                                'time_sw': (['time'], basetime_save_sw), \
                                'time_nw': (['time'], basetime_save_nw), \
                                'y': (['y'], ycoord), \
                                'x': (['x'], xcoord)}, \
                        attrs={'title': 'MCS PF-center composite', \
                            'lon_box':lon_box, \
                            'lat_box':lat_box, \
                            'start_time':start_time, \
                            'end_time':end_time, \
                            'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                            'created_on':time.ctime(time.time())}
                        )

    dsout.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
    dsout.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'

    # dsout['time_ne'].attrs = dsout['time'].attrs
    # dsout['time_se'].attrs = dsout['time'].attrs
    # dsout['time_sw'].attrs = dsout['time'].attrs
    # dsout['time_nw'].attrs = dsout['time'].attrs

    dsout.y.attrs['long_name'] = 'Relative y-distance from MCS PF center'
    dsout.y.attrs['units'] = 'km'

    dsout.x.attrs['long_name'] = 'Relative x-distance from MCS PF center'
    dsout.x.attrs['units'] = 'km'

    dsout.precipitation.attrs['long_name'] = 'MCS precipitation'
    dsout.precipitation.attrs['units'] = 'mm/h'

    dsout.precipitation_ne.attrs['long_name'] = 'MCS precipitation (Northeast moving)'
    dsout.precipitation_ne.attrs['units'] = 'mm/h'

    dsout.precipitation_se.attrs['long_name'] = 'MCS precipitation (Southeast moving)'
    dsout.precipitation_se.attrs['units'] = 'mm/h'

    dsout.precipitation_sw.attrs['long_name'] = 'MCS precipitation (Southwest moving)'
    dsout.precipitation_sw.attrs['units'] = 'mm/h'

    dsout.precipitation_nw.attrs['long_name'] = 'MCS precipitation (Northwest moving)'
    dsout.precipitation_nw.attrs['units'] = 'mm/h'

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netcdf file
    dsout.to_netcdf(path=out_file, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)
    print('Output saved as: ', out_file)

    # import pdb; pdb.set_trace()