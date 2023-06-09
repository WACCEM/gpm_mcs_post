"""
Extracts ERA5 3D variable centered at each MCS saves the output to a netCDF file.
This version processes a subset of tracks provided from input arguments.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "1-June-2022"

import numpy as np
import xarray as xr
import pandas as pd
import sys, os
import time
import yaml
import dask
from dask.distributed import Client, LocalCluster

def extract_era5_prof(filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname):

    # Read ERA5 data
    dse = xr.open_dataset(filename)
    lat_e5 = dse.latitude
    lon_e5 = dse.longitude
    time_e5 = dse.time
    dse = dse.sel(level=slice(min(level_lims), max(level_lims)))
    level = dse.level
    VAR_attrs = {
        'long_name': dse[varname].attrs['long_name'],
        'units': dse[varname].attrs['units'],
    }

    # Number of tracks in the day
    ntracks_day = len(mcs_lat)
    nlevels = len(level)

    # Make array to store output
    VAR_out = np.full((ntracks_day, nlevels, 2*ny+1, 2*nx+1), np.NaN, dtype=float)

    # Loop over each MCS
    for itrack in range(0, ntracks_day):
        # print(f'{itrack}: {mcs_time[itrack]}')
        # Find closest grid point and time index in ERA5 
        lat_idx = np.abs(lat_e5.values - mcs_lat[itrack]).argmin()
        lon_idx = np.abs(lon_e5.values - mcs_lon[itrack]).argmin()
        t5_idx = np.abs(time_e5.values - mcs_time[itrack]).argmin()
        
        # Select the time, level, and region
        # Note index+1 is needed to include the last value in the range
        iVAR = dse[varname].isel(
            time=t5_idx,
            latitude=slice(lat_idx-ny, lat_idx+ny+1),
            longitude=slice(lon_idx-nx, lon_idx+nx+1),
        )
        # Check box size to avoid error (e.g. at date line)
        inx = iVAR.sizes['longitude']
        iny = iVAR.sizes['latitude']
        if (iny > 1) & (inx > 1):
            VAR_out[itrack,:,0:iny,0:inx] = iVAR.values

    # Put output variables to a dictionary for easier acceess
    out_dict = {'VAR':VAR_out, 'VAR_attrs':VAR_attrs, 'level':level}
    print(f'Done processing: {filename}')
    
    return out_dict

if __name__ == "__main__":

    basename = sys.argv[1]
    varname = sys.argv[2]
    config_file = sys.argv[3]
    track_period = sys.argv[4]
    track_start = int(sys.argv[5])
    track_end = int(sys.argv[6])
    digits = int(sys.argv[7])

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    era5_dir = config['era5_dir']
    stats_dir = config['stats_dir']
    output_dir = config['output_dir']
    mcsfile_basename = config['mcsfile_basename']
    nx = config['nx']
    ny = config['ny']
    level_lims = config['level_lims']
    nhours = config['nhours']
    ntimes_max = config['ntimes_max']
    run_parallel = config['run_parallel']
    n_workers = config['n_workers']
    threads_per_worker = config['threads_per_worker']

    # Update output directory by adding a sub-directory year
    start_year = track_period[0:4]
    output_dir = f'{output_dir}{start_year}/'
    mcs_file = f"{stats_dir}{mcsfile_basename}{track_period}.nc"
    os.makedirs(output_dir, exist_ok=True)

    # Set up Dask server
    if run_parallel == 1:
        # Set Dask temporary directory for workers
        dask_tmp_dir = config.get("dask_tmp_dir", "/tmp")
        dask.config.set({'temporary-directory': dask_tmp_dir})
        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

    # Read robust MCS statistics
    dsm = xr.open_dataset(mcs_file)
    # Subset MCS times to reduce array size
    # Most valid MCS data are within 0:ntimes_max
    dsm = dsm.isel(times=slice(0, ntimes_max))
    ntracks_all = dsm.sizes['tracks']

    # Make output filename
    track_start_str = str(track_start).zfill(digits)
    outfilename = f"{output_dir}mcs_era5_{varname}_{track_period}_t{track_start_str}.nc"

    # Select tracks for this part        
    coord_tracks = dsm['tracks'].sel(tracks=slice(track_start, track_end))
    rmcs_lat = dsm['meanlat'].sel(tracks=slice(track_start, track_end))
    rmcs_lon = dsm['meanlon'].sel(tracks=slice(track_start, track_end))
    rmcs_basetime = dsm['base_time'].sel(tracks=slice(track_start, track_end))
    ntracks = len(coord_tracks)

    # Check if longitude is [-180~+180], if so convert it to [0~360] to match ERA5
    if np.nanmin(rmcs_lon) < 0:
        rmcs_lon = rmcs_lon % 360
        print('MCS longitudes are [-180~+180], converted to [0-360] to match ERA5.')

    # Get end times for all tracks
    # Sum over time dimension for valid basetime indices, -1 to get the last valid time index for each track
    # This is the end time index of each track (i.e. +1 equals the lifetime of each track)
    end_time_idx = np.sum(np.isfinite(rmcs_basetime), axis=1)-1
    # Apply fancy indexing to base_time: a tuple that indicates for each track, get the end time index
    end_basetime = rmcs_basetime[(np.arange(0,ntracks), end_time_idx)]
        
    # Get the min/max of all base_times
    min_basetime = rmcs_basetime.sel(times=0).min()
    max_basetime = end_basetime.max()
    # Subtract min base_time by X days to include pre-initiation times from previous days
    offset_days = round(2 * nhours / 24)
    min_basetime_offset = min_basetime.dt.round('H').values - pd.DateOffset(days=offset_days)

    # Make a date list that includes all tracks
    mcs_alldates = pd.date_range(
        start=min_basetime_offset.strftime("%Y-%m-%d"),
        end=max_basetime.dt.strftime("%Y-%m-%d").item(), 
        freq='1D',
    )

    # Select initiation time, and round to the nearest hour
    time0 = rmcs_basetime.isel(times=0).dt.round('H')
    # Get initiation lat/lon    
    rmcs_lon0 = rmcs_lon.isel(times=0).data
    rmcs_lat0 = rmcs_lat.isel(times=0).data

    # Make an array to store the full time series
    nhours_full = nhours + ntimes_max - 1
    full_times = np.ndarray((ntracks, nhours_full), dtype='datetime64[ns]')
    full_lons = np.full((ntracks, nhours_full), np.NaN, dtype=np.float32)
    full_lats = np.full((ntracks, nhours_full), np.NaN, dtype=np.float32)

    # Get MCS track data numpy arrays for better performance
    rmcs_hour0 = time0.data
    rmcs_hours = rmcs_basetime.dt.round('H').data
    rmcs_lons = rmcs_lon.data
    rmcs_lats = rmcs_lat.data

    # Loop over each track
    for itrack in range(0, ntracks):
        # Calculate start/end times prior to initiation
        time0_start = rmcs_hour0[itrack] - pd.offsets.Hour(nhours-1)
        time0_end = rmcs_hour0[itrack] - pd.offsets.Hour(1)
        # Generate hourly time series leading up to -1 h before initiation
        prior_times = np.array(pd.date_range(time0_start, time0_end, freq='1H'))

        # Save full history of times
        full_times[itrack,0:nhours-1] = prior_times
        full_times[itrack,nhours-1:] = rmcs_hours[itrack,:]

        # Repeat initiation lat/lon by X hours (i.e., stay at the initiation location)
        ilon0 = np.repeat(rmcs_lon0[itrack], nhours-1)
        ilat0 = np.repeat(rmcs_lat0[itrack], nhours-1)
        # Save full history of lat/lon
        full_lons[itrack,0:nhours-1] = ilon0
        full_lons[itrack,nhours-1:] = rmcs_lons[itrack,:]
        full_lats[itrack,0:nhours-1] = ilat0
        full_lats[itrack,nhours-1:] = rmcs_lats[itrack,:]

    # Convert to Xarray DataArray
    coord_relativetimes = np.arange(-nhours+1, ntimes_max, 1)
    coord_relativetimes_attrs = {
        'description': 'Relative times for MCS lifecycle',
        'units': 'hour',
    }
    ntimes_full = len(coord_relativetimes)
    
    full_times = xr.DataArray(
        full_times, 
        coords={'tracks':coord_tracks, 'hours':coord_relativetimes}, 
        dims=('tracks','hours'),
    )

    # Convert all MCS times to date strings
    rmcs_dates = full_times.dt.strftime('%Y%m%d')
    # ERA5 data directory is organized by month (yyyymm)
    dirs_month = mcs_alldates.strftime("%Y%m")
    # 3D data is in daily files
    files3d_day = mcs_alldates.strftime("%Y%m%d")
    nfiles = len(mcs_alldates)
    print(f"Total number of ERA5 files: {nfiles}")

    # Create a list to store matchindices for each ERA5 file
    trackindices_all = []
    timeindices_all = []
    results = []

    # Run in serial or parallel
    if run_parallel == 0:

        # Loop over each ERA5 file
        for ifile in range(nfiles):
        # for ifile in range(1):
            idir = dirs_month[ifile]
            iday = files3d_day[ifile]
            filename = f"{era5_dir}{idir}/{basename}{iday}00_{iday}23.nc"
            print(filename)

            # Get all MCS tracks/times indices that are in the same day
            # These tracks use the same ERA5 file
            idx_track, idx_time = np.where(rmcs_dates == files3d_day[ifile])
            # Save track/time indices for the current ERA5 file to the overall list
            trackindices_all.append(idx_track)
            timeindices_all.append(idx_time)

            # Get the track lat/lon/time values
            mcs_lat = full_lats[idx_track, idx_time]
            mcs_lon = full_lons[idx_track, idx_time]
            mcs_time = full_times.values[idx_track, idx_time]

            # Call function to calculate statistics
            result = extract_era5_prof(filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname)
            results.append(result)

        # Serial
        final_result = results

    elif run_parallel == 1:

        # # Set Dask temporary directory for workers
        # dask_tmp_dir = config.get("dask_tmp_dir", "./")
        # dask.config.set({'temporary-directory': dask_tmp_dir})
        # # Initialize dask
        # cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        # client = Client(cluster)

        # Loop over each ERA5 file
        for ifile in range(nfiles):
            idir = dirs_month[ifile]
            iday = files3d_day[ifile]
            filename = f"{era5_dir}{idir}/{basename}{iday}00_{iday}23.nc"
            print(filename)

            # Get all MCS tracks/times indices that are in the same month
            # These tracks use the same ERA5 file
            idx_track, idx_time = np.where(rmcs_dates == files3d_day[ifile])
            # Save track/time indices for the current ERA5 file to the overall list
            trackindices_all.append(idx_track)
            timeindices_all.append(idx_time)

            # Get the track lat/lon/time values
            mcs_lat = full_lats[idx_track, idx_time]
            mcs_lon = full_lons[idx_track, idx_time]
            mcs_time = full_times.data[idx_track, idx_time]

            # Call function to calculate statistics
            result = dask.delayed(extract_era5_prof)(
                filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname
            )
            results.append(result)

        # Trigger dask computation
        final_result = dask.compute(*results)

    VAR_attrs = final_result[0]['VAR_attrs']
    level = final_result[0]['level']
    nlevels = len(level)

    # Create variables for saving output
    fillval = np.NaN
    VAR_out = np.full((ntracks, ntimes_full, nlevels, 2*ny+1, 2*nx+1), fillval, dtype=np.float32)
    # xcoords = np.arange(0, 2*nx+1)
    # ycoords = np.arange(0, 2*ny+1)
    xcoords = np.arange(-nx, nx+1)
    ycoords = np.arange(-ny, ny+1)
    xcoords_attrs = {'long_name':'longitude grids center at MCS', 'units':'0.25deg'}
    ycoords_attrs = {'long_name':'latitude grids center at MCS', 'units':'0.25deg'}

    # Put the results to output track stats variables
    # Loop over each file (parallel return results)
    for ifile in range(nfiles):
    # for ifile in range(1):
        # Get the return results for this pixel file
        iVAR = final_result[ifile]
        if iVAR is not None:
            trackindices = trackindices_all[ifile]
            timeindices = timeindices_all[ifile]
            # Put variable to the MCS track array
            VAR_out[trackindices, timeindices, :, :, :] = iVAR['VAR']

    # Define output dataset
    var_dict = {
        f"{varname}": (['tracks', 'rel_times', 'level', 'y', 'x'], VAR_out, VAR_attrs),
    }
    coord_dict = {
        'tracks': (['tracks'], coord_tracks.data, coord_tracks.attrs),
        'rel_times': (['rel_times'], coord_relativetimes, coord_relativetimes_attrs),
        'level': (['level'], level.data, level.attrs),
        'y': (['y'], ycoords, ycoords_attrs),
        'x': (['x'], xcoords, xcoords_attrs),
    }
    gattr_dict = {
        'Title': 'ERA5 extracted 3D environments for MCS',
        'lon_box_size': nx*2+1,
        'lat_box_size': ny*2+1,
        'Institution': 'Pacific Northwest National Laboratoy',
        'Contact': 'zhe.feng@pnnl.gov',
        'Created_on': time.ctime(time.time()),
    }

    # Define xarray dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    # Update base_time variable dtype as 'double' for better precision
    # bt_dict = {'base_time': {'zlib':True, 'dtype':'float64'}}
    # encoding.update(bt_dict)

    # Write to netcdf file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f'Output saved as: {outfilename}')