"""
Calculate mean ERA5 value from a 3D variable centered at each MCS initiation location  
for a time period prior to initiation and saves the output to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "6-Oct-2021"

import numpy as np
import xarray as xr
import pandas as pd
import sys, os
import time
import yaml
import dask
from dask.distributed import Client, LocalCluster

def calc_era5_prof(filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname):

    # Read ERA5 data
    dse = xr.open_dataset(filename)
    lat_e5 = dse.latitude
    lon_e5 = dse.longitude
    time_e5 = dse.time
    dse = dse.sel(level=slice(min(level_lims), max(level_lims)))
    level = dse.level
    # level = dse.level.sel(level=slice(min(level_lims), max(level_lims)))
    VAR_attrs = {'long_name': dse[varname].attrs['long_name'], 
                 'units': dse[varname].attrs['units'],
                }

    # Number of tracks in the day
    ntracks_day = len(mcs_lat)
    nlevels = len(level)
    # Make array to store output
    VAR_avg = np.full((ntracks_day, nlevels), np.NaN, dtype=float)
    VAR_max = np.full((ntracks_day, nlevels), np.NaN, dtype=float)
    VAR_min = np.full((ntracks_day, nlevels), np.NaN, dtype=float)

    # Loop over each MCS
    for itrack in range(0, ntracks_day):
        # print(f'{itrack}: {mcs_time[itrack]}')
        # Find closest grid point and time index in ERA5 
        lat_idx = np.abs(lat_e5.values - mcs_lat[itrack]).argmin()
        lon_idx = np.abs(lon_e5.values - mcs_lon[itrack]).argmin()
        t5_idx = np.abs(time_e5.values - mcs_time[itrack]).argmin()
        
        # Select the time, level, and region
        # Note index+1 is needed to include the last value in the range
        # iVAR = dse[varname].isel(time=t5_idx, level=slice(lev_idx0, lev_idx1+1), \
        iVAR = dse[varname].isel(time=t5_idx, \
                                latitude=slice(lat_idx-ny, lat_idx+ny+1), \
                                longitude=slice(lon_idx-nx, lon_idx+nx+1))
        # Check box size to avoid error (e.g. at date line)
        if (iVAR.sizes['latitude'] > 1) & (iVAR.sizes['longitude'] > 1):
            VAR_avg[itrack,:] = iVAR.mean(dim=('latitude', 'longitude')).values
            VAR_max[itrack,:] = iVAR.max(dim=('latitude', 'longitude')).values
            VAR_min[itrack,:] = iVAR.max(dim=('latitude', 'longitude')).values

    # Put output variables to a dictionary for easier acceess
    out_dict = {'VAR_avg':VAR_avg, 'VAR_max':VAR_max, 'VAR_min':VAR_min, 
                'VAR_attrs':VAR_attrs, 'level':level}
    print(f'Done processing: {filename}')
    
    return out_dict

    
if __name__ == "__main__":

    basename = sys.argv[1]
    varname = sys.argv[2]
    config_file = sys.argv[3]
    track_period = sys.argv[4]

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
    run_parallel = config['run_parallel']
    n_workers = config['n_workers']
    threads_per_worker = config['threads_per_worker']

    mcs_file = f"{stats_dir}{mcsfile_basename}{track_period}.nc"
    outfilename = f"{output_dir}mcs_preinit_era5_{varname}_{track_period}.nc"
    os.makedirs(output_dir, exist_ok=True)

    # Read robust MCS statistics
    dsm = xr.open_dataset(mcs_file)
    ntracks = dsm.sizes['tracks']
    ntimes = dsm.sizes['times']
    rmcs_lat = dsm['meanlat']
    rmcs_lon = dsm['meanlon']

    # Check if longitude is [-180~+180], if so convert it to [0~360] to match ERA5
    if np.nanmin(rmcs_lon) < 0:
        rmcs_lon = rmcs_lon % 360
        print('MCS longitudes are [-180~+180], converted to [0-360] to match ERA5.')

    # Get end times for all tracks
    rmcs_basetime = dsm.base_time
    # Select initiation time, and round to the nearest hour
    time0 = rmcs_basetime.isel(times=0).dt.round('H')

    # Repeat initiation lat/lon by nhours, then reshape the array to [tracks, nhours]
    rmcs_lon0 = rmcs_lon.isel(times=0).values
    rmcs_lat0 = rmcs_lat.isel(times=0).values
    prior_lons = np.reshape(np.repeat(rmcs_lon0, nhours), (ntracks, nhours))
    prior_lats = np.reshape(np.repeat(rmcs_lat0, nhours), (ntracks, nhours))

    # Make an array to store the times
    prior_times = np.ndarray((ntracks, nhours), dtype='datetime64[ns]')
    # Loop over each track
    for itrack in range(0, ntracks):
        # Calculate time X hours prior to initiation
        time0_1day = time0[itrack].values - pd.offsets.Hour(nhours-1)
        # Generate hourly time series leading up to initiation and save to track array
        prior_times[itrack,:] = np.array(pd.date_range(time0_1day, time0[itrack].values, freq='1H'))
    # Convert to Xarray DataArray
    coord_relativetimes = np.arange(-nhours+1, 1)
    coord_tracks = dsm['tracks']
    prior_times = xr.DataArray(prior_times, 
                                coords={'tracks':coord_tracks, 'hours':coord_relativetimes}, 
                                dims=('tracks','hours'))

    # Make a day list that includes all prior times
    min_basetime = prior_times.min().values
    max_basetime = prior_times.max().values
    mcs_alldates = pd.date_range(start=min_basetime, end=max_basetime, freq='1D')

    # Convert all MCS times to date strings
    rmcs_dates = prior_times.dt.strftime('%Y%m%d')
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
            mcs_lat = prior_lats[idx_track, idx_time]
            mcs_lon = prior_lons[idx_track, idx_time]
            mcs_time = prior_times.values[idx_track, idx_time]

            # Call function to calculate statistics
            result = calc_era5_prof(filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname)
            results.append(result)

        # Serial
        final_result = results
    
    elif run_parallel == 1:

        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

        # Loop over each ERA5 file
        for ifile in range(nfiles):
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
            mcs_lat = prior_lats[idx_track, idx_time]
            mcs_lon = prior_lons[idx_track, idx_time]
            mcs_time = prior_times.values[idx_track, idx_time]

            # Call function to calculate statistics
            result = dask.delayed(calc_era5_prof)(filename, mcs_time, mcs_lat, mcs_lon, level_lims, ny, nx, varname)
            results.append(result)

        # Trigger dask computation
        final_result = dask.compute(*results)

    VAR_attrs = final_result[0]['VAR_attrs']
    level = final_result[0]['level']
    nlevels = len(level)

    # Create variables for saving output
    fillval = np.NaN
    VAR_avg = np.full((ntracks, nhours, nlevels), fillval, dtype=float)
    VAR_max = np.full((ntracks, nhours, nlevels), fillval, dtype=float)
    VAR_min = np.full((ntracks, nhours, nlevels), fillval, dtype=float)

    # Put the results to output track stats variables
    # Loop over each file (parallel return results)
    for ifile in range(nfiles):
        # Get the return results for this pixel file
        iVAR = final_result[ifile]
        if iVAR is not None:
            trackindices = trackindices_all[ifile]
            timeindices = timeindices_all[ifile]
            # Put variable to the MCS track array
            VAR_avg[trackindices, timeindices, :] = iVAR['VAR_avg']
            VAR_max[trackindices, timeindices, :] = iVAR['VAR_max']
            VAR_min[trackindices, timeindices, :] = iVAR['VAR_min']

    
    # Define output dataset
    varlist = {
        f"{varname}_avg": (['tracks', 'rel_times', 'level'], VAR_avg, VAR_attrs),
        f"{varname}_max": (['tracks', 'rel_times', 'level'], VAR_max, VAR_attrs),
        f"{varname}_min": (['tracks', 'rel_times', 'level'], VAR_min, VAR_attrs),
    }
    coordlist = {
        'tracks': (['tracks'], dsm['tracks'], dsm['tracks'].attrs),
        'rel_times': (['rel_times'], coord_relativetimes),
        'level': (['level'], level, level.attrs)
    }
    gattrlist = {
        'Title': 'MCS environments from ERA5',
        'lon_box_size': nx*2+1,
        'lat_box_size': ny*2+1,
        'Institution': 'Pacific Northwest National Laboratoy',
        'Contact': 'zhe.feng@pnnl.gov',
        'Created_on': time.ctime(time.time()),
    }

    # Define xarray dataset
    dsout = xr.Dataset(varlist, coords=coordlist, attrs=gattrlist)
    # Update attributes
    dsout['rel_times'].attrs['description'] = 'Relative times prior to MCS initiation'
    dsout['rel_times'].attrs['units'] = 'hour'

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    # Update base_time variable dtype as 'double' for better precision
    # bt_dict = {'base_time': {'zlib':True, 'dtype':'float64'}}
    # encoding.update(bt_dict)

    # Write to netcdf file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f"Output saved as: {outfilename}")

    import pdb; pdb.set_trace()