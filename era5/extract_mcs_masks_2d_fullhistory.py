"""
Extracts MCS masks centered at each MCS saves the output to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "18-May-2023"

import numpy as np
import xarray as xr
import pandas as pd
# from pandas.tseries.offsets import MonthEnd
import sys, os
import time
import yaml
import dask
from dask.distributed import Client, LocalCluster, wait

#--------------------------------------------------------------------------
def location_to_idx(lat, lon, center):
    """ 
    Convert a latitude and longitude into an index
    
    Args:
        lat: np.array
            Latitude array
        lon: np.array
            Longitufde array
        center: tuple(float)
            location tuple to find (lat, lon)
    
    Returns:
        lat_idx: int
            Index for latitude
        lon_idx: int
            Index for longitude
    """
    # This is for 2D lat/lon
    diff = abs(lat - center[0]) + abs(lon - center[1])
    lat_idx, lon_idx = np.unravel_index(diff.argmin(), diff.shape)
    return lat_idx, lon_idx

#--------------------------------------------------------------------------
def pad_array2d(in_array, lat_idx, lon_idx, ny, nx, ny_d, nx_d):
    """
    Pad 2D array to ny, nx dimensions center at lat_idx, lon_idx.
    
    Args:
        in_array: np.array
            Input 2D array (y, x)
        lat_idx: int
            Center index on latitude (y) dimension
        lon_idx: int
            Center index on longitude (x) dimension
        ny: int
            Number of 1/2 grids to extract data in y dimension
        nx: int
            Number of 1/2 grids to extract data in x dimension
        ny_d: int
            Number of grids in the domain in y dimension
        nx_d: int
            Number of grids in the domain in x dimension

    Returns:
        out_array: np.array
            Padded output array (z, y, x)
    """
    # Constrain lat/lon indices within domain boundary
    iy_min = 0 if (lat_idx-ny < 0) else lat_idx-ny
    iy_max = ny_d if (lat_idx+ny+1 > ny_d) else lat_idx+ny+1
    ix_min = 0 if (lon_idx-nx < 0) else lon_idx-nx
    ix_max = nx_d if (lon_idx+nx+1 > nx_d) else lon_idx+nx+1
    # Number of grids on left, right, bottom, top
    nx_l = lon_idx - ix_min
    nx_r = ix_max - lon_idx - 1
    ny_b = lat_idx - iy_min
    ny_t = iy_max - lat_idx - 1
    # Number of grids to pad on each side
    pnx_l = nx - nx_l
    pnx_r = nx - nx_r
    pny_b = ny - ny_b
    pny_t = ny - ny_t
    # Subset array within the domain
    in_array = in_array[iy_min:iy_max, ix_min:ix_max]
    # Pad array on y & x dimensions
    out_array = np.pad(in_array, ((pny_b,pny_t), (pnx_l,pnx_r)), 'constant', constant_values=np.nan)
    return out_array

#--------------------------------------------------------------------------
def extract_vars(filename, mcs_lat, mcs_lon):

    var_dict = None
    attrs_dict = None

    # Check if file exists
    if os.path.isfile(filename):
        # Read MCS mask
        dsm = xr.open_dataset(filename)
        ny_d = dsm.dims['lat']
        nx_d = dsm.dims['lon']
        lon = dsm['longitude']
        lat = dsm['latitude']
        tb = dsm['tb'].squeeze()
        precipitation = dsm['precipitation'].squeeze()
        cloudtracknumber = dsm['cloudtracknumber'].squeeze()
        cloudnumber = dsm['cloudnumber'].squeeze()

        # Number of tracks in the file
        ntracks_file = len(mcs_lat)
        # Make array to store output
        out_tb = np.full((ntracks_file, 2*ny+1, 2*nx+1), np.NaN, dtype=float)
        out_pcp = np.full((ntracks_file, 2*ny+1, 2*nx+1), np.NaN, dtype=float)
        out_tracknumber = np.full((ntracks_file, 2*ny+1, 2*nx+1), np.NaN, dtype=float)
        out_cloudnumber = np.full((ntracks_file, 2*ny+1, 2*nx+1), np.NaN, dtype=float)

        # Loop over each MCS
        for itrack in range(0, ntracks_file):
            # Find closet lat/lon index to MCS center location
            center = (mcs_lat[itrack], mcs_lon[itrack])
            lat_idx, lon_idx = location_to_idx(lat, lon, center)

            # Extract and pad variables
            _tb = pad_array2d(tb.data, lat_idx, lon_idx, ny, nx, ny_d, nx_d)
            _pcp = pad_array2d(precipitation.data, lat_idx, lon_idx, ny, nx, ny_d, nx_d)
            _tracknumber = pad_array2d(cloudtracknumber.data, lat_idx, lon_idx, ny, nx, ny_d, nx_d)
            _cloudnumber = pad_array2d(cloudnumber.data, lat_idx, lon_idx, ny, nx, ny_d, nx_d)

            iny, inx = _tb.shape
            if (iny == 2*ny+1) & (inx == 2*nx+1):
                out_tb[itrack, :, :] = _tb
                out_pcp[itrack, :, :] = _pcp
                out_tracknumber[itrack, :, :] = _tracknumber
                out_cloudnumber[itrack, :, :] = _cloudnumber

            # Put output variables to a dictionary for easier acceess
            var_dict = {
                'tb': out_tb,
                'precipitation': out_pcp,
                'cloudtracknumber': out_tracknumber,
                'cloudnumber': out_cloudnumber,
            }
            attrs_dict = {
                'tb': tb.attrs,
                'precipitation': precipitation.attrs,
                'cloudtracknumber': cloudtracknumber.attrs,
                'cloudnumber': cloudnumber.attrs,
            }
            # import pdb; pdb.set_trace()

        print(f'Done processing: {filename}')       

    return var_dict, attrs_dict



if __name__ == "__main__":

    config_file = sys.argv[1]
    track_period = sys.argv[2]

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    stats_dir = config['stats_dir']
    pixel_dir = config['pixel_dir'] + track_period + '/'
    output_dir = config['output_dir']
    statsfile_basename = config['statsfile_basename']
    pixelfile_basename = config['pixelfile_basename']
    nx = config['nx']
    ny = config['ny']
    nhours = config['nhours']
    ntimes_max = config['ntimes_max']
    run_parallel = config['run_parallel']
    n_workers = config['n_workers']
    threads_per_worker = config['threads_per_worker']

    stats_file = f"{stats_dir}{statsfile_basename}{track_period}.nc"
    outfilename = f"{output_dir}mcs_2d_mask_{track_period}.nc"
    os.makedirs(output_dir, exist_ok=True)


    # Read robust MCS statistics
    dsm = xr.open_dataset(stats_file)
    # Subset MCS times to reduce array size
    # Most valid MCS data are within 0:ntimes_max
    dsm = dsm.isel(times=slice(0, ntimes_max))

    # dsm = dsm.isel(tracks=slice(1000,1006))

    ntracks = dsm.sizes['tracks']
    rmcs_lat = dsm['meanlat']
    rmcs_lon = dsm['meanlon']

    # # Check if longitude is [-180~+180], if so convert it to [0~360] to match ERA5
    # if np.nanmin(rmcs_lon) < 0:
    #     rmcs_lon = rmcs_lon % 360
    #     print('MCS longitudes are [-180~+180], converted to [0-360] to match ERA5.')

    # Get end times for all tracks
    rmcs_basetime = dsm.base_time
    # Sum over time dimension for valid basetime indices, -1 to get the last valid time index for each track
    # This is the end time index of each track (i.e. +1 equals the lifetime of each track)
    end_time_idx = np.sum(np.isfinite(rmcs_basetime), axis=1)-1
    # Apply fancy indexing to base_time: a tuple that indicates for each track, get the end time index
    end_basetime = rmcs_basetime[(np.arange(0,ntracks), end_time_idx)]
        
    # Get the min/max of all base_times
    min_basetime = rmcs_basetime.sel(times=0).min()
    max_basetime = end_basetime.max()
    # Subtract min base_time by X hours to include pre-initiation times from previous days
    min_basetime_offset = pd.Timestamp(min_basetime.dt.round('min').values) - pd.Timedelta(nhours, unit='h')

    # Make a date/time list that includes all tracks
    mcs_alldates = pd.date_range(
        start=min_basetime_offset.strftime("%Y-%m-%dT%H:%M"),
        end=max_basetime.dt.strftime("%Y-%m-%dT%H:%M").item(), 
        freq='1H',
    )

    # Select initiation time, and round to the nearest minute
    time0 = rmcs_basetime.isel(times=0).dt.round('min')
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
    rmcs_hours = rmcs_basetime.dt.round('min').data
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
    coord_tracks = dsm['tracks']
    full_times = xr.DataArray(
        full_times, 
        coords={'tracks':coord_tracks, 'hours':coord_relativetimes}, 
        dims=('tracks','hours'),
    )

    # Convert all MCS times to date strings
    rmcs_dates = full_times.dt.strftime('%Y%m%d_%H%M')
    # Convert to MCS pixel-level file date/time string format
    files_datetime = mcs_alldates.strftime("%Y%m%d_%H%M")
    nfiles = len(mcs_alldates)
    print(f"Total number of pixel files: {nfiles}")

    # Set up Dask cluster
    if run_parallel >= 1:
        # Set Dask temporary directory for workers
        dask_tmp_dir = config.get("dask_tmp_dir", "/tmp")
        dask.config.set({'temporary-directory': dask_tmp_dir})
        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

    # Create a list to store matchindices for each ERA5 file
    trackindices_all = []
    timeindices_all = []
    results = []

    # Loop over each pixel file
    for ifile in range(nfiles):
    # for ifile in range(1):
        idatetime = files_datetime[ifile]
        filename = f"{pixel_dir}{pixelfile_basename}{idatetime}.nc"
        print(filename)

        # Get all MCS tracks/times indices that are in the same month
        # These tracks use the same ERA5 file
        idx_track, idx_time = np.where(rmcs_dates == idatetime)
        # Save track/time indices for the current ERA5 file to the overall list
        trackindices_all.append(idx_track)
        timeindices_all.append(idx_time)

        # Get the track lat/lon/time values
        mcs_lat = full_lats[idx_track, idx_time]
        mcs_lon = full_lons[idx_track, idx_time]
        # mcs_time = full_times.values[idx_track, idx_time]

        # Call function to calculate statistics
        if run_parallel == 0:
            # Call function to calculate statistics
            result = extract_vars(filename, mcs_lat, mcs_lon)
            results.append(result)
        elif run_parallel >= 1:
            result = dask.delayed(extract_vars)(filename, mcs_lat, mcs_lon)
            results.append(result)
        else:
            print(f'Invalid parallization option run_parallel: {run_parallel}')

    # Trigger dask computation
    if run_parallel == 0:
        final_result = results
    elif run_parallel >= 1:
        final_result = dask.compute(*results)
        wait(final_result)
    else:
        print(f'Invalid parallization option run_parallel: {run_parallel}')


    # Make a variable list and get attributes from one of the returned dictionaries
    # Loop over each return results till one that is not None
    counter = 0
    while counter < nfiles:
        if final_result[counter][0] is not None:
            var_names = list(final_result[counter][0].keys())
            var_attrs = final_result[counter][1]
            break
        counter += 1
    
    # Loop over variable list to create the dictionary entry
    out_dict = {}
    out_dict_attrs = {}
    for ivar in var_names:
        out_dict[ivar] = np.full((ntracks, ntimes_full, 2*ny+1, 2*nx+1), np.NaN, dtype=np.float32)
        out_dict_attrs[ivar] = var_attrs[ivar]

    # Spatial coordinates
    xcoords = np.arange(-nx, nx+1)
    ycoords = np.arange(-ny, ny+1)
    xcoords_attrs = {'long_name':'longitude grids center at MCS', 'units':'0.1deg'}
    ycoords_attrs = {'long_name':'latitude grids center at MCS', 'units':'0.1deg'}

    # Put the results to output track stats variables
    # Loop over each file (parallel return results)
    for ifile in range(nfiles):
    # for ifile in range(1):
        # Get the return results for this pixel file
        iVAR = final_result[ifile][0]
        if iVAR is not None:
            trackindices = trackindices_all[ifile]
            timeindices = timeindices_all[ifile]
            # # Put variable to the MCS track array
            # VAR_out[trackindices, timeindices, :, :] = iVAR['VAR']
            # Loop over each variable and assign values to output dictionary
            for ivar in var_names:
                if iVAR[ivar].ndim == 3:
                    out_dict[ivar][trackindices,timeindices,:,:] = iVAR[ivar]

    # Define a dataset containing all variables
    var_dict = {}
    # Define output variable dictionary
    for key, value in out_dict.items():
        if value.ndim == 4:
            var_dict[key] = (['tracks', 'rel_times', 'y', 'x'], value, out_dict_attrs[key])

    # Define coordinate dictionary
    coord_dict = {
        'tracks': (['tracks'], dsm['tracks'].data, dsm['tracks'].attrs),
        'rel_times': (['rel_times'], coord_relativetimes, coord_relativetimes_attrs),
        'y': (['y'], ycoords, ycoords_attrs),
        'x': (['x'], xcoords, xcoords_attrs),
    }
    # Define global attribute dictionary
    gattr_dict = {
        'Title': 'MCS centered 2D extracted variables',
        'lon_box_size': nx*2+1,
        'lat_box_size': ny*2+1,
        'Institution': 'Pacific Northwest National Laboratory',
        'Contact': 'zhe.feng@pnnl.gov',
        'Created_on': time.ctime(time.time()),
    }

    # Define xarray dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netcdf file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f'Output saved as: {outfilename}')