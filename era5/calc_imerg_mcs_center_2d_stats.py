"""
Calculates domain statistics from extracted MCS center 2D variables and saves to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "19-May-2023"

import os, sys
import time
import numpy as np
import xarray as xr
import dask
from dask.distributed import Client, LocalCluster, wait

#-----------------------------------------------------------------------------------------
def calc_stats_track(in_filename, tracknumber):
    """
    Calculate 2D statistics for a given tracknumber.

    Args:
        in_filename: string
            Input file name.
        tracknumber: int
            Track number to process.

    Returns:
        var_dict: dictionary
            Dictionary containing environmental variables.
        var_attrs: dictionary
            Dictionary containing the attributes of environmental variables.
    """
    print(f'track number: {tracknumber}')

    # Rain rate thresholds to count rain area
    min_rr_thresh = 1.0     # [mm/h]
    heavy_rr_thresh = 10.0  # [mm/h]
    # Approximate pixel area (IMERG: 0.1 degree)
    pixel_area = 100.   # [km^2]

    # Read input data
    # tn = 21141
    # tracknumber = 10005
    # tracknumber = 10000
    ds = xr.open_dataset(in_filename, mask_and_scale=False).sel(tracks=tracknumber).load()
    ntimes = ds.dims['rel_times']
    ny = ds.dims['y']
    nx = ds.dims['x']
    tb = ds['tb']
    precipitation = ds['precipitation']
    cloudtracknumber = ds['cloudtracknumber']
    cloudnumber = ds['cloudnumber']

    # Pixel-file tracknumber need to +1
    px_tracknumber = tracknumber + 1
    # px_tracknumber = tn + 1

    # import matplotlib.pyplot as plt

    # Total domain area (non-missing Tb)
    mask_domain = ~np.isnan(tb)
    # Tracked MCS
    mask_mcs = cloudtracknumber == px_tracknumber
    # Other MCS (excluded the tracked MCS)
    mask_omcs = (cloudtracknumber != px_tracknumber) & (~np.isnan(cloudtracknumber))
    # Non-MCS deep convection
    mask_idc = (cloudnumber > 0) & (np.isnan(cloudtracknumber))

    # Total domain area (number of grid points)
    dim_xy = ('y', 'x')
    npix_domain = (mask_domain == True).sum(dim=dim_xy)

    tb_mcs = tb.where(mask_mcs)
    tb_omcs = tb.where(mask_omcs)
    tb_idc = tb.where(mask_idc)

    pcp_mcs = precipitation.where(mask_mcs)
    pcp_omcs = precipitation.where(mask_omcs)
    pcp_idc = precipitation.where(mask_idc)

    # Tracknumber & cloudnumber
    mcs_tracknumber = cloudtracknumber.where(mask_mcs)
    omcs_tracknumber = cloudtracknumber.where(mask_omcs)
    idc_cloudnumber = cloudnumber.where(mask_idc)

    # Domain mean/max rain rate
    all_pcp_avg = precipitation.where(mask_domain).sum(dim=dim_xy) / npix_domain
    all_pcp_max = precipitation.where(mask_domain).max(dim=dim_xy)

    # Filter weak precipitaiton
    pcp_min = precipitation.where((mask_domain) & (precipitation > min_rr_thresh))
    # Filter non-heavy precipitation
    pcp_heavy = precipitation.where((mask_domain) & (precipitation > heavy_rr_thresh))
    # Total precipitation area (number of pixels)
    pcp_area = (pcp_min > min_rr_thresh).sum(dim=dim_xy) * pixel_area
    pcp_heavy_area = (pcp_heavy > heavy_rr_thresh).sum(dim=dim_xy) * pixel_area

    # Mean rain rate from MCS, other MCS, non-MCS deep convection
    mcs_pcp_avg = pcp_mcs.sum(dim=dim_xy) / npix_domain
    omcs_pcp_avg = pcp_omcs.sum(dim=dim_xy) / npix_domain
    idc_pcp_avg = pcp_idc.sum(dim=dim_xy) / npix_domain

    # Mean/min Tb
    mcs_tb_avg = tb_mcs.mean(dim=dim_xy)
    omcs_tb_avg = tb_omcs.mean(dim=dim_xy)
    idc_tb_avg = tb_idc.mean(dim=dim_xy)
    mcs_tb_min = tb_mcs.min(dim=dim_xy)
    omcs_tb_min = tb_omcs.min(dim=dim_xy)
    idc_tb_min = tb_idc.min(dim=dim_xy)

    # Fractional coverage
    mcs_frac = (mcs_tracknumber > 0).sum(dim=dim_xy) / npix_domain
    omcs_frac = (omcs_tracknumber > 0).sum(dim=dim_xy) / npix_domain
    idc_frac = (idc_cloudnumber > 0).sum(dim=dim_xy) / npix_domain

    # Number of times with valid data
    ntimes_valid = (npix_domain > 0).sum().item()
    # Alternatively, could use: (precipitation.sum(dim=dim_xy) >= 0).sum()
    mcs_number = np.full(ntimes, np.NaN, dtype=np.float32)
    omcs_number = np.full(ntimes, np.NaN, dtype=np.float32)
    idc_number = np.full(ntimes, np.NaN, dtype=np.float32)
    # mcs_area = np.full(ntimes, np.NaN, dtype=np.float32)

    # Loop over each time with valid data
    for tt in range(0, ntimes_valid):
        _mcs_tn = mcs_tracknumber.data[tt,:,:]
        _omcs_tn = omcs_tracknumber.data[tt,:,:]
        _idc_tn = idc_cloudnumber.data[tt,:,:]
        # Count unique values
        _mcs_number, _mcs_counts = np.unique(_mcs_tn, return_counts=True)
        _omcs_number, _omcs_counts = np.unique(_omcs_tn, return_counts=True)
        _idc_number, _idc_counts = np.unique(_idc_tn, return_counts=True)
        # Find indices for not a NaN
        idx_mcs = np.where(~np.isnan(_mcs_number))[0]
        idx_omcs = np.where(~np.isnan(_omcs_number))[0]
        idx_idc = np.where(~np.isnan(_idc_number))[0]
        # Area for each unique clouds
        # _mcs_area = _mcs_counts[idx_mcs]

        # Save the unique number of clouds
        mcs_number[tt] = len(idx_mcs)
        omcs_number[tt] = len(idx_omcs)
        idc_number[tt] = len(idx_idc)


    # Domain mean/max rain rate
    # Total rain area (> 1 mm/h), heavy rain area
    # Mean rain rate from: other MCSs (exclude the tracked one), non-MCS deep convection
    # MCS, other MCS, non-MCS deep convection numbers
    # Mean/min Tb of: MCS, other MCS, non-MCS deep convection
    # Fractional coverage of: all MCS, other MCSs (exclude the tracked one), non-MCS deep convection, shallower clouds/clear-sky

    # Group outputs in dictionaries
    var_dict = {
        'tracknumber': tracknumber,
        'rainrate_mean': all_pcp_avg.data,
        'rainrate_max': all_pcp_max.data,
        'rain_area': pcp_area.data,
        'heavy_rain_area': pcp_heavy_area.data,
        # Rain rate, tb by cloud types
        'mcs_rainrate_mean': mcs_pcp_avg.data,
        'omcs_rainrate_mean': omcs_pcp_avg.data,
        'idc_rainrate_mean': idc_pcp_avg.data,
        'mcs_tb_mean': mcs_tb_avg.data,
        'omcs_tb_mean': omcs_tb_avg.data,
        'idc_tb_mean': idc_tb_avg.data,
        'mcs_tb_min': mcs_tb_min.data,
        'omcs_tb_min': omcs_tb_min.data,
        'idc_tb_min': idc_tb_min.data,
        # Fractional coverage
        'mcs_frac': mcs_frac.data,
        'omcs_frac': omcs_frac.data,
        'idc_frac': idc_frac.data,
        # Number of unique systems
        'mcs_number': mcs_number,
        'omcs_number': omcs_number,
        'idc_number': idc_number,
    }
    # import pdb; pdb.set_trace()
    var_attrs = {
        'tracknumber': {
            'long_name': 'Track number for each MCS'
        },
        'rainrate_mean': {
            'long_name': 'Domain mean rain rate',
            'units': precipitation.attrs['units'],
        },
        'rainrate_max': {
            'long_name': 'Domain max rain rate',
            'units': precipitation.attrs['units'],
        },
        'rain_area': {
            'long_name': 'Domain rain area',
            'units': 'km^2',
            'rr_threshold': min_rr_thresh,
        },
        'heavy_rain_area': {
            'long_name': 'Domain heavy rain area',
            'units': 'km^2',
            'rr_threshold': heavy_rr_thresh,
        },
        # Rain rate, tb by cloud types
        'mcs_rainrate_mean': {
            'long_name': 'MCS domain mean rain rate',
            'units': precipitation.attrs['units'],
        },
        'omcs_rainrate_mean': {
            'long_name': 'Other MCS domain mean rain rate',
            'units': precipitation.attrs['units'],
        },
        'idc_rainrate_mean': {
            'long_name': 'Non-MCS deep convection domain mean rain rate',
            'units': precipitation.attrs['units'],
        },
        'mcs_tb_mean': {
            'long_name': 'MCS mean Tb (in-cloud)',
            'units': tb.attrs['units'],
        },
        'omcs_tb_mean': {
            'long_name': 'Other MCS mean Tb (in-cloud)',
            'units': tb.attrs['units'],
        },
        'idc_tb_mean': {
            'long_name': 'Non-MCS deep convection mean Tb (in-cloud)',
            'units': tb.attrs['units'],
        },
        'mcs_tb_min': {
            'long_name': 'MCS min Tb',
            'units': tb.attrs['units'],
        },
        'omcs_tb_min': {
            'long_name': 'Other MCS min Tb',
            'units': tb.attrs['units'],
        },
        'idc_tb_min': {
            'long_name': 'Non-MCS deep convection min Tb',
            'units': tb.attrs['units'],
        },
        # Fractional coverage
        'mcs_frac': {
            'long_name': 'MCS cloud fraction',
            'units': 'fraction',
        },
        'omcs_frac': {
            'long_name': 'Other MCS cloud fraction',
            'units': 'fraction',
        },
        'idc_frac': {
            'long_name': 'Non-MCS deep convection cloud fraction',
            'units': 'fraction',
        },
        # Number of unique systems
        'mcs_number': {
            'long_name': 'Number of tracked MCS in the domain',
            'units': 'count',
        },
        'omcs_number': {
            'long_name': 'Number of other MCS in the domain',
            'units': 'count',
        },
        'idc_number': {
            'long_name': 'Number of non-MCS deep convection in the domain',
            'units': 'count',
        },
    }
    print(f'Done with track: {tracknumber}')
    # import pdb; pdb.set_trace()
    return var_dict, var_attrs

#-----------------------------------------------------------------------------------------
def work_for_tracks(in_filename, out_filename):
    """
    Calculate 2D statistics for all tracks and saves output to a netCDF file.

    Args:
        in_filename: string
            Input file name.
        out_filename: string
            Output file name.

    Returns:
        out_filename: string
            Output file name.
    """

    # Read input data
    ds = xr.open_dataset(in_filename)    
    # ds = xr.open_dataset(in_filename).sel(tracks=slice(10000,10030))

    tracks = ds['tracks']
    ntracks = ds.dims['tracks']
    ntimes = ds.dims['rel_times']
    ny = ds.dims['y']
    nx = ds.dims['x']

    # Set up Dask cluster
    if run_parallel >= 1:
        # Set Dask temporary directory for workers
        dask_tmp_dir = "/tmp"
        dask.config.set({'temporary-directory': dask_tmp_dir})
        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    results = []

    # Loop over tracks
    for itrack in range(0, ntracks):
        tracknumber = tracks.data[itrack]
        # Serial
        if run_parallel == 0:
            result = calc_stats_track(in_filename, tracknumber)
            results.append(result)
        # Parallel
        elif run_parallel >= 1:
            result = dask.delayed(calc_stats_track)(in_filename, tracknumber)
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

    # import pdb; pdb.set_trace()

    # Make a variable list from one of the returned dictionaries
    var_names = list(final_result[0][0].keys())
    # Get variable attributes from one of the returned dictionaries
    var_attrs = final_result[0][1]

    # Remove tracknumbers from the list
    var_names.remove('tracknumber')
    var_attrs.pop('tracknumber', None)

    # Loop over variable list to create the dictionary entry
    out_dict = {}
    out_dict_attrs = {}
    var2d_dims = (ntracks, ntimes)
    for ivar in var_names:
        out_dict[ivar] = np.full(var2d_dims, np.nan, dtype=np.float32)
        out_dict_attrs[ivar] = var_attrs[ivar]

    # Collect results
    for itrack in range(0, ntracks):
        if final_result[itrack] is not None:
            # Get the return results for this track
            # The result is a tuple: (out_dict, out_dict_attrs)
            # The first entry is the dictionary containing the variables
            iResult = final_result[itrack][0]
            tracknumber = iResult['tracknumber']

            # Double check tracknumber from return to make sure it matches the track
            if tracks.data[itrack] == tracknumber:
                # Loop over each variable and assign values to output dictionary
                for ivar in var_names:
                    if iResult[ivar].ndim == 1:
                        out_dict[ivar][itrack, :] = iResult[ivar]
                        # import pdb; pdb.set_trace()
            else:
                print(f'ERROR: tracknumber does not match: {tracknumber}!')
                sys.exit('Double check results!')

    # Define a dataset containing all PF variables
    var_dict = {}
    # Define output variable dictionary
    for key, value in out_dict.items():
        if value.ndim == 2:
            var_dict[key] = (['tracks', 'rel_times'], value, out_dict_attrs[key])
    coord_dict = {
        'tracks': (['tracks'], tracks.data, tracks.attrs),
        'rel_times': (['rel_times'], ds['rel_times'].data, ds['rel_times'].attrs),
    }
    # Define global attribute dictionary
    gattr_dict = {
        'Title': 'MCS centered 2D variable statistics',
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
    dsout.to_netcdf(path=out_filename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f'Output saved as: {out_filename}')
    # import pdb; pdb.set_trace()

    return out_filename


if __name__ == "__main__":

    # Get track string from input
    # track_string = sys.argv[1]
    filename = sys.argv[1]

    # Get track_string from input filename
    # filename format: mcs_2d_mask_20190101.0000_20200101.0000_txxxxx.nc
    in_basename = 'mcs_2d_mask_'
    # Get file basename
    fn = os.path.basename(filename)
    track_string = fn[len(in_basename):-3]

    # Parallel setup
    run_parallel = 1
    n_workers = 64

    # Get sub-strings
    # track_string format: 20190101.0000_20200101.0000_txxxxx
    track_period = track_string[:-7]
    track_year = track_string[0:4]

    # Input data directory
    in_dir = '/pscratch/sd/f/feng045/waccem/mcs_global/mcs_center_2d/parts/'
    # Output file
    out_dir = '/pscratch/sd/f/feng045/waccem/mcs_global/mcs_center_2d/parts/'
    out_basename = 'mcs_2d_stats_'
    out_filename = f'{out_dir}{out_basename}{track_string}.nc'
    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # Input file
    # in_basename = 'mcs_2d_mask_'
    in_filename = f'{in_dir}{filename}'
    # in_filename = f'{in_dir}{in_basename}{track_string}.nc'

    # Call function to calculate
    result = work_for_tracks(in_filename, out_filename)