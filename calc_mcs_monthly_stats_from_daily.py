"""
Calculates monthly mean MCS statistics from daily mean files and saves output to a netCDF file.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
05/01/2020 - Written.
"""

import numpy as np
import glob, sys, os
import xarray as xr
import time, datetime, calendar, pytz
import yaml

if __name__ == "__main__":

    year = sys.argv[1]
    month = sys.argv[2]
    config_file = sys.argv[3]
    # region = sys.argv[3]
    # year = '2014'
    # month = '01'

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    daily_dir = config['output_daily_dir']
    monthly_dir = config['output_monthly_dir']

    datafiles = sorted(glob.glob(f'{daily_dir}mcs_statsmap_{year}{month}*nc'))
    nfiles = len(datafiles)
    print(year, month)
    print(f'Number of files: {nfiles}')

    outfile = f'{monthly_dir}mcs_statsmap_{year}{month}.nc'
    os.makedirs(monthly_dir, exist_ok=True)

    # Read all data files
    ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested')

    # Simple counts
    mcs_number_ccs = ds.mcs_number_ccs.sum(dim='time')
    mcs_number_pf = ds.mcs_number_pf.sum(dim='time')
    mcs_nhour_ccs = ds.mcs_nhour_ccs.sum(dim='time')
    mcs_nhour_pf = ds.mcs_nhour_pf.sum(dim='time')
    mcs_nhour_speedmcs = ds.mcs_nhour_speedmcs.sum(dim='time')
    initiation_ccs = ds.initiation_ccs.sum(dim='time')

    # Simple averages
    # lifetime_mean = ds.lifetime_mean.mean(dim='time')
    # ccs_area_mean = ds.ccs_area_mean.mean(dim='time')
    # pf_area_mean = ds.pf_area_mean.mean(dim='time')

    # Multiple daily mean by numher of hours per day to get back the daily sum, 
    # then sum over all days, divide by the total number of hours of the month to get monthly mean
    lifetime_mean = (ds.lifetime_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    pf_area_mean = (ds.pf_area_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    ccs_area_mean = (ds.ccs_area_mean * ds.mcs_nhour_ccs).sum(dim='time') / mcs_nhour_ccs

    totalrain_mean = (ds.totalrain_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    totalrainheavy_mean = (ds.totalrainheavy_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    rainrateheavy_mean = (ds.rainrateheavy_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    rainratemax_mean = (ds.rainratemax_mean * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf

    # Propagation Speed
    pf_speed = (ds.pf_speed * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    pf_uspeed = (ds.pf_uspeed * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf
    pf_vspeed = (ds.pf_vspeed * ds.mcs_nhour_pf).sum(dim='time') / mcs_nhour_pf

    mcs_nhour_speedmcs = ds.mcs_nhour_speedmcs.sum(dim='time')
    pf_speed_mcs = (ds.pf_speed_mcs * ds.mcs_nhour_speedmcs).sum(dim='time') / mcs_nhour_speedmcs
    pf_uspeed_mcs = (ds.pf_uspeed_mcs * ds.mcs_nhour_speedmcs).sum(dim='time') / mcs_nhour_speedmcs
    pf_vspeed_mcs = (ds.pf_vspeed_mcs * ds.mcs_nhour_speedmcs).sum(dim='time') / mcs_nhour_speedmcs

    # Compute Epoch Time for the month
    months = np.zeros(1, dtype=int)
    months[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

    ############################################################################
    # Write output file
    var_dict = {
        'mcs_number_ccs': (['time', 'lat', 'lon'], np.expand_dims(mcs_number_ccs, 0)),
        'mcs_number_pf': (['time', 'lat', 'lon'], np.expand_dims(mcs_number_pf, 0)),
        'mcs_nhour_ccs': (['time', 'lat', 'lon'], np.expand_dims(mcs_nhour_ccs, 0)),
        'mcs_nhour_pf': (['time', 'lat', 'lon'], np.expand_dims(mcs_nhour_pf, 0)),
        'mcs_nhour_speedmcs': (['time', 'lat', 'lon'], np.expand_dims(mcs_nhour_speedmcs, 0)),
        'lifetime_mean': (['time', 'lat', 'lon'], np.expand_dims(lifetime_mean, 0)),
        'ccs_area_mean': (['time', 'lat', 'lon'], np.expand_dims(ccs_area_mean, 0)), 
        'pf_area_mean': (['time', 'lat', 'lon'], np.expand_dims(pf_area_mean, 0)),
        'totalrain_mean': (['time', 'lat', 'lon'], np.expand_dims(totalrain_mean, 0)),
        'totalrainheavy_mean': (['time', 'lat', 'lon'], np.expand_dims(totalrainheavy_mean, 0)),
        'rainrateheavy_mean': (['time', 'lat', 'lon'], np.expand_dims(rainrateheavy_mean, 0)),
        'rainratemax_mean': (['time', 'lat', 'lon'], np.expand_dims(rainratemax_mean, 0)),
        'initiation_ccs': (['time', 'lat', 'lon'], np.expand_dims(initiation_ccs, 0)),
        'pf_speed': (['time', 'lat', 'lon'], np.expand_dims(pf_speed, 0)),
        'pf_speed_mcs': (['time', 'lat', 'lon'], np.expand_dims(pf_speed_mcs, 0)),
        'pf_uspeed': (['time', 'lat', 'lon'], np.expand_dims(pf_uspeed, 0)),
        'pf_uspeed_mcs': (['time', 'lat', 'lon'], np.expand_dims(pf_uspeed_mcs, 0)),
        'pf_vspeed': (['time', 'lat', 'lon'], np.expand_dims(pf_vspeed, 0)),
        'pf_vspeed_mcs': (['time', 'lat', 'lon'], np.expand_dims(pf_vspeed_mcs, 0)),
    }
    coord_dict = {
        'lon': (['lon'], ds.lon.data),
        'lat': (['lat'], ds.lat.data),
        'time': (['time'], months),
    }
    gattr_dict = {
        'title': 'MCS monthly statistics map',
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on': time.ctime(time.time()),
    }
    # Define Xarray Dataset
    dsmap = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Define variable attributes
    dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
    dsmap.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'
    dsmap.lat.attrs['long_name'] = 'Latitude'
    dsmap.lat.attrs['units'] = 'degree'
    dsmap.lon.attrs['long_name'] = 'Longitude'
    dsmap.lon.attrs['units'] = 'degree'
    dsmap.mcs_number_ccs.attrs['long_name'] = 'Number of MCS defined by cold cloud shield'
    dsmap.mcs_number_ccs.attrs['units'] = 'count'
    dsmap.mcs_number_pf.attrs['long_name'] = 'Number of MCS defined by precipitation feature'
    dsmap.mcs_number_pf.attrs['units'] = 'count'
    dsmap.mcs_nhour_ccs.attrs['long_name'] = 'Number of MCS hours defined by cold cloud shield'
    dsmap.mcs_nhour_ccs.attrs['units'] = 'hour'
    dsmap.mcs_nhour_pf.attrs['long_name'] = 'Number of MCS hours defined by precipitation feature'
    dsmap.mcs_nhour_pf.attrs['units'] = 'hour'
    dsmap.mcs_nhour_speedmcs.attrs['long_name'] = 'Number of MCS hours for propagation speed (MCS status is met)'
    dsmap.mcs_nhour_speedmcs.attrs['units'] = 'hour'
    dsmap.lifetime_mean.attrs['long_name'] = 'MCS lifetime mean'
    dsmap.lifetime_mean.attrs['units'] = 'hour'
    dsmap.ccs_area_mean.attrs['long_name'] = 'MCS cold cloud shield area mean'
    dsmap.ccs_area_mean.attrs['units'] = 'km2'
    dsmap.pf_area_mean.attrs['long_name'] = 'MCS precipitation feature area mean'
    dsmap.pf_area_mean.attrs['units'] = 'km2'
    dsmap.totalrain_mean.attrs['long_name'] = 'MCS total precipitation mean'
    dsmap.totalrain_mean.attrs['units'] = 'mm'
    dsmap.totalrainheavy_mean.attrs['long_name'] = 'MCS total heavy precipitation (rain rate > 10 mm/h) mean'
    dsmap.totalrainheavy_mean.attrs['units'] = 'mm'
    dsmap.rainrateheavy_mean.attrs['long_name'] = 'MCS heavy rain rate (rain rate > 10 mm/h) mean'
    dsmap.rainrateheavy_mean.attrs['units'] = 'mm/h'
    dsmap.rainratemax_mean.attrs['long_name'] = 'MCS max rain rate mean'
    dsmap.rainratemax_mean.attrs['units'] = 'mm/h'
    dsmap.initiation_ccs.attrs['long_name'] = 'MCS convective initiation hours defined by cold cloud shield'
    dsmap.initiation_ccs.attrs['units'] = 'hour'
    dsmap.initiation_ccs.attrs['long_name'] = 'MCS convective initiation hours defined by cold cloud shield'
    dsmap.initiation_ccs.attrs['units'] = 'hour'
    dsmap.pf_speed.attrs['long_name'] = 'Propagation speed for all times'
    dsmap.pf_speed.attrs['units'] = 'm/s'
    dsmap.pf_uspeed.attrs['long_name'] = 'Propagation speed (x-direction) for all times'
    dsmap.pf_uspeed.attrs['units'] = 'm/s'
    dsmap.pf_uspeed_mcs.attrs['long_name'] = 'Propagation speed (x-direction) when MCS status is met'
    dsmap.pf_uspeed_mcs.attrs['units'] = 'm/s'
    dsmap.pf_vspeed.attrs['long_name'] = 'Propagation speed (y-direction) for all times'
    dsmap.pf_vspeed.attrs['units'] = 'm/s'
    dsmap.pf_vspeed_mcs.attrs['long_name'] = 'Propagation speed (y-direction) when MCS status is met'
    dsmap.pf_vspeed_mcs.attrs['units'] = 'm/s'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsmap.data_vars}

    # Write to netCDF file
    dsmap.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)

    print('Map output saved as: ', outfile)

