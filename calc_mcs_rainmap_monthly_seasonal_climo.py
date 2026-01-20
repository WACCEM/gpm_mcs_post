"""
Script to calculate multi-year seasonal and monthly mean MCS precipitation statistics from monthly rainmap data.

This script processes mcs_rainmap_{year}{month}.nc files to calculate climatological means.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "17-Jan-2026"

import numpy as np
import glob, os, sys
import xarray as xr
import time
import yaml

if __name__ == "__main__":

    config_file = sys.argv[1]

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    datadir = config['output_monthly_dir']
    climo_years = config['climo_years']
    output_climo_dir = config['output_climo_dir']

    pcpfiles = []
    for iyear in climo_years:
        pcpfiles.extend(sorted(glob.glob(f'{datadir}mcs_rainmap_{iyear}*nc')))
    print(f'Number of rainmap files: {len(pcpfiles)}')

    os.makedirs(output_climo_dir, exist_ok=True)

    # Read data
    ds = xr.open_mfdataset(pcpfiles, concat_dim='time', combine='nested')
    nx = ds.sizes['lon']
    ny = ds.sizes['lat']
    lon = ds.lon
    lat = ds.lat

    # Get number of seasons/years
    years = ds.time.dt.year
    nyears = len(np.unique(years))
    print(f'Number of unique years: {nyears}')

    # Get min/max year and create output file names
    min_year = ds.time.min().dt.strftime("%Y").values
    max_year = ds.time.max().dt.strftime("%Y").values
    outfile_season = f'{output_climo_dir}mcs_rainmap_seasonal_climo_{min_year}_{max_year}.nc'
    outfile_month = f'{output_climo_dir}mcs_rainmap_monthly_climo_{min_year}_{max_year}.nc'

    ############################################################################
    # Calculate monthly climatology
    ############################################################################
    
    # Calculate total number of hours for each month
    nhours_month = ds.ntimes.groupby('time.month').sum(dim='time')
    
    # Calculate monthly mean precipitation amounts
    # Convert from mm to mm/day: (total_mm / total_hours * 24 = mm/day)
    precipitation_month = ds.precipitation.groupby('time.month').sum(dim='time').load() / nhours_month * 24
    mcs_precipitation_month = ds.mcs_precipitation.groupby('time.month').sum(dim='time').load() / nhours_month * 24
    
    # Calculate total hours with MCS for each month
    mcs_precipitation_count_month = ds.mcs_precipitation_count.groupby('time.month').sum(dim='time')
    mcs_cloud_count_month = ds.mcs_cloud_count.groupby('time.month').sum(dim='time')
    
    # Calculate MCS precipitation frequency (% of time with MCS precipitation/cloud)
    mcs_precipitation_freq_month = 100. * mcs_precipitation_count_month / nhours_month
    mcs_cloud_freq_month = 100. * mcs_cloud_count_month / nhours_month
    
    # Calculate monthly mean counts (divide by number of years)
    mcs_precipitation_count_month = mcs_precipitation_count_month / nyears
    mcs_cloud_count_month = mcs_cloud_count_month / nyears
    
    # Calculate MCS precipitation fraction (% of total precip from MCS)
    mcs_precipitation_frac_month = 100. * mcs_precipitation_month / precipitation_month
    
    # Calculate mean MCS precipitation intensity (when MCS present)
    mcs_precipitation_intensity_month = mcs_precipitation_month / (mcs_precipitation_count_month / nyears)

    ############################################################################
    # Calculate seasonal climatology
    ############################################################################
    
    # Calculate total number of hours for each season
    nhours_season = ds.ntimes.groupby('time.season').sum(dim='time')
    
    # Calculate seasonal mean precipitation amounts
    # Convert from mm to mm/day: (total_mm / total_hours * 24 = mm/day)
    precipitation_season = ds.precipitation.groupby('time.season').sum(dim='time').load() / nhours_season * 24
    mcs_precipitation_season = ds.mcs_precipitation.groupby('time.season').sum(dim='time').load() / nhours_season * 24
    
    # Calculate total hours with MCS for each season
    mcs_precipitation_count_season = ds.mcs_precipitation_count.groupby('time.season').sum(dim='time')
    mcs_cloud_count_season = ds.mcs_cloud_count.groupby('time.season').sum(dim='time')
    
    # Calculate MCS precipitation frequency (% of time with MCS precipitation/cloud)
    mcs_precipitation_freq_season = 100. * mcs_precipitation_count_season / nhours_season
    mcs_cloud_freq_season = 100. * mcs_cloud_count_season / nhours_season
    
    # Calculate seasonal mean counts (divide by number of years)
    mcs_precipitation_count_season = mcs_precipitation_count_season / nyears
    mcs_cloud_count_season = mcs_cloud_count_season / nyears
    
    # Calculate MCS precipitation fraction (% of total precip from MCS)
    mcs_precipitation_frac_season = 100. * mcs_precipitation_season / precipitation_season
    
    # Calculate mean MCS precipitation intensity (when MCS present)
    mcs_precipitation_intensity_season = mcs_precipitation_season / (mcs_precipitation_count_season / nyears)

    ############################################################################
    # Write seasonal climatology output file
    ############################################################################
    var_dict = {
        'precipitation': (['season', 'lat', 'lon'], precipitation_season.data),
        'mcs_precipitation': (['season', 'lat', 'lon'], mcs_precipitation_season.data),
        'mcs_precipitation_frac': (['season', 'lat', 'lon'], mcs_precipitation_frac_season.data),
        'mcs_precipitation_count': (['season', 'lat', 'lon'], mcs_precipitation_count_season.data),
        'mcs_precipitation_intensity': (['season', 'lat', 'lon'], mcs_precipitation_intensity_season.data),
        'mcs_cloud_count': (['season', 'lat', 'lon'], mcs_cloud_count_season.data),
        'mcs_freq_ccs': (['season', 'lat', 'lon'], mcs_cloud_freq_season.data),
        'mcs_freq_pf': (['season', 'lat', 'lon'], mcs_precipitation_freq_season.data),

    }
    coord_dict = {
        'season': (['season'], precipitation_season.season.data),
        'lat': (['lat'], ds.lat.data),
        'lon': (['lon'], ds.lon.data),
    }
    gattr_dict = {
        'title': 'MCS precipitation climatology',
        'contact':'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on':time.ctime(time.time()),
        'produced_by': __file__,
    }
    dsmap = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    dsmap.season.attrs['long_name'] = 'Season'
    dsmap.season.attrs['units'] = 'none'
    dsmap.lon.attrs['long_name'] = 'Longitude'
    dsmap.lon.attrs['units'] = 'degree'
    dsmap.lat.attrs['long_name'] = 'Latitude'
    dsmap.lat.attrs['units'] = 'degree'
    dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
    dsmap.precipitation.attrs['units'] = 'mm/day'
    dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
    dsmap.mcs_precipitation.attrs['units'] = 'mm/day'
    dsmap.mcs_precipitation_frac.attrs['long_name'] = 'MCS precipitation fraction'
    dsmap.mcs_precipitation_frac.attrs['units'] = '%'
    dsmap.mcs_precipitation_count.attrs['long_name'] = 'Number of hours MCS precipitation is recorded'
    dsmap.mcs_precipitation_count.attrs['units'] = 'hour'
    dsmap.mcs_precipitation_intensity.attrs['long_name'] = 'MCS precipitation intensity'
    dsmap.mcs_precipitation_intensity.attrs['units'] = 'mm/hour'
    dsmap.mcs_cloud_count.attrs['long_name'] = 'Number of hours MCS cloud is recorded'
    dsmap.mcs_cloud_count.attrs['units'] = 'hour'
    dsmap.mcs_freq_ccs.attrs['long_name'] = 'MCS cloud frequency'
    dsmap.mcs_freq_ccs.attrs['units'] = '%'
    dsmap.mcs_freq_pf.attrs['long_name'] = 'MCS precipitation frequency'
    dsmap.mcs_freq_pf.attrs['units'] = '%'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsmap.data_vars}

    dsmap.to_netcdf(path=outfile_season, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
    print('Seasonal climatology output saved as: ', outfile_season)

    ############################################################################
    # Write monthly climatology output file
    ############################################################################
    var_dict = {
        'precipitation': (['month', 'lat', 'lon'], precipitation_month.data),
        'mcs_precipitation': (['month', 'lat', 'lon'], mcs_precipitation_month.data),
        'mcs_precipitation_frac': (['month', 'lat', 'lon'], mcs_precipitation_frac_month.data),
        'mcs_precipitation_count': (['month', 'lat', 'lon'], mcs_precipitation_count_month.data),
        'mcs_precipitation_intensity': (['month', 'lat', 'lon'], mcs_precipitation_intensity_month.data),
        'mcs_cloud_count': (['month', 'lat', 'lon'], mcs_cloud_count_month.data),
        'mcs_freq_ccs': (['month', 'lat', 'lon'], mcs_cloud_freq_month.data),
        'mcs_freq_pf': (['month', 'lat', 'lon'], mcs_precipitation_freq_month.data),

    }
    coord_dict = {
        'month': (['month'], precipitation_month.month.data),
        'lat': (['lat'], ds.lat.data),
        'lon': (['lon'], ds.lon.data),
    }
    gattr_dict = {
        'title': 'MCS precipitation climatology',
        'contact':'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on':time.ctime(time.time()),
        'produced_by': __file__,
    }
    dsmap = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    dsmap.month.attrs['long_name'] = 'Month'
    dsmap.month.attrs['units'] = 'none'
    dsmap.lon.attrs['long_name'] = 'Longitude'
    dsmap.lon.attrs['units'] = 'degree'
    dsmap.lat.attrs['long_name'] = 'Latitude'
    dsmap.lat.attrs['units'] = 'degree'
    dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
    dsmap.precipitation.attrs['units'] = 'mm/day'
    dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
    dsmap.mcs_precipitation.attrs['units'] = 'mm/day'
    dsmap.mcs_precipitation_frac.attrs['long_name'] = 'MCS precipitation fraction'
    dsmap.mcs_precipitation_frac.attrs['units'] = '%'
    dsmap.mcs_precipitation_count.attrs['long_name'] = 'Number of hours MCS precipitation is recorded'
    dsmap.mcs_precipitation_count.attrs['units'] = 'hour'
    dsmap.mcs_precipitation_intensity.attrs['long_name'] = 'MCS precipitation intensity'
    dsmap.mcs_precipitation_intensity.attrs['units'] = 'mm/hour'
    dsmap.mcs_cloud_count.attrs['long_name'] = 'Number of hours MCS cloud is recorded'
    dsmap.mcs_cloud_count.attrs['units'] = 'hour'
    dsmap.mcs_freq_ccs.attrs['long_name'] = 'MCS cloud frequency'
    dsmap.mcs_freq_ccs.attrs['units'] = '%'
    dsmap.mcs_freq_pf.attrs['long_name'] = 'MCS precipitation frequency'
    dsmap.mcs_freq_pf.attrs['units'] = '%'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsmap.data_vars}

    dsmap.to_netcdf(path=outfile_month, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
    print('Monthly climatology output saved as: ', outfile_month)
