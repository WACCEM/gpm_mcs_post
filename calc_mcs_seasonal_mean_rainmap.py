"""
Script to calculate multi- seasonal and monthly mean MCS precipitaiton from monthly data.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "12-Jun-2019"

import numpy as np
import glob, os, sys
import xarray as xr
import pandas as pd
import time, datetime, calendar, pytz

region = sys.argv[1]
#region = 'maritime_continent'

#mcsdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/monthly/'
mcsdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_pt1/monthly/'
mcsfiles = sorted(glob.glob(f'{mcsdir}mcs_rainmap_20????*nc'))
print(f'Number of files: {len(mcsfiles)}')

#outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/climo/'
outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_pt1/climo/'
os.makedirs(outdir, exist_ok=True)
#outfile_season = f'{outdir}mcs_rainmap_seasonal_mean.nc'
#outfile_month = f'{outdir}mcs_rainmap_monthly_mean.nc'

# Read data
ds = xr.open_mfdataset(mcsfiles, concat_dim='time', combine='nested')
nx = ds.dims['lon']
ny = ds.dims['lat']
lon = ds.lon
lat = ds.lat


# Get number of seasons/years
years = pd.date_range(start=ds.time.min().values, end=ds.time.max().values, freq='AS')
nyears = len(years)
#print(nyears)

# Get min/max year and create output file names
min_year = ds.time.min().dt.strftime("%Y").values
max_year = ds.time.max().dt.strftime("%Y").values
outfile_season = f'{outdir}mcs_rainmap_seasonal_mean_{min_year}_{max_year}.nc'
outfile_month = f'{outdir}mcs_rainmap_monthly_mean_{min_year}_{max_year}.nc'

# Calculate 4 season mean
ntimes_season = ds.ntimes.groupby('time.season').sum()

totpcp_season = 24 * ds.precipitation.groupby('time.season').sum(dim='time') / ntimes_season
mcspcp_season = 24 * ds.mcs_precipitation.groupby('time.season').sum(dim='time') / ntimes_season
mcsfrac_season = 100 * mcspcp_season / totpcp_season

mcspcpfreq_season = 100 * ds.mcs_precipitation_count.groupby('time.season').sum(dim='time') / ntimes_season

# Number of hours for MCS precipitation
mcspcpcount_season = ds.mcs_precipitation_count.groupby('time.season').sum(dim='time')
# Mean MCS precipitation intensity
mcspcpintensity_season = ds.mcs_precipitation.groupby('time.season').sum(dim='time') / mcspcpcount_season


# Calculate monthly mean
ntimes_month = ds.ntimes.groupby('time.month').sum()

totpcp_month = 24 * ds.precipitation.groupby('time.month').sum(dim='time') / ntimes_month
mcspcp_month = 24 * ds.mcs_precipitation.groupby('time.month').sum(dim='time') / ntimes_month
mcsfrac_month = 100 * mcspcp_month / totpcp_month

mcspcpfreq_month = 100 * ds.mcs_precipitation_count.groupby('time.month').sum(dim='time') / ntimes_month

# Number of hours for MCS precipitation
mcspcpcount_month = ds.mcs_precipitation_count.groupby('time.month').sum(dim='time')
# Mean MCS precipitation intensity
mcspcpintensity_month = ds.mcs_precipitation.groupby('time.month').sum(dim='time') / mcspcpcount_month




# Write seasonal mean output file
dsmap = xr.Dataset({'precipitation': (['season', 'lat', 'lon'], totpcp_season), \
                    'mcs_precipitation': (['season', 'lat', 'lon'], mcspcp_season), \
                    'mcs_precipitation_frac': (['season', 'lat', 'lon'], mcsfrac_season), \
                    'mcs_precipitation_freq': (['season', 'lat', 'lon'], mcspcpfreq_season), \
                    'mcs_precipitation_intensity': (['season', 'lat', 'lon'], mcspcpintensity_season), \
                    }, \
                    coords={'season': (['season'], totpcp_season.season), \
                            'lat': (['lat'], ds.lat), \
                            'lon': (['lon'], ds.lon)}, \
                attrs={'title': 'MCS precipitation statistics', \
                       'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                       'created_on':time.ctime(time.time())})

dsmap.season.attrs['long_name'] = 'Season'
dsmap.season.attrs['units'] = 'none'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

# dsmap.ntimes.attrs['long_name'] = 'Number of hours in the month'
# dsmap.ntimes.attrs['units'] = 'count'

dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
dsmap.precipitation.attrs['units'] = 'mm/day'

dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
dsmap.mcs_precipitation.attrs['units'] = 'mm/day'

dsmap.mcs_precipitation_frac.attrs['long_name'] = 'MCS precipitation fraction'
dsmap.mcs_precipitation_frac.attrs['units'] = '%'

dsmap.mcs_precipitation_freq.attrs['long_name'] = 'MCS precipitation frequency'
dsmap.mcs_precipitation_freq.attrs['units'] = '%'

dsmap.mcs_precipitation_intensity.attrs['long_name'] = 'MCS precipitation intensity'
dsmap.mcs_precipitation_intensity.attrs['units'] = 'mm/hour'

fillvalue = np.nan
# Set encoding/compression for all variables
comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
encoding = {var: comp for var in dsmap.data_vars}

dsmap.to_netcdf(path=outfile_season, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)

#fillvalue = np.nan
#dsmap.to_netcdf(path=outfile_season, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='season', \
#                encoding={'lon':{'zlib':True, 'dtype':'float32'}, \
#                          'lat':{'zlib':True, 'dtype':'float32'}, \
#                          'precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation_frac':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation_freq':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation_intensity':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          })
print('Seasonal output saved as: ', outfile_season)


# Write monthly mean output file
dsmap = xr.Dataset({'precipitation': (['month', 'lat', 'lon'], totpcp_month), \
                    'mcs_precipitation': (['month', 'lat', 'lon'], mcspcp_month), \
                    'mcs_precipitation_frac': (['month', 'lat', 'lon'], mcsfrac_month), \
                    'mcs_precipitation_freq': (['month', 'lat', 'lon'], mcspcpfreq_month), \
                    'mcs_precipitation_intensity': (['season', 'lat', 'lon'], mcspcpintensity_month), \
                    }, \
                    coords={'month': (['month'], totpcp_month.month), \
                            'lat': (['lat'], ds.lat), \
                            'lon': (['lon'], ds.lon)}, \
                attrs={'title': 'MCS precipitation statistics', \
                       'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                       'created_on':time.ctime(time.time())})

dsmap.month.attrs['long_name'] = 'month'
dsmap.month.attrs['units'] = 'month'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

# dsmap.ntimes.attrs['long_name'] = 'Number of hours in the month'
# dsmap.ntimes.attrs['units'] = 'count'

dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
dsmap.precipitation.attrs['units'] = 'mm/day'

dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
dsmap.mcs_precipitation.attrs['units'] = 'mm/day'

dsmap.mcs_precipitation_frac.attrs['long_name'] = 'MCS precipitation fraction'
dsmap.mcs_precipitation_frac.attrs['units'] = '%'

dsmap.mcs_precipitation_freq.attrs['long_name'] = 'MCS precipitation frequency'
dsmap.mcs_precipitation_freq.attrs['units'] = '%'

dsmap.mcs_precipitation_intensity.attrs['long_name'] = 'MCS precipitation intensity'
dsmap.mcs_precipitation_intensity.attrs['units'] = 'mm/hour'

fillvalue = np.nan
# Set encoding/compression for all variables
comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
encoding = {var: comp for var in dsmap.data_vars}

dsmap.to_netcdf(path=outfile_month, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)

#dsmap.to_netcdf(path=outfile_month, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='month', \
#                encoding={'lon':{'zlib':True, 'dtype':'float32'}, \
#                          'lat':{'zlib':True, 'dtype':'float32'}, \
#                          'precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation_frac':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          'mcs_precipitation_freq':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
#                          })
print('Monthly output saved as: ', outfile_month)
