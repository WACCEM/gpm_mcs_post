"""
Script to calculate multi- seasonal and monthly mean MCS statistics from monthly data.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "19-Jun-2019"

import numpy as np
import glob, os, sys
import xarray as xr
import pandas as pd
import time, datetime, calendar, pytz

region = sys.argv[1]
# region = 'china'

datadir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_pt1/monthly/'
#datadir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/monthly/'
#datadir = f'/global/cscratch1/sd/feng045/waccem/mcs_global/{region}/stats/monthly_satmcs/'
#pcpdir = f'/global/cscratch1/sd/feng045/waccem/mcs_global/{region}/stats/monthly/'
pcpdir = datadir
datafiles = sorted(glob.glob(f'{datadir}mcs_statsmap_20????*nc'))
pcpfiles = sorted(glob.glob(f'{pcpdir}mcs_rainmap_20????*nc'))
print(f'Number of files: {len(datafiles)}')

outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_pt1/climo/'
#outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/climo/'
os.makedirs(outdir, exist_ok=True)
#outfile_season = f'{outdir}mcs_statsmap_seasonal_mean.nc'
#outfile_month = f'{outdir}mcs_statsmap_monthly_mean.nc'

# Read stats data
ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested')
nx = ds.dims['lon']
ny = ds.dims['lat']
lon = ds.lon
lat = ds.lat
#ds

# Read precip data (for total number of hours in each month only)
dsp = xr.open_mfdataset(pcpfiles, concat_dim='time', combine='nested')
nhours = dsp.ntimes

# Get number of seasons/years
years = pd.date_range(start=ds.time.min().values, end=ds.time.max().values, freq='AS')
nyears = len(years)

# Get min/max year and create output file names
min_year = ds.time.min().dt.strftime("%Y").values
max_year = ds.time.max().dt.strftime("%Y").values
outfile_season = f'{outdir}mcs_statsmap_seasonal_mean_{min_year}_{max_year}.nc'
outfile_month = f'{outdir}mcs_statsmap_monthly_mean_{min_year}_{max_year}.nc'

# # Find number of years for month
# uniq_months = np.unique(ds.time.dt.month)
# print(uniq_months)
# for im in range(0, len(uniq_months)):
#     print(im, len(np.where(ds.time.dt.month == uniq_months[0])[0]))


# For number counts, group by month and sum, then divide by total number of years
mcs_number_ccs_month = ds.mcs_number_ccs.groupby('time.month').sum(dim='time').load() / nyears
mcs_number_pf_month = ds.mcs_number_pf.groupby('time.month').sum(dim='time').load() / nyears
mcs_initccs_month = ds.initiation_ccs.groupby('time.month').sum(dim='time').load() / nyears

# Calculate total number of hours each month
nhours_month = nhours.groupby('time.month').sum(dim='time')
# Calculate number of hours for MCS PF over each pixel
mcs_nhour_pf_month = ds.mcs_nhour_pf.groupby('time.month').sum(dim='time')
mcs_nhour_ccs_month = ds.mcs_nhour_ccs.groupby('time.month').sum(dim='time')

# MCS frequency
mcs_freq_ccs_month = 100. * mcs_nhour_ccs_month / nhours_month
mcs_freq_pf_month = 100. * mcs_nhour_pf_month / nhours_month

# MCS size
mcs_pfarea_avg_month = (ds.pf_area_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month
mcs_pfdiam_avg_month = 2 * np.sqrt(mcs_pfarea_avg_month / np.pi)

# MCS lifetime
mcs_lifetime_avg_month = (ds.lifetime_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month
#lifetime_pt_avg_month = ds.lifetime_pt.groupby('time.month').mean(dim='time')

totalrain_avg_month = (ds.totalrain_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month
totalrainheavy_avg_month = (ds.totalrainheavy_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month
rainrateheavy_avg_month = (ds.rainrateheavy_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month
rainratemax_avg_month = (ds.rainratemax_mean * ds.mcs_nhour_pf).groupby('time.month').sum(dim='time') / mcs_nhour_pf_month

# Propagation Speed
nhour_pfspeed_mcs_month = ds.mcs_nhour_speedmcs.groupby('time.month').sum(dim='time')
pfspeed_mcs_avg_month = (ds.pf_speed_mcs * ds.mcs_nhour_speedmcs).groupby('time.month').sum(dim='time') / nhour_pfspeed_mcs_month
pfuspeed_mcs_avg_month = (ds.pf_uspeed_mcs * ds.mcs_nhour_speedmcs).groupby('time.month').sum(dim='time') / nhour_pfspeed_mcs_month
pfvspeed_mcs_avg_month = (ds.pf_vspeed_mcs * ds.mcs_nhour_speedmcs).groupby('time.month').sum(dim='time') / nhour_pfspeed_mcs_month



# For number counts, group by season and sum, then divide by total number of seasons (years)
mcs_number_ccs_season = ds.mcs_number_ccs.groupby('time.season').sum(dim='time').load() / nyears
mcs_number_pf_season = ds.mcs_number_pf.groupby('time.season').sum(dim='time').load() / nyears
mcs_initccs_season = ds.initiation_ccs.groupby('time.season').sum(dim='time').load() / nyears

# Calculate total number of hours each season
nhours_season = nhours.groupby('time.season').sum(dim='time')
# Calculate number of hours for MCS PF over each pixel
mcs_nhour_pf_season = ds.mcs_nhour_pf.groupby('time.season').sum(dim='time')
mcs_nhour_ccs_season = ds.mcs_nhour_ccs.groupby('time.season').sum(dim='time')

# MCS frequency
mcs_freq_ccs_season = 100. * mcs_nhour_ccs_season / nhours_season
mcs_freq_pf_season = 100. * mcs_nhour_pf_season / nhours_season

# MCS size
mcs_pfarea_avg_season = (ds.pf_area_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season
mcs_pfdiam_avg_season = 2 * np.sqrt(mcs_pfarea_avg_season / np.pi)

# MCS lifetime
mcs_lifetime_avg_season = (ds.lifetime_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season
#lifetime_pt_avg_season = ds.lifetime_pt.groupby('time.season').mean(dim='time')

totalrain_avg_season = (ds.totalrain_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season
totalrainheavy_avg_season = (ds.totalrainheavy_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season
rainrateheavy_avg_season = (ds.rainrateheavy_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season
rainratemax_avg_season = (ds.rainratemax_mean * ds.mcs_nhour_pf).groupby('time.season').sum(dim='time') / mcs_nhour_pf_season

# Propagation Speed
nhour_pfspeed_mcs_season = ds.mcs_nhour_speedmcs.groupby('time.season').sum(dim='time')
pfspeed_mcs_avg_season = (ds.pf_speed_mcs * ds.mcs_nhour_speedmcs).groupby('time.season').sum(dim='time') / nhour_pfspeed_mcs_season
pfuspeed_mcs_avg_season = (ds.pf_uspeed_mcs * ds.mcs_nhour_speedmcs).groupby('time.season').sum(dim='time') / nhour_pfspeed_mcs_season
pfvspeed_mcs_avg_season = (ds.pf_vspeed_mcs * ds.mcs_nhour_speedmcs).groupby('time.season').sum(dim='time') / nhour_pfspeed_mcs_season


# Write seasonal mean output file
dsmap = xr.Dataset({'mcs_number': (['season', 'lat', 'lon'], mcs_number_pf_season), \
                    'mcs_number_ccs': (['season', 'lat', 'lon'], mcs_number_ccs_season), \
                    'mcs_freq_ccs': (['season', 'lat', 'lon'], mcs_freq_ccs_season), \
                    'mcs_freq_pf': (['season', 'lat', 'lon'], mcs_freq_pf_season), \
                    'mcs_initiation_ccs': (['season', 'lat', 'lon'], mcs_initccs_season), \
                    'mcs_pfarea': (['season', 'lat', 'lon'], mcs_pfarea_avg_season), \
                    'mcs_pfdiameter': (['season', 'lat', 'lon'], mcs_pfdiam_avg_season), \
                    'mcs_lifetime': (['season', 'lat', 'lon'], mcs_lifetime_avg_season), \
                    'mcs_totalrain': (['season', 'lat', 'lon'], totalrain_avg_season), \
                    'mcs_totalrainheavy': (['season', 'lat', 'lon'], totalrainheavy_avg_season), \
                    'mcs_rainrateheavy': (['season', 'lat', 'lon'], rainrateheavy_avg_season), \
                    'mcs_rainratemax': (['season', 'lat', 'lon'], rainratemax_avg_season), \
                    'mcs_speed': (['season', 'lat', 'lon'], pfspeed_mcs_avg_season), \
                    'mcs_uspeed': (['season', 'lat', 'lon'], pfuspeed_mcs_avg_season), \
                    'mcs_vspeed': (['season', 'lat', 'lon'], pfvspeed_mcs_avg_season), \
                    }, \
                    coords={'season': (['season'], mcs_number_pf_season.season), \
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

dsmap.mcs_number.attrs['long_name'] = 'Number of MCS (defined by PF)'
dsmap.mcs_number.attrs['units'] = 'count'

dsmap.mcs_number_ccs.attrs['long_name'] = 'Number of MCS (defined by CCS)'
dsmap.mcs_number_ccs.attrs['units'] = 'count'

dsmap.mcs_initiation_ccs.attrs['long_name'] = 'MCS initiation count'
dsmap.mcs_initiation_ccs.attrs['units'] = 'count'

dsmap.mcs_freq_ccs.attrs['long_name'] = 'MCS frequency (defined by CCS)'
dsmap.mcs_freq_ccs.attrs['units'] = '%'

dsmap.mcs_freq_pf.attrs['long_name'] = 'MCS frequency (defined by PF)'
dsmap.mcs_freq_pf.attrs['units'] = '%'

dsmap.mcs_pfarea.attrs['long_name'] = 'MCS precipitation feature area'
dsmap.mcs_pfarea.attrs['units'] = 'km2'

dsmap.mcs_pfdiameter.attrs['long_name'] = 'MCS precipitation feature equivalent diameter'
dsmap.mcs_pfdiameter.attrs['units'] = 'km'

dsmap.mcs_lifetime.attrs['long_name'] = 'MCS lifetime'
dsmap.mcs_lifetime.attrs['units'] = 'hour'

dsmap.mcs_totalrain.attrs['long_name'] = 'MCS total precipitation'
dsmap.mcs_totalrain.attrs['units'] = 'mm'

dsmap.mcs_totalrainheavy.attrs['long_name'] = 'MCS total heavy precipitation (rain rate > 10 mm/h)'
dsmap.mcs_totalrainheavy.attrs['units'] = 'mm'

dsmap.mcs_rainrateheavy.attrs['long_name'] = 'MCS heavy rain rate (rain rate > 10 mm/h) mean'
dsmap.mcs_rainrateheavy.attrs['units'] = 'mm/h'

dsmap.mcs_rainratemax.attrs['long_name'] = 'MCS max rain rate'
dsmap.mcs_rainratemax.attrs['units'] = 'mm/h'

dsmap.mcs_speed.attrs['long_name'] = 'MCS precipitation feature propagation speed'
dsmap.mcs_speed.attrs['units'] = 'm/s'

dsmap.mcs_uspeed.attrs['long_name'] = 'MCS precipitation feature propagation speed (zonal direction)'
dsmap.mcs_uspeed.attrs['units'] = 'm/s'

dsmap.mcs_vspeed.attrs['long_name'] = 'MCS precipitation feature propagation speed (meridional direction)'
dsmap.mcs_vspeed.attrs['units'] = 'm/s'

fillvalue = np.nan
# Set encoding/compression for all variables
comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
encoding = {var: comp for var in dsmap.data_vars}

dsmap.to_netcdf(path=outfile_season, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
print('Seasonal output saved as: ', outfile_season)



# Write monthly mean output file
dsmap = xr.Dataset({'mcs_number': (['month', 'lat', 'lon'], mcs_number_pf_month), \
                    'mcs_number_ccs': (['month', 'lat', 'lon'], mcs_number_ccs_month), \
                    'mcs_initiation_ccs': (['month', 'lat', 'lon'], mcs_initccs_month), \
                    'mcs_freq_ccs': (['month', 'lat', 'lon'], mcs_freq_ccs_month), \
                    'mcs_freq_pf': (['month', 'lat', 'lon'], mcs_freq_pf_month), \
                    'mcs_pfarea': (['month', 'lat', 'lon'], mcs_pfarea_avg_month), \
                    'mcs_pfdiameter': (['month', 'lat', 'lon'], mcs_pfdiam_avg_month), \
                    'mcs_lifetime': (['month', 'lat', 'lon'], mcs_lifetime_avg_month), \
                    'mcs_totalrain': (['month', 'lat', 'lon'], totalrain_avg_month), \
                    'mcs_totalrainheavy': (['month', 'lat', 'lon'], totalrainheavy_avg_month), \
                    'mcs_rainrateheavy': (['month', 'lat', 'lon'], rainrateheavy_avg_month), \
                    'mcs_rainratemax': (['month', 'lat', 'lon'], rainratemax_avg_month), \
                    'mcs_speed': (['month', 'lat', 'lon'], pfspeed_mcs_avg_month), \
                    'mcs_uspeed': (['month', 'lat', 'lon'], pfuspeed_mcs_avg_month), \
                    'mcs_vspeed': (['month', 'lat', 'lon'], pfvspeed_mcs_avg_month), \
                    }, \
                    coords={'month': (['month'], mcs_number_pf_month.month), \
                            'lat': (['lat'], ds.lat), \
                            'lon': (['lon'], ds.lon)}, \
                attrs={'title': 'MCS precipitation statistics', \
                       'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                       'created_on':time.ctime(time.time())})

dsmap.month.attrs['long_name'] = 'Month'
dsmap.month.attrs['units'] = 'none'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

dsmap.mcs_number.attrs['long_name'] = 'Number of MCS (defined by PF)'
dsmap.mcs_number.attrs['units'] = 'count'

dsmap.mcs_number_ccs.attrs['long_name'] = 'Number of MCS (defined by CCS)'
dsmap.mcs_number_ccs.attrs['units'] = 'count'

dsmap.mcs_freq_ccs.attrs['long_name'] = 'MCS frequency (defined by CCS)'
dsmap.mcs_freq_ccs.attrs['units'] = '%'

dsmap.mcs_freq_pf.attrs['long_name'] = 'MCS frequency (defined by PF)'
dsmap.mcs_freq_pf.attrs['units'] = '%'

dsmap.mcs_initiation_ccs.attrs['long_name'] = 'MCS initiation count'
dsmap.mcs_initiation_ccs.attrs['units'] = 'count'

dsmap.mcs_pfarea.attrs['long_name'] = 'MCS precipitation feature area'
dsmap.mcs_pfarea.attrs['units'] = 'km2'

dsmap.mcs_pfdiameter.attrs['long_name'] = 'MCS precipitation feature equivalent diameter'
dsmap.mcs_pfdiameter.attrs['units'] = 'km'

dsmap.mcs_lifetime.attrs['long_name'] = 'MCS lifetime'
dsmap.mcs_lifetime.attrs['units'] = 'hour'

dsmap.mcs_totalrain.attrs['long_name'] = 'MCS total precipitation'
dsmap.mcs_totalrain.attrs['units'] = 'mm'

dsmap.mcs_totalrainheavy.attrs['long_name'] = 'MCS total heavy precipitation (rain rate > 10 mm/h)'
dsmap.mcs_totalrainheavy.attrs['units'] = 'mm'

dsmap.mcs_rainrateheavy.attrs['long_name'] = 'MCS heavy rain rate (rain rate > 10 mm/h) mean'
dsmap.mcs_rainrateheavy.attrs['units'] = 'mm/h'

dsmap.mcs_rainratemax.attrs['long_name'] = 'MCS max rain rate'
dsmap.mcs_rainratemax.attrs['units'] = 'mm/h'

dsmap.mcs_speed.attrs['long_name'] = 'MCS precipitation feature propagation speed'
dsmap.mcs_speed.attrs['units'] = 'm/s'

dsmap.mcs_uspeed.attrs['long_name'] = 'MCS precipitation feature propagation speed (zonal direction)'
dsmap.mcs_uspeed.attrs['units'] = 'm/s'

dsmap.mcs_vspeed.attrs['long_name'] = 'MCS precipitation feature propagation speed (meridional direction)'
dsmap.mcs_vspeed.attrs['units'] = 'm/s'

fillvalue = np.nan
# Set encoding/compression for all variables
comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
encoding = {var: comp for var in dsmap.data_vars}

dsmap.to_netcdf(path=outfile_month, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
print('Monthly output saved as: ', outfile_month)

