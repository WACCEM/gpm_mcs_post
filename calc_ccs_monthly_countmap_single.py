"""
Calculates monthly CCS count map and saves output to a netCDF data.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
10/07/2021 - Written.
"""

import numpy as np
import glob, sys, os
import xarray as xr
import time, datetime, calendar, pytz
import yaml

sdate = sys.argv[1]
edate = sys.argv[2]
year = (sys.argv[3])
month = (sys.argv[4]).zfill(2)
config_file = sys.argv[5]

# get inputs from configuration file
stream = open(config_file, 'r')
config = yaml.full_load(stream)
pixel_dir = config['pixelfile_dir']
output_dir = config['output_dir']

mcsfiles = sorted(glob.glob(f'{pixel_dir}/{sdate}_{edate}/mcstrack_{year}{month}??_????.nc'))
print(pixel_dir)
print(year, month)
print('Number of files: ', len(mcsfiles))
os.makedirs(output_dir, exist_ok=True)

out_filename = f'{output_dir}ccs_freqmap_{year}{month}.nc'

# Read data
drop_varlist = ['longitude','latitude','numclouds','pcptracknumber',]
ds = xr.open_mfdataset(mcsfiles, concat_dim='time', combine='nested', 
                        drop_variables=drop_varlist)
print('Finish reading input files.')

ntimes = ds.dims['time']

tb = ds['tb']
tb240_counts = (tb < 240).sum(dim='time')
tb220_counts = (tb < 220).sum(dim='time')
tb200_counts = (tb < 200).sum(dim='time')
rainrate_thres = 1  # [mm/h]
tb240rain_counts = ((tb < 240) * (ds['precipitation'] > rainrate_thres)).sum(dim='time')
tb220rain_counts = ((tb < 220) * (ds['precipitation'] > rainrate_thres)).sum(dim='time')
tb200rain_counts = ((tb < 200) * (ds['precipitation'] > rainrate_thres)).sum(dim='time')

# Compute Epoch Time for the month
months = np.zeros(1, dtype=int)
months[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

# Define xarray dataset for Map
dsout = xr.Dataset({'tb240_counts': (['time', 'lat', 'lon'], tb240_counts.expand_dims('time', axis=0)), \
                    'tb220_counts': (['time', 'lat', 'lon'], tb220_counts.expand_dims('time', axis=0)), \
                    'tb200_counts': (['time', 'lat', 'lon'], tb200_counts.expand_dims('time', axis=0)), \
                    'tb240rain_counts': (['time', 'lat', 'lon'], tb240rain_counts.expand_dims('time', axis=0)), \
                    'tb220rain_counts': (['time', 'lat', 'lon'], tb220rain_counts.expand_dims('time', axis=0)), \
                    'tb200rain_counts': (['time', 'lat', 'lon'], tb200rain_counts.expand_dims('time', axis=0)), \
                    'ntimes': (['time'], xr.DataArray(ntimes).expand_dims('time', axis=0))}, \
                    coords={'time': (['time'], months), \
                            'lat': (['lat'], ds.lat), \
                            'lon': (['lon'], ds.lon)}, \
                attrs={'title': 'CCS occurrence counts', \
                       'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                       'created_on':time.ctime(time.time())})

dsout.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
dsout.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'
dsout.lon.attrs['long_name'] = 'Longitude'
dsout.lon.attrs['units'] = 'degree'
dsout.lat.attrs['long_name'] = 'Latitude'
dsout.lat.attrs['units'] = 'degree'

dsout.ntimes.attrs['long_name'] = 'Number of hours in the month'
dsout.ntimes.attrs['units'] = 'count'

dsout.tb240_counts.attrs['long_name'] = 'Tb < 240K counts'
dsout.tb240_counts.attrs['units'] = 'count'
dsout.tb220_counts.attrs['long_name'] = 'Tb < 220K counts'
dsout.tb220_counts.attrs['units'] = 'count'
dsout.tb220_counts.attrs['long_name'] = 'Tb < 200K counts'
dsout.tb220_counts.attrs['units'] = 'count'

dsout.tb240rain_counts.attrs['long_name'] = f'Tb < 240K & rain rate > {rainrate_thres} mm/h counts'
dsout.tb240rain_counts.attrs['units'] = 'count'
dsout.tb220rain_counts.attrs['long_name'] = f'Tb < 220K & rain rate > {rainrate_thres} mm/h counts'
dsout.tb220rain_counts.attrs['units'] = 'count'
dsout.tb220rain_counts.attrs['long_name'] = f'Tb < 200K & rain rate > {rainrate_thres} mm/h counts'
dsout.tb220rain_counts.attrs['units'] = 'count'

# Set encoding/compression for all variables
comp = dict(zlib=True, dtype='float32')
encoding = {var: comp for var in dsout.data_vars}

fillvalue = np.nan
dsout.to_netcdf(path=out_filename, mode='w', format='NETCDF4', unlimited_dims='time', encoding=encoding)
print('Output saved as: ', out_filename)

import pdb; pdb.set_trace()