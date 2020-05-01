"""
Calculates monthly total and MCS precipitation and saves output to a netCDF data.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
05/01/2020 - Written.
"""

import numpy as np
import glob, sys, os
import xarray as xr
import pandas as pd
import time, datetime, calendar, pytz


sdate = sys.argv[1]
edate = sys.argv[2]
year = (sys.argv[3])
month = (sys.argv[4]).zfill(2)
region = sys.argv[5]
pcpvarname = 'precipitation'

mcsdir = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/mcstracking_ccs4_4h/{sdate}_{edate}/'
outdir = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/stats_ccs4_4h/monthly/'
mcsfiles = sorted(glob.glob(mcsdir + 'mcstrack_' + year + month + '??_????.nc'))
print(mcsdir)
print(year, month)
print('Number of files: ', len(mcsfiles))
os.makedirs(outdir, exist_ok=True)

map_outfile = outdir + 'mcs_rainmap_' + year + month + '.nc'

# Read data
ds = xr.open_mfdataset(mcsfiles, concat_dim='time', combine='nested', drop_variables=['numclouds','tb'])
#ds.load()
print('Finish reading input files.')

ntimes = ds.dims['time']


# Sum MCS precipitation over time, use pcptracknumber > 0 as mask
mcsprecip = ds[pcpvarname].where(ds.pcptracknumber > 0).sum(dim='time')

# Sum total precipitation over time
totprecip = ds[pcpvarname].sum(dim='time')

# Convert all MCS track number to 1 for summation purpose
pcpnumber = ds.pcptracknumber.values
pcpnumber[pcpnumber > 0] = 1

# Convert numpy array to DataArray
mcspcpmask = xr.DataArray(pcpnumber, coords={'time':ds.time, 'lat':ds.lat, 'lon':ds.lon}, dims=['time','lat','lon'])

# Sum MCS PF counts overtime to get number of hours
mcspcpct = mcspcpmask.sum(dim='time')


# Compute Epoch Time for the month
months = np.zeros(1, dtype=int)
months[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

# Define xarray dataset for Map
dsmap = xr.Dataset({'precipitation': (['time', 'lat', 'lon'], totprecip.expand_dims('time', axis=0)), \
                    'mcs_precipitation': (['time', 'lat', 'lon'], mcsprecip.expand_dims('time', axis=0)), \
                    'mcs_precipitation_count': (['time', 'lat', 'lon'], mcspcpct.expand_dims('time', axis=0)), \
                    'ntimes': (['time'], xr.DataArray(ntimes).expand_dims('time', axis=0))}, \
                    coords={'time': (['time'], months), \
                            'lat': (['lat'], ds.lat), \
                            'lon': (['lon'], ds.lon)}, \
                attrs={'title': 'MCS precipitation accumulation', \
                       'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                       'created_on':time.ctime(time.time())})

dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
dsmap.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

dsmap.ntimes.attrs['long_name'] = 'Number of hours in the month'
dsmap.ntimes.attrs['units'] = 'count'

dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
dsmap.precipitation.attrs['units'] = 'mm'

dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
dsmap.mcs_precipitation.attrs['units'] = 'mm'

dsmap.mcs_precipitation_count.attrs['long_name'] = 'Number of hours MCS precipitation is recorded'
dsmap.mcs_precipitation_count.attrs['units'] = 'hour'


fillvalue = np.nan
dsmap.to_netcdf(path=map_outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', \
                encoding={'lon':{'zlib':True, 'dtype':'float32'}, \
                          'lat':{'zlib':True, 'dtype':'float32'}, \
                          'precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
                          'mcs_precipitation':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
                          'mcs_precipitation_count':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
                          })
print('Map output saved as: ', map_outfile)
