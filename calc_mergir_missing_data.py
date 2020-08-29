"""
Calculate number of valid data within a period of time from MERGIR global data.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "26-Feb-2019"

import numpy as np
import glob, os, sys
import xarray as xr
import time, datetime, calendar, pytz

year = sys.argv[1]
month = (sys.argv[2]).zfill(2)

datadir = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/{year}/'
datafiles = sorted(glob.glob(f'{datadir}*merg_{year}{month}????_4km-pixel.nc'))
print(f'Number of files: {len(datafiles)}')

outdir = '/global/cscratch1/sd/feng045/waccem/MERGIR_Global/Regrid/stats/'
map_outfile = f'{outdir}merg_monthly_validcount_{year}{month}.nc'


# Read data
# ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='by_coords')
ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested')

# Total number of times
ntimes = ds.dims['time']

# Get minutes from time
minutes = ds.time.dt.minute
# Find time index for 00 min and 30 min
t00idx = np.where(minutes == 0)[0]
t30idx = np.where(minutes == 30)[0]
ntimes00 = np.count_nonzero(t00idx)
ntimes30 = np.count_nonzero(t30idx)

# Make a mask for valid data (> 0) for all times == 00 min and 30 min, sum over time to get counts
count00 = (ds.Tb.isel(time=t00idx) > 0).sum(dim='time')
count30 = (ds.Tb.isel(time=t30idx) > 0).sum(dim='time')

# Compute Epoch Time for the month
date = np.zeros(1, dtype=int)
date[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

# Define xarray dataset for map
print('Writing map to netCDF file ...')
dsmap = xr.Dataset({'count_00min': (['time', 'lat', 'lon'], count00.expand_dims('time', axis=0)), \
                    'count_30min': (['time', 'lat', 'lon'], count30.expand_dims('time', axis=0)), \
                    'ntimes_00min': (['time'], np.expand_dims(ntimes00, axis=0)), \
                    'ntimes_30min': (['time'], np.expand_dims(ntimes30, axis=0)), \
                    'ntimes': (['time'], np.expand_dims(ntimes, axis=0))}, \
                    coords={'lon': (['lon'], ds.lon), \
                            'lat': (['lat'], ds.lat), \
                            'time': (['time'], date),}, \
                    attrs={'title': 'Valid data counts', \
                           'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                           'created_on':time.ctime(time.time())})

dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
dsmap.time.attrs['units'] = 'seconds since 1970-1-1 0:00:00 0:00'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

dsmap.count_00min.attrs['long_name'] = 'Valid data count at 00 min'
dsmap.count_00min.attrs['units'] = 'count'

dsmap.count_30min.attrs['long_name'] = 'Valid data count at 30 min'
dsmap.count_30min.attrs['units'] = 'count'

dsmap.ntimes_00min.attrs['long_name'] = 'Number of times at 00 min'
dsmap.ntimes_00min.attrs['units'] = 'count'

dsmap.ntimes_30min.attrs['long_name'] = 'Number of times at 30 min'
dsmap.ntimes_30min.attrs['units'] = 'count'

dsmap.ntimes.attrs['long_name'] = 'Number of times in the month'
dsmap.ntimes.attrs['units'] = 'count'

fillvalue = np.nan
dsmap.to_netcdf(path=map_outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', \
                encoding={'count_00min':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}, \
                          'count_30min':{'zlib':True, '_FillValue':fillvalue, 'dtype':'float32'}})

print('Map output saved as: ' + map_outfile)
