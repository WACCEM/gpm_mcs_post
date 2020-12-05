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
startlat = float(sys.argv[6])
endlat = float(sys.argv[7])
startlon = float(sys.argv[8])
endlon = float(sys.argv[9])

mcsdir = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/mcstracking_ccs4_4h/{sdate}_{edate}/'
# outdir = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region{region}/stats_ccs4_4h/monthly/'
# mcsdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/mcstracking_ccs4_4h/{sdate}_{edate}/'
outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/monthly/'

mcsfiles = sorted(glob.glob(mcsdir + 'mcstrack_' + year + month + '??_????.nc'))
print(mcsdir)
print(year, month)
print('Number of files: ', len(mcsfiles))
os.makedirs(outdir, exist_ok=True)

hov_outfile = outdir + 'mcs_rainhov_' + year + month + '.nc'

# Read data
ds = xr.open_mfdataset(mcsfiles, concat_dim='time', combine='nested', drop_variables=['numclouds','tb'])
print('Finish reading input files.')


# Mask out non-MCS precipitation, assign to a tmp array
mcspreciptmp = ds['precipitation'].where(ds.pcptracknumber > 0).values
# Replace NAN with 0, for averaging Hovmoller purpose
mcspreciptmp[np.isnan(mcspreciptmp)] = 0
# Convert numpy array to DataArray
mcspcp = xr.DataArray(mcspreciptmp, coords={'time':ds.time, 'lat':ds.lat, 'lon':ds.lon}, dims=['time','lat','lon'])


# Select a latitude band and time period
mcspreciphov = mcspcp.sel(lat=slice(startlat, endlat), lon=slice(startlon, endlon)).mean(dim='lat')
totpreciphov = ds['precipitation'].sel(lat=slice(startlat, endlat), lon=slice(startlon, endlon)).mean(dim='lat')

# Select time slice matching the chunk
timehov = ds['time']
lonhov = ds.lon.sel(lon=slice(startlon, endlon))

# Convert xarray decoded time back to Epoch Time in seconds
basetime = np.array([tt.tolist()/1e9 for tt in ds.time.values])


# Define xarray dataset for Hovmoller
print('Writing Hovmoller to netCDF file ...')
dshov = xr.Dataset({'precipitation': (['time', 'lon'], totpreciphov), \
                    'mcs_precipitation': (['time', 'lon'], mcspreciphov)}, \
                    coords={'lon': (['lon'], lonhov), \
                            'time': (['time'], basetime)}, \
                    attrs={'title': 'MCS precipitation Hovmoller', \
                           'startlat':startlat, \
                           'endlat':endlat, \
                           'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                           'created_on':time.ctime(time.time())})

dshov.lon.attrs['long_name'] = 'Longitude'
dshov.lon.attrs['units'] = 'degree'

dshov.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)' 
dshov.time.attrs['units'] = 'seconds since 1970-01-01T00:00:00'

dshov.precipitation.attrs['long_name'] = 'All precipitation'
dshov.precipitation.attrs['units'] = 'mm/h'

dshov.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
dshov.mcs_precipitation.attrs['units'] = 'mm/h'

fillvalue = np.nan
dshov.to_netcdf(path=hov_outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', \
                encoding={'precipitation':{'zlib':True, 'dtype':'float32'}, \
                          'mcs_precipitation':{'zlib':True, 'dtype':'float32'}})

print('Hovmoller output saved as: ' + hov_outfile)
