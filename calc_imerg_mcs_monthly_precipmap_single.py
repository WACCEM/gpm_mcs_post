"""
Calculates monthly total and MCS precipitation and saves output to a netCDF data.

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

    sdate = sys.argv[1]
    edate = sys.argv[2]
    year = (sys.argv[3])
    month = (sys.argv[4]).zfill(2)
    # region = sys.argv[5]
    config_file = sys.argv[5]
    # pcpvarname = 'precipitation'

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    pixel_dir = config['pixelfile_dir']
    output_monthly_dir = config['output_monthly_dir']
    pcpvarname = config['pcpvarname']

    mcsfiles = sorted(glob.glob(f'{pixel_dir}/{sdate}_{edate}/mcstrack_{year}{month}??_????.nc'))
    print(pixel_dir)
    print(year, month)
    print('Number of files: ', len(mcsfiles))
    os.makedirs(output_monthly_dir, exist_ok=True)

    map_outfile = f'{output_monthly_dir}mcs_rainmap_{year}{month}.nc'

    # Read data
    ds = xr.open_mfdataset(mcsfiles, concat_dim='time', combine='nested', drop_variables=['numclouds','tb'])
    print('Finish reading input files.')
    ntimes = ds.dims['time']

    # Sum MCS precipitation over time, use cloudtracknumber > 0 as mask
    mcsprecip = ds[pcpvarname].where(ds['cloudtracknumber'] > 0).sum(dim='time')
    #mcsprecip = ds[pcpvarname].where(ds.pcptracknumber > 0).sum(dim='time')

    # Sum total precipitation over time
    totprecip = ds[pcpvarname].sum(dim='time')

    # Convert all MCS track number to 1 for summation purpose
    pcpnumber = ds['pcptracknumber'].values
    pcpnumber[pcpnumber > 0] = 1

    # Convert numpy array to DataArray
    mcspcpmask = xr.DataArray(pcpnumber, coords={'time':ds.time, 'lat':ds.lat, 'lon':ds.lon}, dims=['time','lat','lon'])

    # Sum MCS PF counts overtime to get number of hours
    mcspcpct = mcspcpmask.sum(dim='time')


    # Compute Epoch Time for the month
    months = np.zeros(1, dtype=int)
    months[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

    ############################################################################
    # Write output file
    var_dict = {
        'precipitation': (['time', 'lat', 'lon'], totprecip.expand_dims('time', axis=0)),
        'mcs_precipitation': (['time', 'lat', 'lon'], mcsprecip.expand_dims('time', axis=0)),
        'mcs_precipitation_count': (['time', 'lat', 'lon'], mcspcpct.expand_dims('time', axis=0)),
        'ntimes': (['time'], xr.DataArray(ntimes).expand_dims('time', axis=0)),
    }
    coord_dict = {
        'time': (['time'], months),
        'lat': (['lat'], ds.lat),
        'lon': (['lon'], ds.lon),
    }
    gattr_dict = {
        'title': 'MCS precipitation accumulation',
        'contact':'Zhe Feng, zhe.feng@pnnl.gov',
        'created_on':time.ctime(time.time()),
    }
    dsmap = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
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
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dsmap.data_vars}

    dsmap.to_netcdf(path=map_outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)
    print('Map output saved as: ', map_outfile)
