"""
Calculates monthly total and MCS precipitation Hovmoller diagram and saves output to a netCDF data.

Author: Zhe Feng, zhe.feng@pnnl.gov
"""
import numpy as np
import glob, sys, os
import xarray as xr
import time
import yaml

if __name__ == "__main__":

    sdate = sys.argv[1]
    edate = sys.argv[2]
    year = (sys.argv[3])
    month = (sys.argv[4]).zfill(2)
    config_file = sys.argv[5]

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    pixel_dir = config['pixelfile_dir']
    output_monthly_dir = config['output_monthly_dir']
    pcpvarname = config['pcpvarname']
    startlat = config['startlat']
    endlat = config['endlat']
    startlon = config['startlon']
    endlon = config['endlon']

    mcsfiles = sorted(glob.glob(f'{pixel_dir}/{sdate}_{edate}/mcstrack_{year}{month}??_????.nc'))
    print(pixel_dir)
    print(year, month)
    print('Number of files: ', len(mcsfiles))
    os.makedirs(output_monthly_dir, exist_ok=True)

    hov_outfile = f'{output_monthly_dir}mcs_rainhov_{year}{month}.nc'

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
    var_dict = {
        'precipitation': (['time', 'lon'], totpreciphov.data),
        'mcs_precipitation': (['time', 'lon'], mcspreciphov.data),
    }
    coord_dict = {
        'lon': (['lon'], lonhov.data),
        'time': (['time'], basetime),
    }
    gattr_dict = {
        'title': 'MCS precipitation Hovmoller', \
        'startlat': startlat, \
        'endlat': endlat, \
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
        'created_on': time.ctime(time.time()),
    }
    dshov = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    dshov.lon.attrs['long_name'] = 'Longitude'
    dshov.lon.attrs['units'] = 'degree'
    dshov.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)' 
    dshov.time.attrs['units'] = 'seconds since 1970-01-01T00:00:00'
    dshov.precipitation.attrs['long_name'] = 'All precipitation'
    dshov.precipitation.attrs['units'] = 'mm/h'
    dshov.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
    dshov.mcs_precipitation.attrs['units'] = 'mm/h'

    fillvalue = np.nan
    # Set encoding/compression for all variables
    comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
    encoding = {var: comp for var in dshov.data_vars}

    dshov.to_netcdf(path=hov_outfile, mode='w', format='NETCDF4', unlimited_dims='time', encoding=encoding)
    print('Hovmoller output saved as: ' + hov_outfile)
