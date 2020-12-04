import numpy as np
import glob, os
import xarray as xr
import time, datetime, calendar, pytz


if __name__ == "__main__":

    region = 'asia'
    
    datadir = f'/global/cscratch1/sd/feng045/E3SM/GPM_IMERG/{region}/'
    datafiles = sorted(glob.glob(f'{datadir}mcs_center_composite_20??????_20??????*.nc'))
    print(f'Number of files: {len(datafiles)}')

    outdir = datadir
    out_file = f'{outdir}mcs_center_composite_stats.nc'

    # Read in data
    ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested')
    ntimes = ds.dims['time']

    # Simple mean
    pcp_avg = ds.precipitation.mean(dim='time').load()

    # Rechunk time dimension for dask parallel computation
    pcp = ds.precipitation.chunk({'time': -1})
    quantile = [0.5, 0.9, 0.95, 0.99]
    pcp_percentiles = pcp.quantile(quantile, dim='time', keep_attrs=True).load()

    # Define xarray dataset for output
    dsout = xr.Dataset({'pcp_mean': (['y', 'x'], pcp_avg), \
                        'pcp_quantile': (['quantile','y', 'x'], pcp_percentiles), \
                        }, \
                        coords={'y': (['y'], ds.y), \
                                'x': (['x'], ds.x), \
                                'quantile': (['quantile'], quantile)
                                }, \
                        attrs={'title': 'MCS PF-center composite statistics', \
                                'lon_box':ds.attrs['lon_box'], \
                                'lat_box':ds.attrs['lat_box'], \
                                'ntimes':ntimes, \
                                'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                                'created_on':time.ctime(time.time())}
                        )

    dsout.x.attrs = ds.x.attrs
    dsout.y.attrs = ds.y.attrs

    dsout['pcp_mean'].attrs['long_name'] = 'Mean precipitation'
    dsout['pcp_mean'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile'].attrs['long_name'] = 'Precipitation quantiles'
    dsout['pcp_quantile'].attrs['units'] = 'mm/h'

    # Set encoding/compression for all variables
    comp = dict(zlib=True)
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netcdf file
    dsout.to_netcdf(path=out_file, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
    print('Output saved as: ', out_file)