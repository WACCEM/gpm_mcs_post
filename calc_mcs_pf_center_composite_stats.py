import numpy as np
import glob, os
import xarray as xr
import time


if __name__ == "__main__":

    region = 'asia'
    
#     datadir = f'/global/cscratch1/sd/feng045/E3SM/GPM_IMERG/{region}/'
    datadir = f'/global/cscratch1/sd/feng045/E3SM/GPM_IMERG/{region}/mean_direction/'
    datafiles = sorted(glob.glob(f'{datadir}mcs_center_composite_20??????_20??????*.nc'))
    print(f'Number of files: {len(datafiles)}')

    outdir = datadir
    out_file = f'{outdir}mcs_center_composite_stats.nc'

    # Read in data
    ds = xr.open_mfdataset(datafiles, concat_dim='time', combine='nested')
    ntimes = ds.dims['time']

    # Count the number of valid frames
    ntimes_all = np.count_nonzero(~np.isnan(ds['time']))
    ntimes_ne = np.count_nonzero(~np.isnan(ds['time_ne']))
    ntimes_se = np.count_nonzero(~np.isnan(ds['time_se']))
    ntimes_sw = np.count_nonzero(~np.isnan(ds['time_sw']))
    ntimes_nw = np.count_nonzero(~np.isnan(ds['time_nw']))

    # Simple mean
    pcp_avg = ds['precipitation'].where(~np.isnan(ds['time']), drop=True).mean(dim='time').load()
    pcp_avg_ne = ds['precipitation_ne'].where(~np.isnan(ds['time_ne']), drop=True).mean(dim='time').load()
    pcp_avg_se = ds['precipitation_se'].where(~np.isnan(ds['time_se']), drop=True).mean(dim='time').load()
    pcp_avg_sw = ds['precipitation_sw'].where(~np.isnan(ds['time_sw']), drop=True).mean(dim='time').load()
    pcp_avg_nw = ds['precipitation_nw'].where(~np.isnan(ds['time_nw']), drop=True).mean(dim='time').load()

    # Rechunk time dimension for dask parallel computation
    pcp = ds['precipitation'].where(~np.isnan(ds['time']), drop=True).chunk({'time': -1})
    pcp_ne = ds['precipitation_ne'].where(~np.isnan(ds['time_ne']), drop=True).chunk({'time': -1})
    pcp_se = ds['precipitation_se'].where(~np.isnan(ds['time_se']), drop=True).chunk({'time': -1})
    pcp_sw = ds['precipitation_sw'].where(~np.isnan(ds['time_sw']), drop=True).chunk({'time': -1})
    pcp_nw = ds['precipitation_nw'].where(~np.isnan(ds['time_nw']), drop=True).chunk({'time': -1})
    quantile = [0.5, 0.9, 0.95, 0.99]
    pcp_percentiles = pcp.quantile(quantile, dim='time', keep_attrs=True).load()
    pcp_percentiles_ne = pcp_ne.quantile(quantile, dim='time', keep_attrs=True).load()
    pcp_percentiles_se = pcp_se.quantile(quantile, dim='time', keep_attrs=True).load()
    pcp_percentiles_sw = pcp_sw.quantile(quantile, dim='time', keep_attrs=True).load()
    pcp_percentiles_nw = pcp_nw.quantile(quantile, dim='time', keep_attrs=True).load()

    # Define xarray dataset for output
    dsout = xr.Dataset({'pcp_mean': (['y', 'x'], pcp_avg), \
                        'pcp_mean_ne': (['y', 'x'], pcp_avg_ne), \
                        'pcp_mean_se': (['y', 'x'], pcp_avg_se), \
                        'pcp_mean_sw': (['y', 'x'], pcp_avg_sw), \
                        'pcp_mean_nw': (['y', 'x'], pcp_avg_nw), \
                        'pcp_quantile': (['quantile','y', 'x'], pcp_percentiles), \
                        'pcp_quantile_ne': (['quantile','y', 'x'], pcp_percentiles_ne), \
                        'pcp_quantile_se': (['quantile','y', 'x'], pcp_percentiles_se), \
                        'pcp_quantile_sw': (['quantile','y', 'x'], pcp_percentiles_sw), \
                        'pcp_quantile_nw': (['quantile','y', 'x'], pcp_percentiles_nw), \
                        }, \
                        coords={'y': (['y'], ds.y), \
                                'x': (['x'], ds.x), \
                                'quantile': (['quantile'], quantile)
                                }, \
                        attrs={'title': 'MCS PF-center composite statistics', \
                                'lon_box':ds.attrs['lon_box'], \
                                'lat_box':ds.attrs['lat_box'], \
                                'ntimes_all':ntimes_all, \
                                'ntimes_ne':ntimes_ne, \
                                'ntimes_se':ntimes_se, \
                                'ntimes_sw':ntimes_sw, \
                                'ntimes_nw':ntimes_nw, \
                                'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                                'created_on':time.ctime(time.time())}
                        )

    dsout.x.attrs = ds.x.attrs
    dsout.y.attrs = ds.y.attrs

    dsout['pcp_mean'].attrs['long_name'] = 'Mean precipitation'
    dsout['pcp_mean'].attrs['units'] = 'mm/h'

    dsout['pcp_mean_ne'].attrs['long_name'] = 'Mean precipitation (Northeast moving)'
    dsout['pcp_mean_ne'].attrs['units'] = 'mm/h'

    dsout['pcp_mean_se'].attrs['long_name'] = 'Mean precipitation (Southeast moving)'
    dsout['pcp_mean_se'].attrs['units'] = 'mm/h'

    dsout['pcp_mean_sw'].attrs['long_name'] = 'Mean precipitation (Southwest moving)'
    dsout['pcp_mean_sw'].attrs['units'] = 'mm/h'

    dsout['pcp_mean_nw'].attrs['long_name'] = 'Mean precipitation (Northwest moving)'
    dsout['pcp_mean_nw'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile'].attrs['long_name'] = 'Precipitation quantiles'
    dsout['pcp_quantile'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile_ne'].attrs['long_name'] = 'Precipitation quantiles (Northeast moving)'
    dsout['pcp_quantile_ne'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile_se'].attrs['long_name'] = 'Precipitation quantiles (Southeast moving)'
    dsout['pcp_quantile_se'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile_sw'].attrs['long_name'] = 'Precipitation quantiles (Southwest moving)'
    dsout['pcp_quantile_sw'].attrs['units'] = 'mm/h'

    dsout['pcp_quantile_nw'].attrs['long_name'] = 'Precipitation quantiles (Northwest moving)'
    dsout['pcp_quantile_nw'].attrs['units'] = 'mm/h'

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netcdf file
    dsout.to_netcdf(path=out_file, mode='w', format='NETCDF4_CLASSIC', encoding=encoding)
    print('Output saved as: ', out_file)