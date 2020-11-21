import numpy as np
import glob, os, sys
import xarray as xr
import pandas as pd
import time, datetime, calendar, pytz
import dask
from dask import delayed
from dask.distributed import Client, LocalCluster

def merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR):
    
    # Create a global variable
    var_merge = np.full((1, ny_g, nx_g), np.nan, dtype=float)

    # asia
    var_merge[:, GlatB["asia"]:GlatT["asia"]+1, GlonL["asia"]:GlonR["asia"]+1] = var_list["asia"][:, latB["asia"]:latT["asia"]+1, lonL["asia"]:lonR["asia"]+1]
    # NAM
    var_merge[:, GlatB["nam"]:GlatT["nam"]+1, GlonL["nam"]:GlonR["nam"]+1] = var_list["nam"][:, latB["nam"]:latT["nam"]+1, lonL["nam"]:lonR["nam"]+1]
    # SPAC
    var_merge[:, GlatB["spac"]:GlatT["spac"]+1, GlonL["spac"]:GlonR["spac"]+1] = var_list["spac"][:, latB["spac"]:latT["spac"]+1, lonL["spac"]:lonR["spac"]+1]
    
    return var_merge

def write_netcdf(out_filename, basetime, lat_G, lon_G, ntimes_asia, ntimes_nam, ntimes_spac, \
                    mcs_number_ccs, mcs_number_pf, mcs_nhour_ccs, mcs_nhour_pf, mcs_nhour_speedmcs, \
                    lifetime_mean, ccs_area_mean, pf_area_mean, rainratemax_mean, initiation_ccs, \
                    pf_speed_mcs, pf_uspeed_mcs, pf_vspeed_mcs, \
                    precipitation, mcs_precipitation, mcs_precipitation_count):
    # Write seasonal mean output file
    dsmap = xr.Dataset({'mcs_number_ccs': (['time', 'lat', 'lon'], mcs_number_ccs), \
                        'mcs_number_pf': (['time', 'lat', 'lon'], mcs_number_pf), \
                        'mcs_nhour_ccs': (['time', 'lat', 'lon'], mcs_nhour_ccs), \
                        'mcs_nhour_pf': (['time', 'lat', 'lon'], mcs_nhour_pf), \
                        'mcs_nhour_speedmcs': (['time', 'lat', 'lon'], mcs_nhour_speedmcs), \
                        'lifetime_mean': (['time', 'lat', 'lon'], lifetime_mean), \
                        'ccs_area_mean': (['time', 'lat', 'lon'], ccs_area_mean), \
                        'pf_area_mean': (['time', 'lat', 'lon'], pf_area_mean), \
                        'rainratemax_mean': (['time', 'lat', 'lon'], rainratemax_mean), \
                        'initiation_ccs': (['time', 'lat', 'lon'], initiation_ccs), \
                        'pf_speed_mcs': (['time', 'lat', 'lon'], pf_speed_mcs), \
                        'pf_uspeed_mcs': (['time', 'lat', 'lon'], pf_uspeed_mcs), \
                        'pf_vspeed_mcs': (['time', 'lat', 'lon'], pf_vspeed_mcs), \
                        'precipitation': (['time', 'lat', 'lon'], precipitation), \
                        'mcs_precipitation': (['time', 'lat', 'lon'], mcs_precipitation), \
                        'mcs_precipitation_count': (['time', 'lat', 'lon'], mcs_precipitation_count), \
                        'ntimes_asia': (['time'], ntimes_asia), \
                        'ntimes_nam': (['time'], ntimes_nam), \
                        'ntimes_spac': (['time'], ntimes_spac), \
                        }, \
                        coords={'time': (['time'], basetime), \
                                'lat': (['lat'], lat_G), \
                                'lon': (['lon'], lon_G)}, \
                        attrs={'title': 'MCS monthly mean statistics', \
                            'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                            'created_on':time.ctime(time.time()), \
                            'produced_with':os.path.abspath(sys.argv[0]),\
                            })

    dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
    dsmap.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'

    dsmap.lon.attrs['long_name'] = 'Longitude'
    dsmap.lon.attrs['units'] = 'degree'

    dsmap.lat.attrs['long_name'] = 'Latitude'
    dsmap.lat.attrs['units'] = 'degree'

    # Frequency variables
    dsmap.mcs_number_ccs.attrs['long_name'] = 'Number of MCS defined by cold cloud shield'
    dsmap.mcs_number_ccs.attrs['units'] = 'count'

    dsmap.mcs_number_pf.attrs['long_name'] = 'Number of MCS defined by precipitation feature'
    dsmap.mcs_number_pf.attrs['units'] = 'count'

    dsmap.mcs_nhour_ccs.attrs['long_name'] = 'Number of MCS hours defined by cold cloud shield'
    dsmap.mcs_nhour_ccs.attrs['units'] = 'hour'

    dsmap.mcs_nhour_pf.attrs['long_name'] = 'Number of MCS hours defined by precipitation feature'
    dsmap.mcs_nhour_pf.attrs['units'] = 'hour'

    dsmap.mcs_nhour_speedmcs.attrs['long_name'] = 'Number of MCS hours for propagation speed (MCS status is met)'
    dsmap.mcs_nhour_speedmcs.attrs['units'] = 'hour'

    # Mean variables
    dsmap.lifetime_mean.attrs['long_name'] = 'MCS lifetime mean'
    dsmap.lifetime_mean.attrs['units'] = 'hour'

    dsmap.ccs_area_mean.attrs['long_name'] = 'MCS cold cloud shield area mean'
    dsmap.ccs_area_mean.attrs['units'] = 'km2'

    dsmap.pf_area_mean.attrs['long_name'] = 'MCS precipitation feature area mean'
    dsmap.pf_area_mean.attrs['units'] = 'km2'

    dsmap.rainratemax_mean.attrs['long_name'] = 'MCS max rain rate mean'
    dsmap.rainratemax_mean.attrs['units'] = 'mm/h'

    dsmap.initiation_ccs.attrs['long_name'] = 'MCS convective initiation hours defined by cold cloud shield'
    dsmap.initiation_ccs.attrs['units'] = 'hour'

    dsmap.pf_speed_mcs.attrs['long_name'] = 'Propagation speed when MCS status is met'
    dsmap.pf_speed_mcs.attrs['units'] = 'm/s'

    dsmap.pf_uspeed_mcs.attrs['long_name'] = 'Propagation speed (x-direction) when MCS status is met'
    dsmap.pf_uspeed_mcs.attrs['units'] = 'm/s'

    dsmap.pf_vspeed_mcs.attrs['long_name'] = 'Propagation speed (y-direction) when MCS status is met'
    dsmap.pf_vspeed_mcs.attrs['units'] = 'm/s'

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

    dsmap.to_netcdf(path=out_filename, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', encoding=encoding)
    print('Map output saved as: ', out_filename)

    return True


def combine_regions(sfile_asia, sfile_nam, sfile_spac, pfile_asia, pfile_nam, pfile_spac, mapfile, out_filename):
    # Read landseamask
    ds = xr.open_dataset(mapfile)
    # The landseamask has 2 extra grid point at the beginning and end than the actual IMERG data, 
    # presumable to make it easier to wrap around the 0 deg longitude
    # Remove them so that it matches with the IMERG data
    ds = ds.isel(lon=slice(1,-1))

    # Change longitude from [0,360] to [-180,180], then roll to center over 0 ()
    ds = ds.assign_coords(lon=((ds.lon - 180) % 360) - 180).roll(lon=(ds.dims['lon'] // 2), roll_coords='lon')
    # Subset latitude range within [-60, 60]
    ds = ds.where((ds.lat > -60.) & (ds.lat < 60.), drop=True)
    lon_G = ds.lon
    lat_G = ds.lat

    # Read stats files
    dss_asia = xr.open_dataset(sfile_asia, decode_times=True)
    dss_nam = xr.open_dataset(sfile_nam, decode_times=False)
    dss_spac = xr.open_dataset(sfile_spac, decode_times=False)

    # Read rain files
    dsp_asia = xr.open_dataset(pfile_asia, decode_times=False)
    dsp_nam = xr.open_dataset(pfile_nam, decode_times=False)
    dsp_spac = xr.open_dataset(pfile_spac, decode_times=False)

    # asia (Asia-Pacific)
    Glon_Lidx_asia = np.argmin(np.absolute(lon_G.data - (np.min(dss_asia.lon.data) + 5)))
    Glon_Ridx_asia = np.argmin(np.absolute(lon_G.data - np.max(dss_asia.lon.data)))
    Glat_Bidx_asia = np.argmin(np.absolute(lat_G.data - np.min(dss_asia.lat.data)))
    Glat_Tidx_asia = np.argmin(np.absolute(lat_G.data - np.max(dss_asia.lat.data)))
    nx1 = Glon_Ridx_asia - Glon_Lidx_asia + 1
    ny1 = Glat_Tidx_asia - Glat_Bidx_asia + 1

    # NAM (Europe-North America)
    Glon_Lidx_nam = np.argmin(np.absolute(lon_G.data - (np.min(dss_nam.lon.data))))
    Glon_Ridx_nam = np.argmin(np.absolute(lon_G.data - (np.max(dss_nam.lon.data) - 5)))
    Glat_Bidx_nam = np.argmin(np.absolute(lat_G.data - (np.min(dss_nam.lat.data) + 5)))
    Glat_Tidx_nam = np.argmin(np.absolute(lat_G.data - (np.max(dss_nam.lat.data))))
    nx3 = Glon_Ridx_nam - Glon_Lidx_nam + 1
    ny3 = Glat_Tidx_nam - Glat_Bidx_nam + 1

    # SPAC
    Glon_Lidx_spac = np.argmin(np.absolute(lon_G.data - (np.min(dss_spac.lon.data))))
    Glon_Ridx_spac = np.argmin(np.absolute(lon_G.data - (np.max(dss_spac.lon.data) - 5)))
    Glat_Bidx_spac = np.argmin(np.absolute(lat_G.data - np.min(dss_spac.lat.data)))
    Glat_Tidx_spac = np.argmin(np.absolute(lat_G.data - (np.max(dss_spac.lat.data) - 5)))
    nx6 = Glon_Ridx_spac - Glon_Lidx_spac + 1
    ny6 = Glat_Tidx_spac - Glat_Bidx_spac + 1

    # Find regional index
    # asia
    lon_Lidx_asia = np.argmin(np.absolute(dss_asia.lon.data - (np.min(dss_asia.lon.data) + 5)))
    lon_Ridx_asia = np.argmin(np.absolute(dss_asia.lon.data - (np.max(dss_asia.lon.data))))
    lat_Bidx_asia = np.argmin(np.absolute(dss_asia.lat.data - (np.min(dss_asia.lat.data))))
    lat_Tidx_asia = np.argmin(np.absolute(dss_asia.lat.data - (np.max(dss_asia.lat.data))))
    # NAM 
    lon_Lidx_nam = np.argmin(np.absolute(dss_nam.lon.data - (np.min(dss_nam.lon.data))))
    lon_Ridx_nam = np.argmin(np.absolute(dss_nam.lon.data - (np.max(dss_nam.lon.data) - 5)))
    lat_Bidx_nam = np.argmin(np.absolute(dss_nam.lat.data - (np.min(dss_nam.lat.data) + 5)))
    lat_Tidx_nam = np.argmin(np.absolute(dss_nam.lat.data - (np.max(dss_nam.lat.data))))
    # SPAC
    lon_Lidx_spac = np.argmin(np.absolute(dss_spac.lon.data - (np.min(dss_spac.lon.data))))
    lon_Ridx_spac = np.argmin(np.absolute(dss_spac.lon.data - (np.max(dss_spac.lon.data) - 5)))
    lat_Bidx_spac = np.argmin(np.absolute(dss_spac.lat.data - (np.min(dss_spac.lat.data))))
    lat_Tidx_spac = np.argmin(np.absolute(dss_spac.lat.data - (np.max(dss_spac.lat.data) - 5)))

    # Set up lat/lon indices lists (these are the same for all variables)
    GlatB = {"asia":Glat_Bidx_asia, "nam":Glat_Bidx_nam, "spac":Glat_Bidx_spac}
    GlatT = {"asia":Glat_Tidx_asia, "nam":Glat_Tidx_nam, "spac":Glat_Tidx_spac}
    GlonL = {"asia":Glon_Lidx_asia, "nam":Glon_Lidx_nam, "spac":Glon_Lidx_spac}
    GlonR = {"asia":Glon_Ridx_asia, "nam":Glon_Ridx_nam, "spac":Glon_Ridx_spac}
    latB = {"asia":lat_Bidx_asia, "nam":lat_Bidx_nam, "spac":lat_Bidx_spac}
    latT = {"asia":lat_Tidx_asia, "nam":lat_Tidx_nam, "spac":lat_Tidx_spac}
    lonL = {"asia":lon_Lidx_asia, "nam":lon_Lidx_nam, "spac":lon_Lidx_spac}
    lonR = {"asia":lon_Ridx_asia, "nam":lon_Ridx_nam, "spac":lon_Ridx_spac}
    nx_g = len(lon_G)
    ny_g = len(lat_G)

    # MCS number
    var_list = {"asia":dss_asia.mcs_number_pf.data, "nam":dss_nam.mcs_number_pf.data, "spac":dss_spac.mcs_number_pf.data, }
    mcs_number_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_number_ccs.data, "nam":dss_nam.mcs_number_ccs.data, "spac":dss_spac.mcs_number_ccs.data, }
    mcs_number_ccs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_nhour_ccs.data, "nam":dss_nam.mcs_nhour_ccs.data, "spac":dss_spac.mcs_nhour_ccs.data, }
    mcs_nhour_ccs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_nhour_pf.data, "nam":dss_nam.mcs_nhour_pf.data, "spac":dss_spac.mcs_nhour_pf.data, }
    mcs_nhour_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_nhour_speedmcs.data, "nam":dss_nam.mcs_nhour_speedmcs.data, "spac":dss_spac.mcs_nhour_speedmcs.data, }
    mcs_nhour_speedmcs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.lifetime_mean.data, "nam":dss_nam.lifetime_mean.data, "spac":dss_spac.lifetime_mean.data, }
    lifetime_mean = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.ccs_area_mean.data, "nam":dss_nam.ccs_area_mean.data, "spac":dss_spac.ccs_area_mean.data, }
    ccs_area_mean = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.pf_area_mean.data, "nam":dss_nam.pf_area_mean.data, "spac":dss_spac.pf_area_mean.data, }
    pf_area_mean = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.rainratemax_mean.data, "nam":dss_nam.rainratemax_mean.data, "spac":dss_spac.rainratemax_mean.data, }
    rainratemax_mean = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.initiation_ccs.data, "nam":dss_nam.initiation_ccs.data, "spac":dss_spac.initiation_ccs.data, }
    initiation_ccs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.pf_speed_mcs.data, "nam":dss_nam.pf_speed_mcs.data, "spac":dss_spac.pf_speed_mcs.data, }
    pf_speed_mcs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.pf_uspeed_mcs.data, "nam":dss_nam.pf_uspeed_mcs.data, "spac":dss_spac.pf_uspeed_mcs.data, }
    pf_uspeed_mcs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.pf_vspeed_mcs.data, "nam":dss_nam.pf_vspeed_mcs.data, "spac":dss_spac.pf_vspeed_mcs.data, }
    pf_vspeed_mcs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    # var_list = {"asia":dss_asia.mcs_number_pf.data, "nam":dss_nam.mcs_number_pf.data, "spac":dss_spac.mcs_number_pf.data, }
    # mcs_number_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    # var_list = {"asia":dss_asia.mcs_number_pf.data, "nam":dss_nam.mcs_number_pf.data, "spac":dss_spac.mcs_number_pf.data, }
    # mcs_number_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    # var_list = {"asia":dss_asia.mcs_number_pf.data, "nam":dss_nam.mcs_number_pf.data, "spac":dss_spac.mcs_number_pf.data, }
    # mcs_number_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    # var_list = {"asia":dss_asia.mcs_number_pf.data, "nam":dss_nam.mcs_number_pf.data, "spac":dss_spac.mcs_number_pf.data, }
    # mcs_number_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.precipitation.data, "nam":dsp_nam.precipitation.data, "spac":dsp_spac.precipitation.data, }
    precipitation = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation.data, "nam":dsp_nam.mcs_precipitation.data, "spac":dsp_spac.mcs_precipitation.data, }
    mcs_precipitation = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation_count.data, "nam":dsp_nam.mcs_precipitation_count.data, "spac":dsp_spac.mcs_precipitation_count.data, }
    mcs_precipitation_count = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    ntimes_asia = dsp_asia.ntimes.data
    ntimes_nam = dsp_nam.ntimes.data
    ntimes_spac = dsp_spac.ntimes.data
    # Calculate average number of observed times for the 3 regions
    ntimes_avg = np.nanmean([ntimes_asia, ntimes_nam, ntimes_spac])

    # Write output file
    basetime = dss_nam.time.data
    result = write_netcdf(out_filename, basetime, lat_G, lon_G, ntimes_asia, ntimes_nam, ntimes_spac, \
                            mcs_number_ccs, mcs_number_pf, mcs_nhour_ccs, mcs_nhour_pf, mcs_nhour_speedmcs, \
                            lifetime_mean, ccs_area_mean, pf_area_mean, rainratemax_mean, initiation_ccs, \
                            pf_speed_mcs, pf_uspeed_mcs, pf_vspeed_mcs, \
                            precipitation, mcs_precipitation, mcs_precipitation_count)
    # import pdb; pdb.set_trace()
    return result


def main():
    # Parallel options: 0-run serially, 1-run parallel (uses Dask)
    run_parallel = 1
    # Set up Dask
    n_workers = 12
    threads_per_worker = 1

    start_year = 2000
    end_year = 2019

    rootdir = '/global/cscratch1/sd/feng045/waccem/mcs_region/'
    dir_asia = f'{rootdir}asia/stats_ccs4_4h/monthly/'
    dir_nam = f'{rootdir}nam/stats_ccs4_4h/monthly/'
    dir_spac = f'{rootdir}spac/stats_ccs4_4h/monthly/'
    mapfile = '/global/project/projectdirs/m1867/zfeng/gpm/map_data/IMERG_land_sea_mask.nc'

    outdir = '/global/cscratch1/sd/feng045/waccem/mcs_region/global/monthly/'
    os.makedirs(outdir, exist_ok=True)

    print(os.path.abspath(sys.argv[0]))

    # Create a list of dates in the format 'yearmonths'
    dates = []
    for iy in range(start_year, end_year+1, 1):
        dates = dates + [d.strftime('%Y%m') for d in pd.date_range(start=f'{iy}-01', end=f'{iy}-12', freq='MS')]

    
    # Initialize dask
    if run_parallel == 1:
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker)
        client = Client(cluster)

    results = []
    for idate in dates:
        print(idate)
        sfile_asia = f'{dir_asia}mcs_statsmap_{idate}.nc'
        sfile_nam = f'{dir_nam}mcs_statsmap_{idate}.nc'
        sfile_spac = f'{dir_spac}mcs_statsmap_{idate}.nc'
        pfile_asia = f'{dir_asia}mcs_rainmap_{idate}.nc'
        pfile_nam = f'{dir_nam}mcs_rainmap_{idate}.nc'
        pfile_spac = f'{dir_spac}mcs_rainmap_{idate}.nc'

        sfile_asia_ok = os.path.isfile(sfile_asia)
        sfile_nam_ok = os.path.isfile(sfile_nam)
        sfile_spac_ok = os.path.isfile(sfile_spac)
        pfile_asia_ok = os.path.isfile(pfile_asia)
        pfile_nam_ok = os.path.isfile(pfile_nam)
        pfile_spac_ok = os.path.isfile(pfile_spac)

        if sfile_asia_ok & sfile_nam_ok & sfile_spac_ok & \
            pfile_asia_ok & pfile_nam_ok & pfile_spac_ok:
            # print('Found all files')
            out_filename = f'{outdir}mcs_statsmap_{idate}.nc'

            if run_parallel == 0:
                # Call function to merge regions
                result = combine_regions(sfile_asia, sfile_nam, sfile_spac, pfile_asia, pfile_nam, pfile_spac, mapfile, out_filename)
            elif run_parallel == 1:
                status = delayed(combine_regions)(sfile_asia, sfile_nam, sfile_spac, pfile_asia, pfile_nam, pfile_spac, mapfile, out_filename)
                results.append(status)
        else:
            print('Not all region files are found!')
            if sfile_asia_ok == False:
                print(f'Missing: {sfile_asia}')
            if sfile_nam_ok == False:
                print(f'Missing: {sfile_nam}')
            if sfile_spac_ok == False:
                print(f'Missing: {sfile_spac}')
            if pfile_asia_ok == False:
                print(f'Missing: {pfile_asia}')
            if pfile_nam_ok == False:
                print(f'Missing: {pfile_nam}')
            if pfile_spac_ok == False:
                print(f'Missing: {pfile_spac}')
        # import pdb; pdb.set_trace()

    if run_parallel == 1:
        # Collect results from Dask
        results = dask.compute(*results)


if __name__ == "__main__":
    main()