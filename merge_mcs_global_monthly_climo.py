import numpy as np
import glob, os, sys
import xarray as xr
import time

def merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR):
    
    # Create a global variable
    var_merge = np.full((12, ny_g, nx_g), np.nan, dtype=float)

    # asia
    var_merge[:, GlatB["asia"]:GlatT["asia"]+1, GlonL["asia"]:GlonR["asia"]+1] = var_list["asia"][:, latB["asia"]:latT["asia"]+1, lonL["asia"]:lonR["asia"]+1]
    # NAM
    var_merge[:, GlatB["nam"]:GlatT["nam"]+1, GlonL["nam"]:GlonR["nam"]+1] = var_list["nam"][:, latB["nam"]:latT["nam"]+1, lonL["nam"]:lonR["nam"]+1]
    # SPAC
    var_merge[:, GlatB["spac"]:GlatT["spac"]+1, GlonL["spac"]:GlonR["spac"]+1] = var_list["spac"][:, latB["spac"]:latT["spac"]+1, lonL["spac"]:lonR["spac"]+1]
    
    return var_merge

def write_netcdf(out_filename,
                    lat_G,
                    lon_G,
                    months,
                    mcs_number,
                    mcs_lifetime,
                    mcs_totalrain,
                    mcs_totalrainheavy,
                    mcs_rainrateheavy,
                    mcs_rainratemax,
                    mcs_freq_ccs,
                    mcs_freq_pf,
                    mcs_initiation_ccs,
                    mcs_pfdiameter,
                    mcs_speed,
                    mcs_uspeed,
                    mcs_vspeed,
                    precipitation,
                    mcs_precipitation,
                    mcs_precipitation_freq,
                    mcs_precipitation_frac,
                    mcs_precipitation_intensity,
                    ):
    # Write monthly mean output file
    dsmap = xr.Dataset({'mcs_number': (['month', 'lat', 'lon'], mcs_number), \
                        'mcs_freq_ccs': (['month', 'lat', 'lon'], mcs_freq_ccs), \
                        'mcs_freq_pf': (['month', 'lat', 'lon'], mcs_freq_pf), \
                        'mcs_initiation_ccs': (['month', 'lat', 'lon'], mcs_initiation_ccs), \
                        'mcs_pfdiameter': (['month', 'lat', 'lon'], mcs_pfdiameter), \
                        'mcs_lifetime': (['month', 'lat', 'lon'], mcs_lifetime), \
                        'mcs_totalrain': (['month', 'lat', 'lon'], mcs_totalrain), \
                        'mcs_totalrainheavy': (['month', 'lat', 'lon'], mcs_totalrainheavy), \
                        'mcs_rainrateheavy': (['month', 'lat', 'lon'], mcs_rainrateheavy), \
                        'mcs_rainratemax': (['month', 'lat', 'lon'], mcs_rainratemax), \
                        'mcs_speed': (['month', 'lat', 'lon'], mcs_speed), \
                        'mcs_uspeed': (['month', 'lat', 'lon'], mcs_uspeed), \
                        'mcs_vspeed': (['month', 'lat', 'lon'], mcs_vspeed), \
                        'precipitation': (['month', 'lat', 'lon'], precipitation), \
                        'mcs_precipitation': (['month', 'lat', 'lon'], mcs_precipitation), \
                        'mcs_precipitation_frac': (['month', 'lat', 'lon'], mcs_precipitation_frac), \
                        'mcs_precipitation_freq': (['month', 'lat', 'lon'], mcs_precipitation_freq), \
                        'mcs_precipitation_intensity': (['month', 'lat', 'lon'], mcs_precipitation_intensity), \
                        }, \
                        coords={'month': (['month'], months), \
                                'lat': (['lat'], lat_G), \
                                'lon': (['lon'], lon_G)}, \
                        attrs={'title': 'MCS monthly mean climatology', \
                                'contact':'Zhe Feng, zhe.feng@pnnl.gov', \
                                'created_on':time.ctime(time.time()), \
                                'produced_with':os.path.abspath(sys.argv[0]),\
                                })

    dsmap.month.attrs['long_name'] = 'Month'
    dsmap.month.attrs['units'] = 'none'

    dsmap.lon.attrs['long_name'] = 'Longitude'
    dsmap.lon.attrs['units'] = 'degree'

    dsmap.lat.attrs['long_name'] = 'Latitude'
    dsmap.lat.attrs['units'] = 'degree'

    dsmap.mcs_number.attrs['long_name'] = 'Number of MCS'
    dsmap.mcs_number.attrs['units'] = 'count'

    dsmap.mcs_freq_ccs.attrs['long_name'] = 'MCS frequency (defined by CCS)'
    dsmap.mcs_freq_ccs.attrs['units'] = '%'

    dsmap.mcs_freq_pf.attrs['long_name'] = 'MCS frequency (defined by PF)'
    dsmap.mcs_freq_pf.attrs['units'] = '%'

    dsmap.mcs_initiation_ccs.attrs['long_name'] = 'MCS initiation count'
    dsmap.mcs_initiation_ccs.attrs['units'] = 'count'

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

    dsmap.precipitation.attrs['long_name'] = 'Total precipitation'
    dsmap.precipitation.attrs['units'] = 'mm/day'

    dsmap.mcs_precipitation.attrs['long_name'] = 'MCS precipitation'
    dsmap.mcs_precipitation.attrs['units'] = 'mm/day'

    dsmap.mcs_precipitation_frac.attrs['long_name'] = 'MCS precipitation fraction'
    dsmap.mcs_precipitation_frac.attrs['units'] = '%'

    dsmap.mcs_precipitation_freq.attrs['long_name'] = 'MCS precipitation frequency'
    dsmap.mcs_precipitation_freq.attrs['units'] = '%'

    dsmap.mcs_precipitation_intensity.attrs['long_name'] = 'MCS precipitation intensity'
    dsmap.mcs_precipitation_intensity.attrs['units'] = 'mm/hour'

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
    print(f'Reading input data ...')
    dss_asia = xr.open_dataset(sfile_asia, decode_times=True)
    dss_nam = xr.open_dataset(sfile_nam, decode_times=False)
    dss_spac = xr.open_dataset(sfile_spac, decode_times=False)

    # Read rain files
    dsp_asia = xr.open_dataset(pfile_asia, decode_times=False)
    dsp_nam = xr.open_dataset(pfile_nam, decode_times=False)
    dsp_spac = xr.open_dataset(pfile_spac, decode_times=False)

    print(f'Finished reading input data. Merging tiles ... ')
    
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

    # Merge tiles for each variable
    var_list = {"asia":dss_asia.mcs_number.data, "nam":dss_nam.mcs_number.data, "spac":dss_spac.mcs_number.data, }
    mcs_number = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_lifetime.data, "nam":dss_nam.mcs_lifetime.data, "spac":dss_spac.mcs_lifetime.data, }
    mcs_lifetime = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_totalrain.data, "nam":dss_nam.mcs_totalrain.data, "spac":dss_spac.mcs_totalrain.data, }
    mcs_totalrain = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_totalrainheavy.data, "nam":dss_nam.mcs_totalrainheavy.data, "spac":dss_spac.mcs_totalrainheavy.data, }
    mcs_totalrainheavy = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)
    
    var_list = {"asia":dss_asia.mcs_rainrateheavy.data, "nam":dss_nam.mcs_rainrateheavy.data, "spac":dss_spac.mcs_rainrateheavy.data, }
    mcs_rainrateheavy = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_rainratemax.data, "nam":dss_nam.mcs_rainratemax.data, "spac":dss_spac.mcs_rainratemax.data, }
    mcs_rainratemax = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)
    
    var_list = {"asia":dss_asia.mcs_freq_ccs.data, "nam":dss_nam.mcs_freq_ccs.data, "spac":dss_spac.mcs_freq_ccs.data, }
    mcs_freq_ccs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_freq_pf.data, "nam":dss_nam.mcs_freq_pf.data, "spac":dss_spac.mcs_freq_pf.data, }
    mcs_freq_pf = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_initiation_ccs.data, "nam":dss_nam.mcs_initiation_ccs.data, "spac":dss_spac.mcs_initiation_ccs.data, }
    mcs_initiation_ccs = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_pfdiameter.data, "nam":dss_nam.mcs_pfdiameter.data, "spac":dss_spac.mcs_pfdiameter.data, }
    mcs_pfdiameter = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_speed.data, "nam":dss_nam.mcs_speed.data, "spac":dss_spac.mcs_speed.data, }
    mcs_speed = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_uspeed.data, "nam":dss_nam.mcs_uspeed.data, "spac":dss_spac.mcs_uspeed.data, }
    mcs_uspeed = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dss_asia.mcs_vspeed.data, "nam":dss_nam.mcs_vspeed.data, "spac":dss_spac.mcs_vspeed.data, }
    mcs_vspeed = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.precipitation.data, "nam":dsp_nam.precipitation.data, "spac":dsp_spac.precipitation.data, }
    precipitation = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation.data, "nam":dsp_nam.mcs_precipitation.data, "spac":dsp_spac.mcs_precipitation.data, }
    mcs_precipitation = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation_freq.data, "nam":dsp_nam.mcs_precipitation_freq.data, "spac":dsp_spac.mcs_precipitation_freq.data, }
    mcs_precipitation_freq = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation_frac.data, "nam":dsp_nam.mcs_precipitation_frac.data, "spac":dsp_spac.mcs_precipitation_frac.data, }
    mcs_precipitation_frac = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    var_list = {"asia":dsp_asia.mcs_precipitation_intensity.data, "nam":dsp_nam.mcs_precipitation_intensity.data, "spac":dsp_spac.mcs_precipitation_intensity.data, }
    mcs_precipitation_intensity = merge_tiles(nx_g, ny_g, var_list, GlatB, GlatT, GlonL, GlonR, latB, latT, lonL, lonR)

    months = dss_spac.month

    # Write output file
    print(f'Writing output data ... ')
    result = write_netcdf(out_filename, 
                            lat_G,
                            lon_G,
                            months,
                            mcs_number,
                            mcs_lifetime,
                            mcs_totalrain,
                            mcs_totalrainheavy,
                            mcs_rainrateheavy,
                            mcs_rainratemax,
                            mcs_freq_ccs,
                            mcs_freq_pf,
                            mcs_initiation_ccs,
                            mcs_pfdiameter,
                            mcs_speed,
                            mcs_uspeed,
                            mcs_vspeed,
                            precipitation,
                            mcs_precipitation,
                            mcs_precipitation_freq,
                            mcs_precipitation_frac,
                            mcs_precipitation_intensity,
                            )
    # import pdb; pdb.set_trace()
    return result


def main():
    # period = '2000_2019'
    period = '2001_2019'
    # period = '2014_2019'

    rootdir = '/global/cscratch1/sd/feng045/waccem/mcs_region/'
    dir_asia = f'{rootdir}asia/stats_ccs4_4h/climo/'
    # dir_asia = f'{rootdir}apac/stats_ccs4_4h/climo/'
    dir_nam = f'{rootdir}nam/stats_ccs4_4h/climo/'
    dir_spac = f'{rootdir}spac/stats_ccs4_4h/climo/'
    mapfile = '/global/project/projectdirs/m1867/zfeng/gpm/map_data/IMERG_land_sea_mask.nc'

    outdir = '/global/cscratch1/sd/feng045/waccem/mcs_region/global/'
    out_filename = f'{outdir}mcs_monthly_mean_map_{period}.nc'
    os.makedirs(outdir, exist_ok=True)

    # Input file names
    sfile_asia = f'{dir_asia}mcs_statsmap_monthly_mean_{period}.nc'
    sfile_nam = f'{dir_nam}mcs_statsmap_monthly_mean_{period}.nc'
    sfile_spac = f'{dir_spac}mcs_statsmap_monthly_mean_{period}.nc'
    pfile_asia = f'{dir_asia}mcs_rainmap_monthly_mean_{period}.nc'
    pfile_nam = f'{dir_nam}mcs_rainmap_monthly_mean_{period}.nc'
    pfile_spac = f'{dir_spac}mcs_rainmap_monthly_mean_{period}.nc'

    # import pdb; pdb.set_trace()
    result = combine_regions(sfile_asia, sfile_nam, sfile_spac, pfile_asia, pfile_nam, pfile_spac, mapfile, out_filename)

if __name__ == "__main__":
    main()