"""
Maps IR+IMERG robust MCS tracking statistics back to pixel-level map and saves output as daily netCDF data.
This version does not calculate propagation speed variables.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
01/30/2020 - Written.
"""
__author__ = "Zhe Feng, zhe.feng@pnnl.gov"

import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset
import sys, glob, os.path
import time, datetime, calendar, pytz
from pytz import utc

def make_base_time(filename):
    """
    Create basetime from a filename.
    This is a much faster way, if basetime in the file is exactly the same as filename.
    """
    timestr = os.path.basename(filename)[-16:-3]
    year = int(timestr[0:4])
    month = int(timestr[4:6])
    day = int(timestr[6:8])
    hour = int(timestr[9:11])
    minute = int(timestr[11:13])
    second = 0
    bt = calendar.timegm(datetime.datetime(year, month, day, hour, minute, second, tzinfo=utc).timetuple())
    return bt, filename

def datetime_to_timestamp(dtime):
    """
    Convert datetime64 (when reading base_time from netCDF) to epoch time
    """
    if dtime.tolist() == None:
        return np.nan
    else:
        return dtime.tolist()/1e9


sdate = sys.argv[1]
edate = sys.argv[2]
year = int(sys.argv[3])
month = int(sys.argv[4])
day = int(sys.argv[5])
region = sys.argv[6]
#sdate = '20181001'
#edate = '20190310'
#year = 2018
#month = 10

stats_path = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/stats_ccs4_4h/'
pixel_path = f'/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/{region}/mcstracking_ccs4_4h/'

yearstr = str(year)
monthstr = str(month).zfill(2)
daystr = str(day).zfill(2)
# stats_file = f'{stats_path}robust_mcs_tracks_{sdate}_{edate}.nc'
stats_file = f'{stats_path}robust_mcs_tracks_extc_{sdate}_{edate}.nc'
# Open stats file
stats = xr.open_dataset(stats_file)
#stats.load()

pixel_files = sorted(glob.glob(f'{pixel_path}{sdate}_{edate}/mcstrack_*nc'))
print(f'Found {len(pixel_files)} mcstrack files')

output_dir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/daily/'
output_file = f'{output_dir}mcs_statsmap_{yearstr}{monthstr}{daystr}.nc'
os.makedirs(output_dir, exist_ok=True)

# Get map dimensions
dspix = xr.open_dataset(pixel_files[0])
nx = dspix.dims['lon']
ny = dspix.dims['lat']
lon = dspix.lon
lat = dspix.lat

# Set up a dictionary that keyed by base_time, valued by filename.
# This is a clever way to move between time and pixel files
time_to_pixfilename = {key:value for key, value in map(make_base_time, pixel_files)}


# Get initial time for each track
inittime = stats.base_time.isel(times=0).load()
# Get the month
initmonth = inittime.dt.month
initday = inittime.dt.day
# Get max track length
ntimes = stats.dims['times']

# Select track indices with initial time within the current month
trackidx = stats.tracks.where((initmonth == month) & (initday == day), drop=True)
nmcs = len(trackidx)
print('Number of MCS: ', nmcs)
# Convert track indices to integer
trackidx = np.array(trackidx.values).astype(int)

# Subset variables to tracks in the current month
base_time = stats.base_time.sel(tracks=trackidx).load()
lifetime = (stats.length.sel(tracks=trackidx) * stats.time_resolution_hour).load()
tracks = stats.tracks.sel(tracks=trackidx).load()
ccs_area = stats.ccs_area.sel(tracks=trackidx).load()
pf_area = stats.pf_area.sel(tracks=trackidx, nmaxpf=0).load()
# mcs_status = stats.pf_mcsstatus.sel(tracks=trackidx).load()
starttrackresult = stats.starttrackresult.sel(tracks=trackidx).load()
pf_maxrainrate = stats.pf_maxrainrate.sel(tracks=trackidx).load()
total_rain = stats.total_rain.sel(tracks=trackidx).load()
total_heavyrain = stats.total_heavyrain.sel(tracks=trackidx).load()
rainrate_heavyrain = stats.rainrate_heavyrain.sel(tracks=trackidx).load()
# cloudnumber = stats.cloudnumber.sel(tracks=trackidx)
# speed = stats.movement_r_meters_per_second.sel(tracks=trackidx).load()
# # movement_x and movement_y are in [km]/timestep, divide by 3.6 (1000m/3600s) to convert to [m/s]
# uspeed = (stats.movement_storm_x.sel(tracks=trackidx)/3.6).load()
# vspeed = (stats.movement_storm_y.sel(tracks=trackidx)/3.6).load()


# Convert datetime64 (when reading base_time from netCDF) to epoch time
base_times = np.array([datetime_to_timestamp(t) for t in base_time.data.ravel()])
base_times = np.reshape(base_times, base_time.shape)

# Double check if datetime matches the filename keyed from the dictionary
#base_times.data[5,1], time_to_pixfilename[base_times[5,1]]


# Create variables for maps
nt_uniq = len(np.unique(base_times))
# map_lifetime = np.zeros((1, ny, nx))
map_lifetime_all = np.zeros((nmcs, ny, nx))*np.NAN

map_ccsarea = np.zeros((ny, nx))
map_pfarea = np.zeros((ny, nx))
map_rainrateheavy = np.zeros((ny, nx))
map_rainratemax = np.zeros((ny, nx))
map_totalrainheavy = np.zeros((ny, nx))
map_totalrain = np.zeros((ny, nx))
map_nhour_ccs = np.zeros((ny, nx))
map_nhour_pf = np.zeros((ny, nx))
map_nmcs_ccs = np.zeros((ny, nx))
map_nmcs_pf = np.zeros((ny, nx))
map_init_ccs = np.zeros((ny, nx))
# map_init_pf = np.zeros((ny, nx))
# map_genesis_ccs = np.zeros((ny, nx))
# map_genesis_pf = np.zeros((ny, nx))


nframes = 0

# Loop over each MCS track
for imcs in range(nmcs):
# for imcs in range(0, 2):
    itrack = tracks.values[imcs] + 1
    ilifetime = lifetime.values[imcs]
#     idxconvinit = lifecycle_index.values[imcs,0]
#     idxgenesis = lifecycle_index.values[imcs,1]
    istartstatus = starttrackresult.values[imcs]
    print(f'Track number: {itrack}')
    
    temp_track = np.zeros((ny, nx))*np.nan
    temp_nmcs_ccs = np.zeros((ny, nx))
    temp_nmcs_pf = np.zeros((ny, nx))
    
    # Loop over each time
    for it in range(ntimes):
#        imcsstatus = mcs_status.values[imcs,it]
#         ilifestage = lifecycle_stage.values[imcs,it]
#         icn = cloudnumber.values[imcs,it]
        
        if np.isnan(base_times[imcs, it]):
            continue
        else:
            pixfname = time_to_pixfilename.get(base_times[imcs, it], 'None')
        # Check to make sure the basetime key exist in the dictionary before proceeding
        if (pixfname != 'None'):
            # print(pixfname)
            
            # Read pixel data
            ds = Dataset(pixfname)
            cloudtracknumber = ds.variables['cloudtracknumber'][0,:,:]
            pcptracknumber = ds.variables['pcptracknumber'][0,:,:]
#             fcloudnumber = ds.variables['cloudnumber'][0,:,:]
            # Find cloudnumber matching the current one
#             idx_c = np.where(fcloudnumber == icn)
            idx_c = np.where(cloudtracknumber == itrack)
            idx_p = np.where(pcptracknumber == itrack)
            # Check the number of pixels found for this track
            if (len(idx_c[0]) == 0):
                print(f'{pixfname}')
                print(f'Warning: track #{itrack} has no matching pixel found! Something is not right.')
            
            # Set each frame with a number
            # This is suitable for variables that change with time. e.g. PF area
            # Get values from stats variables
            iccsarea = ccs_area.values[imcs, it]
            ipfarea = pf_area.values[imcs, it]
            itotalrain = total_rain.values[imcs, it]
            iheavyrain = total_heavyrain.values[imcs, it]
            irainrateheavy = rainrate_heavyrain.values[imcs, it]
            # Get max rain rate from all PFs
            imaxrainrate = np.nanmax(pf_maxrainrate.values[imcs, it, :])

            # Put stats value onto the map
            temp_c = np.zeros((ny, nx))*np.NAN
            temp_p = np.zeros((ny, nx))*np.NAN
            temp_totalrain = np.full((ny, nx), np.NAN)
            temp_totalrainheavy = np.full((ny, nx), np.NAN)
            temp_rainrateheavy = np.full((ny, nx), np.NAN)
            temp_rainratemax = np.full((ny, nx), np.NAN)
            temp_pfspeed = np.full((ny, nx), np.NAN)
            temp_uspeed = np.full((ny, nx), np.NAN)
            temp_vspeed = np.full((ny, nx), np.NAN)

            temp_c[idx_c] = iccsarea
            temp_p[idx_p] = ipfarea
            temp_totalrain[idx_p] = itotalrain
            temp_totalrainheavy[idx_p] = iheavyrain
            temp_rainrateheavy[idx_p] = irainrateheavy
            temp_rainratemax[idx_p] = imaxrainrate

            # If MCS start status is not a split, check for initiation/genesis
            if (istartstatus != 13):
                # If time step is == 0 (initiation)
                if (it == 0):
                    map_init_ccs[idx_c] += 1

            # Make a temporary copy of the output array
            temp_map_ccsarea = np.copy(map_ccsarea)
            temp_map_pfarea = np.copy(map_pfarea)

            temp_map_totalrain = np.copy(map_totalrain)
            temp_map_totalrainheavy = np.copy(map_totalrainheavy)
            temp_map_rainrateheavy = np.copy(map_rainrateheavy)
            temp_map_rainratemax = np.copy(map_rainratemax)

            # Stack two arrays and use nansum, then save it back to output array
            map_ccsarea[:,:] = np.nansum(np.dstack((temp_map_ccsarea, temp_c)), axis=2)
            map_pfarea[:,:] = np.nansum(np.dstack((temp_map_pfarea, temp_p)), axis=2)

            map_totalrain[:,:] = np.nansum(np.dstack((temp_map_totalrain, temp_totalrain)), axis=2)
            map_totalrainheavy[:,:] = np.nansum(np.dstack((temp_map_totalrainheavy, temp_totalrainheavy)), axis=2)
            map_rainrateheavy[:,:] = np.nansum(np.dstack((temp_map_rainrateheavy, temp_rainrateheavy)), axis=2)
            map_rainratemax[:,:] = np.nansum(np.dstack((temp_map_rainratemax, temp_rainratemax)), axis=2)

            map_nhour_ccs[idx_c] += 1
            map_nhour_pf[idx_p] += 1
            
            # Set all pixels to the track number, this is a mask for the entire duration of the MCS
            temp_track[idx_c] = itrack
            temp_nmcs_ccs[idx_c] = 1
            temp_nmcs_pf[idx_p] = 1

            nframes+=1
            
            ds.close()
    
    # Set the entire track map to a single value
    # This is suitable for mapping single value variables such as lifetime
    temp_track[np.where(temp_track == itrack)] = ilifetime
    
    # Save value to output map
    map_lifetime_all[imcs,:,:] = temp_track[:,:]
    # Add the count for the number of MCS 
    map_nmcs_ccs[:,:] += temp_nmcs_ccs[:,:]
    map_nmcs_pf[:,:] += temp_nmcs_pf[:,:]

percentiles = [50,75,90,95]
map_lifetime_pts = np.nanpercentile(map_lifetime_all, percentiles, axis=0)

# Calculate conditional mean (divide sum by total number of hours at each pixel)
map_ccsarea_avg = map_ccsarea / map_nhour_ccs
map_pfarea_avg = map_pfarea / map_nhour_pf

map_totalrain_avg = map_totalrain / map_nhour_pf
map_totalrainheavy_avg = map_totalrainheavy / map_nhour_pf
map_rainrateheavy_avg = map_rainrateheavy / map_nhour_pf
map_rainratemax_avg = map_rainratemax / map_nhour_pf

map_lifetime_avg = np.nanmean(map_lifetime_all, axis=0)


# Compute Epoch Time for the month
months = np.zeros(1, dtype=int)
months[0] = calendar.timegm(datetime.datetime(int(year), int(month), int(day), 0, 0, 0, tzinfo=pytz.UTC).timetuple())

# Define Xarray Dataset
dsmap = xr.Dataset({'mcs_number_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_nmcs_ccs, 0)), \
                    'mcs_number_pf': (['time', 'lat', 'lon'], np.expand_dims(map_nmcs_pf, 0)), \
                    'mcs_nhour_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_ccs, 0)), \
                    'mcs_nhour_pf': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_pf, 0)), \
                    'lifetime_mean': (['time', 'lat', 'lon'], np.expand_dims(map_lifetime_avg, 0)), \
                    'ccs_area_mean': (['time', 'lat', 'lon'], np.expand_dims(map_ccsarea_avg, 0)),  \
                    'pf_area_mean': (['time', 'lat', 'lon'], np.expand_dims(map_pfarea_avg, 0)), \

                    'totalrain_mean': (['time', 'lat', 'lon'], np.expand_dims(map_totalrain_avg, 0)), \
                    'totalrainheavy_mean': (['time', 'lat', 'lon'], np.expand_dims(map_totalrainheavy_avg, 0)), \
                    'rainrateheavy_mean': (['time', 'lat', 'lon'], np.expand_dims(map_rainrateheavy_avg, 0)), \
                    'rainratemax_mean': (['time', 'lat', 'lon'], np.expand_dims(map_rainratemax_avg, 0)), \

                    'initiation_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_init_ccs, 0)), \
#                     'initiation_pf': (['time', 'lat', 'lon'], np.expand_dims(map_init_pf, 0)), \
#                     'genesis_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_genesis_ccs, 0)), \
#                     'genesis_pf': (['time', 'lat', 'lon'], np.expand_dims(map_genesis_pf, 0)), \
                    'lifetime_pt': (['time', 'percentiles', 'lat', 'lon'], np.expand_dims(map_lifetime_pts, 0)), \
                    'ntimes': (['time'], xr.DataArray(nframes).expand_dims('time', axis=0)), \
                   }, \
                    coords={'lon': (['lon'], lon), \
                            'lat': (['lat'], lat), \
                            'percentiles': (['percentiles'], percentiles), \
                            'time': (['time'], months), \
                            }, \
                    attrs={'title': 'MCS daily statistics map', \
                           'total_number_of_times': nframes, \
                           'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
                           'created_on': time.ctime(time.time())})

dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
dsmap.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'

dsmap.lat.attrs['long_name'] = 'Latitude'
dsmap.lat.attrs['units'] = 'degree'

dsmap.lon.attrs['long_name'] = 'Longitude'
dsmap.lon.attrs['units'] = 'degree'

dsmap.ntimes.attrs['long_name'] = 'Number of times with MCS'
dsmap.ntimes.attrs['units'] = 'count'

# Frequency variables
dsmap.mcs_number_ccs.attrs['long_name'] = 'Number of MCS defined by cold cloud shield'
dsmap.mcs_number_ccs.attrs['units'] = 'count'

dsmap.mcs_number_pf.attrs['long_name'] = 'Number of MCS defined by precipitation feature'
dsmap.mcs_number_pf.attrs['units'] = 'count'

dsmap.mcs_nhour_ccs.attrs['long_name'] = 'Number of MCS hours defined by cold cloud shield'
dsmap.mcs_nhour_ccs.attrs['units'] = 'hour'

dsmap.mcs_nhour_pf.attrs['long_name'] = 'Number of MCS hours defined by precipitation feature'
dsmap.mcs_nhour_pf.attrs['units'] = 'hour'

# Mean variables
dsmap.lifetime_mean.attrs['long_name'] = 'MCS lifetime mean'
dsmap.lifetime_mean.attrs['units'] = 'hour'

dsmap.ccs_area_mean.attrs['long_name'] = 'MCS cold cloud shield area mean'
dsmap.ccs_area_mean.attrs['units'] = 'km2'

dsmap.pf_area_mean.attrs['long_name'] = 'MCS precipitation feature area mean'
dsmap.pf_area_mean.attrs['units'] = 'km2'

dsmap.totalrain_mean.attrs['long_name'] = 'MCS total precipitation mean'
dsmap.totalrain_mean.attrs['units'] = 'mm'

dsmap.totalrainheavy_mean.attrs['long_name'] = 'MCS total heavy precipitation (rain rate > 10 mm/h) mean'
dsmap.totalrainheavy_mean.attrs['units'] = 'mm'

dsmap.rainrateheavy_mean.attrs['long_name'] = 'MCS heavy rain rate (rain rate > 10 mm/h) mean'
dsmap.rainrateheavy_mean.attrs['units'] = 'mm/h'

dsmap.rainratemax_mean.attrs['long_name'] = 'MCS max rain rate mean'
dsmap.rainratemax_mean.attrs['units'] = 'mm/h'

dsmap.initiation_ccs.attrs['long_name'] = 'MCS convective initiation hours defined by cold cloud shield'
dsmap.initiation_ccs.attrs['units'] = 'hour'

# dsmap.initiation_pf.attrs['long_name'] = 'MCS convective initiation hours defined by precipitation feature'
# dsmap.initiation_pf.attrs['units'] = 'hour'

# dsmap.genesis_ccs.attrs['long_name'] = 'MCS genesis hours defined by cold cloud shield'
# dsmap.genesis_ccs.attrs['units'] = 'hour'

# dsmap.genesis_pf.attrs['long_name'] = 'MCS genesis hours defined by precipitation feature'
# dsmap.genesis_pf.attrs['units'] = 'hour'

# Percentile variables
dsmap.lifetime_pt.attrs['long_name'] = 'MCS lifetime percentiles'
dsmap.lifetime_pt.attrs['units'] = 'hour'

fillvalue = np.nan
# Set encoding/compression for all variables
comp = dict(zlib=True, _FillValue=fillvalue, dtype='float32')
encoding = {var: comp for var in dsmap.data_vars}

dsmap.to_netcdf(path=output_file, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='time', \
                encoding=encoding)

print('Map output saved as: ', output_file)

