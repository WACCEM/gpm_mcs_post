"""
Maps IR+IMERG robust MCS tracking statistics back to pixel-level map and saves output as monthly netCDF data.

Author: Zhe Feng, zhe.feng@pnnl.gov
History:
04/04/2019 - Written.
06/18/2019 - Added saving lifetime_pt.
"""

import numpy as np
import xarray as xr
from netCDF4 import Dataset
import sys, glob, os.path
import time, datetime, calendar, pytz
from pytz import utc
import yaml

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
        # return dtime.tolist()/1e9
        return int(dtime.tolist()/1e9)


sdate = sys.argv[1]
edate = sys.argv[2]
year = int(sys.argv[3])
month = int(sys.argv[4])
# region = sys.argv[5]
config_file = sys.argv[5]
#sdate = '20181001'
#edate = '20190310'
#year = 2018
#month = 10

# get inputs from configuration file
stream = open(config_file, 'r')
config = yaml.full_load(stream)
# stats_file = config['stats_file']
stats_dir = config['stats_dir']
pixel_dir = config['pixelfile_dir']
output_monthly_dir = config['output_monthly_dir']
varname_speed = config['varname_speed']
varname_mov_x = config['varname_mov_x']
varname_mov_y = config['varname_mov_y']

yearstr = str(year)
monthstr = str(month).zfill(2)
stats_file = f'{stats_dir}mcs_tracks_final_{sdate}_{edate}.nc'
# stats_file = f'{stats_dir}robust_mcs_tracks_extc_{sdate}_{edate}.nc'
print(f'{stats_file}')
# Open stats file
stats = xr.open_dataset(stats_file)
stats.load()

pixel_files = sorted(glob.glob(f'{pixel_dir}/{sdate}_{edate}/mcstrack_*nc'))
print(f'Found {len(pixel_files)} mcstrack files')

output_file = f'{output_monthly_dir}mcs_statsmap_{yearstr}{monthstr}.nc'
os.makedirs(output_monthly_dir, exist_ok=True)

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
inittime = stats['base_time'].isel(times=0)
# Get the month
initmonth = inittime.dt.month
# Get max track length
ntimes = stats.dims['times']

# Select track indices with initial time within the current month
trackidx = stats['tracks'].where((initmonth == month), drop=True)
nmcs = len(trackidx)
print('Number of MCS: ', nmcs)
# Convert track indices to integer
trackidx = np.array(trackidx.values).astype(int)

# Subset variables to tracks in the current month
start_split_cloudnumber = stats['start_split_cloudnumber'].sel(tracks=trackidx).load()
base_time = stats['base_time'].sel(tracks=trackidx).load()
# lifetime = (stats['length'].sel(tracks=trackidx) * stats.time_resolution_hour).load()
lifetime = (stats['track_duration'].sel(tracks=trackidx) * stats.time_resolution_hour).load()
tracks = stats['tracks'].sel(tracks=trackidx).load()
ccs_area = stats['ccs_area'].sel(tracks=trackidx).load()
pf_area = stats['pf_area'].sel(tracks=trackidx, nmaxpf=0).load()
mcs_status = stats['pf_mcsstatus'].sel(tracks=trackidx).load()
# start_status = stats['starttrackresult'].sel(tracks=trackidx).load()
start_status = stats['start_status'].sel(tracks=trackidx).load()
pf_maxrainrate = stats['pf_maxrainrate'].sel(tracks=trackidx).load()
total_rain = stats['total_rain'].sel(tracks=trackidx).load()
total_heavyrain = stats['total_heavyrain'].sel(tracks=trackidx).load()
rainrate_heavyrain = stats['rainrate_heavyrain'].sel(tracks=trackidx).load()

# Check if speed variable exist in the track stats dataset
ds_varnames = list(stats.keys())
if varname_speed in ds_varnames:
    speed_exist = True
    # cloudnumber = stats.cloudnumber.sel(tracks=trackidx)
    speed = stats[varname_speed].sel(tracks=trackidx)
    # movement_x and movement_y are in [km]/timestep, divide by 3.6 (1000m/3600s) to convert to [m/s]
    uspeed = stats[varname_mov_x].sel(tracks=trackidx)/3.6
    vspeed = stats[varname_mov_y].sel(tracks=trackidx)/3.6
else:
    speed_exist = False
# import pdb; pdb.set_trace()


# Convert datetime64 (when reading base_time from netCDF) to epoch time
base_times = np.array([datetime_to_timestamp(t) for t in base_time.data.ravel()])
base_times = np.reshape(base_times, base_time.shape)
# import pdb; pdb.set_trace()

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
# map_nhour_ccdbz20 = np.zeros((ny, nx))
map_nmcs_ccs = np.zeros((ny, nx))
map_nmcs_pf = np.zeros((ny, nx))
map_init_ccs = np.zeros((ny, nx))
# map_init_pf = np.zeros((ny, nx))
# map_genesis_ccs = np.zeros((ny, nx))
# map_genesis_pf = np.zeros((ny, nx))
# map_corearea = np.zeros((ny, nx))
# map_coremaxis = np.zeros((ny, nx))
map_pfspeed = np.zeros((ny, nx))
map_uspeed = np.zeros((ny, nx))
map_vspeed = np.zeros((ny, nx))
map_pfspeed_mcs = np.zeros((ny, nx))
map_uspeed_mcs = np.zeros((ny, nx))
map_vspeed_mcs = np.zeros((ny, nx))
# map_pfspeed_mcs_ne = np.zeros((ny, nx))
# map_uspeed_mcs_ne = np.zeros((ny, nx))
# map_vspeed_mcs_ne = np.zeros((ny, nx))
# map_pfspeed_mcs_se = np.zeros((ny, nx))
# map_uspeed_mcs_se = np.zeros((ny, nx))
# map_vspeed_mcs_se = np.zeros((ny, nx))
# map_pfspeed_mcs_sw = np.zeros((ny, nx))
# map_uspeed_mcs_sw = np.zeros((ny, nx))
# map_vspeed_mcs_sw = np.zeros((ny, nx))
# map_pfspeed_mcs_nw = np.zeros((ny, nx))
# map_uspeed_mcs_nw = np.zeros((ny, nx))
# map_vspeed_mcs_nw = np.zeros((ny, nx))
map_nhour_speedmcs = np.zeros((ny, nx))
# map_nhour_speedmcs_ne = np.zeros((ny, nx))
# map_nhour_speedmcs_se = np.zeros((ny, nx))
# map_nhour_speedmcs_sw = np.zeros((ny, nx))
# map_nhour_speedmcs_nw = np.zeros((ny, nx))



nframes = 0

# Loop over each MCS track
for imcs in range(nmcs):

    itrack = tracks.values[imcs] + 1
    ilifetime = lifetime.values[imcs]
    istartstatus = start_status.values[imcs]
    istart_splitcloudnumber = start_split_cloudnumber[imcs]
    print(f'Track number: {itrack}')
    
    temp_track = np.zeros((ny, nx))*np.nan
    temp_nmcs_ccs = np.zeros((ny, nx))
    temp_nmcs_pf = np.zeros((ny, nx))
    
    # Loop over each time
    for it in range(ntimes):
        imcsstatus = mcs_status.values[imcs,it]
        
        if np.isnan(base_times[imcs, it]):
            continue
        else:
            pixfname = time_to_pixfilename.get(base_times[imcs, it], 'None')
        # Check to make sure the basetime key exist in the dictionary before proceeding
        if (pixfname != 'None'):
#             print(pixfname)
            
            # Read pixel data
            ds = Dataset(pixfname)
            cloudtracknumber = ds.variables['cloudtracknumber'][0,:,:]
            pcptracknumber = ds.variables['pcptracknumber'][0,:,:]
#             fcloudnumber = ds.variables['cloudnumber'][0,:,:]
            # Find cloudnumber matching the current one
#             idx_c = np.where(fcloudnumber == icn)
            # idx_c = np.where(cloudtracknumber == itrack)
            # idx_p = np.where(pcptracknumber == itrack)
            idx_c = cloudtracknumber == itrack
            idx_p = pcptracknumber == itrack
            npix_c = np.count_nonzero(idx_c)
            npix_p = np.count_nonzero(idx_p)
            
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

#             icorearea = core_area.values[imcs, it]
#             icoremaxis = core_maxis.values[imcs, it]


            # Put stats value onto the map
            temp_c = np.full((ny, nx), np.NAN)
            temp_p = np.full((ny, nx), np.NAN)
            temp_totalrain = np.full((ny, nx), np.NAN)
            temp_totalrainheavy = np.full((ny, nx), np.NAN)
            temp_rainrateheavy = np.full((ny, nx), np.NAN)
            temp_rainratemax = np.full((ny, nx), np.NAN)
            temp_pfspeed = np.full((ny, nx), np.NAN)
            temp_uspeed = np.full((ny, nx), np.NAN)
            temp_vspeed = np.full((ny, nx), np.NAN)

            temp_pfspeed = np.full((ny, nx), np.NAN)
            temp_uspeed = np.full((ny, nx), np.NAN)
            temp_vspeed = np.full((ny, nx), np.NAN)
            # temp_pfspeed_ne = np.full((ny, nx), np.NAN)
            # temp_uspeed_ne = np.full((ny, nx), np.NAN)
            # temp_vspeed_ne = np.full((ny, nx), np.NAN)
            # temp_pfspeed_se = np.full((ny, nx), np.NAN)
            # temp_uspeed_se = np.full((ny, nx), np.NAN)
            # temp_vspeed_se = np.full((ny, nx), np.NAN)
            # temp_pfspeed_sw = np.full((ny, nx), np.NAN)
            # temp_uspeed_sw = np.full((ny, nx), np.NAN)
            # temp_vspeed_sw = np.full((ny, nx), np.NAN)
            # temp_pfspeed_nw = np.full((ny, nx), np.NAN)
            # temp_uspeed_nw = np.full((ny, nx), np.NAN)
            # temp_vspeed_nw = np.full((ny, nx), np.NAN)

            temp_c[idx_c] = iccsarea
            temp_p[idx_p] = ipfarea
            
            # Assign single PF value to the entire PF mask
            # This is to increase area for certain extreme statistics (e.g. max rain rate)
            # to make results less noisy. Absolute location wise this is just an approximation.
            temp_totalrain[idx_p] = itotalrain
            temp_totalrainheavy[idx_p] = iheavyrain
            temp_rainrateheavy[idx_p] = irainrateheavy
            temp_rainratemax[idx_p] = imaxrainrate
            # The following gives more precise location mapping
            # Use PF mask as proxy location to map core values
            if speed_exist:
                ispeed = speed.values[imcs, it]
                iuspeed = uspeed.values[imcs, it]
                ivspeed = vspeed.values[imcs, it]
                temp_pfspeed[idx_p] = ispeed
                temp_uspeed[idx_p] = iuspeed
                temp_vspeed[idx_p] = ivspeed

            # # MCS propagating to NE
            # if ((iuspeed > 0) & (ivspeed > 0)):
            #     temp_pfspeed_ne[idx_p] = ispeed
            #     temp_uspeed_ne[idx_p] = iuspeed
            #     temp_vspeed_ne[idx_p] = ivspeed
            # # MCS propagating to SE
            # if ((iuspeed > 0) & (ivspeed < 0)):
            #     temp_pfspeed_se[idx_p] = ispeed
            #     temp_uspeed_se[idx_p] = iuspeed
            #     temp_vspeed_se[idx_p] = ivspeed
            # # MCS propagating to SW
            # if ((iuspeed < 0) & (ivspeed < 0)):
            #     temp_pfspeed_sw[idx_p] = ispeed
            #     temp_uspeed_sw[idx_p] = iuspeed
            #     temp_vspeed_sw[idx_p] = ivspeed
            # # MCS propagating to NW
            # if ((iuspeed < 0) & (ivspeed > 0)):
            #     temp_pfspeed_nw[idx_p] = ispeed
            #     temp_uspeed_nw[idx_p] = iuspeed
            #     temp_vspeed_nw[idx_p] = ivspeed
            
            # If MCS start status is not a split, check for initiation/genesis
            # if (istartstatus != 13):
            if (np.isnan(istart_splitcloudnumber)):
                # If time step is == 0 (initiation)
                if (it == 0):
                    map_init_ccs[idx_c] += 1

            # Keep speed if MCS status is met
            if (speed_exist == True) & (imcsstatus == 1):
                temp_map_pfspeed_mcs = np.copy(map_pfspeed_mcs)
                temp_map_uspeed_mcs = np.copy(map_uspeed_mcs)
                temp_map_vspeed_mcs = np.copy(map_vspeed_mcs)
                map_pfspeed_mcs[:,:] = np.nansum(np.dstack((temp_map_pfspeed_mcs, temp_pfspeed)), axis=2)
                map_uspeed_mcs[:,:] = np.nansum(np.dstack((temp_map_uspeed_mcs, temp_uspeed)), axis=2)
                map_vspeed_mcs[:,:] = np.nansum(np.dstack((temp_map_vspeed_mcs, temp_vspeed)), axis=2)
                map_nhour_speedmcs[idx_p] += 1

                # # MCS propagating to NE
                # if ((iuspeed > 0) & (ivspeed > 0)):
                #     temp_map_pfspeed_mcs = np.copy(map_pfspeed_mcs_ne)
                #     temp_map_uspeed_mcs = np.copy(map_uspeed_mcs_ne)
                #     temp_map_vspeed_mcs = np.copy(map_vspeed_mcs_ne)
                #     map_pfspeed_mcs_ne[:,:] = np.nansum(np.dstack((temp_map_pfspeed_mcs, temp_pfspeed_ne)), axis=2)
                #     map_uspeed_mcs_ne[:,:] = np.nansum(np.dstack((temp_map_uspeed_mcs, temp_uspeed_ne)), axis=2)
                #     map_vspeed_mcs_ne[:,:] = np.nansum(np.dstack((temp_map_vspeed_mcs, temp_vspeed_ne)), axis=2)
                #     map_nhour_speedmcs_ne[idx_p] += 1

                # # MCS propagating to SE
                # if ((iuspeed > 0) & (ivspeed < 0)):
                #     temp_map_pfspeed_mcs = np.copy(map_pfspeed_mcs_se)
                #     temp_map_uspeed_mcs = np.copy(map_uspeed_mcs_se)
                #     temp_map_vspeed_mcs = np.copy(map_vspeed_mcs_se)
                #     map_pfspeed_mcs_se[:,:] = np.nansum(np.dstack((temp_map_pfspeed_mcs, temp_pfspeed_se)), axis=2)
                #     map_uspeed_mcs_se[:,:] = np.nansum(np.dstack((temp_map_uspeed_mcs, temp_uspeed_se)), axis=2)
                #     map_vspeed_mcs_se[:,:] = np.nansum(np.dstack((temp_map_vspeed_mcs, temp_vspeed_se)), axis=2)
                #     map_nhour_speedmcs_se[idx_p] += 1

                # # MCS propagating to SW
                # if ((iuspeed < 0) & (ivspeed < 0)):
                #     temp_map_pfspeed_mcs = np.copy(map_pfspeed_mcs_sw)
                #     temp_map_uspeed_mcs = np.copy(map_uspeed_mcs_sw)
                #     temp_map_vspeed_mcs = np.copy(map_vspeed_mcs_sw)
                #     map_pfspeed_mcs_sw[:,:] = np.nansum(np.dstack((temp_map_pfspeed_mcs, temp_pfspeed_sw)), axis=2)
                #     map_uspeed_mcs_sw[:,:] = np.nansum(np.dstack((temp_map_uspeed_mcs, temp_uspeed_sw)), axis=2)
                #     map_vspeed_mcs_sw[:,:] = np.nansum(np.dstack((temp_map_vspeed_mcs, temp_vspeed_sw)), axis=2)
                #     map_nhour_speedmcs_sw[idx_p] += 1

                # # MCS propagating to NW
                # if ((iuspeed < 0) & (ivspeed > 0)):
                #     temp_map_pfspeed_mcs = np.copy(map_pfspeed_mcs_nw)
                #     temp_map_uspeed_mcs = np.copy(map_uspeed_mcs_nw)
                #     temp_map_vspeed_mcs = np.copy(map_vspeed_mcs_nw)
                #     map_pfspeed_mcs_nw[:,:] = np.nansum(np.dstack((temp_map_pfspeed_mcs, temp_pfspeed_nw)), axis=2)
                #     map_uspeed_mcs_nw[:,:] = np.nansum(np.dstack((temp_map_uspeed_mcs, temp_uspeed_nw)), axis=2)
                #     map_vspeed_mcs_nw[:,:] = np.nansum(np.dstack((temp_map_vspeed_mcs, temp_vspeed_nw)), axis=2)
                #     map_nhour_speedmcs_nw[idx_p] += 1
                    
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
            
            if speed_exist:
                temp_map_pfspeed = np.copy(map_pfspeed)
                temp_map_uspeed = np.copy(map_uspeed)
                temp_map_vspeed = np.copy(map_vspeed)
                map_pfspeed[:,:] = np.nansum(np.dstack((temp_map_pfspeed, temp_pfspeed)), axis=2)
                map_uspeed[:,:] = np.nansum(np.dstack((temp_map_uspeed, temp_uspeed)), axis=2)
                map_vspeed[:,:] = np.nansum(np.dstack((temp_map_vspeed, temp_vspeed)), axis=2)

            map_nhour_ccs[idx_c] += 1
            map_nhour_pf[idx_p] += 1
            # if (it > 2):
            #     import pdb; pdb.set_trace()
#             map_nhour_ccdbz20[idx_ccdbz20] += 1
            
            # Set all pixels to the track number, this is a mask for the entire duration of the MCS
#             temp_track[np.where(cloudtracknumber == itrack)] = itrack
#             temp_nmcs[np.where(cloudtracknumber == itrack)] = 1
            temp_track[idx_c] = itrack
            temp_nmcs_ccs[idx_c] = 1
            temp_nmcs_pf[idx_p] = 1

            nframes+=1
            
            ds.close()
        else:
            print(f'No pixel-file found: {base_times[imcs, it]}')
            # import pdb; pdb.set_trace()
    
    # Set the entire track map to a single value
    # This is suitable for mapping single value variables such as lifetime
    temp_track[np.where(temp_track == itrack)] = ilifetime
    
    # Save value to output map
#     map_lifetime[0,:,:] += temp_track[:,:]
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

map_pfspeed_avg = map_pfspeed / map_nhour_pf
map_pfspeed_mcs_avg = map_pfspeed_mcs / map_nhour_speedmcs
map_uspeed_avg = map_uspeed / map_nhour_pf
map_uspeed_mcs_avg = map_uspeed_mcs / map_nhour_speedmcs
map_vspeed_avg = map_vspeed / map_nhour_pf
map_vspeed_mcs_avg = map_vspeed_mcs / map_nhour_speedmcs

# map_uspeed_mcs_avg_ne = map_uspeed_mcs_ne / map_nhour_speedmcs_ne
# map_vspeed_mcs_avg_ne = map_vspeed_mcs_ne / map_nhour_speedmcs_ne
# map_pfspeed_mcs_avg_ne = map_pfspeed_mcs_ne / map_nhour_speedmcs_ne

# map_uspeed_mcs_avg_se = map_uspeed_mcs_se / map_nhour_speedmcs_se
# map_vspeed_mcs_avg_se = map_vspeed_mcs_se / map_nhour_speedmcs_se
# map_pfspeed_mcs_avg_se = map_pfspeed_mcs_se / map_nhour_speedmcs_se

# map_uspeed_mcs_avg_sw = map_uspeed_mcs_sw / map_nhour_speedmcs_sw
# map_vspeed_mcs_avg_sw = map_vspeed_mcs_sw / map_nhour_speedmcs_sw
# map_pfspeed_mcs_avg_sw = map_pfspeed_mcs_sw / map_nhour_speedmcs_sw

# map_uspeed_mcs_avg_nw = map_uspeed_mcs_nw / map_nhour_speedmcs_nw
# map_vspeed_mcs_avg_nw = map_vspeed_mcs_nw / map_nhour_speedmcs_nw
# map_pfspeed_mcs_avg_nw = map_pfspeed_mcs_nw / map_nhour_speedmcs_nw

map_lifetime_avg = np.nanmean(map_lifetime_all, axis=0)


# Compute Epoch Time for the month
months = np.zeros(1, dtype=int)
months[0] = calendar.timegm(datetime.datetime(int(year), int(month), 1, 0, 0, 0, tzinfo=pytz.UTC).timetuple())

# Define Xarray Dataset
dsmap = xr.Dataset({'mcs_number_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_nmcs_ccs, 0)), \
                    'mcs_number_pf': (['time', 'lat', 'lon'], np.expand_dims(map_nmcs_pf, 0)), \
                    'mcs_nhour_ccs': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_ccs, 0)), \
                    'mcs_nhour_pf': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_pf, 0)), \
                    'mcs_nhour_speedmcs': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_speedmcs, 0)), \
                    # 'mcs_nhour_speedmcs_ne': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_speedmcs_ne, 0)), \
                    # 'mcs_nhour_speedmcs_se': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_speedmcs_se, 0)), \
                    # 'mcs_nhour_speedmcs_sw': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_speedmcs_sw, 0)), \
                    # 'mcs_nhour_speedmcs_nw': (['time', 'lat', 'lon'], np.expand_dims(map_nhour_speedmcs_nw, 0)), \
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
                    'pf_speed': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_avg, 0)), \
                    'pf_speed_mcs': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_mcs_avg, 0)), \
                    'pf_uspeed': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_avg, 0)), \
                    'pf_uspeed_mcs': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg, 0)), \
                    'pf_vspeed': (['time', 'lat', 'lon'], np.expand_dims(map_vspeed_avg, 0)), \
                    'pf_vspeed_mcs': (['time', 'lat', 'lon'], np.expand_dims(map_vspeed_mcs_avg, 0)), \
                    # 'pf_speed_mcs_ne': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_mcs_avg_ne, 0)), \
                    # 'pf_uspeed_mcs_ne': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_ne, 0)), \
                    # 'pf_vspeed_mcs_ne': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_ne, 0)), \
                    # 'pf_speed_mcs_se': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_mcs_avg_se, 0)), \
                    # 'pf_uspeed_mcs_se': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_se, 0)), \
                    # 'pf_vspeed_mcs_se': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_se, 0)), \
                    # 'pf_speed_mcs_sw': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_mcs_avg_sw, 0)), \
                    # 'pf_uspeed_mcs_sw': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_sw, 0)), \
                    # 'pf_vspeed_mcs_sw': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_sw, 0)), \
                    # 'pf_speed_mcs_nw': (['time', 'lat', 'lon'], np.expand_dims(map_pfspeed_mcs_avg_nw, 0)), \
                    # 'pf_uspeed_mcs_nw': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_nw, 0)), \
                    # 'pf_vspeed_mcs_nw': (['time', 'lat', 'lon'], np.expand_dims(map_uspeed_mcs_avg_nw, 0)), \

                    'lifetime_pt': (['time', 'percentiles', 'lat', 'lon'], np.expand_dims(map_lifetime_pts, 0)), \
                    'ntimes': (['time'], xr.DataArray(nframes).expand_dims('time', axis=0)), \
                   }, \
                    coords={'lon': (['lon'], lon), \
                            'lat': (['lat'], lat), \
                            'percentiles': (['percentiles'], percentiles), \
                            'time': (['time'], months), \
                            }, \
                    attrs={'title': 'MCS monthly statistics map', \
                           'total_number_of_times': nframes, \
                           'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
                           'created_on': time.ctime(time.time())})

dsmap.time.attrs['long_name'] = 'Epoch Time (since 1970-01-01T00:00:00)'
dsmap.time.attrs['units'] = 'Seconds since 1970-1-1 0:00:00 0:00'

# dsmap.percentiles.attrs['long_name'] = 'Percentile values'
# dsmap.percentiles.attrs['units'] = 'unitless'

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

dsmap.mcs_nhour_speedmcs.attrs['long_name'] = 'Number of MCS hours for propagation speed (MCS status is met)'
dsmap.mcs_nhour_speedmcs.attrs['units'] = 'hour'

# dsmap.mcs_nhour_speedmcs_ne.attrs['long_name'] = 'Number of NE MCS hours for propagation speed (MCS status is met)'
# dsmap.mcs_nhour_speedmcs_ne.attrs['units'] = 'hour'

# dsmap.mcs_nhour_speedmcs_se.attrs['long_name'] = 'Number of SE MCS hours for propagation speed (MCS status is met)'
# dsmap.mcs_nhour_speedmcs_se.attrs['units'] = 'hour'

# dsmap.mcs_nhour_speedmcs_sw.attrs['long_name'] = 'Number of SW MCS hours for propagation speed (MCS status is met)'
# dsmap.mcs_nhour_speedmcs_sw.attrs['units'] = 'hour'

# dsmap.mcs_nhour_speedmcs_nw.attrs['long_name'] = 'Number of NW MCS hours for propagation speed (MCS status is met)'
# dsmap.mcs_nhour_speedmcs_nw.attrs['units'] = 'hour'

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

dsmap.pf_speed.attrs['long_name'] = 'Propagation speed for all times'
dsmap.pf_speed.attrs['units'] = 'm/s'

dsmap.pf_speed_mcs.attrs['long_name'] = 'Propagation speed when MCS status is met'
dsmap.pf_speed_mcs.attrs['units'] = 'm/s'

# dsmap.pf_speed_mcs_ne.attrs['long_name'] = 'NE Propagation speed when MCS status is met'
# dsmap.pf_speed_mcs_ne.attrs['units'] = 'm/s'

# dsmap.pf_speed_mcs_se.attrs['long_name'] = 'SE Propagation speed when MCS status is met'
# dsmap.pf_speed_mcs_se.attrs['units'] = 'm/s'

# dsmap.pf_speed_mcs_sw.attrs['long_name'] = 'SW Propagation speed when MCS status is met'
# dsmap.pf_speed_mcs_sw.attrs['units'] = 'm/s'

# dsmap.pf_speed_mcs_nw.attrs['long_name'] = 'NW Propagation speed when MCS status is met'
# dsmap.pf_speed_mcs_nw.attrs['units'] = 'm/s'

dsmap.pf_uspeed.attrs['long_name'] = 'Propagation speed (x-direction) for all times'
dsmap.pf_uspeed.attrs['units'] = 'm/s'

dsmap.pf_uspeed_mcs.attrs['long_name'] = 'Propagation speed (x-direction) when MCS status is met'
dsmap.pf_uspeed_mcs.attrs['units'] = 'm/s'

# dsmap.pf_uspeed_mcs_ne.attrs['long_name'] = 'NE Propagation speed (x-direction) when MCS status is met'
# dsmap.pf_uspeed_mcs_ne.attrs['units'] = 'm/s'

# dsmap.pf_uspeed_mcs_se.attrs['long_name'] = 'SE Propagation speed (x-direction) when MCS status is met'
# dsmap.pf_uspeed_mcs_se.attrs['units'] = 'm/s'

# dsmap.pf_uspeed_mcs_sw.attrs['long_name'] = 'SW Propagation speed (x-direction) when MCS status is met'
# dsmap.pf_uspeed_mcs_sw.attrs['units'] = 'm/s'

# dsmap.pf_uspeed_mcs_nw.attrs['long_name'] = 'NW Propagation speed (x-direction) when MCS status is met'
# dsmap.pf_uspeed_mcs_nw.attrs['units'] = 'm/s'

dsmap.pf_vspeed.attrs['long_name'] = 'Propagation speed (y-direction) for all times'
dsmap.pf_vspeed.attrs['units'] = 'm/s'

dsmap.pf_vspeed_mcs.attrs['long_name'] = 'Propagation speed (y-direction) when MCS status is met'
dsmap.pf_vspeed_mcs.attrs['units'] = 'm/s'

# dsmap.pf_vspeed_mcs_ne.attrs['long_name'] = 'NE Propagation speed (y-direction) when MCS status is met'
# dsmap.pf_vspeed_mcs_ne.attrs['units'] = 'm/s'

# dsmap.pf_vspeed_mcs_se.attrs['long_name'] = 'SE Propagation speed (y-direction) when MCS status is met'
# dsmap.pf_vspeed_mcs_se.attrs['units'] = 'm/s'

# dsmap.pf_vspeed_mcs_sw.attrs['long_name'] = 'SW Propagation speed (y-direction) when MCS status is met'
# dsmap.pf_vspeed_mcs_sw.attrs['units'] = 'm/s'

# dsmap.pf_vspeed_mcs_nw.attrs['long_name'] = 'NW Propagation speed (y-direction) when MCS status is met'
# dsmap.pf_vspeed_mcs_nw.attrs['units'] = 'm/s'


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

