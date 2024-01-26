"""
Calculate representative ERA5 environments within each tracked object mask and 
saves the output matching track statistics to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "01-25-2024"

import numpy as np
import xarray as xr
import pandas as pd
import sys, os
import copy
import time
import yaml
import warnings
import diag_functions as afwa
from joblib import Parallel, delayed

def calc_era5_prof(
        filenameZ, 
        filenameT, 
        filenameRH, 
        filenameU, 
        filenameV, 
        idx_track, 
        track_time, 
        level_lims,
    ):


    # Find unique track_times
    unique_track_time = np.unique(track_time)
    n_unique_maskfiles = len(unique_track_time)
    # Make a list for all tracking mask files for this day
    mask_filenames = []
    for ifile in range(n_unique_maskfiles):
        # Convert current track time to a datetime string matching the mask file format
        track_timestring = pd.to_datetime(unique_track_time[ifile]).strftime('%Y%m%d_%H%M%S')
        # Tracking mask filename
        mask_filename = f'{mask_dir}{trackmask_basename}{track_timestring}.nc'
        maskfile_exist = os.path.isfile(mask_filename)
        if maskfile_exist:
            mask_filenames.append(mask_filename)
        else:
            print(f'WARNING: tracking mask file is missing: {mask_filename}')

    # Read and combine all tracking mask files
    dsm = xr.open_mfdataset(mask_filenames, concat_dim='time', combine='nested')
    # Drop the 2D latitude/longitude variables
    dsm = dsm.drop_vars(['longitude','latitude'])
    # Rename the lat/lon dimensions to latitude/longitude to be consistent with ERA5 dimension names
    # This is needed in Xarray for masking between different DataSets
    dsm = dsm.rename({'lon':'longitude', 'lat':'latitude'})
    # Get the min/max lat/lon from the mask file
    lon_m = dsm['longitude']
    lat_m = dsm['latitude']
    min_lon, max_lon = lon_m.min().item(), lon_m.max().item()
    min_lat, max_lat = lat_m.min().item(), lat_m.max().item()
    # Get variables from mask file
    time_mask = dsm['time']
    mask_tracknumbers = dsm['cloudtracknumber']
    PW = dsm['PW']

    #print(f" track_time in calc   {track_time} ")

    # Read ERA5 data, subset domain to match mask file
    dseZ = xr.open_dataset(filenameZ).sel(latitude=slice(max_lat, min_lat), longitude=slice(min_lon, max_lon))
    dseT = xr.open_dataset(filenameT).sel(latitude=slice(max_lat, min_lat), longitude=slice(min_lon, max_lon))
    dseR = xr.open_dataset(filenameRH).sel(latitude=slice(max_lat, min_lat), longitude=slice(min_lon, max_lon))
    dseZsfc = xr.open_dataset(era5_terrain_file).sel(latitude=slice(max_lat, min_lat), longitude=slice(min_lon, max_lon))
    #dseU = xr.open_dataset(filenameU)
    #dseV = xr.open_dataset(filenameV)

    # Reverse vertical coordinate (era5 puts top of atmos at first k index)
    dseZ = dseZ.reindex(level=list(reversed(dseZ.level)))
    dseT = dseT.reindex(level=list(reversed(dseT.level)))
    dseR = dseR.reindex(level=list(reversed(dseR.level)))
    #dseU = dseU.reindex(level=list(reversed(dseU.level)))
    #dseV = dseV.reindex(level=list(reversed(dseV.level)))

    # Subset pressure levels (reduce unncessary levels in the stratosphere)
    dseZ = dseZ.sel(level=slice(max(level_lims), min(level_lims)))
    dseT = dseT.sel(level=slice(max(level_lims), min(level_lims)))
    dseR = dseR.sel(level=slice(max(level_lims), min(level_lims)))
    #dseU = dseU.sel(level=slice(max(level_lims), min(level_lims)))
    #dseV = dseV.sel(level=slice(max(level_lims), min(level_lims)))    

    # Get surface elevation [m]
    z_sfc = dseZsfc['Z'].squeeze() / 9.80665

    # grab coordinate frames
    # lat_e5 = dseZ.latitude
    # lon_e5 = dseZ.longitude
    time_e5 = dseZ.time
    #level = dseZ.level  #pressure level
    pressure = dseZ.level
    # Convert pressure unit to Pa
    pressure_pa = pressure.data * 100

    #U_attrs = {'long_name': dseU.U.attrs['long_name'], 'units': dseU.U.attrs['units'],}
    #V_attrs = {'long_name': dseV.V.attrs['long_name'], 'units': dseV.V.attrs['units'],}
    #T_attrs = {'long_name': dseT.T.attrs['long_name'], 'units': dseT.T.attrs['units'],}
    #R_attrs = {'long_name': dseR.R.attrs['long_name'], 'units': dseR.R.attrs['units'],}
    #Z_attrs = {'long_name': dseZ.Z.attrs['long_name'], 'units': dseZ.Z.attrs['units'],}

    # Number of tracks in the day
    # ntracks_day = 5
    ntracks_day = len(track_time)
    nlevels = len(pressure)
    
    # # Find level indices matching the limits
    # lev_idx = np.nonzero(np.in1d(dse.level.values, level.values))[0]
    # lev_idx0 = min(lev_idx)
    # lev_idx1 = max(lev_idx)

    # Make arrays to store output
    meanMUCAPE = np.full((ntracks_day), np.NaN, dtype=np.float32) 
    meanMUCIN = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #meanMULFC = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #meanMULCL = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #meanMUEL = np.full((ntracks_day), np.NaN, dtype=np.float32)        

    #medianMUCAPE = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #medianMUCIN = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #medianMULFC = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #medianMULCL = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #medianMUEL = np.full((ntracks_day), np.NaN, dtype=np.float32)

    maxMUCAPE = np.full((ntracks_day), np.NaN, dtype=np.float32)
    maxMUCIN = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #maxMULFC = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #maxMULCL = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #maxMUEL = np.full((ntracks_day), np.NaN, dtype=np.float32)   

    minMUCAPE = np.full((ntracks_day), np.NaN, dtype=np.float32)
    minMUCIN = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #minMULFC = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #minMULCL = np.full((ntracks_day), np.NaN, dtype=np.float32)
    #minMUEL = np.full((ntracks_day), np.NaN, dtype=np.float32)

    meanPW = np.full((ntracks_day), np.NaN, dtype=np.float32)
    maxPW = np.full((ntracks_day), np.NaN, dtype=np.float32)
    minPW = np.full((ntracks_day), np.NaN, dtype=np.float32)


    print(f' ntracks_day {ntracks_day} ')
    
    # Loop over each track
    for itrack in range(0, ntracks_day):
        print(f'{itrack}: {track_time[itrack]}')
        
        # Find closest grid point and time index in ERA5 
        # lat_idx = np.abs(lat_e5.values - mcs_lat[itrack]).argmin()
        # lon_idx = np.abs(lon_e5.values - mcs_lon[itrack]).argmin()
        t5_idx = np.abs(time_e5.values - track_time[itrack]).argmin()

        # Find closest time index in mask file
        tmask_idx = np.abs(time_mask.values - track_time[itrack]).argmin()

        # Current track number (+1 to match mask file)
        itrack_num = idx_track[itrack] + 1
        # Select matching time
        mask_time = mask_tracknumbers.isel(time=tmask_idx)
        # Get current mask
        _mask = mask_time.where(mask_time == itrack_num, other=0).load()

        # Select matching time from ERA5 DataSet
        iT = dseT.isel(time=t5_idx)
        iR = dseR.isel(time=t5_idx)
        iZ = dseZ.isel(time=t5_idx)
        # Subset ERA5 variables within the mask, 
        # setting drop=True crops the data outside the mask
        _tk = iT['T'].where(_mask == itrack_num, drop=True)
        _rh = iR['R'].where(_mask == itrack_num, drop=True)
        _z = iZ['Z'].where(_mask == itrack_num, drop=True) / 9.80665
        _z_sfc = z_sfc.where(_mask == itrack_num, drop=True)

        # Select matching time from mask DataSet
        iPW = PW.isel(time=tmask_idx)
        # Subset mask variables within the mask
        _PW = iPW.where(_mask == itrack_num, drop=True)
        # import pdb; pdb.set_trace()

        # Filter 3D variables below surface (replaced with -999)
        # They will not be included in searching for max ThetaE below 500 hPa for the most unstable level
        # Vertical coordinates (Z, P) are not filtered
        fillval = -999.
        _tk_f = _tk.where(_z > _z_sfc, other=fillval)
        _rh_f = _rh.where(_z > _z_sfc, other=fillval)
       
        # Make a 3D array of pressure by repeating the profile ny, nx times
        _nx = _tk.sizes['longitude']
        _ny = _tk.sizes['latitude']
        p3d = np.tile(pressure_pa, (_ny * _nx)).reshape(_nx, _ny, nlevels)
        # Swap dimension from 2 to 0 such that: [z, y, x]
        # p3d = np.swapaxes(p3d, 2, 0)
        # Transpose 3D pressure array to match the other variables: [z, y, x]
        p3d = np.transpose(p3d, (2,1,0))
        
        # Call AFWA diagnostics on data filtered below surface
        ostat, afwaMUCAPE, afwaMUCIN, afwaMULCL, afwaMULFC, afwaMUEL, afwaMULPL = afwa.diag_functions.diag_map(_tk_f, _rh_f, p3d, _z, 1, 1)
        if ostat == 1:
            # Use a context manager to temporarily suppress the warning
            with warnings.catch_warnings():
                # Filter out the specific warning
                warnings.filterwarnings("ignore", category=RuntimeWarning)
                # Replace undefined values with NaN
                afwaMUCAPE[afwaMUCAPE < 0] = np.NaN
                afwaMUCIN[afwaMUCIN < -999] = np.NaN
                afwaMULCL[afwaMULCL < 0] = np.NaN
                afwaMULFC[afwaMULFC < 0] = np.NaN
                afwaMUEL[afwaMUEL < 0] = np.NaN
                afwaMULPL[afwaMULPL < 0] = np.NaN
                
                # Calculate stats of afwa vars:
                meanMUCAPE[itrack] = np.nanmean(afwaMUCAPE)
                meanMUCIN[itrack] = np.nanmean(afwaMUCIN)

                maxMUCAPE[itrack] = np.nanmax(afwaMUCAPE)
                maxMUCIN[itrack]  = np.nanmax(afwaMUCIN)
                #maxMULFC[itrack]  = np.nanmax(afwaMULFC)
                #maxMULCL[itrack]  = np.nanmax(afwaMULCL)
                #maxMUEL[itrack]   = np.nanmax(afwaMUEL)  

                minMUCAPE[itrack] = np.nanmin(afwaMUCAPE)
                minMUCIN[itrack] = np.nanmin(afwaMUCIN)

                meanPW[itrack] = np.nanmean(_PW)
                maxPW[itrack] = np.nanmax(_PW)
                minPW[itrack] = np.nanmin(_PW)


        ## define attributes for output vars
        MUCAPE_attrs = {'long_name': 'MU CAPE in area', 'units': 'J/kg' }
        MUCIN_attrs  = {'long_name': 'MU CIN in area', 'units': 'J/kg' }
        MULCL_attrs  = {'long_name': 'MU LCL in area', 'units': 'm ASL'}
        MULFC_attrs  = {'long_name': 'MU LFC in area', 'units': 'm ASL'}
        MUEL_attrs   = {'long_name': 'MU EL in area', 'units': 'm ASL'}
        PW_attrs = PW.attrs

        #maxMUCAPE_attrs = {'long_name': 'max MU CAPE in area', 'units': 'J/kg' }
        #maxMUCIN_attrs  = {'long_name': 'max MU CIN (i.e., least oppressive CIN) in area', 'units': 'J/kg' }
        #maxMULCL_attrs  = {'long_name': 'max MU LCL in area', 'units': 'm ASL'}
        #maxMULFC_attrs  = {'long_name': 'max MU LFC in area', 'units': 'm ASL'}
        #maxMUEL_attrs   = {'long_name': 'max MU EL in area', 'units': 'm ASL'}

        #minMUCAPE_attrs = {'long_name': 'min MU CAPE in area', 'units': 'J/kg' }
        #minMUCIN_attrs  = {'long_name': 'min MU CIN (i.e., most oppressive CIN) in area', 'units': 'J/kg' }
        #minMULCL_attrs  = {'long_name': 'min MU LCL in area', 'units': 'm ASL'}
        #minMULFC_attrs  = {'long_name': 'min MU LFC in area', 'units': 'm ASL'}
        #minMUEL_attrs   = {'long_name': 'min MU EL in area', 'units': 'm ASL'}

        #meanMUCAPE_attrs = {'long_name': 'mean MU CAPE in area', 'units': 'J/kg' }
        #meanMUCIN_attrs  = {'long_name': 'mean MU CIN in area', 'units': 'J/kg' }
        #meanMULCL_attrs  = {'long_name': 'mean MU LCL in area', 'units': 'm ASL'}
        #meanMULFC_attrs  = {'long_name': 'mean MU LFC in area', 'units': 'm ASL'}
        #meanMUEL_attrs   = {'long_name': 'mean MU EL in area', 'units': 'm ASL'}

        #medianMUCAPE_attrs = {'long_name': 'median MU CAPE in area', 'units': 'J/kg' }
        #medianMUCIN_attrs  = {'long_name': 'median MU CIN in area', 'units': 'J/kg' }
        #medianMULCL_attrs  = {'long_name': 'median MU LCL in area', 'units': 'm ASL'}
        #medianMULFC_attrs  = {'long_name': 'median MU LFC in area', 'units': 'm ASL'}
        #medianMUEL_attrs   = {'long_name': 'median MU EL in area', 'units': 'm ASL'}

  
    # Put output variables to a dictionary for easier acceess
    var_dict = {
        'meanMUCAPE':meanMUCAPE,
        'maxMUCAPE':maxMUCAPE,
        'minMUCAPE':minMUCAPE,
        'meanMUCIN':meanMUCIN, 
        'maxMUCIN':maxMUCIN, 
        'minMUCIN':minMUCIN, 

        'meanPW':meanPW,
        'maxPW':maxPW,
        'minPW':minPW,
    }
    # Output variable attribute dictionary, the keys must match those in var_dict
    var_attrs = {
        'meanMUCAPE':MUCAPE_attrs,
        'maxMUCAPE':MUCAPE_attrs,
        'minMUCAPE':MUCAPE_attrs,
        'meanMUCIN':MUCIN_attrs,
        'maxMUCIN':MUCIN_attrs,
        'minMUCIN':MUCIN_attrs,

        'meanPW':PW_attrs,
        'maxPW':PW_attrs,
        'minPW':PW_attrs,
    } 

    print(f'Done processing: {filenameT}')
    
    # import pdb; pdb.set_trace()
    return var_dict, var_attrs



if __name__ == "__main__":

    ##basename = sys.argv[1]
    ##varname = sys.argv[2]
    ##config_file = sys.argv[3]
    ##track_period = sys.argv[4]

    # hardcode set inputs:
    
    ##stream = open(config_file, 'r')
    ##config = yaml.full_load(stream)
    era5_dir = '/pscratch/sd/j/jmarquis/era5data/'
    era5_terrain_file = '/pscratch/sd/j/jmarquis/era5data/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc'
    
    #for MCSs
    # stats_dir = '/pscratch/sd/j/jmarquis/ERA5_waccem/MCStracking/allstats/'
    # output_dir = '/pscratch/sd/j/jmarquis/ERA5_waccem/Bandpassed/MCSenvs/'
    # mcsfile_basename = 'mcs_tracks_final_'
    
    # for MPs
    stats_dir = '/pscratch/sd/j/jmarquis/ERA5_waccem/Bandpassed/vortstats/'
    mask_dir = '/pscratch/sd/j/jmarquis/ERA5_waccem/Bandpassed/vortracking/'
    output_dir = '/pscratch/sd/j/jmarquis/ERA5_waccem/Bandpassed/MPenvs/'
    trackstats_basename = 'trackstats_final_'
    trackmask_basename = 'vorbpf_tracks_morevar'

    track_period = '20040501.0000_20040831.2300'
    level_lims = [1000,100]  # 100mb, 1000mb
    # nlevels = config['nlevels']

    # Parallel setup
    run_parallel = 0
    n_workers = 64

    stats_file = f"{stats_dir}{trackstats_basename}{track_period}.nc"
    outfilename = f"{output_dir}mp_tracks_era5_envs_{track_period}.nc"
    os.makedirs(output_dir, exist_ok=True)

    print(f'stats_file:     {stats_file}')
    print(f'outfile:      {output_dir}{outfilename}')


    # Read track statistics file
    dsm = xr.open_dataset(stats_file)
    ntracks = dsm.sizes['tracks']
    ntimes = dsm.sizes['times']
    print(f'start ntimes:  {ntimes}')

    # rmcs_lat = dsm['meanlat']
    # rmcs_lon = dsm['meanlon']
    # # Check if longitude is [-180~+180], if so convert it to [0~360] to match ERA5
    # if np.nanmin(rmcs_lon) < 0:
    #     rmcs_lon = rmcs_lon % 360
    #     print('MCS longitudes are [-180~+180], converted to [0-360] to match ERA5.')

    # Get end times for all tracks
    tracks_basetime = dsm.base_time
    # Sum over time dimension for valid basetime indices, -1 to get the last valid time index for each track
    # This is the end time index of each track (i.e. +1 equals the lifetime of each track)
    end_time_idx = np.sum(np.isfinite(tracks_basetime), axis=1)-1
    # Apply fancy indexing to base_time: a tuple that indicates for each track, get the end time index
    end_basetime = tracks_basetime[(np.arange(0,ntracks), end_time_idx)]

    # Get the min/max of all base_times
    min_basetime = tracks_basetime.sel(times=0).min()
    max_basetime = end_basetime.max()
    
    # Make a date list that includes all tracks
    tracks_alldates = pd.date_range(start=min_basetime.values, end=max_basetime.values, freq='1D')
    # Convert all track times and ERA5 times to date strings
    tracks_dates = tracks_basetime.dt.strftime('%Y%m%d')
    # ERA5 data directory is organized by month (yyyymm)
    dirs_month = tracks_alldates.strftime("%Y%m")
    # 3D data is in daily files
    files3d_day = tracks_alldates.strftime("%Y%m%d")
    nfiles = len(tracks_alldates)

    # TODO: setting nfiles to a smaller number for testing
    nfiles = 2
    print(f"Total number of ERA5 files: {nfiles}")


    # Create a list to store matchindices for each ERA5 file
    trackindices_all = []
    timeindices_all = []
    results = []

    p_filenameZ = []
    p_filenameT = []
    p_filenameRH = []
    p_filenameU = []
    p_filenameV = []
    p_idx_track = []
    p_track_time = []
    p_level_lims = []

    # Loop over dates that have MCSs in them
    for ifile in range(nfiles):
        idir = dirs_month[ifile]
        iday = files3d_day[ifile]
        print(iday)
                
        # 3d era5 var files for the day: 
        filenameT = f"{era5_dir}{idir}/e5.oper.an.pl.128_130_t.ll025sc.{iday}00_{iday}23.nc"
        filenameRH = f"{era5_dir}{idir}/e5.oper.an.pl.128_157_r.ll025sc.{iday}00_{iday}23.nc"
        filenameU = f"{era5_dir}{idir}/e5.oper.an.pl.128_131_u.ll025uv.{iday}00_{iday}23.nc"
        filenameV = f"{era5_dir}{idir}/e5.oper.an.pl.128_132_v.ll025uv.{iday}00_{iday}23.nc"
        filenameZ = f"{era5_dir}{idir}/e5.oper.an.pl.128_129_z.ll025sc.{iday}00_{iday}23.nc"
        
        # These tracks use the same ERA5 file
        idx_track, idx_time = np.where(tracks_dates == files3d_day[ifile])

        # Save track/time indices for the current ERA5 file to the overall list
        trackindices_all.append(idx_track)
        timeindices_all.append(idx_time)

        # Get the track lat/lon/time values in the same day
        track_time = tracks_basetime.values[idx_track, idx_time]

        # check if there are MCSs at current era5 time:
        tracks_present =  len(track_time) 
        print(f'tracks_present   {tracks_present} ')

        # Put inputs in lists (for parallelization)
        p_filenameZ.append(filenameZ)
        p_filenameT.append(filenameT)
        p_filenameRH.append(filenameRH)
        p_filenameU.append(filenameU)
        p_filenameV.append(filenameV)
        p_idx_track.append(idx_track)
        p_track_time.append(track_time)
        p_level_lims.append(level_lims)

        # if mcs present at era5 time 
        if tracks_present > 0:
            # Serial
            if run_parallel == 0:
                # Call function to calculate statistics
                result = calc_era5_prof(
                    p_filenameZ[ifile], 
                    p_filenameT[ifile], 
                    p_filenameRH[ifile], 
                    p_filenameU[ifile], 
                    p_filenameV[ifile], 
                    p_idx_track[ifile], 
                    p_track_time[ifile], 
                    p_level_lims[ifile],
                )
                results.append(result)          
    # end for ifile in range(nfiles):

    # Parallel
    if run_parallel == 1:
        results = Parallel(n_jobs=n_workers)(
            delayed(calc_era5_prof)(
                p_filenameZ[ifile], 
                p_filenameT[ifile], 
                p_filenameRH[ifile], 
                p_filenameU[ifile], 
                p_filenameV[ifile], 
                p_idx_track[ifile], 
                p_track_time[ifile], 
                p_level_lims[ifile],
            ) for ifile in range(0, nfiles)
        )
    # Final returned outputs
    final_result = results

    #-------------------------------------------------------------------------
    # Collect returned data and write to output
    #-------------------------------------------------------------------------

    # Make a variable list and get attributes from one of the returned dictionaries
    # Loop over each return results till one that is not None
    counter = 0
    while counter < nfiles:
        if final_result[counter] is not None:
            var2d_names = list(final_result[counter][0].keys())
            var_attrs = final_result[counter][1]
            break
        counter += 1

    # Loop over variable list to create the dictionary entry
    out_dict = {}
    out_dict_attrs = {}
    for ivar in var2d_names:
        out_dict[ivar] = np.full((ntracks, ntimes), np.NaN, dtype=np.float32)
        out_dict_attrs[ivar] = var_attrs[ivar]

    # Put the results to output track stats variables
    # Loop over each file (parallel return results)
    for ifile in range(nfiles):
        # Get the return results for this pixel file
        iVAR2d = final_result[ifile][0]
        # iVAR3d = final_result[ifile][1]
        if iVAR2d is not None:
            trackindices = trackindices_all[ifile]
            timeindices = timeindices_all[ifile]
            # Loop over each variable and assign values to output dictionary
            for ivar in var2d_names:
                if iVAR2d[ivar].ndim == 1:
                    out_dict[ivar][trackindices,timeindices] = iVAR2d[ivar]

    # Define a dataset containing all variables
    var_dict = {}
    # Define output variable dictionary
    for key, value in out_dict.items():
        if value.ndim == 2:
            var_dict[key] = (['tracks', 'times'], value, out_dict_attrs[key])

    # Define coordinate list
    coord_dict = {
        'tracks': (['tracks'], np.arange(0, ntracks)),
        'times': (['times'], np.arange(0, ntimes)),
    }

    # Define global attributes
    gattr_dict = {
        'Title': 'Track environments',
        'Created_on': time.ctime(time.time()),
    }

    # Define output Xarray dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netcdf file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f'Output saved as: {outfilename}')
    # import pdb; pdb.set_trace()


   