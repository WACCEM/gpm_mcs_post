"""
This script finds MCS track numbers that overlaps with TCs and saves them to a netCDF file.
The IBTrACS data need to be updated to include all years between 2000-2019.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "26-Mar-2020"

import numpy as np
import glob, os, sys
import xarray as xr
import time, datetime, calendar, pytz
import yaml
import dask
from dask import delayed
from dask.distributed import Client, LocalCluster

def get_mcs_tc_tracknumber(mcsfile, dstc):
    """
    Find MCS track numbers that overlap with TC.
    """

    # Read MCS pixel-file
    dsmcs = xr.open_dataset(mcsfile, decode_times=False, mask_and_scale=False)
    lon = dsmcs['lon']
    lat = dsmcs['lat']
    # Get MCS cloud mask
    cloudtracknumber = dsmcs.cloudtracknumber.squeeze().data
    # Get MCS time string (yyyymmddhhmm)
    # mcs_filename = mcsfile.split('/')[-1]
    # mcs_filename = os.path.basename(mcsfile)
    # mcstime = mcs_filename.split('_')[1] + mcs_filename.split('_')[2][0:4]
    mcstime = dsmcs['time'].data.item()
    tctime = dstc['base_time'].data

    # Find the closest time between TC and MCS
    tidx = np.argmin(np.abs(tctime - mcstime))

    # Create a 2D mask for TC
    mask = np.zeros((cloudtracknumber.shape), dtype='int')

    # If only 1 matched TC time is found
    if np.size(tidx) == 1:
        # Calculate the distance between MCS grid points to the TC center 
        # (Note only lon/lat DataArray works here, numpy ndarray does not work)
        d = (lat - dstc.lat.data[tidx])**2 + (lon - dstc.lon.data[tidx])**2
        # Create a circular mask within TC ROCI
        mask[d <= (dstc.roci.data[tidx]/110)**2] = 1
    # IF more than 1 matched TC time is found
    else:
        # Loop over each time
        for k in np.arange(np.size(tidx)):
            a = np.zeros((cloudtracknumber.shape), dtype='int')
            a = ((lat - dstc.lat.data[tidx][k])**2 + (lon - dstc.lon.data[tidx][k])**2) <= (dstc.roci.data[tidx][k]/110)**2
            # Add to the mask
            mask = a + mask
    
    # Find MCS mask pixels that overlaps with TC mask
    overlap = cloudtracknumber[(mask > 0)]
    # Get the unique MCS tracknumbers
    # mcs_tracknumber_tc = np.unique(overlap[~np.isnan(overlap)]).astype(np.int32)
    mcs_tracknumber_tc = np.unique(overlap[(overlap > 0)])

    return mcs_tracknumber_tc.tolist()


def write_netcdf(outfile, mcs_tracknumber, mcs_nhours):
    """
    Write output netCDF file.
    """
    # Get number of tracks
    ntracks = len(mcs_tracknumber)

    # Define output var lists
    varlist = {'mcs_tracknumber': (['tracks'], mcs_tracknumber), \
               'mcs_nhours': (['tracks'], mcs_nhours),            }
    coordlist = {'tracks': (['tracks'], np.arange(0, ntracks, 1))}
    attrlist = {'title': 'MCS tracknumbers near Tropical Cyclones', \
                'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
                'created_on': time.ctime(time.time())}

    # Define Xarray Dataset
    dsout = xr.Dataset(varlist, coords=coordlist, attrs=attrlist)

    # Define variable attributes
    dsout.mcs_tracknumber.attrs['long_name'] = 'MCS tracknumbers that overlap with TC'
    dsout.mcs_tracknumber.attrs['comments'] = 'Subtract tracknumbers by 1 to match with MCS statistics file track indices'
    dsout.mcs_tracknumber.attrs['units'] = 'unitless'

    dsout.mcs_nhours.attrs['long_name'] = 'Number of hours each MCS overlap with TC'
    dsout.mcs_nhours.attrs['units'] = 'hours'

    # Write to netCDF file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')

    return


if __name__ == '__main__':


    indates = sys.argv[1]
    config_file = sys.argv[2]
    
    inyear = indates[0:4]

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    pixel_dir = config['pixelfile_dir']
    output_dir = config['output_dir']
    tc_file = config['tc_file']
    run_parallel = config['run_parallel']
    n_workers = config['n_workers']

    # Output file
    outfile = f'{output_dir}mcs_tc_tracknumbers_{indates}.nc'
    os.makedirs(output_dir, exist_ok=True)

    begin_time = datetime.datetime.now()

    # Find all pixel-level files
    pixelfiles = sorted(glob.glob(f'{pixel_dir}{indates}/mcstrack_{inyear}????_????.nc'))
    print(f'Number of pixel files: {len(pixelfiles)}')

    # # Get the region lat/lon boundary
    # ds = xr.open_dataset(pixelfiles[0])
    # ny = ds.dims['lat']; nx = ds.dims['lon']
    # lon = ds.lon; lat = ds.lat
    # lonmin, lonmax = lon.min().data, lon.max().data
    # latmin, latmax = lat.min().data, lat.max().data
    # import pdb; pdb.set_trace()

    # Read TC data
    dstc = xr.open_dataset(tc_file, decode_times=False)
    num_tc = dstc.sizes['storms']
    storms = dstc['storms'].data
    tctime = dstc['base_time'].data

    # Subset the TC dataset to keep those within the region and in the same year
    # buffer_zone = 5  # [degree]
    # dstc = dstc.where((dstc.lon >= (lonmin - buffer_zone)) & (dstc.lon <= (lonmax + buffer_zone)) & 
    #                     (dstc.lat >= (latmin - buffer_zone)) & (dstc.lat <= (latmax + buffer_zone)) &
    #                     (dstc.year == int(inyear)), drop=True)

    if run_parallel == 1:
        # Set Dask temporary directory for workers
        dask_tmp_dir = '/tmp'
        dask.config.set({'temporary-directory': dask_tmp_dir})
        # Initialize dask
        cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
        client = Client(cluster)

    saved_tracknumber = []
    # Loop over each MCS pixel file
    for ifile in range(0, len(pixelfiles)):

        ifilename = pixelfiles[ifile]
        print(os.path.basename(ifilename))

        # Call function to find MCS tracknumber within TCs
        if run_parallel == 0:
            tracknumber = get_mcs_tc_tracknumber(ifilename, dstc)
        elif run_parallel == 1:
            tracknumber = delayed(get_mcs_tc_tracknumber)(ifilename, dstc)
        
        # Append to the tracknumber list
        saved_tracknumber.append(tracknumber)

    # Collect results from Dask
    if run_parallel == 1:
        results = dask.compute(*saved_tracknumber)
        # print("Done processing {} files.".format(len(results))
        # Flatten the append list
        tracknumber_list = [item for sublist in results for item in sublist]
    elif run_parallel == 0:
        tracknumber_list = [item for sublist in saved_tracknumber for item in sublist]

    # Get unique track numbers, counts. Counts mean number of hours each MCS overlap with AR
    mcs_tracknumber, mcs_nhours = np.unique(tracknumber_list, return_counts=True)
    ntracks = len(mcs_tracknumber)
    print(f'Found MCS near TC: {ntracks}')

    # Write output to file
    write_netcdf(outfile, mcs_tracknumber, mcs_nhours)

    # Print processing time
    print(datetime.datetime.now() - begin_time)
