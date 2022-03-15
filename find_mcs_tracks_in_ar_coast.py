"""
This script finds MCS track numbers that overlap with Atmospheric Rivers that hit the west coast of 
Europe, North and South America and save them to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "01-May-2020"

import numpy as np
import glob, os, sys
import xarray as xr
import time, datetime, calendar, pytz
import yaml
import dask
from dask import delayed
from dask.distributed import Client, LocalCluster
from scipy.ndimage import binary_dilation, generate_binary_structure

def make_timestring(filename):
    """
    Create timestring from a filename.
    """
    timestr = os.path.basename(filename)[-16:-3]
    year = (timestr[0:4])
    month = (timestr[4:6])
    day = (timestr[6:8])
    hour = (timestr[9:11])
    minute = (timestr[11:13])
    filetime = f'{year}-{month}-{day}T{hour}:{minute}'
    return filename, filetime

def get_mcs_ar_tracknumber(mcs_filename, ar_mask, landmask, bounds):

    # Get the domain subset bounds
    min_y = bounds['min_y']
    max_y = bounds['max_y']
    min_x = bounds['min_x']
    max_x = bounds['max_x']

    # Read MCS pixel file
    dspx = xr.open_dataset(mcs_filename, mask_and_scale=False)
    cloudtracknumber = dspx['cloudtracknumber'].squeeze().data[min_y:max_y, min_x:max_x]

    # Find MCS mask pixels that overlaps with AR mask
    overlap = cloudtracknumber[((ar_mask == 1) & (landmask == 1))]
    # Get the unique MCS tracknumbers
    mcs_tracknumber_ar = np.unique(overlap[(overlap > 0)])
    
    # Output as a list
    return mcs_tracknumber_ar.tolist()

def write_netcdf(outfile, mcs_tracknumber, mcs_nhours):
    # Get number of tracks
    ntracks = len(mcs_tracknumber)

    # Define output var lists
    varlist = {'mcs_tracknumber': (['tracks'], mcs_tracknumber), \
               'mcs_nhours': (['tracks'], mcs_nhours),            }
    coordlist = {'tracks': (['tracks'], np.arange(0, ntracks, 1))}
    attrlist = {'title': 'MCS tracknumbers near Atmospheric Rivers', \
                'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
                'created_on': time.ctime(time.time())}

    # Define Xarray Dataset
    dsout = xr.Dataset(varlist, coords=coordlist, attrs=attrlist)

    # Define variable attributes
    dsout.mcs_tracknumber.attrs['long_name'] = 'MCS tracknumbers that overlap with AR'
    dsout.mcs_tracknumber.attrs['comments'] = 'Subtract tracknumbers by 1 to match with MCS statistics file track indices'
    dsout.mcs_tracknumber.attrs['units'] = 'unitless'

    dsout.mcs_nhours.attrs['long_name'] = 'Number of hours each MCS overlap with AR'
    dsout.mcs_nhours.attrs['units'] = 'hours'

    # Write to netCDF file
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4_CLASSIC', unlimited_dims='tracks')
    print(f'Output saved: {outfile}')



if __name__ == '__main__':

    indates = sys.argv[1]
    config_file = sys.argv[2]
    
    inyear = indates[0:4]

    # get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    pixel_dir = config['pixelfile_dir']
    output_dir = config['output_dir']
    ar_dir = config['ar_dir']
    landmask_file = config['landmask_file']
    n_workers = config['n_workers']

    # AR tracking file
    arfile = f'{ar_dir}tag_MERRA2_hourlyIVT_{inyear}.nc'

    outfile = f'{output_dir}mcs_ar_tracknumbers_{indates}.nc'
    os.makedirs(output_dir, exist_ok=True)

    begin_time = datetime.datetime.now()

    # Find all pixel-level files
    pixelfiles = sorted(glob.glob(f'{pixel_dir}{indates}/mcstrack_{inyear}????_????.nc'))
    print(f'Number of pixel files: {len(pixelfiles)}')

    # Define regions to remove AR. 
    # The Europe region is the same with Rutz et al. (2019) JGR Fig. 4.
    # West cost of North and South America
    lon_boxes = [[-15,10], [-140,-115], [-76,-70]]
    lat_boxes = [[35,62], [35,62], [-62,-30]]

    # Get landmask
    dslm = xr.open_dataset(landmask_file)
    lonlm = dslm.lon
    latlm = dslm.lat
    lon2d, lat2d = np.meshgrid(lonlm, latlm)
    # Define coast mask (land fraction between 20%-90%)
    coastmask = (dslm.landseamask < 90) & (dslm.landseamask > 20)

    # Defines shape of growth. This grows one pixel as a cross
    dilationstructure = generate_binary_structure(2,1)
    # Dilate coast mask to include surrounding area of the coastline
    coastmaskdilate = binary_dilation(coastmask.data, structure=dilationstructure, iterations=5).astype(coastmask.dtype)

    # Make a copy of the dilated coast mask
    # coastmask_west = np.copy(coastmaskdilate)
    coastmask_keep = np.zeros(coastmaskdilate.shape, dtype=np.int8)

    # Loop over each boxes, mark coastal regions
    for ii in range(len(lon_boxes)):
        mask = (lon2d >= min(lon_boxes[ii])) & (lon2d <= max(lon_boxes[ii])) & \
               (lat2d >= min(lat_boxes[ii])) & (lat2d <= max(lat_boxes[ii])) & \
               (coastmaskdilate == 1)
        coastmask_keep[mask] = 1

    # Get the region lat/lon boundary
    ds = xr.open_dataset(pixelfiles[0])
    ny = ds.dims['lat']
    nx = ds.dims['lon']
    lon = ds.lon
    lat = ds.lat
    lonmin, lonmax = lon.min().data, lon.max().data
    latmin, latmax = lat.min().data, lat.max().data

    # Read AR data (1 file per year)
    dsar = xr.open_dataset(arfile)

    # For some reason, the first/last point from AR was not included even though lonmin/lonmax appears to be exactly matching the AR lat/lon
    # Add a very small buffer value so that the subset AR region is exactly the same as the MCS region
    dsar = dsar.sel(lon=slice(lonmin-0.01, lonmax+0.01), lat=slice(latmin-0.01, latmax+0.01))
    nyar = dsar.dims['lat']
    nxar = dsar.dims['lon']

    # Check to make sure the AR subset domain size is the same with MCS domain size
    if (nx != nxar) | (ny != nyar):
        print("Error: AR subset domain size does not match MCS domain size!")
        sys.exit()

    # Set up a dictionary that keyed by filename, valued by filetime.
    # This way the pixelfile name can be used to get a time string to locate the corresponding AR time
    time_to_pixfilename = {key:value for key, value in map(make_timestring, pixelfiles)}


    # Find boundary of the coastal region
    yid, xid = np.where(coastmask_keep == 1)
    min_x, max_x = np.min(xid), np.max(xid)
    min_y, max_y = np.min(yid), np.max(yid)
    bounds = {'min_y':min_y, 'max_y':max_y, 'min_x':min_x, 'max_x':max_x}
    # Subset the coastal region
    coastmask_keep = coastmask_keep[min_y:max_y, min_x:max_x]

    # import matplotlib.pyplot as plt
    # import pdb; pdb.set_trace()

    # Initialize dask
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
    client = Client(cluster)

    saved_tracknumber = []
    # for ifile in range(0, 10):
    for ifile in range(0, len(pixelfiles)):

        ifilename = pixelfiles[ifile]
        # Get the time string 
        itime = time_to_pixfilename.get(ifilename)
        print(itime)

        # Find the neartest AR data time, and subset the domain
        ar_mask = dsar.ar_binary_tag.sel(time=itime, method='nearest').data[min_y:max_y, min_x:max_x]

        # Call function to find MCS tracknumber within ARs
        # tracknumber = get_mcs_ar_tracknumber(ifilename, ar_mask, coastmask_keep, bounds)
        tracknumber = delayed(get_mcs_ar_tracknumber)(ifilename, ar_mask, coastmask_keep, bounds)
        
        # Extend the tracknumber list
        saved_tracknumber.append(tracknumber)

    # Collect results from Dask
    results = dask.compute(*saved_tracknumber)
    # print("Done processing {} files.".format(len(results))

    # Flatten the append list
    # tracknumber_list = [item for sublist in saved_tracknumber for item in sublist]
    tracknumber_list = [item for sublist in results for item in sublist]

    # Get unique track numbers, counts. Counts mean number of hours each MCS overlap with AR
    mcs_tracknumber, mcs_nhours = np.unique(tracknumber_list, return_counts=True)
    ntracks = len(mcs_tracknumber)

    # Write output to file
    write_netcdf(outfile, mcs_tracknumber, mcs_nhours)

    # Print processing time
    print(datetime.datetime.now() - begin_time)


