import numpy as np
import glob, os, sys
import xarray as xr
import time, datetime, calendar, pytz
import dask
from dask import delayed
from dask.distributed import Client, LocalCluster

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
#     bt = calendar.timegm(datetime.datetime(year, month, day, hour, minute, second, tzinfo=utc).timetuple())
    return filename, filetime

def get_mcs_ar_tracknumber(mcs_filename, ar_mask, landmask):
    # Read MCS pixel file
    dspx = xr.open_dataset(mcs_filename)
    # cloudtracknumber = dspx.cloudtracknumber.astype(np.int32).squeeze()
    cloudtracknumber = dspx.cloudtracknumber.squeeze().data
    
    # Find MCS mask pixels that overlaps with AR mask
    # overlap = cloudtracknumber[np.where(ar_mask == 1)]
    overlap = cloudtracknumber[np.where((ar_mask == 1) & (landmask == 1))]
    # Get the unique MCS tracknumbers
    mcs_tracknumber_ar = np.unique(overlap[~np.isnan(overlap)]).astype(np.int32)
    
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

    region = sys.argv[1]
    indates = sys.argv[2]
    # region = 'npac'
    # indates = '20150101_20151231'
    inyear = indates[0:4]

    # Number of workers for Dask
    n_workers = 20

    statsdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_pt1/'
    pixeldir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/mcstracking_ccs4_pt1/{indates}/'
 #   statsdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/stats_ccs4_4h/'
 #   pixeldir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/{region}/mcstracking_ccs4_4h/{indates}/'
    #ardir = '/global/cscratch1/sd/feng045/waccem/AR_Tempest/'
    ardir = '/global/cscratch1/sd/feng045/waccem/AR_Tempest_hourly/'

    pixelfiles = sorted(glob.glob(f'{pixeldir}mcstrack_{inyear}????_????.nc'))
    #arfile = f'{ardir}MERRA2.ar_tag.Tempest_v1.3hourly.{inyear}.nc'
    arfile = f'{ardir}tag_MERRA2_hourlyIVT_{inyear}.nc'
    print(f'Number of pixel files: {len(pixelfiles)}')

    landmaskfile = f'/global/project/projectdirs/m1867/zfeng/gpm/map_data/IMERG_landmask_{region}.nc'

    outdir = statsdir
    # outfile = f'{outdir}mcs_ar_tracknumbers_{indates}.nc'
    outfile = f'{outdir}mcs_ar_tracknumbers_{indates}.nc'

    # Get landmask
    dslm = xr.open_dataset(landmaskfile)
    landmask = (dslm.landseamask < 90).data

    # Get the region lat/lon boundary
    ds = xr.open_dataset(pixelfiles[0])
    ny = ds.dims['lat']
    nx = ds.dims['lon']
    lon = ds.lon
    lat = ds.lat
    lonmin, lonmax = lon.min().data, lon.max().data
    latmin, latmax = lat.min().data, lat.max().data
    # print(lonmin, lonmax, latmin, latmax)
    # print(ds.dims)

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

    # Initialize dask
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
    client = Client(cluster)

    saved_tracknumber = []
    # for ifile in range(0, 60):
    for ifile in range(0, len(pixelfiles)):

        ifilename = pixelfiles[ifile]
        # Get the time string 
        itime = time_to_pixfilename.get(ifilename)
        print(itime)

        # Find the neartest AR data time
        ar_mask = dsar.ar_binary_tag.sel(time=itime, method='nearest').data

        # Call function to find MCS tracknumber within ARs
        # tracknumber = get_mcs_ar_tracknumber(ifilename, ar_mask, landmask)
        tracknumber = delayed(get_mcs_ar_tracknumber)(ifilename, ar_mask, landmask)
        
        # Extend the tracknumber list
        # saved_tracknumber.extend(tracknumber)
        saved_tracknumber.append(tracknumber)

    # Collect results from Dask
    results = dask.compute(*saved_tracknumber)
    # print("Done processing {} files.".format(len(results))

    # Flatten the append list
    tracknumber_list = [item for sublist in results for item in sublist]

    # Get unique track numbers, counts. Counts mean number of hours each MCS overlap with AR
    mcs_tracknumber, mcs_nhours = np.unique(tracknumber_list, return_counts=True)
    ntracks = len(mcs_tracknumber)

    # import pdb; pdb.set_trace()

    # Get unique track numbers, counts. Counts mean number of hours each MCS overlap with AR
    # mcs_tracknumber, mcs_nhours = np.unique(saved_tracknumber, return_counts=True)

    # Write output to file
    write_netcdf(outfile, mcs_tracknumber, mcs_nhours)



