import numpy as np
import os, glob
import time, datetime, calendar
from netCDF4 import Dataset, num2date, chartostring
from scipy.ndimage import label, binary_dilation, generate_binary_structure
from skimage.measure import regionprops
from math import pi
from scipy.stats import skew
import xarray as xr
import pandas as pd

######################################################
def location_to_idx(lat, lon, center_lat, center_lon):
    """ 
    Convert a latitude and longitude into an index
    
    Parameters
    ----------
    lat: latitude array 
        latitude values of grid
    lon: longitude array 
        latitude values of grid
    center_lat: float
        latitude location to find
    center_lon: float
        longitude location to find
    """
    # This is for 2D lat/lon
    diff = abs(lat - center_lat) + abs(lon - center_lon)
    lat_idx, lon_idx = np.unravel_index(diff.argmin(), diff.shape)
    return lat_idx, lon_idx

# def pad_array(s_array, ny, nx):
#     """ Pad s_array out to span dimensions. Centering Array, offset distributed right. """
#     s_array_dims = s_array.shape
#     d1 = int(np.ceil(2*ny - s_array_dims[0] ))
#     d2 = int(np.ceil(2*nx - s_array_dims[1] ))

#     d_array = np.pad(s_array, ((d1,d1), (d2,d2)), 'constant', constant_values=np.nan)
#     return d_array[0:ny*2, 0:nx*2]

######################################################
def composite_mcs_center(
    pixel_filename,
    idx_track,
    pf_lon,
    pf_lat,
    nx,
    ny,
    ):
    """ 
    Composite MCS PF-center for each MCS track from a pixel-level file
    
    Parameters
    ----------
    pixel_filename: string
        File name of the pixel-level data
    idx_track: int array 
        Track indices for each MCS in the pixel-level file
    pf_lon: float array
        Longitudes of MCS PF center
    pf_lat: float array
        Latitudes of MCS PF center
    nx: int
        Number of grids in x-direction for the MCS window (x size: nx * 2 + 1)
    ny: int
        Number of grids in y-direction for the MCS window (y size: ny * 2 + 1)
    """

    # Create arrays for composite
    nmatchcloud = len(idx_track)
    # npf_out = 0
    rain_out = np.full((nmatchcloud, ny*2+1, nx*2+1), 0, dtype=np.float)

    # Read pixel data file
    if os.path.isfile(pixel_filename):
        print(pixel_filename)

        # Read pixel-level data
        ds = xr.open_dataset(pixel_filename, decode_times=False)
        lon = ds.longitude.values
        lat = ds.latitude.values
        cloudid_basetime = ds.base_time.values
        cloudnumbermap = ds.cloudtracknumber.squeeze().values
        rawrainratemap = ds.precipitation.squeeze().values
        xdim = ds.dims['lon']
        ydim = ds.dims['lat']
        ds.close()

        if (nmatchcloud > 0):

            for imatchcloud in range(nmatchcloud):

                # Track number needs to add 1
                itracknum = idx_track[imatchcloud] + 1
                # print(itracknum)

                # Find closest pixel location to the MCS PF center
                liy, lix = location_to_idx(lat, lon, pf_lat[imatchcloud], pf_lon[imatchcloud])
                # Take a window center at the MCS PF
                # rain_cut = pad_array(rawrainratemap[liy-ny:liy+ny, lix-nx:lix+nx], ny, nx)
                # cloudnumber_cut = pad_array(cloudnumbermap[liy-ny:liy+ny, lix-nx:lix+nx], ny, nx)
                rain_cut = rawrainratemap[liy-ny:liy+ny+1, lix-nx:lix+nx+1]
                cloudnumber_cut = cloudnumbermap[liy-ny:liy+ny+1, lix-nx:lix+nx+1]

                # Mask out area not matching the MCS track number
                imask = cloudnumber_cut == itracknum
                rain_cut[imask == False] = 0

                # Check number of pixels matching the MCS
                inpix_cloud = np.count_nonzero(imask)
                if inpix_cloud > 0:
                    # Save the MCS rain map
                    rain_out[imatchcloud, :, :] = rain_cut
                    # npf_out += 1
                else:
                    print(f'Warning: no cloud matching track # {itracknum}')
                    # import pdb; pdb.set_trace()

    # Sum over space for each MCS
    rain_sum = np.sum(rain_out, axis=(1,2))
    # Find the ones with non-zero rainfall
    idx_valid = np.where(rain_sum > 0)[0]
    # This is the number of valid MCS PF
    npf_out = len(idx_valid)
    rain_out = rain_out[idx_valid, :, :]

    # If no valid MCS, replace output base time as NaN
    # This frame will be excluded before final output in the main function
    if npf_out == 0:
        cloudid_basetime = np.array(np.NaN)

    return (rain_out, npf_out, cloudid_basetime)