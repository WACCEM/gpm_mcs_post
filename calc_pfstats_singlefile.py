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


def sort_renumber(labelcell_number2d, min_cellpix):
    """
    Sorts 2D labeled cells by size, and removes cells smaller than min_cellpix.
    
    Args:
        labelcell_number2d: np.ndarray()
            Labeled cell number array in 2D.
        min_cellpix: float
            Minimum number of pixel to count as a cell.

    Returns:
        sortedlabelcell_number2d: np.ndarray(int)
            Sorted labeled cell number array in 2D.
        sortedcell_npix: np.ndarray(int)
            Number of pixels for each labeled cell in 1D.
    """

    # Create output arrays
    sortedlabelcell_number2d = np.zeros(np.shape(labelcell_number2d), dtype=int) 

    # Get number of labeled cells
    nlabelcells = np.nanmax(labelcell_number2d)

    # Check if there is any cells identified
    if (nlabelcells > 0):

        labelcell_npix = np.full(nlabelcells, -999, dtype=int)

        # Loop over each labeled cell
        for ilabelcell in range(1, nlabelcells+1):
            # Count number of pixels for the cell
            ilabelcell_npix = np.count_nonzero(labelcell_number2d == ilabelcell)
            # Check if cell satisfies size threshold
            if (ilabelcell_npix > min_cellpix):
                labelcell_npix[ilabelcell-1] = ilabelcell_npix

        # Check if any of the cells passes the size threshold test
        ivalidcells = np.array(np.where(labelcell_npix > 0))[0, :]
        ncells = len(ivalidcells)

        if (ncells > 0):
            # Isolate cells that satisfy size threshold
            # Add one since label numbers start at 1 and indices, which validcells reports starts at 0
            labelcell_number1d = np.copy(ivalidcells) + 1 
            labelcell_npix = labelcell_npix[ivalidcells]

            # Sort cells from largest to smallest and get the sorted index
            order = np.argsort(labelcell_npix)
            order = order[::-1]  # Reverses the order

            # Sort the cells by size
            sortedcell_npix = np.copy(labelcell_npix[order])
            sortedcell_number1d = np.copy(labelcell_number1d[order])

            # Loop over the 2D cell number to re-number them by size
            cellstep = 0
            for icell in range(0, ncells):
                # Find 2D indices that match the cell number
                sortedcell_indices = np.where(labelcell_number2d == sortedcell_number1d[icell])
                # Get one of the dimensions from the 2D indices to count the size
    #             nsortedcellindices = np.shape(sortedcell_indices)[1]
                nsortedcellindices = len(sortedcell_indices[1])
                # Check if the size matches the sorted cell size
                if (nsortedcellindices == sortedcell_npix[icell]):
                    # Renumber the cell in 2D
                    cellstep += 1
                    sortedlabelcell_number2d[sortedcell_indices] = np.copy(cellstep)

        else:
            # Return an empty array
            sortedcell_npix = np.zeros(0)
    else:
        # Return an empty array
        sortedcell_npix = np.zeros(0)
    
    return sortedlabelcell_number2d, sortedcell_npix


def calc_pfstats_singlefile(pixel_filename, pixel_filebase, rmcs_basetime, rr_min, pf_min_npix, pixel_radius, nmaxpf):

    pixelfile_path = os.path.dirname(pixel_filename)
    prelength = len(pixelfile_path)+1 + len(pixel_filebase)
    ittdatetimestring = pixel_filename[prelength:(prelength+13)]

    

    # Read pixel data file
    if os.path.isfile(pixel_filename):
        print(pixel_filename)
        # cloudiddata = Dataset(pixel_filename, 'r')
        # cloudid_basetime = cloudiddata['base_time'][:]
        # cloudnumbermap = cloudiddata['cloudtracknumber'][:]
        # rawrainratemap = cloudiddata['precipitation'][:]
        # cloudiddata.close()

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

        # Find all matching time indices from robust MCS stats file to the pixel file
        matchindices = np.array(np.where(np.abs(rmcs_basetime - cloudid_basetime) < 1))
        # The returned match indices are for [tracks, times] dimensions respectively
        idx_track = matchindices[0]
        idx_time = matchindices[1]

        nmatchcloud = len(idx_track)

        if (nmatchcloud > 0):
            # Create arrays for PF statistics
            pf_npf = np.full(nmatchcloud, -999, dtype=int)
            pf_pfnpix = np.full((nmatchcloud, nmaxpf), -999, dtype=int)
            pf_pflon = np.full((nmatchcloud, nmaxpf), np.nan, dtype=float)
            pf_pflat = np.full((nmatchcloud, nmaxpf), np.nan, dtype=float)
            pf_maxrate = np.full((nmatchcloud, nmaxpf), np.nan, dtype=float)

            # Generate a dilation structure. This grows one pixel as a cross
            dilationstructure = generate_binary_structure(2,1)  

            for imatchcloud in range(nmatchcloud):
                # Intialize array for only MCS data
                filteredrainratemap = np.full((ydim, xdim), np.nan, dtype=float)

                # Track number needs to add 1
                itracknum = idx_track[imatchcloud] + 1
                # itime = idx_time[imatchcloud]

                # Get a mask for the current cloud track
                # icloudmask = cloudnumbermap == itracknum
                # Get location indices of the cloud
                icloudlocationy, icloudlocationx = np.where(cloudnumbermap == itracknum)
                inpix_cloud = len(icloudlocationy)

                # Fill array with MCS data
                filteredrainratemap[icloudlocationy, icloudlocationx] = np.copy(rawrainratemap[icloudlocationy,icloudlocationx])

                ## Set edges of boundary
                miny = np.nanmin(icloudlocationy)
                if miny <= 10:
                    miny = 0
                else:
                    miny = miny - 10

                maxy = np.nanmax(icloudlocationy)
                if maxy >= ydim - 10:
                    maxy = ydim
                else:
                    maxy = maxy + 11

                minx = np.nanmin(icloudlocationx)
                if minx <= 10:
                    minx = 0
                else:
                    minx = minx - 10

                maxx = np.nanmax(icloudlocationx)
                if maxx >= xdim - 10:
                    maxx = xdim
                else:
                    maxx = maxx + 11

                ## Isolate smaller region around cloud shield
                subrainratemap = np.copy(filteredrainratemap[miny:maxy, minx:maxx])

                # Get dimensions of subsetted region
                subdimy, subdimx = np.shape(subrainratemap)

                # Derive precipitation feature statistics
                # print('Calculating precipitation statistics')
                ipfy, ipfx = np.where(subrainratemap > rr_min)
                nrainpix = len(ipfy)

                if nrainpix > 0:
                    # Create binary map
                    binarypfmap = subrainratemap > rr_min

                    # Dilate (aka smooth)
                    # dilatedbinarypfmap = binary_dilation(binarypfmap, structure=dilationstructure, iterations=1).astype(filteredrainratemap.dtype)

                    # Label precipitation features
                    # pfnumberlabelmap, numpf = label(dilatedbinarypfmap)
                    pfnumberlabelmap, numpf = label(binarypfmap)

                    # Sort and renumber PFs, and remove small PFs
                    pf_number, pf_npix = sort_renumber(pfnumberlabelmap, pf_min_npix)
                    # Update number of PFs after sorting and renumbering
                    npf_new = np.nanmax(pf_number)
                    numpf = npf_new
                    pfnumberlabelmap = pf_number
                    del pf_number, npf_new

                    if numpf > 0:
                        npf_save = np.nanmin([nmaxpf, numpf])
                        # print('PFs present, initializing matrices')

                        # Initialize matrices
                        pfnpix = np.zeros(npf_save, dtype=float)
                        pflon = np.full(npf_save, np.nan, dtype=float)
                        pflat = np.full(npf_save, np.nan, dtype=float)
                        pfmaxrate = np.full(npf_save, np.nan, dtype=float)

                        # print('Looping through each feature to calculate statistics')
                        print(('Number of PFs ' + str(numpf)))
                        ###############################################
                        # Loop through each feature
                        for ipf in range(1, npf_save+1):
                            # Find associated indices
                            iipfy, iipfx = np.where(pfnumberlabelmap == ipf)
                            niipfpix = len(iipfy)

                            # Compare the size of the PF to make sure it matches with the sorted PF
                            if (niipfpix != pf_npix[ipf-1]):
                                sys.exit("pixel count not match")
                            else:
                                # Compute statistics

                                # Basic statistics
                                pfnpix[ipf-1] = np.copy(niipfpix)
                                pflon[ipf-1] = np.nanmean(lat[iipfy[:]+miny, iipfx[:]+minx])
                                pflat[ipf-1] = np.nanmean(lon[iipfy[:]+miny, iipfx[:]+minx])

                                pfmaxrate[ipf-1] = np.nanmax(subrainratemap[iipfy[:], iipfx[:]])
                        
                        # Save precipitation feature statisitcs
                        pf_npf[imatchcloud] = np.copy(numpf)
                        pf_pflon[imatchcloud, 0:npf_save]= pflon[0:npf_save]
                        pf_pflat[imatchcloud, 0:npf_save] = pflat[0:npf_save]
                        pf_pfnpix[imatchcloud, 0:npf_save] = pfnpix[0:npf_save]
                        pf_maxrate[imatchcloud, 0:npf_save] = pfmaxrate[0:npf_save]

                        # print(imatchcloud, inpix_cloud)

            # if (nmatchcloud > 1) & (inpix_cloud > 2000):
            import pdb; pdb.set_trace()

        # import pdb; pdb.set_trace()

    return nmatchcloud, matchindices, pf_npf, pf_pflon, pf_pflat, pf_pfnpix, 