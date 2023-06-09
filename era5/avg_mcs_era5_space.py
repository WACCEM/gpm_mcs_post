import os, sys
import yaml
import xarray as xr
import numpy as np

if __name__ == '__main__':

    infile = sys.argv[1]
    # config_file = sys.argv[2]

    # # Read configuration from yaml file
    # stream = open(config_file, "r")
    # config = yaml.full_load(stream)

    # in_dir3d = config["in_dir3d"]
    # in_dir2d = config["in_dir2d"]
    # in_dir2d_derived = config["in_dir2d_derived"]
    # output_dir = config["output_dir"]

    # Define average spatial box radius in [degrees]
    # The total average width is 2 * box_radius (+/- box_radius)
    box_radius = 1.5

    # Output directory
    indir = os.path.dirname(infile)
    outdir = f'{indir}/avg_space/'
    os.makedirs(outdir, exist_ok=True)

    # Output filename
    infn = os.path.basename(infile)
    outfile = f'{outdir}{infn}'

    # Read input data
    ds = xr.open_dataset(infile)
    # Get coordinates in degree
    x = ds['x'] * 0.25
    y = ds['y'] * 0.25
    # Find indices for x, y dimensions
    xid = np.where((x >= -1*box_radius) & (x <= box_radius))[0]
    yid = np.where((y >= -1*box_radius) & (y <= box_radius))[0]

    # Subset spatial area, then average
    dsout = ds.isel(x=slice(min(xid), max(xid)), y=slice(min(yid), max(yid))).mean(dim=('y', 'x'), keep_attrs=True)
    # Add global attributes
    dsout.attrs['spatial_avg_width_degree'] = box_radius * 2

    # Set encoding/compression for all variables
    comp = dict(zlib=True)
    encoding = {var: comp for var in dsout.data_vars}
    dsout.to_netcdf(path=outfile, mode='w', format='NETCDF4', unlimited_dims='tracks', encoding=encoding)
    print('Output saved as: ', outfile)