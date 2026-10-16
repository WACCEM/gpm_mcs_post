"""
Calculate the maximum Tropical Cyclone ROCI (radius of the last closed isobar) from the IBTrACS dataset.
"""
import numpy as np
import pandas as pd
import xarray as xr
import calendar, datetime, pytz, time

if __name__ == '__main__':
        
    # Input IBTrACS CSV file (Global)
    # https://doi.org/10.25921/82ty-9e16
    # filename = f'/global/cfs/cdirs/m1867/zfeng/gpm/IBTrACS/ibtracs.since1980.list.v04r00.csv'
    # IBTrACS CSV file (North Atlantic)
    filename = f'/global/cfs/cdirs/m1867/zfeng/gpm/IBTrACS/ibtracs.NA.list.v04r00.csv'
    
    start_year = 2004
    end_year = 2022

    # Output file name
    # outfilename = f'/global/cfs/cdirs/m1867/zfeng/gpm/IBTrACS/ibtracs_{start_year}_{end_year}.nc'
    outfilename = f'/global/cfs/cdirs/m1867/zfeng/gpm/IBTrACS/ibtracs_NA_{start_year}_{end_year}.nc'

    # Only read several columns
    usecols=['SEASON','SID','LON','LAT','ISO_TIME','IFLAG','WMO_AGENCY','USA_ROCI','TOKYO_R30_LONG',\
            'REUNION_R34_NE','REUNION_R34_SE','REUNION_R34_NW','REUNION_R34_SW','BOM_ROCI']
    data = pd.read_csv(filename, usecols=usecols, skiprows=[1])

    n_data = np.size(data['SEASON'])
    year = np.zeros(n_data, dtype=np.int16)
    month = np.zeros(n_data, dtype=np.int16)
    day = np.zeros(n_data, dtype=np.int16)
    hr = np.zeros(n_data, dtype=np.int16)
    min = np.zeros(n_data, dtype=np.int16)
    sec = np.zeros(n_data, dtype=np.int16)
    roci = np.zeros(n_data, dtype=np.float32)
    base_time = np.full(n_data, np.nan, dtype=float)

    # Save date and ROCI
    for i in np.arange(n_data):
        # if i%5000 == 0:
        #     print(i)
        tmp = np.asarray(data['ISO_TIME'][i:i+1])
        year[i] = int((tmp[0].split()[0]).split('-')[0])
        month[i] = int((tmp[0].split()[0]).split('-')[1])
        day[i] = int((tmp[0].split()[0]).split('-')[2])
        hour = ((tmp[0].split()[1]).split())[0].split(':')
        # hr[i] = float(float(hour[0]) + float(hour[1])/60 + float(hour[2])/3600)
        hr[i] = int(hour[0])
        min[i] = int(hour[1])
        sec[i] = int(hour[2])

        if (year[i] >= start_year) & (year[i] <= end_year):
            print(tmp)
            # Comput Epoch time
            base_time[i] = calendar.timegm(
                datetime.datetime(year[i], month[i], day[i], hr[i], min[i], sec[i], tzinfo=pytz.UTC).timetuple()
            )

            if len(np.asarray(data['USA_ROCI'][i:i+1])[0].strip()) !=0:
                roci[i] = np.float(np.asarray(data['USA_ROCI'][i:i+1]))
            elif len(np.asarray(data['TOKYO_R30_LONG'][i:i+1])[0].strip()) !=0:
                roci[i] = np.float(np.asarray(data['TOKYO_R30_LONG'][i:i+1]))
            elif len(np.asarray(data['BOM_ROCI'][i:i+1])[0].strip()) !=0:
                roci[i] = np.float(np.asarray(data['BOM_ROCI'][i:i+1]))
            elif len(np.asarray(data['REUNION_R34_NE'][i:i+1])[0].strip()) !=0 or \
                len(np.asarray(data['REUNION_R34_SE'][i:i+1])[0].strip()) !=0 or \
                len(np.asarray(data['REUNION_R34_NW'][i:i+1])[0].strip()) !=0 or \
                len(np.asarray(data['REUNION_R34_SW'][i:i+1])[0].strip()) !=0:
                a = [
                    np.float(data['REUNION_R34_NE'][i:i+1]), 
                    np.float(data['REUNION_R34_SE'][i:i+1]), 
                    np.float(data['REUNION_R34_NW'][i:i+1]),
                    np.float(data['REUNION_R34_SW'][i:i+1]),
                ]
                roci[i] = np.float(np.max(a))
            else:
                roci[i] = 0.0


    # Convert DOI unit from [mile] to [km]
    roci = roci * 1.60934

    # Subset time period and ROCI > 0 & ROCI < 2000 km
    # There are some data with ROCI up to 4000-8000 km, clearly wrong
    id = np.where(
        (year >= start_year) & (year <= end_year) & 
        (roci > 0) & (roci < 2000)
    )[0]    
    lon = np.asarray(data['LON'])
    lat = np.asarray(data['LAT'])

    # Number of storm records
    nstorms = len(id)

    # Define output dictionaries
    var_dict = {
        'lon': (['storms'], lon[id]),
        'lat': (['storms'], lat[id]),
        'roci': (['storms'], roci[id]),
        'base_time': (['storms'], base_time[id]),
        'year': (['storms'], year[id]),
        'month': (['storms'], month[id]),
        'day': (['storms'], day[id]),
        'hour': (['storms'], hr[id]),
        'minute': (['storms'], min[id]),
        'second': (['storms'], sec[id]),
    }
    coord_dict = {
        'storms': (['storms'], np.arange(0, nstorms, 1))
    }
    gattr_dict = {
        'title': 'IBTrACS post-processed ROCI', \
        'contact': 'Zhe Feng, zhe.feng@pnnl.gov', \
        'created_on': time.ctime(time.time())
    }
    # Define Xarray Dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)

    # Set encoding/compression for all variables
    comp = dict(zlib=True)
    encoding = {var: comp for var in dsout.data_vars}

    # Write to netCDF file
    dsout.to_netcdf(path=outfilename, mode='w', format='NETCDF4', unlimited_dims='storms', encoding=encoding)
    print(f'Output saved: {outfilename}')

    import pdb; pdb.set_trace()


