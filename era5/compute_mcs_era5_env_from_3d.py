"""
Calculates ERA5 2D environmental variables from 3D profiles and saves to a netCDF file.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "14-June-2022"

import os, sys
import time
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import diag_functions as afwa
from wrf import interplevel
from itertools import repeat
from multiprocessing import Pool

#-----------------------------------------------------------------------------------------
def calc_shear(level_shear, u, v, z_agl, u10, v10):
    """
    Calculates vertical wind shear and direction at specified AGL levels

    Args:
        level_shear: np.array(float)
            Vertical height level AGL (km).
        u: xr.DataArray
            U components of wind.
        v: xr.DataArray
            V components of wind.
        z_agl: xr.DataArray
            Height above ground level.
        u10: xr.DataArray
            10m U components of wind.
        v10: xr.DataArray
            10m V components of wind.    

    Returns:
        shear_mag: xr.DataArray
            Bulk wind shear magnitude.
        shear_dir: xr.DataArray
            Bulk wind shear direction.
    """
    # Interpolate U, V to specific AGL levels
    u_agl = interplevel(u, z_agl, level_shear)
    v_agl = interplevel(v, z_agl, level_shear)
    # Calculate wind speed
    wspd_agl = np.sqrt(u_agl**2 + v_agl**2)
    wspd10 = np.sqrt(u10**2 + v10**2)

    # Wind shear magnitude
    shear_mag = wspd_agl - wspd10

    # Wind shear direction
    ushear = u_agl - u10
    vshear = v_agl - v10
    # .metpy.dequantify() converts data back to a unit-naive array
    # https://unidata.github.io/MetPy/latest/tutorials/xarray_tutorial.html?highlight=dequantify#
    shear_dir = mpcalc.wind_direction(ushear * units('m/s'), vshear * units('m/s')).metpy.dequantify()

    # Assign attributes
    shear_mag_attrs = {'long_name':'Bulk wind shear magnitude', 'units':'m/s'}
    shear_dir_attrs = {'long_name':'Bulk wind shear direction', 'units':'degree'}
    shear_mag = shear_mag.assign_attrs(shear_mag_attrs)
    shear_dir = shear_dir.assign_attrs(shear_dir_attrs)

    return shear_mag, shear_dir

#-----------------------------------------------------------------------------------------
def calc_shear_2level(u, v, z, level_lower, level_upper):
    """
    Calculates vertical wind shear and direction between two levels

    Args:
        u: xr.DataArray
            U components of wind.
        v: xr.DataArray
            V components of wind.
        z: xr.DataArray
            Height above ground level.
        level_lower: float
            Lower height level.
        level_upper: float
            Upper height level.    

    Returns:
        shear_mag: xr.DataArray
            Bulk wind shear magnitude.
        shear_dir: xr.DataArray
            Bulk wind shear direction.
    """
    # Interpolate U, V at lower & upper levels (on HAMSL coordinate)
    u_lower = interplevel(u, z, level_lower)
    v_lower = interplevel(v, z, level_lower)
    u_upper = interplevel(u, z, level_upper)
    v_upper = interplevel(v, z, level_upper)
    # Calculate wind speed
    wspd_lower = np.sqrt(u_lower**2 + v_lower**2)
    wspd_upper = np.sqrt(u_upper**2 + v_upper**2)
    # Wind shear magnitude
    shear_mag = wspd_upper - wspd_lower
    # Wind shear direction
    ushear = u_upper - u_lower
    vshear = v_upper - v_lower
    # .metpy.dequantify() converts data back to a unit-naive array
    # https://unidata.github.io/MetPy/latest/tutorials/xarray_tutorial.html?highlight=dequantify#
    shear_dir = mpcalc.wind_direction(ushear * units('m/s'), vshear * units('m/s')).metpy.dequantify()

    # Assign attributes
    shear_mag_attrs = {'long_name':f'Bulk wind shear magnitude between {level_lower}m and {level_upper}m', 'units':'m/s'}
    shear_dir_attrs = {'long_name':f'Bulk wind shear direction between {level_lower}m and {level_upper}m', 'units':'degree'}
    shear_mag = shear_mag.assign_attrs(shear_mag_attrs)
    shear_dir = shear_dir.assign_attrs(shear_dir_attrs)

    return shear_mag, shear_dir

#-----------------------------------------------------------------------------------------
def calc_envs_track(file_names, tracknumber):
    """
    Calculate environments for a given tracknumber.

    Args:
        file_names: dictionary
            Dictionary containing input file names.
        tracknumber: int
            Track number to process.

    Returns:
        var_dict: dictionary
            Dictionary containing environmental variables.
        var_attrs: dictionary
            Dictionary containing the attributes of environmental variables.
    """
    print(f'track number: {tracknumber}')

    # Get file names from dictionary
    file_T = file_names['file_T'] 
    file_Q = file_names['file_Q']
    file_R = file_names['file_R']
    file_Z = file_names['file_Z']
    file_U = file_names['file_U']
    file_V = file_names['file_V']
    file_W = file_names['file_W']
    file_10u = file_names['file_10u']
    file_10v = file_names['file_10v']
    file_2t = file_names['file_2t']
    file_2d = file_names['file_2d']
    file_sp = file_names['file_sp']
    file_z_sfc = file_names['file_z_sfc']

    # Define levels (HAGL) to calculate surface wind shear
    level_shear = [2000, 4000, 6000, 8000, 10000]   # [m AGL]
    level_theta_e = [3000, 5000, 7000]   # [m AGL]

    # Read 3D variables
    # Temperature
    ds_t = xr.open_dataset(file_T).sel(tracks=tracknumber).load()
    # tracks = ds_t['tracks']
    # ntracks = ds_t.dims['tracks']
    ntimes = ds_t.dims['rel_times']
    nz = ds_t.dims['level']
    ny = ds_t.dims['y']
    nx = ds_t.dims['x']
    # Reverse vertical coordinate
    ds_t = ds_t.reindex(level=list(reversed(ds_t.level))).load()
    tk = ds_t['T']

    # Water vapor
    ds_q = xr.open_dataset(file_Q).sel(tracks=tracknumber).load()
    ds_q = ds_q.reindex(level=list(reversed(ds_q.level)))
    # Convert unit to g/kg
    q = ds_q['Q'] * 1000
    q_attrs = ds_q['Q'].attrs
    q_attrs['units'] = 'g kg**-1'

    # Relative humidity
    ds_rh = xr.open_dataset(file_R).sel(tracks=tracknumber).load()
    ds_rh = ds_rh.reindex(level=list(reversed(ds_rh.level)))
    rh = ds_rh['R']
    rh_attrs = rh.attrs

    # Geopotential height
    ds_z = xr.open_dataset(file_Z).sel(tracks=tracknumber).load()
    ds_z = ds_z.reindex(level=list(reversed(ds_z.level)))
    # Convert unit to m
    z = ds_z['Z'] / 9.80665

    # U wind
    ds_u = xr.open_dataset(file_U).sel(tracks=tracknumber).load()
    ds_u = ds_u.reindex(level=list(reversed(ds_u.level)))
    u = ds_u['U']

    # V wind
    ds_v = xr.open_dataset(file_V).sel(tracks=tracknumber).load()
    ds_v = ds_v.reindex(level=list(reversed(ds_v.level)))
    v = ds_v['V']

    # W wind
    ds_w = xr.open_dataset(file_W).sel(tracks=tracknumber).load()
    ds_w = ds_w.reindex(level=list(reversed(ds_w.level)))
    w = ds_w['W']
    w_attrs = ds_w['W'].attrs

    # Read 2D variables, subset the tracks
    ds_10u = xr.open_dataset(file_10u).sel(tracks=tracknumber).load()
    u10 = ds_10u['VAR_10U']
    ds_10v = xr.open_dataset(file_10v).sel(tracks=tracknumber).load()
    v10 = ds_10v['VAR_10V']
    ds_2t = xr.open_dataset(file_2t).sel(tracks=tracknumber).load()
    t2 = ds_2t['VAR_2T']
    ds_2d = xr.open_dataset(file_2d).sel(tracks=tracknumber).load()
    td2 = ds_2d['VAR_2D']
    ds_sp = xr.open_dataset(file_sp).sel(tracks=tracknumber).load()
    sp = ds_sp['SP']
    ds_z_sfc = xr.open_dataset(file_z_sfc).sel(tracks=tracknumber).load()
    z_sfc = ds_z_sfc['Z'] / 9.80665

    # Get pressure level, convert unit to Pa
    p = ds_t.level.data * 100
    # Make a 3D array of pressure by repeating the profile ny, nx times
    p3d = np.tile(p, (ny * nx)).reshape(ny, nx, nz)
    # Swap dimension from 2 to 0 such that: [z, y, x]
    p3d = np.swapaxes(p3d, 2, 0)

    # Calculate dew point temperature
    td = mpcalc.dewpoint_from_relative_humidity(tk * units('K'), rh * units('%'))
    # Calculate equivalent potential temperature
    theta_e = mpcalc.equivalent_potential_temperature(p3d * units('Pa'), tk * units('K'), td).metpy.dequantify()
    # theta = mpcalc.potential_temperature(p3d * units('Pa'), tk * units('K')).metpy.dequantify()
    # Calculate equivalent potential temperature at 2m
    theta_e2m = mpcalc.equivalent_potential_temperature(sp * units('Pa'), t2 * units('K'), td2 * units('K')).metpy.dequantify()
    # Filter values below surface
    theta_e_f = theta_e.where(z > z_sfc, other=np.NaN)
    # theta_f = theta.where(z > z_sfc, other=np.NaN)

    # Convert HAMSL to HAGL by subtracting surface height
    z_agl = z - z_sfc
    z_agl = z_agl.where(z_agl >= 0)

    # Filter 3D variables below surface (replaced with -999)
    # They will not be included in searching for max ThetaE below 500 hPa for the most unstable level
    # Vertical coordinates (Z, P) are not filtered
    fillval = -999.
    tk_f = tk.where(z > z_sfc, other=fillval)
    rh_f = rh.where(z > z_sfc, other=fillval)
    q_f = q.where(z > z_sfc, other=fillval)
    u_f = u.where(z > z_sfc, other=np.NaN)
    v_f = v.where(z > z_sfc, other=np.NaN)
    w_f = w.where(z > z_sfc, other=np.NaN)

    # Create arrays to store outputs
    var2d_dims = (ntimes, ny, nx)
    mucape = np.full(var2d_dims, np.NaN, dtype=np.float32)
    mucin = np.full(var2d_dims, np.NaN, dtype=np.float32)
    lcl = np.full(var2d_dims, np.NaN, dtype=np.float32)
    lfc = np.full(var2d_dims, np.NaN, dtype=np.float32)
    el = np.full(var2d_dims, np.NaN, dtype=np.float32)
    lpl = np.full(var2d_dims, np.NaN, dtype=np.float32)

    theta_e_diff_0to3km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    theta_e_diff_0to5km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    theta_e_diff_0to7km = np.full(var2d_dims, np.NaN, dtype=np.float32)

    shear_mag_0to2km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_0to4km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_0to6km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_0to10km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_4to8km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_6to10km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_mag_lcl_el = np.full(var2d_dims, np.NaN, dtype=np.float32)

    shear_dir_0to2km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_0to4km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_0to6km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_0to10km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_4to8km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_6to10km = np.full(var2d_dims, np.NaN, dtype=np.float32)
    shear_dir_lcl_el = np.full(var2d_dims, np.NaN, dtype=np.float32)

    # Loop over times
    for itime in range(0, ntimes):
        # Proceed if this time has valid data (MCS exists)
        # if np.nanmax(z_sfc[itime, :, :]) > 0:
        if np.count_nonzero(~np.isnan(u10[itime, :, :])) > 0:
            _tk = tk_f[itime, :, :, :]
            _q = q_f[itime, :, :, :]
            _rh = rh_f[itime, :, :, :]
            _z = z[itime, :, :, :]
            _z_agl = z_agl[itime, :, :, :]
            _u = u_f[itime, :, :, :]
            _v = v_f[itime, :, :, :]
            _theta_e = theta_e_f[itime, :, :, :]
            # _z_sfc = z_sfc[itime, :, :]
            _u10 = u10[itime, :, :]
            _v10 = v10[itime, :, :]
            _theta_e2m = theta_e2m[itime, :, :]

            # Call AFWA diagnostics on data filtered below surface
            ostat, _mucape, _mucin, _lcl, _lfc, _el, _lpl = afwa.diag_functions.diag_map(_tk, _rh, p3d, _z, 1, 1)
            if ostat == 1:
                # Calculate shear between LCL and EL
                shear_mag_ll, shear_dir_ll = calc_shear_2level(_u, _v, _z, _lcl, _el)
                shear_mag_ll = shear_mag_ll.data
                shear_dir_ll = shear_dir_ll.data
                # Replace values where LCL or EL is not defined (-9.999e30)
                mask_ll = np.logical_or(_lcl < 0, _el < 0)
                shear_mag_ll[mask_ll] = np.NaN
                shear_dir_ll[mask_ll] = np.NaN
                shear_mag_lcl_el[itime, :, :] = shear_mag_ll
                shear_dir_lcl_el[itime, :, :] = shear_dir_ll

                # Replace undefined values with NaN
                _mucape[_mucape < 0] = np.NaN
                _mucin[_mucin < -999] = np.NaN
                _lcl[_lcl < 0] = np.NaN
                _lfc[_lfc < 0] = np.NaN
                _el[_el < 0] = np.NaN
                _lpl[_lpl < 0] = np.NaN
                # Save to output arrays
                mucape[itime, :, :] = _mucape
                mucin[itime, :, :] = _mucin
                lcl[itime, :, :] = _lcl
                lfc[itime, :, :] = _lfc
                el[itime, :, :] = _el
                lpl[itime, :, :] = _lpl

            # Calculate ThetaE difference
            # Interpolate ThetaE to specific AGL levels
            theta_e_agl = interplevel(_theta_e, _z_agl, level_theta_e)
            # Calculate ThetaE difference
            theta_e_diff = _theta_e2m - theta_e_agl
            theta_e_diff_0to3km[itime, :, :] = theta_e_diff.sel(level=3000)
            theta_e_diff_0to5km[itime, :, :] = theta_e_diff.sel(level=5000)
            theta_e_diff_0to7km[itime, :, :] = theta_e_diff.sel(level=7000)

            # Calculate shear between two levels
            shear_mag_4_8, shear_dir_4_8 = calc_shear_2level(_u, _v, _z_agl, 4000, 8000)
            shear_mag_4to8km[itime, :, :] = shear_mag_4_8
            shear_dir_4to8km[itime, :, :] = shear_dir_4_8
            shear_mag_6_10, shear_dir_6_10 = calc_shear_2level(_u, _v, _z_agl, 6000, 10000)
            shear_mag_6to10km[itime, :, :] = shear_mag_6_10
            shear_dir_6to10km[itime, :, :] = shear_dir_6_10

            # Calculate shear between surface and upper levels
            shear_mag, shear_dir = calc_shear(level_shear, _u, _v, _z_agl, _u10, _v10)
            shear_mag_0to2km[itime, :, :] = shear_mag.sel(level=2000)
            shear_mag_0to4km[itime, :, :] = shear_mag.sel(level=4000)
            shear_mag_0to6km[itime, :, :] = shear_mag.sel(level=6000)
            shear_mag_0to10km[itime, :, :] = shear_mag.sel(level=10000)

            shear_dir_0to2km[itime, :, :] = shear_dir.sel(level=2000)
            shear_dir_0to4km[itime, :, :] = shear_dir.sel(level=4000)
            shear_dir_0to6km[itime, :, :] = shear_dir.sel(level=6000)
            shear_dir_0to10km[itime, :, :] = shear_dir.sel(level=10000)

    # Get select level RH, Q, W
    rh_925mb = rh_f.sel(level=925).data
    rh_850mb = rh_f.sel(level=850).data
    rh_700mb = rh_f.sel(level=700).data
    rh_600mb = rh_f.sel(level=600).data
    rh_500mb = rh_f.sel(level=500).data
    q_925mb = q_f.sel(level=925).data
    q_850mb = q_f.sel(level=850).data
    q_700mb = q_f.sel(level=700).data
    q_600mb = q_f.sel(level=600).data
    q_500mb = q_f.sel(level=500).data
    w_850mb = w_f.sel(level=850).data
    w_700mb = w_f.sel(level=700).data
    w_500mb = w_f.sel(level=500).data

    # Group outputs in dictionaries
    var_dict = {
        'tracknumber': tracknumber,
        'MUCAPE': mucape,
        'MUCIN': mucin,
        'LCL': lcl,
        'LFC': lfc,
        'EL': el,
        'LPL': lpl,
        'rh_925mb': rh_925mb,
        'rh_850mb': rh_850mb,
        'rh_700mb': rh_700mb,
        'rh_600mb': rh_600mb,
        'rh_500mb': rh_500mb,
        'q_925mb': q_925mb,
        'q_850mb': q_850mb,
        'q_700mb': q_700mb,
        'q_600mb': q_600mb,
        'q_500mb': q_500mb,
        'w_850mb': w_850mb,
        'w_700mb': w_700mb,
        'w_500mb': w_500mb,
        'theta_e_diff_0to3km': theta_e_diff_0to3km,
        'theta_e_diff_0to5km': theta_e_diff_0to5km,
        'theta_e_diff_0to7km': theta_e_diff_0to7km,
        'shear_mag_0to2km': shear_mag_0to2km,
        'shear_mag_0to4km': shear_mag_0to4km,
        'shear_mag_0to6km': shear_mag_0to6km,
        'shear_mag_0to10km': shear_mag_0to10km,
        'shear_mag_4to8km': shear_mag_4to8km,
        'shear_mag_6to10km': shear_mag_6to10km,
        'shear_mag_lcl_el': shear_mag_lcl_el,
        'shear_dir_0to2km': shear_dir_0to2km,
        'shear_dir_0to4km': shear_dir_0to4km,
        'shear_dir_0to6km': shear_dir_0to6km,
        'shear_dir_0to10km': shear_dir_0to10km,
        'shear_dir_4to8km': shear_dir_4to8km,
        'shear_dir_6to10km': shear_dir_6to10km,
        'shear_dir_lcl_el': shear_dir_lcl_el,
    }
    shear_mag_attrs = {'long_name':'Bulk wind shear magnitude', 'units':'m/s'}
    shear_dir_attrs = {'long_name':'Bulk wind shear direction', 'units':'degree'}
    theta_e_diff_attrs = {'long_name':'Equivalent potential temperature difference between layers', 'units':'K'}
    var_attrs = {
        'tracknumber': {
            'long_name': 'Track number for each MCS'
        },
        'MUCAPE': {
            'long_name': 'Most unstable convective available potential energy',
            'units': 'J/kg',
            # '_FillValue': fillval,
        },
        'MUCIN': {
            'long_name': 'Most unstable convective inhibition',
            'units': 'J/kg',
            # '_FillValue': fillval,
        },
        'LCL': {
            'long_name': 'Lifted condensation level',
            'units': 'm',
            # '_FillValue': fillval,
        },
        'LFC': {
            'long_name': 'Level of free convection',
            'units': 'm',
            # '_FillValue': fillval,
        },
        'EL': {
            'long_name': 'Equilibrium level',
            'units': 'm',
            # '_FillValue': fillval,
        },
        'LPL': {
            'long_name': 'Most unstable lifted parcel level',
            'units': 'm',
            # '_FillValue': fillval,
        },
        'rh_925mb': rh_attrs,
        'rh_850mb': rh_attrs,
        'rh_700mb': rh_attrs,
        'rh_600mb': rh_attrs,
        'rh_500mb': rh_attrs,
        'q_925mb': q_attrs,
        'q_850mb': q_attrs,
        'q_700mb': q_attrs,
        'q_600mb': q_attrs,
        'q_500mb': q_attrs,
        'w_850mb': w_attrs,
        'w_700mb': w_attrs,
        'w_500mb': w_attrs,
        'theta_e_diff_0to3km': theta_e_diff_attrs,
        'theta_e_diff_0to5km': theta_e_diff_attrs,
        'theta_e_diff_0to7km': theta_e_diff_attrs,
        'shear_mag_0to2km': shear_mag_attrs,
        'shear_mag_0to4km': shear_mag_attrs,
        'shear_mag_0to6km': shear_mag_attrs,
        'shear_mag_0to10km': shear_mag_attrs,
        'shear_mag_4to8km': shear_mag_attrs,
        'shear_mag_6to10km': shear_mag_attrs,
        'shear_mag_lcl_el': {
            'long_name': 'Bulk wind shear magnitude between LCL and EL',
            'units': 'm/s',
        },
        'shear_dir_0to2km': shear_dir_attrs,
        'shear_dir_0to4km': shear_dir_attrs,
        'shear_dir_0to6km': shear_dir_attrs,
        'shear_dir_0to10km': shear_dir_attrs,
        'shear_dir_4to8km': shear_dir_attrs,
        'shear_dir_6to10km': shear_dir_attrs,
        'shear_dir_lcl_el': {
            'long_name': 'Bulk wind shear direction between LCL and EL',
            'units': 'degree',
        },
    }
    # import pdb; pdb.set_trace()
    return var_dict, var_attrs

#-----------------------------------------------------------------------------------------
def work_for_tracks(file_names, out_filename):
    """
    Calculate environments for all tracks and saves output to a netCDF file.

    Args:
        file_names: dictionary
            Dictionary containing input file names.
        out_filename: string
            Output file name.

    Returns:
        out_filename: string
            Output file name.
    """

    file_T = file_names['file_T']

    # Read 3D variables
    ds_t = xr.open_dataset(file_T)
    tracks = ds_t['tracks']
    ntracks = ds_t.dims['tracks']
    ntimes = ds_t.dims['rel_times']
    ny = ds_t.dims['y']
    nx = ds_t.dims['x']

    final_result = []
    # Serial
    if run_parallel == 0:
        for itrack in range(0, ntracks):
        # for itrack in range(91, 92):
            tracknumber = tracks.data[itrack]
            result = calc_envs_track(file_names, tracknumber)
            final_result.append(result)
    # Parallel
    elif run_parallel >= 1:
        pool = Pool(n_workers)
        final_result = pool.starmap(calc_envs_track, zip(repeat(file_names), tracks.data))
        pool.close()
    else:
        sys.exit('Valid parallelization flag not set.')

    # Make a variable list from one of the returned dictionaries
    var_names = list(final_result[0][0].keys())
    # Get variable attributes from one of the returned dictionaries
    var_attrs = final_result[0][1]

    # Remove tracknumbers from the list
    var_names.remove('tracknumber')
    var_attrs.pop('tracknumber', None)

    # Loop over variable list to create the dictionary entry
    out_dict = {}
    out_dict_attrs = {}
    var2d_dims = (ntracks, ntimes, ny, nx)
    for ivar in var_names:
        out_dict[ivar] = np.full(var2d_dims, np.nan, dtype=np.float32)
        out_dict_attrs[ivar] = var_attrs[ivar]

    # Collect results
    for itrack in range(0, ntracks):
    # for itrack in range(791, 792):
        if final_result[itrack] is not None:
            # Get the return results for this track
            # The result is a tuple: (out_dict, out_dict_attrs)
            # The first entry is the dictionary containing the variables
            iResult = final_result[itrack][0]
            tracknumber = iResult['tracknumber']

            # Double check tracknumber from return to make sure it matches the track
            if tracks.data[itrack] == tracknumber:
                # Loop over each variable and assign values to output dictionary
                for ivar in var_names:
                    if iResult[ivar].ndim == 3:
                        out_dict[ivar][itrack, :, :, :] = iResult[ivar]
                        # import pdb; pdb.set_trace()
            else:
                print(f'ERROR: tracknumber does not match: {tracknumber}!')
                sys.exit('Double check results!')

    # Define a dataset containing all PF variables
    var_dict = {}
    # Define output variable dictionary
    for key, value in out_dict.items():
        if value.ndim == 4:
            var_dict[key] = (['tracks', 'rel_times', 'y', 'x'], value, out_dict_attrs[key])
    coord_dict = {
        'tracks': (['tracks'], tracks.data, tracks.attrs),
        'rel_times': (['rel_times'], ds_t['rel_times'].data, ds_t['rel_times'].attrs),
        'y': (['y'], ds_t['y'].data, ds_t['y'].attrs),
        'x': (['x'], ds_t['x'].data, ds_t['x'].attrs),
    }
    gattr_dict = {
        'Title': 'ERA5 computed 2D environment data for MCS',
        'Institution': 'Pacific Northwest National Laboratory',
        'Contact': 'zhe.feng@pnnl.gov',
        'Created_on': time.ctime(time.time()),
    }
    # Define xarray dataset
    dsout = xr.Dataset(var_dict, coords=coord_dict, attrs=gattr_dict)
    # Set encoding/compression for all variables
    comp = dict(zlib=True, dtype='float32')
    encoding = {var: comp for var in dsout.data_vars}
    # Write to netcdf file
    dsout.to_netcdf(path=out_filename, mode='w', format='NETCDF4', 
                    unlimited_dims='tracks', encoding=encoding)
    print(f'Output: {out_filename}')
    # import pdb; pdb.set_trace()
    return out_filename


if __name__ == "__main__":

    # Get track string from input
    track_string = sys.argv[1]

    # Parallel setup
    run_parallel = 1
    n_workers = 20

    # Get sub-strings
    # track_string format: 20190101.0000_20200101.0000_t10100
    track_period = track_string[:-7]
    track_year = track_string[0:4]

    # ERA5 input data directory
    in_dir3d = '/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_3d/'
    in_dir2d = '/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_2d/'
    # Output file directory
    out_dir = f'/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_2d_derived/{track_year}/'
    out_basename = 'mcs_era5_2D_ENVS_'
    out_filename = f'{out_dir}{out_basename}{track_string}.nc'
    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # 2D variable files
    file_10u = f'{in_dir2d}mcs_era5_VAR_10U_{track_period}.nc'
    file_10v = f'{in_dir2d}mcs_era5_VAR_10V_{track_period}.nc'
    file_2t = f'{in_dir2d}mcs_era5_VAR_2T_{track_period}.nc'
    file_2d = f'{in_dir2d}mcs_era5_VAR_2D_{track_period}.nc'
    file_sp = f'{in_dir2d}mcs_era5_SP_{track_period}.nc'
    file_z_sfc = f'{in_dir2d}mcs_era5_Z_{track_period}.nc'

    # 3D variable files
    file_T = f'{in_dir3d}{track_year}/mcs_era5_T_{track_string}.nc'
    file_Q = f'{in_dir3d}{track_year}/mcs_era5_Q_{track_string}.nc'
    file_R = f'{in_dir3d}{track_year}/mcs_era5_R_{track_string}.nc'
    file_Z = f'{in_dir3d}{track_year}/mcs_era5_Z_{track_string}.nc'
    file_U = f'{in_dir3d}{track_year}/mcs_era5_U_{track_string}.nc'
    file_V = f'{in_dir3d}{track_year}/mcs_era5_V_{track_string}.nc'
    file_W = f'{in_dir3d}{track_year}/mcs_era5_W_{track_string}.nc'

    # Check input files existance
    fc_10u = os.path.isfile(file_10u)
    fc_10v = os.path.isfile(file_10v)
    fc_2t = os.path.isfile(file_2t)
    fc_2d = os.path.isfile(file_2d)
    fc_sp = os.path.isfile(file_sp)
    fc_z_sfc = os.path.isfile(file_z_sfc)
    fc_T = os.path.isfile(file_T)
    fc_Q = os.path.isfile(file_Q)
    fc_R = os.path.isfile(file_R)
    fc_Z = os.path.isfile(file_Z)
    fc_U = os.path.isfile(file_U)
    fc_V = os.path.isfile(file_V)
    fc_W = os.path.isfile(file_W)
    # Print missing file name
    if fc_10u == False:
        print(f'Missing: {file_10u}')
    if fc_10v == False:
        print(f'Missing: {file_10v}')
    if fc_2t == False:
        print(f'Missing: {file_2t}')
    if fc_2d == False:
        print(f'Missing: {file_2d}')
    if fc_sp == False:
        print(f'Missing: {file_sp}')
    if fc_z_sfc == False:
        print(f'Missing: {file_z_sfc}')
    if fc_T == False:
        print(f'Missing: {file_T}')
    if fc_Q == False:
        print(f'Missing: {file_Q}')
    if fc_R == False:
        print(f'Missing: {file_R}')
    if fc_Z == False:
        print(f'Missing: {file_Z}')
    if fc_U == False:
        print(f'Missing: {file_U}')
    if fc_V == False:
        print(f'Missing: {file_V}')
    if fc_W == False:
        print(f'Missing: {file_W}')

    # Check if all input files exist
    if fc_10u & fc_10v & fc_2t & fc_2d & fc_sp & fc_z_sfc & \
        fc_T & fc_Q & fc_R & fc_Z & fc_U & fc_V & fc_W:

        # Group input file names in a dictionary
        file_names = {
            'file_T': file_T, 
            'file_Q': file_Q, 
            'file_R': file_R, 
            'file_Z': file_Z, 
            'file_U': file_U, 
            'file_V': file_V, 
            'file_W': file_W, 
            'file_10u': file_10u, 
            'file_10v': file_10v, 
            'file_2t': file_2t,
            'file_2d': file_2d,
            'file_sp': file_sp,
            'file_z_sfc': file_z_sfc,
        }
        # Call function to calculate
        result = work_for_tracks(file_names, out_filename)
    else:
        print('ERROR: missing input file(s), code exists.')
        sys.exit()
    # import pdb; pdb.set_trace()