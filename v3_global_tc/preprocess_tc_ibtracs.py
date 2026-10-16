"""
Extract global Tropical Cyclone ROCI from IBTrACS v04r01.

Input:
    ./data/ibtracs.since1980.list.v04r01.csv

Output:
    ./data/ibtracs_global_1998_2025.nc

Notes:
    - ROCI-related variables in IBTrACS v04r01 are in nautical miles (nmile).
    - This script converts nmile to km using 1 nmile = 1.852 km.
    - base_time is written as Unix epoch seconds.
    - The script includes sanity checks to prevent accidentally writing kiloseconds,
      milliseconds, microseconds, or nanoseconds.
"""

import os
import time
import numpy as np
import pandas as pd
import xarray as xr


# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
filename = "./data/ibtracs.since1980.list.v04r01.csv"

start_year = 1998
end_year = 2025

outfilename = f"./data/ibtracs_global_{start_year}_{end_year}.nc"

# IBTrACS v04r01 radius units are nautical miles.
NMILE_TO_KM = 1.852

# Maximum physically plausible ROCI threshold after unit conversion.
max_roci_km = 2000.0


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
def to_numeric(series):
    """Convert a pandas Series to numeric values, coercing invalid strings to NaN."""
    return pd.to_numeric(series, errors="coerce")


def unix_seconds_from_datetime64(iso_time):
    """
    Convert a pandas datetime Series to Unix epoch seconds robustly.

    Important:
    Do NOT assume that ``iso_time.astype("int64")`` is always in nanoseconds.
    Depending on pandas / numpy versions, timezone-aware datetimes can be
    backed by microsecond or other resolutions. Using ``Timestamp.timestamp()``
    avoids the resolution-dependent 1e3 error that produces kiloseconds
    instead of seconds.

    Returns
    -------
    base_time : pandas Series, float64
        Unix epoch seconds. NaN where ISO_TIME is invalid.
    """
    base_time = pd.Series(np.nan, index=iso_time.index, dtype="float64")
    valid = iso_time.notna()

    # Robust conversion to Unix seconds.
    # For ~1e5 IBTrACS rows this is fast enough and avoids unit ambiguity.
    base_time.loc[valid] = (
        iso_time.loc[valid]
        .apply(lambda t: t.timestamp())
        .astype("float64")
    )

    return base_time

def check_base_time(base_time, year, start_year, end_year):
    """
    Sanity check base_time.

    For 1998-2025, Unix seconds should be O(1e9), not O(1e6).
    O(1e6) indicates kiloseconds and will break TC-MCS time matching.
    """
    bt = np.asarray(base_time, dtype="float64")
    yy = np.asarray(year, dtype="float64")

    ok = np.isfinite(bt) & np.isfinite(yy) & (yy >= start_year) & (yy <= end_year)
    if not np.any(ok):
        raise RuntimeError("No finite base_time values available for sanity check.")

    bt_min = np.nanmin(bt[ok])
    bt_max = np.nanmax(bt[ok])

    dt_min = pd.to_datetime(bt_min, unit="s", utc=True)
    dt_max = pd.to_datetime(bt_max, unit="s", utc=True)

    print("------------------------------------------------------------")
    print("base_time sanity check")
    print("------------------------------------------------------------")
    print(f"base_time range: {bt_min:.3f} to {bt_max:.3f}")
    print(f"as UTC datetime: {dt_min} to {dt_max}")

    # Expected Unix seconds range for requested period, with a 2-year loose margin.
    expected_min = pd.Timestamp(f"{start_year - 2}-01-01", tz="UTC").timestamp()
    expected_max = pd.Timestamp(f"{end_year + 2}-12-31", tz="UTC").timestamp()

    if not (expected_min <= bt_min <= expected_max and expected_min <= bt_max <= expected_max):
        raise RuntimeError(
            "base_time sanity check failed. Values are not plausible Unix seconds.\n"
            f"Expected roughly between {expected_min:.0f} and {expected_max:.0f}; "
            f"got {bt_min:.3f} to {bt_max:.3f}.\n"
            "If values are around 1e6 for 1998-2025, they are kiloseconds. "
            "Check that the conversion uses / 1.0e9, not / 1.0e12 or // 10**12."
        )


def fill_roci_from_column(data, roci_nmile, roci_source, colname):
    """Fill missing ROCI values from one source column."""
    vals = to_numeric(data[colname])
    use = roci_nmile.isna() & vals.notna()
    roci_nmile.loc[use] = vals.loc[use]
    roci_source.loc[use] = colname
    return roci_nmile, roci_source


if __name__ == "__main__":

    # ------------------------------------------------------------
    # Columns to read
    # ------------------------------------------------------------
    usecols = [
        "SEASON",
        "SID",
        "NAME",
        "BASIN",
        "TRACK_TYPE",
        "LON",
        "LAT",
        "ISO_TIME",
        "IFLAG",
        "WMO_AGENCY",
        "USA_ROCI",
        "TOKYO_R30_LONG",
        "REUNION_R34_NE",
        "REUNION_R34_SE",
        "REUNION_R34_SW",
        "REUNION_R34_NW",
        "BOM_ROCI",
    ]

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Input file not found: {filename}")

    # IBTrACS CSV has a second header-like row containing units.
    # skiprows=[1] skips that row.
    data = pd.read_csv(
        filename,
        usecols=usecols,
        skiprows=[1],
        dtype=str,
        low_memory=False,
    )

    # ------------------------------------------------------------
    # Parse time and check available years
    # ------------------------------------------------------------
    iso_time = pd.to_datetime(data["ISO_TIME"], errors="coerce", utc=True)
    iso_year = iso_time.dt.year
    season = to_numeric(data["SEASON"])

    print("------------------------------------------------------------")
    print("IBTrACS file year check")
    print("------------------------------------------------------------")
    print(f"Input file: {filename}")
    print(f"Requested period: {start_year}-{end_year}")

    if season.notna().any():
        print(f"Available SEASON range: {int(season.min())}-{int(season.max())}")
        print(f"Rows with SEASON == {end_year}: {int((season == end_year).sum())}")
    else:
        print("WARNING: SEASON could not be parsed.")

    if iso_year.notna().any():
        print(f"Available ISO_TIME year range: {int(iso_year.min())}-{int(iso_year.max())}")
        print(f"Rows with ISO_TIME year == {end_year}: {int((iso_year == end_year).sum())}")
    else:
        raise RuntimeError("ISO_TIME could not be parsed. Cannot continue.")

    print("\nLast 10 ISO_TIME year counts:")
    print(iso_year.value_counts().sort_index().tail(10))

    if not ((season == end_year).any() or (iso_year == end_year).any()):
        print(
            f"\nWARNING: No {end_year} records found in this local file. "
            f"The script will still process available records within {start_year}-{end_year}."
        )

    # ------------------------------------------------------------
    # Convert basic variables
    # ------------------------------------------------------------
    lon = to_numeric(data["LON"])
    lat = to_numeric(data["LAT"])

    year = iso_time.dt.year.astype("float64")
    month = iso_time.dt.month.astype("float64")
    day = iso_time.dt.day.astype("float64")
    hour = iso_time.dt.hour.astype("float64")
    minute = iso_time.dt.minute.astype("float64")
    second = iso_time.dt.second.astype("float64")

    # Correct Unix epoch time in seconds.
    # This is the key line that prevents the previous /1000 error.
    base_time = unix_seconds_from_datetime64(iso_time)

    # ------------------------------------------------------------
    # Extract ROCI using priority logic
    # ------------------------------------------------------------
    roci_nmile = pd.Series(np.nan, index=data.index, dtype="float64")
    roci_source = pd.Series("", index=data.index, dtype="object")

    roci_nmile, roci_source = fill_roci_from_column(data, roci_nmile, roci_source, "USA_ROCI")
    roci_nmile, roci_source = fill_roci_from_column(data, roci_nmile, roci_source, "TOKYO_R30_LONG")
    roci_nmile, roci_source = fill_roci_from_column(data, roci_nmile, roci_source, "BOM_ROCI")

    reunion_cols = [
        "REUNION_R34_NE",
        "REUNION_R34_SE",
        "REUNION_R34_SW",
        "REUNION_R34_NW",
    ]
    reunion_r34 = data[reunion_cols].apply(pd.to_numeric, errors="coerce")
    reunion_r34_max = reunion_r34.max(axis=1, skipna=True)

    use_reunion = roci_nmile.isna() & reunion_r34_max.notna()
    roci_nmile.loc[use_reunion] = reunion_r34_max.loc[use_reunion]
    roci_source.loc[use_reunion] = "MAX_REUNION_R34"

    # Convert nautical miles to km.
    roci_km = roci_nmile * NMILE_TO_KM

    # ------------------------------------------------------------
    # Subset target period and quality control
    # ------------------------------------------------------------
    mask = (
        (year >= start_year)
        & (year <= end_year)
        & lon.notna()
        & lat.notna()
        & base_time.notna()
        & roci_km.notna()
        & (roci_km > 0)
        & (roci_km < max_roci_km)
    )

    idx = np.where(mask.to_numpy())[0]
    nrecords = len(idx)

    print("------------------------------------------------------------")
    print("ROCI extraction summary")
    print("------------------------------------------------------------")
    print(f"Valid records after filtering: {nrecords}")

    if nrecords == 0:
        raise RuntimeError("No valid ROCI records found after filtering.")

    print("\nROCI source counts:")
    print(roci_source.loc[mask].value_counts())

    imax = roci_km.loc[mask].idxmax()
    print("\nMaximum extracted ROCI record:")
    print(f"SID        : {data.loc[imax, 'SID']}")
    print(f"NAME       : {data.loc[imax, 'NAME']}")
    print(f"BASIN      : {data.loc[imax, 'BASIN']}")
    print(f"ISO_TIME   : {data.loc[imax, 'ISO_TIME']}")
    print(f"LAT/LON    : {lat.loc[imax]}, {lon.loc[imax]}")
    print(f"ROCI       : {roci_km.loc[imax]:.2f} km")
    print(f"Source     : {roci_source.loc[imax]}")

    # Sanity check before writing.
    check_base_time(base_time.loc[mask].to_numpy(), year.loc[mask].to_numpy(), start_year, end_year)

    # ------------------------------------------------------------
    # Build xarray Dataset
    # ------------------------------------------------------------
    storms = np.arange(nrecords, dtype=np.int32)

    dsout = xr.Dataset(
        data_vars={
            "sid": (["storms"], data["SID"].fillna("").to_numpy(dtype=str)[idx]),
            "name": (["storms"], data["NAME"].fillna("").to_numpy(dtype=str)[idx]),
            "basin": (["storms"], data["BASIN"].fillna("").to_numpy(dtype=str)[idx]),
            "track_type": (["storms"], data["TRACK_TYPE"].fillna("").to_numpy(dtype=str)[idx]),
            "lon": (["storms"], lon.to_numpy(dtype=np.float32)[idx]),
            "lat": (["storms"], lat.to_numpy(dtype=np.float32)[idx]),
            "roci": (["storms"], roci_km.to_numpy(dtype=np.float32)[idx]),
            "roci_source": (["storms"], roci_source.to_numpy(dtype=str)[idx]),
            "base_time": (["storms"], base_time.to_numpy(dtype=np.float64)[idx]),
            "year": (["storms"], year.to_numpy(dtype=np.int16)[idx]),
            "month": (["storms"], month.to_numpy(dtype=np.int16)[idx]),
            "day": (["storms"], day.to_numpy(dtype=np.int16)[idx]),
            "hour": (["storms"], hour.to_numpy(dtype=np.int16)[idx]),
            "minute": (["storms"], minute.to_numpy(dtype=np.int16)[idx]),
            "second": (["storms"], second.to_numpy(dtype=np.int16)[idx]),
        },
        coords={
            "storms": (["storms"], storms),
        },
        attrs={
            "title": "IBTrACS v04r01 global post-processed ROCI",
            "source_file": filename,
            "period_requested": f"{start_year}-{end_year}",
            "time_filter": "ISO_TIME calendar year",
            "roci_priority": "USA_ROCI > TOKYO_R30_LONG > BOM_ROCI > max(REUNION_R34 quadrants)",
            "roci_original_units": "nautical miles",
            "roci_output_units": "km",
            "roci_conversion": "1 nautical mile = 1.852 km",
            "max_roci_km_threshold": max_roci_km,
            "created_on": time.ctime(time.time()),
        },
    )

    # Variable attributes
    dsout["sid"].attrs = {
        "long_name": "IBTrACS storm identifier",
    }
    dsout["name"].attrs = {
        "long_name": "storm name",
    }
    dsout["basin"].attrs = {
        "long_name": "IBTrACS basin",
    }
    dsout["track_type"].attrs = {
        "long_name": "IBTrACS track type",
    }
    dsout["lon"].attrs = {
        "long_name": "storm center longitude",
        "units": "degrees_east",
    }
    dsout["lat"].attrs = {
        "long_name": "storm center latitude",
        "units": "degrees_north",
    }
    dsout["roci"].attrs = {
        "long_name": "radius of outermost closed isobar or proxy radius",
        "units": "km",
        "note": "Extracted using priority logic; TOKYO_R30_LONG and REUNION_R34 are wind-radius proxies.",
    }
    dsout["roci_source"].attrs = {
        "long_name": "source variable used to derive roci",
    }
    dsout["base_time"].attrs = {
        "long_name": "Unix epoch time",
        "standard_name": "time",
        "units": "seconds since 1970-01-01 00:00:00 UTC",
        "calendar": "standard",
        "note": "Stored as Unix seconds, not kiloseconds.",
    }
    dsout["year"].attrs = {"long_name": "UTC year from ISO_TIME"}
    dsout["month"].attrs = {"long_name": "UTC month from ISO_TIME"}
    dsout["day"].attrs = {"long_name": "UTC day from ISO_TIME"}
    dsout["hour"].attrs = {"long_name": "UTC hour from ISO_TIME"}
    dsout["minute"].attrs = {"long_name": "UTC minute from ISO_TIME"}
    dsout["second"].attrs = {"long_name": "UTC second from ISO_TIME"}

    # Compression.
    # Avoid applying zlib encoding to string variables.
    encoding = {}
    for var in dsout.data_vars:
        if np.issubdtype(dsout[var].dtype, np.number):
            encoding[var] = {"zlib": True, "complevel": 4}

    dsout.to_netcdf(
        path=outfilename,
        mode="w",
        format="NETCDF4",
        unlimited_dims="storms",
        encoding=encoding,
    )

    print("------------------------------------------------------------")
    print(f"Output saved: {outfilename}")
    print("------------------------------------------------------------")

