**Procedures to find MCSs in AR**


*  Find MCS tracks within AR over land/coast (uses Dask parallelization)

	For Europe, MCSs over coastal region in Western Europe (15W-10E, 35N-62N) within ARs are identified. The region is the same with that defined in Rutz et al. (2019) JGR (ARTMIP project).
	
	- *find_mcs_tracks_in_ar_coast.py*
	
    Shell script to run multiple years serially:
	- *loop_find_mcs_tracks_ar.sh*

**Procedures to find MCSs in TC:**

*  Find MCS tracks within TC (uses Dask parallelization) :
	
    - *find_mcs_tracks_in_tc.py*
	
    Shell script to run multiple years and regions (sometimes doesn't work when running multiple year/regions, Dask complains memory issue):
	*loop_find_mcs_tracks_tc.sh*

**Procedures to filter AR and TC**


*  Filter MCSs identified in AR and TC together (for regions that has AR: NAM, EUROPE):
	- *filter_mcs_tracks_ar_tc.py*
	
    Shell script to run multiple years serially
	- *loop_filter_mcs_tracks_ar_tc.sh*

*  Filter MCSs identified in TC (for regions with TC only):
	- *filter_mcs_tracks_tc.py*
	
    Shell script to run multiple years serially
	- *loop_filter_mcs_tracks_tc.sh*

*  Move original pixel-level MCS track files to mcstracking_orig directory to preserve them (important, otherwise they will be overwritten in the next step)

*  Rerun labeling pixel-level MCS tracks step using the filtered robust MCS statistics files (robust_mcs_tracks_extc_yyyymmdd_yyyymmdd.nc)
	
    IDL script to run multiple years for a region (should be run in interactive node): (FLEXTRKR directory)
	- *loop_run_gpm_irpf_mcs.pro*
	
    Or just submit a job to run the labeling pixel-level MCS step


**Procedures to calculate daily/monthly statistics on pixel grid**

*  Calculate daily mean MCS statistics:
	- *map_imerg_mcs_stats_byday.py*
	- *map_imerg_mcs_stats_byday_nospeed.py*

	Calculate daily mean MCS statistics from robust MCS statistics files and maps onto the native pixel grid.
	The "_nospeed" version excludes propagation speed variables, useful for data files before advection speed step is run.
	For large regions with lots of MCS (e.g. Asia, SPac) , directly calculating monthly mean statistics would take too long to run.
	This approach calculates daily mean statistics and saves to daily netCDF files, which are much faster to run with more parallization. The daily output netCDF files can be further averaged to monthly files.

*  Calculate monthly mean MCS statistics:
	- *map_imerg_mcs_stats_bymonth.py*
	- *map_imerg_mcs_stats_bymonth_nospeed.py*


	Calculate monthly mean MCS statistics from robust MCS statistics files and maps onto the native pixel grid.
	The "_nospeed" version excludes propagation speed variables, useful for data files before advection speed step is run.
	For smaller region with fewer MCS (e.g., NAM), monthly statistics can be run directly.


	- *calc_mcs_monthly_stats_from_daily.py*
	- *calc_mcs_monthly_stats_from_daily_nospeed.py*

	
	Calculates monthly mean MCS statistics from daily mean files.

	- *calc_imerg_mcs_monthly_precipmap_single.py*

	Calculates monthly mean total, MCS precipitation amount, MCS precipitation hours and saves to a netCDF file.

	- *calc_imerg_mcs_monthly_preciphov_single.py*
	
	Calculates monthly Hovmoller diagram of total, MCS precipitation and saves to a netCDF file.

*  Calculate seasonal mean MCS statistics:
	- *calc_mcs_seasonal_mean_rainmap.py*

	Calculates seasonal mean total, MCS precipitation amount and saves to a netCDF file.

	
	- *calc_mcs_seasonal_mean_statsmap.py*
	- *calc_mcs_seasonal_mean_statsmap_nospeed.py*

	Calculate seasonal mean MCS statistics (e.g., MCS number, lifetime, area, max rain rate, etc.) and saves to a netCDF file.

*  Job submission scripts using NERSC TaskFarmer for parallization are in the /scripts directory.