**Procedures to find MCSs in AR**


*  Find MCS tracks within AR over land/coast (uses Dask parallelization)
	○ For Europe, MCSs over coastal region in Western Europe (15W-10E, 35N-62N) within ARs are identified. The region is the same with that defined in Rutz et al. (2019) JGR (ARTMIP project).
	
	*find_mcs_tracks_in_ar_coast.py*
	
    Shell script to run multiple years serially:
	*loop_find_mcs_tracks_ar.sh*

**Procedures to find MCSs in TC:**

*  Find MCS tracks within TC (uses Dask parallelization) :
	
    Python code: *find_mcs_tracks_in_tc.py*
	
    Shell script to run multiple years and regions (sometimes doesn't work when running multiple year/regions, Dask complains memory issue):
	*loop_find_mcs_tracks_tc.sh*

**Procedures to filter AR and TC**


*  Filter MCSs identified in AR and TC together (for regions that has AR: NAM, EUROPE):
	*filter_mcs_tracks_ar_tc.py*
	
    Shell script to run multiple years serially
	*loop_filter_mcs_tracks_ar_tc.sh*

*  Filter MCSs identified in TC (for regions with TC only):
	*filter_mcs_tracks_tc.py*
	
    Shell script to run multiple years serially
	*loop_filter_mcs_tracks_tc.sh*

*  Move original pixel-level MCS track files to mcstracking_orig directory to preserve them (important, otherwise they will be overwritten in the next step)

*  Rerun labeling pixel-level MCS tracks step using the filtered robust MCS statistics files (robust_mcs_tracks_extc_yyyymmdd_yyyymmdd.nc)
	
    IDL script to run multiple years for a region (should be run in interactive node): (FLEXTRKR directory)
	*loop_run_gpm_irpf_mcs.pro*
	
    Or just submit a job to run the labeling pixel-level MCS step
