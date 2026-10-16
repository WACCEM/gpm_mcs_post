# TC exclusion for the global MCS dataset V3 (IMERG V07B, 1998–2025)

Scripts used to remove MCS tracks that overlap tropical cyclones (TCs) from the V3 global MCS tracking dataset (the `extc` files). They were written and modified by Sihan Zhang for the V3 dataset.

They supersede the scripts of the same names in the top-level folder of this repository, which are not used for V3. The top-level versions use only the single TC record nearest in time at each hour and convert ROCI with a statute-mile factor; they are kept for the earlier regional datasets (the top-level `preprocess_tc_ibtracs.py` currently reads the North Atlantic IBTrACS file).

## Workflow

1. `preprocess_tc_ibtracs.py`
   Reads `./data/ibtracs.since1980.list.v04r01.csv` (IBTrACS v04r01, https://doi.org/10.25921/82ty-9e16) and writes `./data/ibtracs_global_1998_2025.nc`, one record per best-track point (`lon`, `lat`, `roci` in km, `base_time` in Unix seconds). The radius is the first available of `USA_ROCI`, `TOKYO_R30_LONG`, `BOM_ROCI`, or the largest `REUNION_R34_*` quadrant, converted from nautical miles with 1 nautical mile = 1.852 km; only records with 0 < ROCI < 2000 km are kept.

2. `find_mcs_tracks_in_tc.py <yyyymmdd.hhmm_yyyymmdd.hhmm> <config.yml>`
   For each hourly pixel file of the tracking period, finds the MCS track numbers (`cloudtracknumber`) that overlap a circle of radius ROCI around **any** TC record within ±`tc_time_window_hours` (default 3) of the MCS time, using great-circle distance and handling the date line. Writes `mcs_tc_tracknumbers_<dates>.nc` to `output_dir`.
   Config keys: `pixelfile_dir`, `output_dir`, `run_parallel`, `n_workers`, and optionally `tc_time_window_hours`.
   **Note:** the path of the TC file is hard-coded in the script (variable `tc_file`) and has to be edited.

3. `filter_mcs_tracks_ar_tc.py` and `mcsmip_filter_mcs_tracks_ar_tc.py` `<dates> <config.yml>`
   Remove the flagged tracks from `mcs_tracks_final_<dates>.nc` (MCSMIP: `mcsmip_mcs_tracks_final_<dates>.nc`) and write `..._final_extc_<dates>.nc` with renumbered tracks. Config key: `stats_dir`. The AR part of the original script is commented out, so only TCs are filtered.

4. Re-run the PyFLEXTRKR pixel-mapping step (`mapfeature_driver`) with the `*_final_extc_` statistics to write the TC-excluded pixel files (`mcstracking_extc_*`).
