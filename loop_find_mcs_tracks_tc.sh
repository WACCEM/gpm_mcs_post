#!/bin/bash

# This script runs the Python code find_mcs_tracks_in_tc.py for a list of regions and a number of years
# The Python code has Dask parallel built-in and should be run in an interactive node

conda activate /global/homes/f/feng045/envs/p37

# declare -a REGIONS=("nam" "asia" "spac")
declare -a REGIONS=("apac")
STARTYEAR=2014
ENDYEAR=2019

# Loop over each region
for region in "${REGIONS[@]}"; do

    # Loop over each year
    for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
        echo ${region} ${iyear}0101_${iyear}1231
        python find_mcs_tracks_in_tc.py ${region} ${iyear}0101_${iyear}1231
    done

done
