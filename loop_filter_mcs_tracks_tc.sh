#!/bin/bash

# This script runs the Python code filter_mcs_tracks_tc.py for a list of regions and a number of years
# For unknown reason, writing the output netCDF file in newer version of Xarray (0.15.0+) result in segmentation fault.
# For now, run this in Xarray (0.14.0) and Python 3.7 or lower
conda activate /global/homes/f/feng045/envs/py37

declare -a REGIONS=("apac")
STARTYEAR=2016
ENDYEAR=2016

# Loop over each region
for region in "${REGIONS[@]}"; do

    # Loop over each year
    for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
        echo ${iyear}
        python filter_mcs_tracks_tc.py $region ${iyear}0101_${iyear}1231
    done

done
