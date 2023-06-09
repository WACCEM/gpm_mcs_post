#!/bin/bash

# This script runs the Python code find_mcs_tracks_in_tc.py for a list of regions and a number of years
# The Python code has Dask parallel built-in and should be run in an interactive node

# conda activate /global/homes/f/feng045/envs/p37
module load python
source activate /global/common/software/m1867/python/testflex

# declare -a REGIONS=("nam" "asia" "spac")
# declare -a REGIONS=("asia")
STARTYEAR=2002
ENDYEAR=2020
config_file='config_gpm_global.yml'

# Loop over each year
for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
    # echo ${region} ${iyear}0101_${iyear}1231
    # python find_mcs_tracks_in_tc.py ${region} ${iyear}0101_${iyear}1231
    year1=$((${iyear}+1))
    idates=${iyear}0101.0000_${year1}0101.0000
    echo ${idates}
    python find_mcs_tracks_in_tc.py ${idates} ${config_file}
done

# # Loop over each region
# for region in "${REGIONS[@]}"; do

#     # Loop over each year
#     for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
#         # echo ${region} ${iyear}0101_${iyear}1231
#         # python find_mcs_tracks_in_tc.py ${region} ${iyear}0101_${iyear}1231
#         year1=$((${iyear}+1))
#         idates=${iyear}0101.0000_${year1}0101.0000
#         echo ${idates}
#         # python find_mcs_tracks_in_tc.py ${idates} ${config_file}
#     done

# done
