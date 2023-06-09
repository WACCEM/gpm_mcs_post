#!/bin/bash

# This script runs the Python code filter_mcs_tracks_ar_tc.py for a list of years
# conda activate /global/homes/f/feng045/envs/py37
module load python
source activate /global/common/software/m1867/python/testflex

# declare -a REGIONS=("nam")
STARTYEAR=$1
ENDYEAR=$2
config_file='config_gpm_global.yml'

# Loop over each year
for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
    year1=$((${iyear}+1))
    idates=${iyear}0101.0000_${year1}0101.0000
    echo ${idates}
    python filter_mcs_tracks_ar_tc.py ${idates} ${config_file}
done

# # Loop over each region
# for region in "${REGIONS[@]}"; do

#     # Loop over each year
#     for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
#         echo ${iyear}
#         python filter_mcs_tracks_ar_tc.py $region ${iyear}0101_${iyear}1231
#     done

# done
