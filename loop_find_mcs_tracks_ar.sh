#!/bin/bash

# This script runs the Python code find_mcs_tracks_in_ar.py for a number of years
# The Python code has Dask parallel built-in and should be run in an interactive node

# conda activate /global/homes/f/feng045/envs/p37
module load python
source activate /global/common/software/m1867/python/testflex

START=2020
END=2020
config_file='config_gpm_global.yml'

# Loop over each year
for iyear in $(seq $START $END); do
#   echo ${iyear}0101_${iyear}1231
    year1=$((${iyear}+1))
    idates=${iyear}0101.0000_${year1}0101.0000
    python find_mcs_tracks_in_ar_coast.py ${idates} ${config_file}
    # EU-NAM
    # python find_mcs_tracks_in_ar_coast.py nam ${iyear}0101_${iyear}1231
done
