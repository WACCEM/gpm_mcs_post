#!/bin/bash

# This script runs the Python code calc_mcs_center_composite_rainfall.py for a number of years
# The Python code uses Dask parallelization, so this script should be run in an exclusive interactive node
conda activate /global/homes/f/feng045/envs/p37

STARTYEAR=$1
ENDYEAR=$2

# Loop over each year
for iyear in $(seq ${STARTYEAR} ${ENDYEAR}); do
    echo ${iyear}
    python calc_mcs_center_composite_rainfall.py ${iyear}0601T00 ${iyear}0901T00
done

