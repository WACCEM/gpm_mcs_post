#!/bin/bash

# This script runs calc_mcs_seasonal_mean_rainmap.py for a list of regions
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 03/14/2020

# Activate conda environment
conda activate /global/homes/f/feng045/envs/p37

## declare an array containing region names
#declare -a regions=("asia" "sio" "europe" "africa" "nam" "sam" "npac" "spac")
#declare -a regions=("apac" "afcsam" "eunam" "npac" "spac")
#declare -a regions=("apac" "afcsam" "eunam" "npac" "spac")
declare -a regions=("asia")

## now loop through the list of regions
for ireg in "${regions[@]}"; do
   echo "$ireg"
   python calc_mcs_seasonal_mean_rainmap.py ${ireg}
   python calc_mcs_seasonal_mean_statsmap.py ${ireg}
done
