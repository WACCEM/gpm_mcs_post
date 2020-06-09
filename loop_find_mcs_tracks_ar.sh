#!/bin/bash

# This script runs the Python code find_mcs_tracks_in_ar.py for a number of years
# The Python code has Dask parallel built-in and should be run in an interactive node

conda activate /global/homes/f/feng045/envs/p37

START=2014
END=2019

# Loop over each year
for iyear in $(seq $START $END); do
  echo ${iyear}0101_${iyear}1231
  # EU-NAM
  python find_mcs_tracks_in_ar_coast.py nam ${iyear}0101_${iyear}1231
done
