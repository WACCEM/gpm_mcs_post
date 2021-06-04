#!/bin/bash

# This script calls subset robust MCS tracks Python code for a specified region and each year
# Zhe Feng, 12/08/2020

if [[ $# -ne 8 ]] ; then
  echo 'Usage: run_subset_robust_mcs_tracks.sh lon_min lon_max lat_min lat_max start_year end_year input_dir output_dir'
  exit 1
fi

conda activate /global/cfs/cdirs/m1867/zfeng/envs/py38

lonmin=$1
lonmax=$2
latmin=$3
latmax=$4
syear=$5
eyear=$6
indir=$7
outdir=$8
# indir="/global/cscratch1/sd/feng045/waccem/mcs_region/spac/stats_ccs4_4h/robust/filtered/"
# outdir = f'/global/cscratch1/sd/feng045/waccem/mcs_region/spac/subset_amazon/'

# Loop over year
for iyear in $(seq $syear $eyear); do
  infile=$(ls -1 ${indir}robust_mcs_tracks_extc_${iyear}*nc)
  echo $infile
  python subset_robust_mcs_tracks.py $lonmin $lonmax $latmin $latmax $infile $outdir
done