#!/bin/bash

# This script tars files in each monthly directory to a file in WACCEM project directory.
# Zhe Feng, 07/27/2020

if [[ $# -ne 6 ]] ; then
  echo 'Usage: tar_mcstracking_pixel_bymonth.sh start_year end_year start_month end_month input_dir output_dir'
  exit 1
fi

# Example input
# syear=2014
# eyear=2014
# smonth=01
# emonth=02
# indir="/global/cscratch1/sd/feng045/waccem/mcs_region/apac/subset_darwin/"
# outdir="/global/project/projectdirs/m1867/www/mcs_global/darwin/mcstracking/"

syear=$1
eyear=$2
smonth=$3
emonth=$4
indir=$5
outdir=$6

# Loop over year
for iyear in $(seq $syear $eyear); do

  # Make an output directory for this year
  ioutdir="${outdir}${iyear}/"
  mkdir -p ${ioutdir}

  # Loop over month
  for imon in $(seq -f "%02g" $smonth $emonth); do
    # Go in each monthly directory gets rid of the input dirctory info in the tar file
    cd ${indir}${iyear}/${imon}

#    echo "tar -cvzf ${ioutdir}/${imon}.tar.gz *.nc"
    tar -cvzf ${ioutdir}/${imon}.tar.gz *.nc

    cd -
  done

done
