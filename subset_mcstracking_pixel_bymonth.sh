#!/bin/bash

# This script subsets the MCS tracking pixel-level files to a sub-region and put each month in a separate directory
# Zhe Feng, 07/27/2020

if [[ $# -ne 11 ]] ; then
  echo 'Usage: subset_mcstracking_pixel_bymonth.sh region lon_min lon_max lat_min lat_max start_year end_year start_month end_month input_dir output_dir'
  exit 1
fi

module load nco

# Example input
# region="apac"
# lonmin=110.
# lonmax=160.
# latmin=-30.
# latmax=0.
# syear=2014
# eyear=2014
# smonth=01
# emonth=12
# indir="/global/cscratch1/sd/feng045/waccem/mcs_region/${region}/mcstracking_ccs4_4h/"
# outdir="/global/cscratch1/sd/feng045/waccem/mcs_region/${region}/subset_darwin/"

# Get input arguments
region=$1
lonmin=$2
lonmax=$3
latmin=$4
latmax=$5
syear=$6
eyear=$7
smonth=$8
emonth=$9
indir=${10}
outdir=${11}

# Loop over year
for iyear in $(seq $syear $eyear); do

  # Loop over month
  for imon in $(seq -f "%02g" $smonth $emonth); do

    # Create output directory for each month
    ioutdir="${outdir}/${iyear}/${imon}/"
    mkdir -p ${ioutdir}

    echo $iyear-$imon

    # Find all files in this month
    infiles=$(ls -1 ${indir}/${iyear}0101_${iyear}1231/mcstrack_${iyear}${imon}*nc)
    
    # #echo $infiles
    
    # Loop over list input file
    for file in ${infiles}; do
      # Strip file path to get file name
      fname=$(basename "$file")
      # Add output path to filename
      outname="${ioutdir}${fname}"
      echo "$outname"
      #echo "ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname"
      # Run ncks to subset a region
      ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname
    done

  done
  
done