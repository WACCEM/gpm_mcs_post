#!/bin/bash

# Subset the MCS tracking pixel-level files to a sub-region and compress the files in each month to a tar file.
# Zhe Feng, 06/04/2021

if [[ $# -ne 10 ]] ; then
  echo 'Usage: subset_tar_mcstracking_pixel_bymonth.sh lon_min lon_max lat_min lat_max start_year end_year start_month end_month input_dir output_dir'
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
lonmin=$1
lonmax=$2
latmin=$3
latmax=$4
syear=$5
eyear=$6
smonth=$7
emonth=$8
indir=${09}
outdir=${10}

# Scratch directory to store the temporary subset files
scratchdir="/global/cscratch1/sd/feng045/tmp/"

# Loop over year
for iyear in $(seq $syear $eyear); do

  # Loop over month
  for imon in $(seq -f "%02g" $smonth $emonth); do

    # Create output directory for each month in the scratch directory
    iscratchdir="${scratchdir}/${iyear}/${imon}/"
    mkdir -p ${iscratchdir}

    echo $iyear-$imon

    # Find all files in this month
    infiles=$(ls -1 ${indir}/${iyear}0101_${iyear}1231/mcstrack_${iyear}${imon}*nc)
    
    # echo $infiles
    
    # Loop over list input file
    for file in ${infiles}; do
      # Strip file path to get file name
      fname=$(basename "$file")
      # Add output path to filename
      outname="${iscratchdir}${fname}"
      echo "$outname"
      #echo "ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname"
      # Run ncks to subset a region
      ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname
    done

    # Make an output directory for this year
    ioutdir="${outdir}${iyear}/"
    mkdir -p ${ioutdir}

    # Compress the files for this month
    # echo "tar -cvzf ${ioutdir}${imon}.tar.gz ${iscratchdir}*.nc"
    tar -cvzf ${ioutdir}${imon}.tar.gz ${iscratchdir}*.nc

  done
  
done