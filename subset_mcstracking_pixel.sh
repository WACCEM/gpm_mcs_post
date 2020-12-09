#!/bin/bash

# This script subsets the MCS tracking pixel-level files to a sub-region
# Zhe Feng, 04/13/2020

module load nco

inyear=2019
inmonth=07
inday=0?

#lonmin=100.
#lonmax=130.
#latmin=-15.
#latmax=5.

lonmin=-15.
lonmax=40.
latmin=-15.
latmax=25.

region="spac"

datadir="/global/cscratch1/sd/feng045/waccem/mcs_region/${region}/mcstracking_ccs4_4h/"
infiles=$(ls -1 ${datadir}${inyear}0101_${inyear}1231/mcstrack_${inyear}${inmonth}${inday}*nc)

#outpath="/global/cscratch1/sd/feng045/waccem/mcs_region/apac/subset4latos/"
outpath="/global/cscratch1/sd/feng045/waccem/mcs_region/${region}/subset_africa/"
mkdir -p "$outpath"

#echo $infiles

# Loop over list input file
for file in ${infiles}; do
  # Strip file path to get file name
  fname=$(basename "$file")
  # Add output path to filename
  #outname="${outpath}IMERG_${fname}"
  outname="${outpath}${fname}"
  echo "$outname"
  # echo "ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname"
  # Run ncks to subset a region
  ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $file -o $outname
done
