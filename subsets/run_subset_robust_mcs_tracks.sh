#!/bin/bash

# This script calls subset robust MCS tracks Python code for a specified region and each year
# The region and period are read from an ASCII file "input.txt"
# To run this script:
# >run_subset_robust_mcs_tracks.sh input.txt

conda activate /global/common/software/m1867/python/py38

# Read each line into var1, var2, ...
# The variables must be in the exact order to work
eval $(awk '{print "var"NR"="$2}' $1)
# Assign to specific variables
lonmin=$var1
lonmax=$var2
latmin=$var3
latmax=$var4
syear=$var5
eyear=$var6
smonth=$var7
emonth=$var8
# indir=$var9
# outdir=$var10
indir=$var11
outdir=$var12

# Loop over year
for iyear in $(seq $syear $eyear); do
  infile=$(ls -1 ${indir}/robust_mcs_tracks_extc_${iyear}*nc)
  echo $infile
  python subset_robust_mcs_tracks.py $lonmin $lonmax $latmin $latmax $infile $outdir
done

# Tar all track stats files into a single file
cd ${outdir}
tar -cvzf robust_mcs_track_stats.tar.gz robust_mcs_tracks_*.nc
cd -