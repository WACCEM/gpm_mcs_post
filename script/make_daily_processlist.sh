#!/bin/bash

# This script generates a list of year/month/day for processing MCS daily statistics for TaskFarmer
# The pattern on each line looks like this
#run_map_mcs_stats.sh 20160101 20161231 2016 12 01 asia

#region="apac"
#region="spac"
#region="africa"
#region="afcsam"
region=$1
# Start/end year
START=2016
END=2019

# Output list file name
outfile="processlist_map_mcs_stats_daily_${region}_${START}_${END}"
#outfile="processlist_map_mcs_stats_daily_${region}"
# Create an empty file
> $outfile

# Loop over year
for year in $(seq $START $END); do
  # Loop over 365 days
  for d in {0..364}; do 
    echo run_map_mcs_stats_byday.sh ${year}0101 ${year}1231 $(date -d "${year}-01-01 + $d days" +'%Y %m %d';) ${region} >> ${outfile}
  done
done

echo ${outfile}
