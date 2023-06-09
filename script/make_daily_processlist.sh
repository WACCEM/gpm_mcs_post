#!/bin/bash

# This script generates a list of year/month/day for processing MCS daily statistics for TaskFarmer
# The pattern on each line looks like this
#run_map_mcs_stats.sh 20160101 20161231 2016 12 01 asia

#region="apac"
#region="spac"
#region="africa"
#region="afcsam"
# region=$1
config='~/program/gpm_mcs_post/config_gpm_global.yml'
# Start/end year
START=2020
END=2020

# Output list file name
# outfile="processlist_map_mcs_stats_daily_${region}_${START}_${END}"
outfile="processlist_map_mcs_stats_daily_global_${START}_${END}"
#outfile="processlist_map_mcs_stats_daily_${region}"
# Create an empty file
> $outfile

# Loop over year
for year in $(seq $START $END); do
    year1=$((year+1))
    sdate=${year}0101.0000
    edate=${year1}0101.0000
    echo ${sdate}
    # Loop over 365 days
    for d in {0..364}; do 
        echo run_map_mcs_stats_byday.sh ${sdate} ${edate} $(date -d "${year}-01-01 + $d days" +'%Y %m %d';) ${config} >> ${outfile}
    done
done

echo ${outfile}
