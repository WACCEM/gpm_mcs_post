#!/bin/bash

#run_mcs_monthly_precip.sh 20140101 20141231 2014 3 eus

#region="asia"
#region="sio"
#region="africa"
#region="europe"
#region="eunam"
#region="afcsam"
#region="sam"
#region="nam"
#region="npac"
#region="spac"
#region="apac"
region=$1
config=$2
START=2020
END=2020
#runscript="run_map_mcs_stats.sh"
runscript="run_mcs_monthly_precip.sh"
#runscript="run_mcs_monthly_rainhov.sh"
# runscript="run_map_mcs_monthly_stats.sh"
# runscript="run_calc_mcs_monthly_stats.sh"

listname="processlist_mcs_monthly_rain_"${region}
#listname="processlist_mcs_monthly_rainhov_"${region}
# listname="processlist_map_mcs_stats_"${region}
# listname="processlist_calc_mcs_monthly_stats_"${region}

# Create an empty file, will overwrite if exists
:> ${listname}

for iyear in $(seq $START $END); do
  for imon in {01..12}; do
    echo ${runscript} ${iyear}0120 ${iyear}0228 ${iyear} ${imon} ${config} >> ${listname}
    # echo ${runscript} ${iyear}0101 ${iyear}1231 ${iyear} $imon $region >> ${listname}
    # echo ${runscript} ${iyear} $imon $region >> ${listname}
  done
done

echo ${listname}
