#!/bin/bash

#run_mcs_monthly_precip.sh 20140101 20141231 2014 3 eus

config='~/program/gpm_mcs_post/config_wrf_saag.yml'
# config='~/program/gpm_mcs_post/config_gpm_saag.yml'
START=2018
END=2019
#runscript="run_map_mcs_stats.sh"
#runscript="run_mcs_monthly_precip.sh"
runscript="run_mcs_monthly_rainhov.sh"
# runscript="run_map_mcs_monthly_stats.sh"
# runscript="run_calc_mcs_monthly_stats.sh"

listname="processlist_mcs_monthly_rain_wrf_saag"
# listname="processlist_mcs_monthly_rain_"${region}
#listname="processlist_mcs_monthly_rainhov_"${region}
# listname="processlist_map_mcs_stats_"${region}
# listname="processlist_calc_mcs_monthly_stats_"${region}

# Create an empty file, will overwrite if exists
> ${listname}

for iyear in $(seq $START $END); do
  for imon in {01..12}; do
    # echo ${runscript} ${iyear}0120 ${iyear}0228 ${iyear} ${imon} ${config} >> ${listname}
    iyear1=$((iyear+1))
    sdate=${iyear}0601.0000
    edate=${iyear1}0601.0000
    # echo ${sdate}_${edate} ${imon}
    echo ${runscript} ${sdate} ${edate} ${iyear} ${imon} ${config} >> ${listname}
    # echo ${runscript} ${iyear} $imon $region >> ${listname}
  done
done

echo ${listname}
