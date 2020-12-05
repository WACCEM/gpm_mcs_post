#!/bin/bash

#run_mcs_monthly_precip.sh 20140101 20141231 2014 3 eus

region=$1
START=2019
END=2019
startlat=-5
endlat=5
startlon=60
endlon=180

runscript="run_mcs_monthly_rainhov.sh"

listname="processlist_mcs_monthly_rainhov_"${region}

# Create an empty file, will overwrite if exists
:> ${listname}

for iyear in $(seq $START $END); do
  for imon in {01..12}; do
    echo ${runscript} ${iyear}0101 ${iyear}1231 ${iyear} $imon $region $startlat $endlat $startlon $endlon >> ${listname}
  done
done

echo ${listname}
