#!/bin/bash

#region="europe"
region=$1
startyear=2014
endyear=2018

for iyear in $(seq ${startyear} ${endyear}); do
    echo ${iyear}
    python filter_mcs_tracks_ar.py $region ${iyear}0101_${iyear}1231
done
