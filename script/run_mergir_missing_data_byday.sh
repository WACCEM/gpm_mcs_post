#!/bin/bash

conda activate /global/cfs/cdirs/m1867/zfeng/envs/py38
# conda activate /global/homes/f/feng045/envs/p37

for iyear in {2020..2020}; do
  for imon in {01..01}; do
    for iday in {01..31}; do
      echo ${iyear} ${imon} ${iday}
      python /global/homes/f/feng045/program/gpm_mcs_post/calc_mergir_missing_data_byday.py ${iyear} ${imon} ${iday}
    done
  done
done
