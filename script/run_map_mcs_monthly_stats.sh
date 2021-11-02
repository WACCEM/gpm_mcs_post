#!/bin/bash

conda activate /global/cfs/cdirs/m1867/zfeng/envs/py38

#python /global/homes/f/feng045/program/gpm_mcs_post/map_imerg_mcs_stats_bymonth_nospeed.py $1 $2 $3 $4 $5
python /global/homes/f/feng045/program/gpm_mcs_post/map_imerg_mcs_stats_bymonth.py $1 $2 $3 $4 $5
