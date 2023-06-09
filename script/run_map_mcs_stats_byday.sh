#!/bin/bash

# source activate /global/common/software/m1867/python/py38
source activate /global/common/software/m1867/python/testflex

python /global/homes/f/feng045/program/gpm_mcs_post/map_imerg_mcs_stats_byday.py $1 $2 $3 $4 $5 $6
#python /global/homes/f/feng045/program/gpm_mcs_post/map_imerg_mcs_stats_byday_nospeed.py $1 $2 $3 $4 $5 $6
