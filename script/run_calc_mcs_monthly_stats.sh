#!/bin/bash

conda activate /global/homes/f/feng045/envs/p37

# python /global/homes/f/feng045/program/gpm_mcs_post/calc_mcs_monthly_stats_from_daily_nospeed.py $1 $2 $3

python /global/homes/f/feng045/program/gpm_mcs_post/calc_mcs_monthly_stats_from_daily.py $1 $2 $3
