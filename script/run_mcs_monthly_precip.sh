#!/bin/bash

conda activate /global/homes/f/feng045/envs/py37

python /global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_monthly_precipmap_single.py $1 $2 $3 $4 $5
