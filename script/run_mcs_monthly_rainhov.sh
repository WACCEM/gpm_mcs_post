#!/bin/bash

source activate /global/common/software/m1867/python/py38

python /global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_monthly_preciphov_single.py $1 $2 $3 $4 $5
