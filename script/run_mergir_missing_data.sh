#!/bin/bash

conda activate /global/cfs/cdirs/m1867/zfeng/envs/py38
# conda activate /global/homes/f/feng045/envs/p37

python /global/homes/f/feng045/program/gpm_mcs_post/calc_mergir_missing_data.py $1 $2
