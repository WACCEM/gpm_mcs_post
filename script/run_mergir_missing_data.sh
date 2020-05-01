#!/bin/bash

#source /global/common/cori/software/python/3.6-anaconda-4.4/bin/activate /global/homes/f/feng045/envs/py36
conda activate /global/homes/f/feng045/envs/py37

python /global/homes/f/feng045/program/waccem/gpm/calc_mergir_missing_data.py $1 $2
