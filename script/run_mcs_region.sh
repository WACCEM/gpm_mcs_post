#!/bin/bash

source /global/common/cori/software/python/3.6-anaconda-4.4/bin/activate /global/homes/f/feng045/envs/cori-env

python /global/homes/f/feng045/program/waccem/gpm/define_mcs_track_region.py $1 $2
