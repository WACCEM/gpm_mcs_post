#!/bin/bash

conda activate /global/common/software/m1867/python/py38

# $1: BASE_NAME (e5.oper.an.sfc.128_137_tcwv.ll025sc.)
# $2: VAR_NAME (TCWV)
# $3: DIR_NAME (sfc)
# $4: CONFIG_FILE (config_gpm_era5_2d_spac.yml)
# $5: TRACK_DATES (20190101_20191231)
python /global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_2d.py $1 $2 $3 $4 $5