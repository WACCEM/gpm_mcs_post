#!/bin/bash

# This script calculates multi-year MCS climatology from monthly mean data
# Calls scripts:
#   - calc_mcs_seasonal_mean_rainmap.py
#   - calc_mcs_seasonal_mean_statsmap.py
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 04/05/2022

# Activate conda environment
module load python
source activate /global/common/software/m1867/python/testflex

config_file='config_gpm_global.yml'

python calc_mcs_seasonal_mean_rainmap.py ${config_file}
python calc_mcs_seasonal_mean_statsmap.py ${config_file}

