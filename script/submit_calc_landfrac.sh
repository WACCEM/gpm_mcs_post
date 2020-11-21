#!/bin/bash

startyear=2008
endyear=2012
# region="asia"
region="spac"
# region="nam"

# Loop over each year
for iyear in $(seq $startyear $endyear); do

  sed "s/YEAR/"${iyear}"/g;s/REGION/"${region}"/g" slurm_calc_landfrac_template > slurm_landfrac_${region}_${iyear}.sbat

  echo slurm_landfrac_${region}_${iyear}.sbat
  sbatch slurm_landfrac_${region}_${iyear}.sbat

done