#!/bin/bash

# This script makes Job Array tasks to calculate MCS center 2D stats for all years
# The Python code: calc_imerg_mcs_center_2d_stats.py
# The slurm.submit_calc_mcs_center_stats_template contains the Job Array setup
# For calculating a single file containing 1000 tracks, it takes ~1 min on a Perlmutter CPU node (running 64 CPU)
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 05/24/2023

submit_job="no"

# Input data
indir='/pscratch/sd/f/feng045/waccem/mcs_global/mcs_center_2d/parts/'
basename='mcs_2d_mask_'

# Python code name
codename='calc_imerg_mcs_center_2d_stats.py'

# Slurm script
slurm_template='slurm.submit_calc_mcs_center_stats_template'

# Create a emtpy task file
exedir=$PWD
taskname=${exedir}/tasks_calc_mcs_center_2d_stats.txt
> ${taskname}
echo ${taskname}

# Pass file lists to the task file
# infiles=$(ls -1 ${indir}${basename}*.nc)
infiles=$(ls -1 ${indir}${basename}2000*.nc)
# Remove directory, pass to awk to make the commands
ls ${infiles} | xargs -n 1 basename | awk -v codename="$codename" '{print "python", codename, $1}' >> ${taskname}

# Get the number of files (tasks)
nfiles=$(ls -1 ${infiles} | wc -l)
echo "Number of tasks: " ${nfiles}

# Replace NTASKS with nfiles in the slurm script
slurm_fname=${slurm_template}_2d.sh
sed "s/NTASKS/"${nfiles}"/g" ${slurm_template} > ${slurm_fname}
echo ${exedir}/${slurm_fname}

# Submit slurm script
if [[ "${submit_job}" == "yes" ]]; then
    sbatch ${slurm_fname}
fi