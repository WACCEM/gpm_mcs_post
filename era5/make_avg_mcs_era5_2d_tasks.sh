#!/bin/bash

# This script makes taskfarmer lists to average ERA5 2D data for all years
# The runscript (run_avg_mcs_era5_space.sh) calls the Python code (avg_mcs_era5_space.py)
# The slurm.submit_avg_era5_template contains taskfarmer setup
# For averaging 2D derived data: --nodes=3, --ntasks-per-node=128, export THREADS=64, < 10 min for all years
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 03/10/2023

submit_job="no"

indir="/pscratch/sd/f/feng045/waccem/mcs_global/era5_2d/"

basename="mcs_era5_"

runscript="$PWD/run_avg_mcs_era5_space.sh"
slurm_template="$PWD/slurm.submit_avg_era5_template"

# Create a emtpy task file
exedir=$PWD
taskname=${exedir}/tasks_avg_mcs_era5_2d
> ${taskname}
echo ${taskname}

# Pass file lists to the task file
infiles=$(ls -1 ${indir}${basename}*.nc)
ls ${infiles} | awk -v runscript="$runscript" '{print runscript, $1}' >> ${taskname}

# Create slurm script
slurm_fname=${slurm_template}_2d.sh
sed "s/YEAR/2d/g" ${slurm_template} > ${slurm_fname}
echo ${slurm_fname}

# Submit slurm script
if [[ "${submit_job}" == "yes" ]]; then
    sbatch ${slurm_fname}
fi