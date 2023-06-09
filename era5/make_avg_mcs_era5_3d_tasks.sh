#!/bin/bash

# This script makes taskfarmer lists to average ERA5 3D/2D derived data spatially for each year
# The runscript (run_avg_mcs_era5_space.sh) calls the Python code (avg_mcs_era5_space.py)
# The slurm.submit_avg_era5_template contains taskfarmer setup
# For averaging 3D data: using --nodes=5, --ntasks-per-node=128, < 10 min each year
# For averaging 2D derived data: --nodes=2, --ntasks-per-node=128, < 5 min each year
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 03/09/2023

syear=2000
eyear=2001

submit_job="no"

indir="/pscratch/sd/f/feng045/waccem/mcs_global/era5_3d/"
# indir="/pscratch/sd/f/feng045/waccem/mcs_global/era5_2d_derived/"

basename="mcs_era5_"

runscript="$PWD/run_avg_mcs_era5_space.sh"
slurm_template="$PWD/slurm.submit_avg_era5_template"

# # Create a emtpy task file
# exedir=$PWD
# taskname=${exedir}/tasks_avg_mcs_era5
# > ${taskname}
# echo ${taskname}

# Loop over year
for iyear in $(seq $syear $eyear); do

    # Create a emtpy task file
    exedir=$PWD
    taskname=${exedir}/tasks_avg_mcs_era5_${iyear}
    > ${taskname}
    echo ${taskname}

    # Find all input files for the year
    idir="${indir}${iyear}/"
    infiles=$(ls -1 ${idir}${basename}*.nc)

    # Check if there are files found
    if [[ ${infiles} ]]; then
        # Pass file lists to the task file
        ls ${infiles} | awk -v runscript="$runscript" '{print runscript, $1}' >> ${taskname}

        # Create slurm script
        slurm_fname=${slurm_template}_${iyear}.sh
        sed "s/YEAR/"${iyear}"/g" ${slurm_template} > ${slurm_fname}
        echo ${slurm_fname}

        # Submit slurm script
        if [[ "${submit_job}" == "yes" ]]; then
            sbatch ${slurm_fname}
        fi
    fi
done