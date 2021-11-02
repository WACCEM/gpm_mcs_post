#!/bin/bash
#SBATCH -A m2645
#SBATCH --job-name=subset
#SBATCH -p debug
#SBATCH -N 2
#SBATCH -c 64
#SBATCH --time=00:30:00
#SBATCH -L SCRATCH,project
#SBATCH -C haswell
#SBATCH --output=log_%A_%a.log

# To submit the job for a specific range of parts:
# sbatch --array <indexlist>[%<limit>] slurm_subset_jobarray.sh

# For example, to submit all 20 jobs, but limit to running 2 at a time:
# sbatch --array [0-19]%2 slurm_subset_jobarray.sh
# For running a number of specific jobs:
# sbatch --array [3,5,12,30-33]%6 slurm_subset_jobarray.sh

date
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=64
runcommands.sh tasks_$SLURM_ARRAY_TASK_ID
date
