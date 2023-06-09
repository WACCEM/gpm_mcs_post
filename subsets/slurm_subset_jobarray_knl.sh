#!/bin/bash
#SBATCH -A m2637
#SBATCH --job-name=subset
#SBATCH -C knl
#SBATCH -q regular
#SBATCH -N 2
#SBATCH -c 100
#SBATCH --time=12:00:00
#SBATCH -L SCRATCH,project
#SBATCH --output=log_%A_%a.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov

# To submit the job for a specific range of parts:
# sbatch --array=<indexlist>[%<limit>] slurm_subset_jobarray_knl.sh

# For example, to submit all 20 jobs, but limit to running 2 at a time:
# sbatch --array=0-19%2 slurm_subset_jobarray_knl.sh
# For running a number of specific jobs:
# sbatch --array=3,5,12,30-33%6 slurm_subset_jobarray_knl.sh

date
cd /global/homes/f/feng045/program/gpm_mcs_post/subsets

export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=100
runcommands.sh tasks_$SLURM_ARRAY_TASK_ID
date
