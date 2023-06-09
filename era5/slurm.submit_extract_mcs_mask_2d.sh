#!/bin/bash
#SBATCH -A m1867
#SBATCH -J mcs_mask
#SBATCH -t 02:00:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --output=log_mcs_mask_2d_%A_%a.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov

# To submit the job for a specific range of parts:
# sbatch --array=<indexlist>[%<limit>] slurm_jobarray.sh
# For example, to submit all 20 jobs, but limit to running 2 at a time:
# sbatch --array=1-21%2 slurm_jobarray.sh
# For running a number of specific jobs:
# sbatch --array=3,5,12,15-18%6 slurm_jobarray.sh

date
source activate /global/common/software/m1867/python/py310

# Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p tasks_extract_mcs_masks_2d)
echo $LINE
# Run the line as a command
$LINE

date