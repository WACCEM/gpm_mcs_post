#!/bin/bash
#SBATCH -A m1867
#SBATCH -J 2008
#SBATCH -t 00:10:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --output=log_mcs_mask_2d_2008_%A_%a.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --array=1-35

date
source activate /global/common/software/m1867/python/py310

# Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p tasks_extract_mcs_masks_2d_2008.txt)
echo $LINE
# Run the line as a command
$LINE

date
