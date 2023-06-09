#!/bin/bash
#SBATCH -A m1867
#SBATCH -J stats2d
#SBATCH -t 00:03:00
#SBATCH -q regular
#SBATCH -C cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --output=log_mcs_center_2d_stats.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
#SBATCH --array=1-22

date
source activate /global/common/software/m1867/python/py310

# Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p tasks_calc_mcs_center_2d_stats.txt)
echo $LINE
# Run the line as a command
$LINE

date