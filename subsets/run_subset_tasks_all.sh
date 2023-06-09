#!/bin/bash
# This script runs subsetting tasks for MCS tracking data from an input ASCII file
# To run
# ./run_subset_tasks_all.sh subset_input_region.txt

# Get subset filename from input
subset_file=$1

# Go to subset script directory
cd /global/homes/f/feng045/program/gpm_mcs_post/subsets

# 1. Remove pre-existing temporary directories
echo "Removing existing tmp data in: ${SCRATCH}/tmp/20??"
rm -fr $SCRATCH/tmp/20??

# 2. Create task lists (also make output sub-directories)
echo "Making subset tasks ..."
./make_subset_mcstracking_pixel_tasks.sh ${subset_file}

# 3. Submit the slurm job using jobarray:
sbatch --array=0-19 slurm_subset_jobarray_knl.sh

# 4. Run subset MCS track stats files
echo "Running subset for MCS track stats ..."
./run_subset_robust_mcs_tracks.sh ${subset_file}

# 5. Submit the slurm job to tar the subset data
# sbatch slurm_tarfiles_knl.sh