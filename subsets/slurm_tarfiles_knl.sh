#!/bin/bash
#SBATCH -A m2637
#SBATCH --job-name=tarfiles
##SBATCH -p regular
#SBATCH -p debug
#SBATCH -N 2
#SBATCH -c 64
#SBATCH --time=00:30:00
#SBATCH -L SCRATCH,project
#SBATCH -C knl
#SBATCH --output=log_tar_files.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
date
cd /global/homes/f/feng045/program/gpm_mcs_post/subsets
export PATH=$PATH:/usr/common/tig/taskfarmer/1.5/bin:$(pwd)
export THREADS=64
runcommands.sh tasks_tar_files
date
