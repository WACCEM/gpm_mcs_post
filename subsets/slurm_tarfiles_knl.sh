#!/bin/bash
#SBATCH -A m2645
#SBATCH --job-name=tarfiles
#SBATCH -p regular
#SBATCH -N 2
#SBATCH -c 64
#SBATCH --time=03:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C knl
#SBATCH --output=log_tar_files.log
date
export THREADS=64
runcommands.sh tasks_tar_files
date
