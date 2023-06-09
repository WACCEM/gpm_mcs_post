#!/bin/bash
#SBATCH -A m1867
#SBATCH -J ERA_2d
#SBATCH -t 00:10:00
#SBATCH -q debug
#SBATCH -C cpu
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=128
#SBATCH --exclusive
#SBATCH --output=log_avg_era5_3d_2d.log
#SBATCH --mail-type=END
#SBATCH --mail-user=zhe.feng@pnnl.gov
date
module load taskfarmer
export THREADS=64
cd /global/homes/f/feng045/program/gpm_mcs_post/era5
runcommands.sh tasks_avg_mcs_era5_2d
date