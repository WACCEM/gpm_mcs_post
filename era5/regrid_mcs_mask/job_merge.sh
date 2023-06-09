#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 05:30:00

module load parallel
source activate regrid
srun parallel --jobs 21 bash ./allstep_nco.sh {} < subdirlist.txt
