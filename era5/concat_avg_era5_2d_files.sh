#!/bin/bash

# This script concatenates ERA5 spatial mean 2D derived variable files from multiple parts into a single file
# Author: Zhe Feng, zhe.feng@pnnl.gov
# Date: 03/09/2023

# Load E3SM environment, which contains NCO
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

syear=2000
eyear=2020

# Input/output directory
indir="/pscratch/sd/f/feng045/waccem/mcs_global/era5_2d_derived/"
outdir="/pscratch/sd/f/feng045/waccem/mcs_global/era5_2d_derived/mean_300km/"
# Variable names
varnames=(
    "mcs_era5_2D_ENVS_"
)

# Loop over year
for iyear in $(seq $syear $eyear); do
    # Input directory
    echo ${iyear}
    idir=${indir}${iyear}/avg_space/

    # Loop over variables
    for ivar in "${varnames[@]}"; do
        # echo ${ivar}
        # Find all files and put into an array
        readarray -t filenames <<< "$(find ${idir}${ivar}*nc)"
        # Get the first filename
        fn=$(basename ${filenames[0]})
        # Filename format: mcs_era5_Z_20000601.0000_20010101.0000_t20200.nc
        substr=${fn#*${ivar}}   # Remove everything up to and including ${ivar}
        idatetime=${substr%%_t*nc} # Remove everything after the next "_t*nc"
        # Make output filename
        outfn=${outdir}${ivar}${idatetime}.nc
        # echo ${filenames}

        # Concatenante the files
        ncrcat -h -O ${idir}${ivar}*nc ${outfn}
        echo ${outfn}
    done
done