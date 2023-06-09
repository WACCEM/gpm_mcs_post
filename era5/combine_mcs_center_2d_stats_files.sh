#!/bin/bash

# Combines MCS center 2D stats files (split to multiple parts) for the same year into a single file

# Load E3SM environment first (nco)
# loade3smpm

syear=2000
eyear=2020

indir="/pscratch/sd/f/feng045/waccem/mcs_global/mcs_center_2d/parts/"
outdir="/pscratch/sd/f/feng045/waccem/mcs_global/mcs_center_stats/"
filebasename="mcs_2d_stats_"

# Loop over year
for year in $(seq $syear $eyear); do
    echo "Working on ... " ${year}
    # Substitute base_name, vars, date, region in the slurm template
    # Get the next year
    year1=$((year+1))
    # Special treatment for 2000 (track starts on 20000601)
    if [ ${year} -eq 2000 ]
    then
        date="${year}0601.0000_${year1}0101.0000"
    else
        date="${year}0101.0000_${year1}0101.0000"
    fi
    outfile=${outdir}/${filebasename}${date}.nc
    # echo "ncrcat -h " ${indir}/${filebasename}${date}_t*.nc ${outdir}/${filebasename}${date}.nc
    ncrcat -h -O ${indir}/${filebasename}${date}_t*.nc ${outfile}
    echo "Output saved: " ${outfile}
done