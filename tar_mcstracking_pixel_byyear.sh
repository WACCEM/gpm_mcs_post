#!/bin/bash

# This script tars files in each year directory to a file in WACCEM project directory.
# Zhe Feng, 07/27/2020

region=$1
#region="spac"
#region="asia"

indir="/global/cscratch1/sd/liunana/IR_IMERG_Combined/mcs_region/${region}/mcstracking_ccs4_4h/"
outdir="/global/project/projectdirs/m1867/zfeng/gpm/mcs_region/${region}/mcstracking_ccs_4h/"
mkdir -p ${outdir}

# Loop over year
for iyear in {2000..2019}; do

    # Go in each monthly directory gets rid of the input dirctory info in the tar file
    cd ${indir}${iyear}0101_${iyear}1231/

#    echo "tar -cvzf ${outdir}/${iyear}0101_${iyear}1231.tar.gz *.nc"
    tar -cvzf ${outdir}/${iyear}0101_${iyear}1231.tar.gz *.nc

    cd -

done
