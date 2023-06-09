#!/bin/bash
# This script makes slurm job array scripts to subset MCS pixel files within a region and period
# The region and period are read from an ASCII file "input.txt"
# To run this script:
# 1. Remove pre-existing temporary directories
# rm -fr $scratchdir/20??
# 2. Create task lists (also make output sub-directories)
# make_subset_mcstracking_pixel_tasks.sh input.txt
# 3. Submit the slurm job using jobarray:
# sbatch --array=0-19 slurm_subset_jobarray_knl.sh
# 4. Submit the slurm job to tar the subset data
# sbatch slurm_tarfiles_knl.sh

# Read each line into var1, var2, ...
# The variables must be in the exact order to work
eval $(awk '{print "var"NR"="$2}' $1)
# Assign to specific variables
lonmin=$var1
lonmax=$var2
latmin=$var3
latmax=$var4
syear=$var5
eyear=$var6
smonth=$var7
emonth=$var8
indir=$var9
outdir=$var10

# Scratch directory to store the temporary subset files
scratchdir="/global/cscratch1/sd/feng045/tmp/"
exedir=$PWD

# Create a tar script file
tar_scriptname=${exedir}/tasks_tar_files
> ${tar_scriptname}
# Create the directory for the final output tar file
mkdir -p ${outdir}

jobnum=0
# Loop over year
for iyear in $(seq $syear $eyear); do

    # Put tar job in tar script file
    echo "tar -cvzf ${outdir}/${iyear}.tar.gz ${scratchdir}/${iyear}" >> ${tar_scriptname}

    # Create a emtpy task file
    taskname=${exedir}/tasks_${jobnum}
    > ${taskname}

    echo $iyear ${taskname}

    # Make start/end date for each year
    iyear1=$((iyear+1))
    sdate=${iyear}0101.0000
    edate=${iyear1}0101.0000

    # Loop over month
    for imon in $(seq -f "%02g" $smonth $emonth); do

        # Create output directory for each month in the scratch directory
        ioutdir="${scratchdir}/${iyear}/${imon}/"
        mkdir -p ${ioutdir}

        # Find all files in this month
        # infiles=$(ls -1 ${indir}/${iyear}0101_${iyear}1231/mcstrack_${iyear}${imon}*nc)
        infiles=$(ls -1 ${indir}/${sdate}_${edate}/mcstrack_${iyear}${imon}*nc)

        # Check if there are files found
        if [[ ${infiles} ]]; then
            # # Create a emtpy task file
            # taskname=${exedir}/tasks_${jobnum}
            # > ${taskname}

            # echo $iyear-$imon ${taskname}

            # Pass file lists to the task file
            ls ${infiles} | awk -v dir="$ioutdir" -v lonmin="$lonmin" -v lonmax="$lonmax" -v latmin="$latmin" -v latmax="$latmax" \
                '{print "ncks_pixel_file.sh", $1, dir, lonmin, lonmax, latmin, latmax}' >> ${taskname}
            
        # else
        #     echo "No files found for ${iyear}-${imon}"
        fi
    done
    ((jobnum=jobnum+1))
done

echo "Tar script file: ${tar_scriptname}"