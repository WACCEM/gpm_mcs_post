#!/bin/bash
# Make a job array script to extract MCS center 2D masks

code_name="/global/homes/f/feng045/program/gpm_mcs_post/era5/extract_mcs_masks_2d_fullhistory.py"
config_name="/global/homes/f/feng045/program/gpm_mcs_post/era5/config_gpm_mask_2d_global.yml"

template="template_slurm_mcs_masks_2d"
start_year=2000
end_year=2020

# Create a emtpy task file
task_name="tasks_extract_mcs_masks_2d"
> ${task_name}
# Slurm job name
# slurm_name="slurm.submit_extract_mcs_mask_2d.sh"

# Loop over year
for year in $(seq ${start_year} ${end_year}); do  

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

    echo "python ${code_name} ${config_name} ${date}" >> ${task_name}

    # let "tasknum+=1"
    # echo ${year}: ${task_name}
done

# sed "s/REGION/"${region}"/g" ${template} > ${slurm_name}
echo ${task_name}
# echo ${slurm_name}
