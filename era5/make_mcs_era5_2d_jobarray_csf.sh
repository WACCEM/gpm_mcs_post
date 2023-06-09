#!/bin/bash
# Make a job array script to calculate ERA5 2D variables for MCS.
# __author__ = "Zhe.Feng@pnnl.gov"
# __date__ = "3-May-2022"

regions=("asia" "spac" "nam")
base_name=("e5.csf.")
vars=("CSF")
dirs=("sfc")
# regions=("4scream")
code_name="/global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_2d_preinitiation.py"
config_filebase="/global/u1/f/feng045/program/gpm_mcs_post/config_gpm_era5_2d_"
template="template_slurm_mcs_era5_2d"
start_year=2000
end_year=2019

# tasknum=0
# Loop over region
for region in "${regions[@]}"; do
    # Create a emtpy task file
    task_name="tasks_era5_${region}_2d"
    > ${task_name}
    # Slurm job name
    slurm_name="slurm_era5_${region}_2d"

    # Loop over variables
    for ((i = 0; i < ${#base_name[@]}; ++i)); do
        # Loop over year
        for year in $(seq ${start_year} ${end_year}); do  

            # # Substitute base_name, vars, date, region in the slurm template
            # # Get the next year
            # year1=$((year+1))
            # # Special treatment for 2000 (track starts on 20000601)
            # if [ ${year} -eq 2000 ]
            # then
            #     # date="${year}0601.0000_${year1}0101.0000"
            #     date="${year}0601_${year}1231"
            # else
            #     # date="${year}0101.0000_${year1}0101.0000"
            #     date="${year}0101_${year}1231"
            # fi
            date="${year}0101_${year}1231"
            # date="${year}0120_${year}0228"

            # echo "run_gpm_era5_2d.sh ${base_name[$i]} ${vars[$i]} ${dirs[$i]} ${config_filebase}${region}.yml ${date}" >> ${task_name}
            echo "python ${code_name} ${base_name[$i]} ${vars[$i]} ${dirs[$i]} ${config_filebase}${region}.yml ${date}" >> ${task_name}

            # let "tasknum+=1"
            # echo ${year}: ${task_name}
        done
    done

    sed "s/REGION/"${region}"/g" ${template} > ${slurm_name}
    echo ${task_name}
    echo ${slurm_name}
done