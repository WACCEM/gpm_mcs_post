#!/bin/bash
# Make a job array script to calculate ERA5 3D variables for MCS.
# __author__ = "Zhe.Feng@pnnl.gov"
# __date__ = "31-Aug-2021"

# base_name=("e5.oper.an.pl.128_130_t.ll025sc." "e5.oper.an.pl.128_133_q.ll025sc." "e5.oper.an.pl.128_129_z.ll025sc."
#            "e5.oper.an.pl.128_135_w.ll025sc." "e5.oper.an.pl.128_157_r.ll025sc."
#            "e5.oper.an.pl.128_131_u.ll025uv." "e5.oper.an.pl.128_132_v.ll025uv." 
#            "e5.oper.an.pl.128_155_d.ll025sc." "e5.oper.an.pl.128_060_pv.ll025sc.")
# vars=("T" "Q" "Z"
#       "W" "R" 
#       "U" "V"
#       "D" "PV")

base_name=("e5.oper.an.pl.128_130_t.ll025sc." "e5.oper.an.pl.128_133_q.ll025sc." "e5.oper.an.pl.128_129_z.ll025sc."
           "e5.oper.an.pl.128_131_u.ll025uv." "e5.oper.an.pl.128_132_v.ll025uv." 
           "e5.oper.an.pl.128_135_w.ll025sc." "e5.oper.an.pl.128_157_r.ll025sc.")
vars=("T" "Q" "Z"
      "U" "V"
      "W" "R")
# regions=("asia" "spac" "nam")
regions=("4scream")
# code_name="/global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_3d.py"
# code_name="/global/homes/f/feng045/program/gpm_mcs_post/extract_imerg_mcs_era5_3d.py"
code_name="/global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_3d_preinitiation.py"
config_filebase="/global/u1/f/feng045/program/gpm_mcs_post/config_gpm_era5_3d_"
template="template_slurm_mcs_era5_3d"
# start_year=2000
# end_year=2019
start_year=2020
end_year=2020

# Loop over region
for region in "${regions[@]}"; do
    # Create a emtpy task file
    task_name="tasks_era5_${region}_3d"
    > ${task_name}
    # Slurm job name
    slurm_name="slurm_era5_${region}_3d"

    # Loop over variables
    for ((i = 0; i < ${#base_name[@]}; ++i)); do
        # Loop over year
        for year in $(seq ${start_year} ${end_year}); do  

            # Substitute base_name, vars, date, region in the slurm template
            # date="${year}0101_${year}1231"
            date="${year}0120_${year}0228"

            echo "python ${code_name} ${base_name[$i]} ${vars[$i]} ${config_filebase}${region}.yml ${date}" >> ${task_name}
        done
    done

    sed "s/REGION/"${region}"/g" ${template} > ${slurm_name}
    echo ${task_name}
    echo ${slurm_name}
done