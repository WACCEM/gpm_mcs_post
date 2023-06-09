#!/bin/bash

base_name=("e5.oper.an.pl.128_130_t.ll025sc." "e5.oper.an.pl.128_133_q.ll025sc." "e5.oper.an.pl.128_157_r.ll025sc."
           "e5.oper.an.pl.128_131_u.ll025uv." "e5.oper.an.pl.128_132_v.ll025uv." "e5.oper.an.pl.128_135_w.ll025sc."
           "e5.oper.an.pl.128_155_d.ll025sc." "e5.oper.an.pl.128_060_pv.ll025sc.")
vars=("T" "Q" "R"
      "U" "V" "W" 
      "D" "PV")
# base_name=("e5.oper.an.pl.128_131_u.ll025uv." "e5.oper.an.pl.128_132_v.ll025uv.")
# vars=("U" "V")
# regions=("asia")
regions=("4scream")

# Loop over region
for region in "${regions[@]}"; do
    # Loop over year
    # for year in {2000..2019}; do
    for year in {2020..2020}; do
        # Loop over variables
        for ((i = 0; i < ${#base_name[@]}; ++i)); do
            # Substitute base_name, vars, date, region in the slurm template
            # date="${year}0101_${year}1231"
            date="${year}0120_${year}0228"
            template="template_slurm_mcs_era5_3d_${region}"
            slurm_name="slurm_era5_${vars[$i]}_${date}_${region}"

            sed "s/BASENAME/"${base_name[$i]}"/g;s/VAR/"${vars[$i]}"/g;s/DATES/"${date}"/g" ${template} > ${slurm_name}

            # Submit the job
            echo ${slurm_name}
            # sbatch ${slurm_name}
        done
    done
done
