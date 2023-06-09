#!/bin/bash

base_name=("e5.oper.an.sfc.128_137_tcwv.ll025sc." "e5.oper.an.sfc.128_167_2t.ll025sc." "e5.oper.an.sfc.128_168_2d.ll025sc."
           "e5.oper.an.sfc.128_134_sp.ll025sc." "e5.oper.an.sfc.128_165_10u.ll025sc." "e5.oper.an.sfc.128_166_10v.ll025sc."
           "e5.oper.an.sfc.128_151_msl.ll025sc." "e5.oper.an.sfc.128_059_cape.ll025sc." "e5.oper.an.sfc.128_034_sstk.ll025sc."
           "e5.oper.an.sfc.128_232_ie.ll025sc." "e5.oper.an.sfc.128_231_ishf.ll025sc." "e5.oper.an.sfc.128_235_skt.ll025sc."
           "e5.oper.an.vinteg.162_084_viwvd.ll025sc.")
vars=("TCWV" "VAR_2T" "VAR_2D"
      "SP" "VAR_10U" "VAR_10V"
      "MSL" "CAPE" "SSTK"
      "IE" "ISHF" "SKT"
      "VIWVD")
dirs=("sfc" "sfc" "sfc" 
      "sfc" "sfc" "sfc" 
      "sfc" "sfc" "sfc"
      "sfc" "sfc" "sfc"
      "vinteg")
regions=("asia")
# regions=("4scream")
# base_name=("e5.oper.an.sfc.128_059_cape.ll025sc.")
# vars=("CAPE")
# dirs=("sfc")

# Loop over region
for region in "${regions[@]}"; do
    # Loop over year
    for year in {2019..2019}; do
        # Loop over variables
        for ((i = 0; i < ${#base_name[@]}; ++i)); do
            # Substitute base_name, vars, date, region in the slurm template
            date="${year}0101_${year}1231"
            # date="${year}0120_${year}0228"
            template="template_slurm_mcs_era5_2d"
            slurm_name="slurm_era5_${vars[$i]}_${date}_${region}"

            sed "s/BASENAME/"${base_name[$i]}"/g;s/VAR/"${vars[$i]}"/g;s/DIR/"${dirs[$i]}"/g;s/DATES/"${date}"/g;s/REGION/"${region}"/g" ${template} > ${slurm_name}

            # Submit the job
            echo ${slurm_name}
            # sbatch ${slurm_name}
        done
    done
done
