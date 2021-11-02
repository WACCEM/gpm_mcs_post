#!/bin/bash
# Make a job array script to calculate ERA5 2D variables for MCS.
# __author__ = "Zhe.Feng@pnnl.gov"
# __date__ = "31-Aug-2021"

base_name=("e5.csf." "e5.oper.an.sfc.128_137_tcwv.ll025sc."  "e5.oper.an.sfc.128_059_cape.ll025sc." 
           "e5.oper.an.sfc.128_034_sstk.ll025sc." "e5.oper.an.sfc.128_235_skt.ll025sc."
           "e5.oper.an.sfc.128_167_2t.ll025sc." "e5.oper.an.sfc.128_168_2d.ll025sc."
           "e5.oper.an.sfc.128_134_sp.ll025sc." "e5.oper.an.sfc.128_151_msl.ll025sc." 
           "e5.oper.an.sfc.128_165_10u.ll025sc." "e5.oper.an.sfc.128_166_10v.ll025sc."           
           "e5.oper.an.sfc.128_232_ie.ll025sc." "e5.oper.an.sfc.128_231_ishf.ll025sc." 
           "e5.oper.an.vinteg.162_084_viwvd.ll025sc.")
vars=("CSF" "TCWV" "CAPE" 
      "SSTK" "SKT"
      "VAR_2T" "VAR_2D"
      "SP" "MSL" 
      "VAR_10U" "VAR_10V" 
      "IE" "ISHF" 
      "VIWVD")
dirs=("vinteg" "sfc" "sfc" 
      "sfc" "sfc" 
      "sfc" "sfc" 
      "sfc" "sfc" 
      "sfc" "sfc" 
      "sfc" "sfc"
      "vinteg")
regions=("asia" "spac" "nam")
# base_name=("e5.oper.an.sfc.128_137_tcwv.ll025sc." "e5.oper.an.sfc.128_059_cape.ll025sc." 
#            "e5.oper.an.sfc.128_034_sstk.ll025sc." "e5.oper.an.sfc.128_235_skt.ll025sc."
#            "e5.oper.an.sfc.128_167_2t.ll025sc." "e5.oper.an.sfc.128_168_2d.ll025sc."
#            "e5.oper.an.sfc.128_134_sp.ll025sc." "e5.oper.an.sfc.128_151_msl.ll025sc." 
#            "e5.oper.an.sfc.128_165_10u.ll025sc." "e5.oper.an.sfc.128_166_10v.ll025sc."           
#            "e5.oper.an.sfc.128_232_ie.ll025sc." "e5.oper.an.sfc.128_231_ishf.ll025sc."
#            "e5.oper.an.vinteg.162_084_viwvd.ll025sc.")
# vars=("TCWV" "CAPE" 
#       "SSTK" "SKT"
#       "VAR_2T" "VAR_2D"
#       "SP" "MSL" 
#       "VAR_10U" "VAR_10V" 
#       "IE" "ISHF"
#       "VIWVD")
# dirs=("sfc" "sfc" 
#       "sfc" "sfc" 
#       "sfc" "sfc" 
#       "sfc" "sfc" 
#       "sfc" "sfc" 
#       "sfc" "sfc"
#       "vinteg")
# regions=("4scream")
# code_name="/global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_2d.py"
code_name="/global/homes/f/feng045/program/gpm_mcs_post/calc_imerg_mcs_era5_2d_preinitiation.py"
config_filebase="/global/u1/f/feng045/program/gpm_mcs_post/config_gpm_era5_2d_"
template="template_slurm_mcs_era5_2d"
start_year=2000
end_year=2019
# start_year=2020
# end_year=2020
# base_name=("e5.oper.an.sfc.128_059_cape.ll025sc.")
# vars=("CAPE")
# dirs=("sfc")

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

            # Substitute base_name, vars, date, region in the slurm template
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