#!/bin/bash
source activate regrid
echo "start time: " $(date)

drc_in_root="/pscratch/sd/w/wu59/pyflextrkr/mcstracking"
drc_out_root="/pscratch/sd/w/wu59/pyflextrkr/regrid/nco-out/mcstracking_regrid025"
grid_source="/global/homes/w/wu59/Flex/regrid/regrid_nco/grids_481x1440.nc"

# loop for each year-each subdirectory passing through
for sub_dir in $1 ; do
    drc_in="${drc_in_root}/${sub_dir}"
    drc_pr="${drc_out_root}/pr_out/${sub_dir}"
    drc_tb="${drc_out_root}/tb_out/${sub_dir}"
    drc_mcs="${drc_out_root}/mcs_out/${sub_dir}"

    # conservative remap for precip
    ncremap --no_stdin -v precipitation --alg_typ=conserve -g $grid_source \
    -I $drc_in -O $drc_pr &

    # biliner remap for tb
    ncremap --no_stdin -v tb -a bilinear -g $grid_source \
    -I $drc_in -O $drc_tb &

    # nearest interpolation for all other variables-most time-timeconsuming step
    ncremap --mem_mb=0 --no_stdin -v cloudtracknumber,pcptracknumber,cloudtracknumber_nomergesplit,merge_tracknumbers,split_tracknumbers,cloudnumber,cloudtype \
    -a neareststod -g $grid_source \
    -I $drc_in -O $drc_mcs &
    # wait for all the regrid
    wait

    # loop for every single hour
    for f in ${drc_pr}/*.nc ;do

          pr_file=$f
          tb_file=${f/pr_out/tb_out}
          mcs_file=${f/pr_out/mcs_out}
          #merge pr and tb file to rewritten tb file
          ncks -A $pr_file $tb_file
          #merge file to mcs file
          ncks -A $tb_file $mcs_file
          #remove some extra variables- final output
          ncks -O -x -v lat_bnds,lon_bnds,gw,area $mcs_file $mcs_file
    done
    rm -rf $drc_pr  $drc_tb
done
echo "end time: " $(date)
