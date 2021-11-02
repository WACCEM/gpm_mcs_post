#!/bin/bash

#module load cray-netcdf
#module load nco
#module load nco/4.7.9-intel

filename=${1##*/}
outfile="$2"/${filename}

lonmin=$3
lonmax=$4
latmin=$5
latmax=$6

ncks -O -Q -d lat,$latmin,$latmax -d lon,$lonmin,$lonmax $1 ${outfile}