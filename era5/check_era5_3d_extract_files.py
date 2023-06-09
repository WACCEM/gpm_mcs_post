"""
Check the extracted ERA5 3D files, make a new task list and slurm script for 
the missing files, or files too small, or files with incomplete tracks.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "25-Jun-2022"
import numpy as np
import xarray as xr
import pandas as pd
import yaml
import enum
import os, sys
import textwrap
import subprocess

#----------------------------------------------------------------------------
class SIZE_UNIT(enum.Enum):
    # Enum for size units
    BYTES = 1
    KB = 2
    MB = 3
    GB = 4

#----------------------------------------------------------------------------
def convert_unit(size_in_bytes, unit):
    """ Convert the size from bytes to other units like KB, MB or GB"""
    if unit == SIZE_UNIT.KB:
        return size_in_bytes/1024
    elif unit == SIZE_UNIT.MB:
        return size_in_bytes/(1024*1024)
    elif unit == SIZE_UNIT.GB:
        return size_in_bytes/(1024*1024*1024)
    else:
        return size_in_bytes

#----------------------------------------------------------------------------
def get_file_size(file_name, size_type = SIZE_UNIT.BYTES):
    """ Get file in size in given unit like KB, MB or GB"""
    size = os.path.getsize(file_name)
    return convert_unit(size, size_type)


if __name__ == "__main__":

    # Start/end years to process
    start_year = 2002
    end_year = 2002

    code_dir = '/global/homes/f/feng045/program/gpm_mcs_post/script/'
    code_name = f'extract_mcs_era5_3d_fullhistory_bytracks.py'
    config_file = f'config_gpm_era5_3d_global.yml'
    slurm_basename = f'slurm.submit_era5_3d_'
    data_dir = '/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_3d/'

    # Number of jobs allowed to run simultaneously in job array
    # njobs_run = 60
    # Submit slurm job
    submit_job = False

    # Minimum expected file size [MB]
    min_filesize_MB = 70

    # Number of tracks to process per part
    ntracks_part = 100
    # Set the number of digits for 0 padding
    # This should be set to the digit for the maximum number of tracks
    # e.g., ntracks = 32138, digits = 5
    digits = 5

    # ERA5 variable name and file base names
    var_dict = {
        "T": "e5.oper.an.pl.128_130_t.ll025sc.",
        "Q": "e5.oper.an.pl.128_133_q.ll025sc.", 
        "Z": "e5.oper.an.pl.128_129_z.ll025sc.",
        "U": "e5.oper.an.pl.128_131_u.ll025uv.", 
        "V": "e5.oper.an.pl.128_132_v.ll025uv.",
        "W": "e5.oper.an.pl.128_135_w.ll025sc.", 
        "R": "e5.oper.an.pl.128_157_r.ll025sc.",
    }

    # Make start/end dates for tracking periods
    # start_dates = ['20000601.0000']
    start_dates = pd.date_range(f'{start_year}-01-01', f'{end_year}-01-01', freq='1YS').strftime('%Y%m%d.%H%M')
    end_dates = pd.date_range(f'{start_year+1}-01-01', f'{end_year+1}-01-01', freq='1YS').strftime('%Y%m%d.%H%M')

    # Loop over each period
    for ii in range(len(start_dates)):
        track_period = f'{start_dates[ii]}_{end_dates[ii]}'

        # Get inputs from configuration file
        stream = open(config_file, 'r')
        config = yaml.full_load(stream)
        stats_dir = config['stats_dir']
        mcsfile_basename = config['mcsfile_basename']
        ntimes_max = config['ntimes_max']

        # Read robust MCS statistics
        mcs_file = f"{stats_dir}{mcsfile_basename}{track_period}.nc"
        dsm = xr.open_dataset(mcs_file)
        # Subset MCS times to reduce array size
        # Most valid MCS data are within 0:ntimes_max
        dsm = dsm.isel(times=slice(0, ntimes_max))
        ntracks_all = dsm.sizes['tracks']

        # Number of parts
        nparts = np.floor(ntracks_all / ntracks_part).astype(int)
        # Make a list for track start/end 
        track_start = []
        track_end = []
        for ii in range(0, nparts): 
            track_start.append(ii*ntracks_part)
            track_end.append((ii+1)*ntracks_part-1)
        # Add the last part to the list
        track_start.append(track_end[-1]+1)
        track_end.append(ntracks_all)
        # Update total number of parts
        nparts = len(track_start)

        # Create the list of job tasks needed by SLURM...
        syear = track_period[0:4]
        task_filename = f'tasks_era5_3d_{syear}_rerun.txt'
        task_file = open(task_filename, "w")
        ntasks = 0
        for varname in var_dict:
            basename = var_dict[varname]
            # print(varname, '->', basename)
            for ipart in range(0, nparts): 
                # Make file name
                fname = f'{data_dir}{syear}/mcs_era5_{varname}_{track_period}_t{track_start[ipart]:05}.nc'
                # Check if file exists
                if os.path.isfile(fname) is False:
                    # Print missing file name
                    print(f'{fname} Missing')
                    # Create command to produce the file
                    cmd = f'python {code_name} ' \
                        f'{basename} {varname} {config_file} {track_period} ' \
                        f'{track_start[ipart]} {track_end[ipart]} {digits}'
                    task_file.write(f"{cmd}\n")
                    ntasks += 1
                else:
                    # Check if file size is too small
                    if (get_file_size(fname, SIZE_UNIT.MB) < min_filesize_MB):
                        # If file is at least 1 MB (not empty)
                        if (get_file_size(fname, SIZE_UNIT.MB) > 1.0):
                            # Get number of tracks in file
                            ds = xr.open_dataset(fname)
                            nt_file = ds.sizes['tracks']
                            # Get number of tracks expected
                            nt_expect = track_end[ipart] - track_start[ipart]
                            # If the two does not equal, the file is not correct
                            if nt_file != nt_expect:
                                # Print file name
                                print(f'{fname} Incomplete track')
                                # Create command to produce the file
                                cmd = f'python {code_name} ' \
                                    f'{basename} {varname} {config_file} {track_period} ' \
                                    f'{track_start[ipart]} {track_end[ipart]} {digits}'
                                task_file.write(f"{cmd}\n")
                                ntasks += 1
                        else:
                            # Print file name
                            print(f'{fname} Too small')
                            # Create command to produce the file
                            cmd = f'python {code_name} ' \
                                f'{basename} {varname} {config_file} {track_period} ' \
                                f'{track_start[ipart]} {track_end[ipart]} {digits}'
                            task_file.write(f"{cmd}\n")
                            ntasks += 1
        task_file.close()

        # Check ntasks, non-zero means there are files missing
        if ntasks > 0:
            print(task_filename)

            # Create a SLURM submission script for the above task list...
            slurm_filename = f"{slurm_basename}{syear}_rerun"
            slurm_file = open(slurm_filename, "w")
            text = f"""\
                #!/bin/bash
                #SBATCH --job-name={syear}
                #SBATCH -A m1867
                #SBATCH --time=00:45:00
                #SBATCH -p regular
                #SBATCH -N 1
                #SBATCH -c 64
                #SBATCH -C haswell
                #SBATCH --exclusive
                #SBATCH --output=log_era5_3d_{syear}_%A_%a.log
                #SBATCH --mail-type=END
                #SBATCH --mail-user=zhe.feng@pnnl.gov
                #SBATCH --array=1-{ntasks}

                date
                source activate /global/common/software/m1867/python/testflex
                cd {code_dir}

                # Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
                LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {task_filename})
                echo $LINE
                # Run the line as a command
                $LINE

                date
                """
            slurm_file.writelines(textwrap.dedent(text))
            slurm_file.close()
            print(slurm_filename)

            # Run command
            if submit_job == True:
                # cmd = f'sbatch --array=1-{ntasks}%{njobs_run} {slurm_filename}'
                cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
                print(cmd)
                # subprocess.run(f'{cmd}', shell=True)

        else:
            print(f'No missing files for {syear}!')
            # Remove task file
            os.remove(task_filename)
