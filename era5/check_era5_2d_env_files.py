"""
Check the ERA5 2D environment files, make a new task list and slurm script for 
the missing files, or files too small, or files with incomplete tracks.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "4-Aug-2022"
import numpy as np
import os, sys
import glob
import enum
import textwrap
import subprocess
import xarray as xr

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
    start_year = 2011
    end_year = 2021

    # Submit slurm job
    submit_job = False

    code_dir = '/global/homes/f/feng045/program/gpm_mcs_post/'
    code_name = f'compute_mcs_era5_env_from_3d.py'
    slurm_basename = f'slurm.submit_tasks_'
    in_dir3d = '/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_3d/'
    in_basename3d = 'mcs_era5_T_'
    in_basename2d = 'mcs_era5_2D_ENVS_'

    # Minimum expected 2D ENV file size [MB]
    min_filesize_MB = 100

    years = np.arange(start_year, end_year, 1)
    # Loop over years
    for iyear in years:
        # Find 3D files
        idir = f'{in_dir3d}{iyear}/'
        files = sorted(glob.glob(f'{idir}{in_basename3d}*nc'))

        # Create the list of job tasks needed by SLURM...
        task_filename = f'slurm_tasks_2d_env_{iyear}_rerun.txt'
        task_file = open(task_filename, "w")
        ntasks = 0
        # Loop over files
        for ifile in files:
            # Get filename without the path
            fname3d = os.path.basename(ifile)
            # Get the date/time string
            # Filename format: mcs_era5_T_20190101.0000_20200101.0000_t10000.nc
            nlead_char = len(in_basename3d)
            date_str = fname3d[nlead_char:-3]

            # Make 2D ENV file name
            fname = f'{idir}{in_basename2d}{date_str}.nc'
            # Check if file exists
            if (os.path.isfile(fname) is False):
                # Print missing file name
                print(f'{fname} Missing')
                # Create command to produce the file
                cmd = f'python {code_name} {date_str}'
                task_file.write(f"{cmd}\n")
                ntasks += 1
            else:
                # Check if file size is too small
                if (get_file_size(fname, SIZE_UNIT.MB) < min_filesize_MB):    
                    # If file is at least 1 MB (not empty)
                    if (get_file_size(fname, SIZE_UNIT.MB) > 1.0):
                        # Get number of tracks in 2D file
                        ds = xr.open_dataset(fname)
                        nt_file = ds.sizes['tracks']
                        # Get number of tracks in 3D file
                        ds3d = xr.open_dataset(ifile)
                        nt_expect = ds3d.sizes['tracks']
                        # If the two does not equal, the file is not correct
                        if nt_file != nt_expect:
                            # Print file name
                            print(f'{fname} Incomplete track')
                            # Create command to produce the file
                            cmd = f'python {code_name} {date_str}'
                            task_file.write(f"{cmd}\n")
                            ntasks += 1
                    else:            
                        # Print file name
                        print(f'{fname} Too small')
                        # Create command to produce the file
                        cmd = f'python {code_name} {date_str}'
                        task_file.write(f"{cmd}\n")
                        ntasks += 1
        task_file.close()

        # Check ntasks, non-zero means there are files missing
        if ntasks > 0:
            print(task_filename)

            # Create a SLURM submission script for the above task list...
            slurm_filename = f"{slurm_basename}{iyear}_rerun"
            slurm_file = open(slurm_filename, "w")
            text = f"""\
                #!/bin/bash
                #SBATCH --job-name=2D_{iyear}
                #SBATCH -A m1867
                #SBATCH -p regular
                #SBATCH --time=00:10:00
                #SBATCH -p regular
                #SBATCH -N 1
                #SBATCH -c 64
                #SBATCH -C haswell
                #SBATCH --exclusive
                #SBATCH --output=log_era5_{iyear}_%A_%a.log
                #SBATCH --mail-type=END
                #SBATCH --mail-user=zhe.feng@pnnl.gov
                #SBATCH --array=1-{ntasks}
                
                date
                source activate /global/common/software/m1867/python/py38
                cd {code_dir}

                LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {task_filename}) 
                echo $LINE
                $LINE
                """
            slurm_file.writelines(textwrap.dedent(text))
            slurm_file.close()
            print(f'{slurm_filename}')

            # Run command
            if submit_job == True:
                cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
                print(cmd)
                subprocess.run(f'{cmd}', shell=True)
        
        else:
            print(f'No missing files for {iyear}!')
            # Remove task file
            os.remove(task_filename)

    # import pdb; pdb.set_trace()