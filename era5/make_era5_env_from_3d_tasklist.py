"""
Make slurm scripts for calculating ERA5 2D environmental variables from 3D profiles for each year.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "4-Aug-2022"

import numpy as np
import os
import glob
import textwrap
import subprocess

if __name__ == "__main__":

    start_year = 2002
    end_year = 2002

    # Submit slurm job
    submit_job = True
  
    # ERA5 input data directory
    in_dir3d = '/global/cscratch1/sd/feng045/waccem/mcs_global/global/era5_3d/'
    in_basename = 'mcs_era5_T_'
    code_dir = '/global/homes/f/feng045/program/gpm_mcs_post/'
    code_name = 'compute_mcs_era5_env_from_3d.py'

    years = np.arange(start_year, end_year+1, 1)
    # Loop over years
    for iyear in years:
        # Find files
        idir = f'{in_dir3d}{iyear}/'
        files = sorted(glob.glob(f'{idir}{in_basename}*nc'))

        ntasks = 0
        if len(files) > 0:
            # Create the list of job tasks needed by SLURM...
            task_filename = f"slurm_tasks_2d_env_{iyear}.txt"
            task_file = open(task_filename, "w")

            # Loop over files
            for ifile in files:
                # Get filename without the path
                fname = os.path.basename(ifile)
                # Get the date/time string
                # Filename format: mcs_era5_T_20190101.0000_20200101.0000_t10000.nc
                nlead_char = len(in_basename)
                date_str = fname[nlead_char:-3]
                
                # Generate the command and add it to the task list...
                cmd = f'python {code_name} {date_str}'
                task_file.write(f"{cmd}\n")
                ntasks += 1
            task_file.close()
            print(f'{task_filename}')
                
            # Create a SLURM submission script for the above task list...
            # slurm_file = open("slurm.submit_tasks", "w")
            slurm_filename = f"slurm.submit_tasks_{iyear}"
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
                
                date
                """
            slurm_file.writelines(textwrap.dedent(text))
            slurm_file.close()
            print(f'{slurm_filename}')

        # Run command
        if submit_job == True:
            cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
            print(cmd)
            subprocess.run(f'{cmd}', shell=True)
        # import pdb; pdb.set_trace()
