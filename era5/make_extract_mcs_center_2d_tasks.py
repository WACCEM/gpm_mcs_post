"""
Make task list and slurm scripts for extracting MCS center 2D variables.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "19-May-2023"
import numpy as np
import xarray as xr
import pandas as pd
import yaml
import textwrap
import subprocess

if __name__ == "__main__":
    code_dir = '/global/homes/f/feng045/program/gpm_mcs_post/era5/'
    code_name = f'extract_mcs_masks_2d_fullhistory_bytracks.py'
    config_file = f'config_gpm_mask_2d_global.yml'
    slurm_basename = f'slurm.submit_extract_mcs_mask_2d_'

    # Number of jobs allowed to run simultaneously in job array
    # njobs_run = 60
    # Submit slurm job
    submit_job = True
    
    # Start/end years to process
    start_year = 2000
    end_year = 2000

    # Number of tracks to process per part
    ntracks_part = 1000
    # Set the number of digits for 0 padding
    # This should be set to the digit for the maximum number of tracks
    # e.g., ntracks = 32138, digits = 5
    digits = 5

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
        statsfile_basename = config['statsfile_basename']
        ntimes_max = config['ntimes_max']

        # Read robust MCS statistics
        mcs_file = f"{stats_dir}{statsfile_basename}{track_period}.nc"
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
        task_filename = f'tasks_extract_mcs_masks_2d_{syear}.txt'
        task_file = open(task_filename, "w")
        ntasks = 0
        for ipart in range(0, nparts): 
            cmd = f'python {code_name} ' \
                f'{config_file} {track_period} ' \
                f'{track_start[ipart]} {track_end[ipart]} {digits}'
            # print(cmd)
            task_file.write(f"{cmd}\n")
            ntasks += 1
        task_file.close()
        print(task_filename)

        # Create a SLURM submission script for the above task list...
        slurm_filename = f"{slurm_basename}{syear}.sh"
        slurm_file = open(slurm_filename, "w")
        text = f"""\
            #!/bin/bash
            #SBATCH -A m1867
            #SBATCH -J {syear}
            #SBATCH -t 00:10:00
            #SBATCH -q regular
            #SBATCH -C cpu
            #SBATCH --nodes=1
            #SBATCH --ntasks-per-node=128
            #SBATCH --exclusive
            #SBATCH --output=log_mcs_mask_2d_{syear}_%A_%a.log
            #SBATCH --mail-type=END
            #SBATCH --mail-user=zhe.feng@pnnl.gov
            #SBATCH --array=1-{ntasks}

            date
            source activate /global/common/software/m1867/python/py310

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
            cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
            print(cmd)
            subprocess.run(f'{cmd}', shell=True)

        # import pdb; pdb.set_trace()