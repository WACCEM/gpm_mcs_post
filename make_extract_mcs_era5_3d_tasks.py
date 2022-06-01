import numpy as np
import xarray as xr
import pandas as pd
import yaml
import textwrap

if __name__ == "__main__":

    code_dir = '/global/homes/f/feng045/program/gpm_mcs_post/'
    code_name = f'extract_mcs_era5_3d_fullhistory_bytracks.py'
    config_file = f'config_gpm_era5_3d_global.yml'

    # Start/end years to process
    start_year = 2019
    end_year = 2020

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
        task_filename = f'tasks_era5_3d_{syear}.txt'
        task_file = open(task_filename, "w")
        ntasks = 0
        for varname in var_dict:
            basename = var_dict[varname]
            # print(varname, '->', basename)
            for ipart in range(0, nparts): 
                cmd = f'python {code_name} ' \
                    f'{basename} {varname} {config_file} {track_period} ' \
                    f'{track_start[ipart]} {track_end[ipart]} {digits}'
                # print(cmd)
                task_file.write(f"{cmd}\n")
                ntasks += 1
        task_file.close()
        print(task_filename)

        # Create a SLURM submission script for the above task list...
        slurm_filename = f"slurm.submit_era5_3d_{syear}"
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
            ##SBATCH --array=1-{ntasks}

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