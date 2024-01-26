import os
import numpy as np
import pandas as pd
import textwrap
import subprocess

if __name__ == "__main__":

    submit_job = True
    
    # Start/end dates
    # start_date = '2000-06-01'
    start_date = '2001-01-01'
    end_date = '2020-12-31'

    task_filename = f'tasks_combine_pixelfiles_daily.txt'
    slurm_filename = f'slurm.combine_pixelfiles_daily.sh'

    # Input/output directory
    in_dir = '/pscratch/sd/f/feng045/waccem/mcs_global/mcstracking/'
    out_dir = '/pscratch/sd/f/feng045/waccem/mcs_global/mcstracking_daily/'
    basename = 'mcstrack_'

    os.makedirs(out_dir, exist_ok=True)

    # Make a list of dates
    date_range = pd.date_range(start=start_date, end=end_date, freq='D')
    # Get unique years & months
    unique_years = np.unique(date_range.year)
    # unique_months = np.unique(date_range.month)

    # Create a job task needed by SLURM...
    task_file = open(task_filename, "w")
    ntasks = 0

    # Loop over each date
    for idate in date_range:
        iday = idate.strftime('%Y%m%d')
        iyear_str = idate.strftime('%Y')
        # Find sub-directory starting with the year
        idir = [d for d in os.listdir(in_dir) if os.path.isdir(os.path.join(in_dir, d)) and d.startswith(iyear_str)]
        dir_year = f'{in_dir}{idir[0]}/'
        # Input/output file names
        infiles = f'{dir_year}{basename}{iday}*nc'
        outfile = f'{out_dir}{basename}{iday}.nc'
        # Task command
        cmd = f'ncrcat -O -h {infiles} {outfile}'
        # import pdb; pdb.set_trace()
        # Add to task file
        task_file.write(f"{cmd}\n")
        ntasks += 1
    
    task_file.close()
    print(task_filename)

    # Create a SLURM submission script for the above task list...
    slurm_file = open(slurm_filename, "w")
    text = f"""\
        #!/bin/bash
        #SBATCH -A m1867
        #SBATCH -J catfiles
        #SBATCH -t 00:30:00
        #SBATCH -q debug
        #SBATCH -C cpu
        #SBATCH --nodes=8
        #SBATCH --ntasks-per-node=256
        #SBATCH --exclusive
        #SBATCH --output=log_combine_pixelfiles_daily.log
        #SBATCH --mail-type=END
        #SBATCH --mail-user=zhe.feng@pnnl.gov
        date
        module load taskfarmer
        source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
        export THREADS=256
        cd {os.getcwd()}
        runcommands.sh {task_filename}
        date
    """
    slurm_file.writelines(textwrap.dedent(text))
    slurm_file.close()
    print(slurm_filename)

    # Run command
    if submit_job == True:
        # cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
        cmd = f'sbatch {slurm_filename}'
        print(cmd)
        subprocess.run(f'{cmd}', shell=True)

    # import pdb; pdb.set_trace()