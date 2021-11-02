"""
Finds missing MCS ERA5 files and provide the job array command to resubmit the failed jobs.
"""
__author__ = "Zhe.Feng@pnnl.gov"
__date__ = "31-Aug-2021"

import sys, os
import yaml
import pandas as pd

if __name__ == "__main__":
    config_file = sys.argv[1]
    start_year = sys.argv[2]
    end_year = sys.argv[3]

    # Get inputs from configuration file
    stream = open(config_file, 'r')
    config = yaml.full_load(stream)
    data_dir = config['output_dir']
    basename = 'mcs_tracks_era5_'
    # Get region name from config list, e.g., 'config_gpm_era5_2d_spac.yml'
    region = config_file.split('.')[0].split('_')[-1]
    # slurm_name = f'slurm_era5_{region}_2d'
    slurm_name = f'slurm_era5_{region}_3d'

    # Variable list, must be in the same order as the job array task list
    # var_list = ['CSF', 'TCWV', 'VAR_2T', 'VAR_2D', 'SP', 'VAR_10U', 'VAR_10V',
    #             'MSL', 'CAPE', 'SSTK', 'IE', 'ISHF', 'SKT', 'VIWVD']
    # var_list = ["CSF", "TCWV", "CAPE", "SSTK", "SKT", "VAR_2T", "VAR_2D", 
    #             "SP", "MSL", "VAR_10U", "VAR_10V", "IE", "ISHF", "VIWVD"]
    # var_list = ["W", "R", "Q", "D", "U", "V"]
    var_list = ["T", "Q", "Z", "U", "V", "W", "R"]

    # Generate time marks within the start/end datetime
    syear_list = pd.date_range(start=start_year, end=end_year, freq='YS', closed='left').strftime('%Y%m%d')
    eyear_list = pd.date_range(start=start_year, end=end_year, freq='Y', closed='left').strftime('%Y%m%d')

    job_num = []
    line_num = 1
    # Loop over each variable in the list
    for var in var_list:
        # Loop over each year
        for ii in range(len(syear_list)):
            # Check if this file exists
            fname = f"{data_dir}{basename}{var}_{syear_list[ii]}_{eyear_list[ii]}.nc"
            # Record the job number (line number)
            if os.path.isfile(fname) == False:
                job_num.append(line_num)
                print(f'{line_num}: {os.path.basename(fname)}')
            # Increment line number
            line_num = line_num + 1

    # Print all jobs with missing files
    separator = ","
    if len(job_num) > 0:
        print(f'Resubmit missing jobs:')
        print('sbatch --array='+separator.join(map(str, job_num))+f' {slurm_name}')
    else:
        print(f'All files present!')

    # import pdb; pdb.set_trace()