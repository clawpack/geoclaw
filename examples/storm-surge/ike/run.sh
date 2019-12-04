#!/bin/sh
#
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=apam                # The account name for the job.
#SBATCH --job-name= Run Storm Jobs    # The job name.
#SBATCH -c 10                         # The number of cpu cores to use.
#SBATCH --time=1:30:00                # The time the job will take to run.
#SBATCH --mem-per-cpu=15gb            # The memory the job will use per cpu core.
#SBATCH --nodes=1                     # The memory the job will use per cpu core.


make new
make all
# End of script
