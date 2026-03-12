#!/bin/bash
 
#SBATCH --job-name=CONFIRMsumm
#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --partition long.cpu
#SBATCH --mail-type=END,FAIL
#SBATCH -o ./slurmoutputs/slurm-CONFIRM_summary.out

# Load the R module
module load r/4.1.0-gcc-11.1.0

# Set the R_LIBS_USER environment variable
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.1:$R_LIBS_USER
 
# Execute the job
time Rscript summary.R 

exit 0