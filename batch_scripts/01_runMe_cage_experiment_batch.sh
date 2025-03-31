#!/bin/bash
#SBATCH -N 1               	                                # number of nodes (no MPI, so we only use a single node)
#SBATCH -n 20            	                                # number of cores
#SBATCH --time=24:00:00    	                                # walltime allocation, which has the format (D-HH:MM:SS), here set to 1 hour
#SBATCH --mem=5GB         	                                # memory required per node (here set to 4 GB)
#SBATCH --output=slurm_outputs/01_runMe_cage/slurm-%A_%a.out

# Notification configuration
#SBATCH --array=1-36
#SBATCH --mail-type=END					    	# Send a notification email when the job is done (=END)

#loading modules
module load R/4.3.1

# Execute the program
MIN_VALUES=FALSE
MAX_VALUES=TRUE

Rscript code/HPC_versions/01_runME_cage_experiment_hpc.R
