#!/bin/bash
#SBATCH --output=CVp0A1.out
#SBATCH -J CVp0A1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_p0_A1_old_20sims.R
