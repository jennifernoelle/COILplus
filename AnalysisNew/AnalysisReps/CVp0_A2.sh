#!/bin/bash
#SBATCH --output=CVp0A2.out
#SBATCH -J CVp0A2
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_p0_A2_old_100sims.R
