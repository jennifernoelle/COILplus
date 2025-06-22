#!/bin/bash
#SBATCH --output=CVp75A3.out
#SBATCH -J CVp75A3
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_p75_A3_old_200sims.R
