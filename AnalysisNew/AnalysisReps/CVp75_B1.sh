#!/bin/bash
#SBATCH --output=CVp75B1.out
#SBATCH -J CVp75B1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_p75_B1_b_20sims.R
