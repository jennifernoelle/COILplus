#!/bin/bash
#SBATCH --output=CVp75B3.out
#SBATCH -J CVp75B3
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_p75_B3_nb_200sims.R
