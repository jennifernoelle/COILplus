#!/bin/bash
#SBATCH --output=CVb.out
#SBATCH -J CVb
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4b_cv_COIL_75_100_500sims.R
