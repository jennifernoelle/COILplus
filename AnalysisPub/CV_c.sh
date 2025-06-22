#!/bin/bash
#SBATCH --output=CVc.out
#SBATCH -J CVc
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4c_cv_COILp_b_75_100_500sims.R
