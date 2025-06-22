#!/bin/bash
#SBATCH --output=CVpdC4retry.out
#SBATCH -J CVpdC4retry
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=12GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_C4_b_400sims.R
