#!/bin/bash
#SBATCH --output=CVpdA1.out
#SBATCH -J CVpdA1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_A1_old_20sims.R
