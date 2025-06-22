#!/bin/bash
#SBATCH --output=CVpdC1.out
#SBATCH -J CVpdC1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_C1_b_20sims.R
