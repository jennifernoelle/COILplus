#!/bin/bash
#SBATCH --output=CVpdC2.out
#SBATCH -J CVpdC2
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_C2_b_100sims.R
