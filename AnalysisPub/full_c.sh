#!/bin/bash
#SBATCH --output=full_c.out
#SBATCH -J full_c
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3c_analysis_COILp_b_75_100_500sims.R
