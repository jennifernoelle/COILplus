#!/bin/bash
#SBATCH --output=full_b.out
#SBATCH -J full_b
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3b_analysis_COIL_75_100_500sims.R
