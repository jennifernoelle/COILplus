#!/bin/bash
#SBATCH --output=full_a.out
#SBATCH -J full_a
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3a_analysis_COIL_0_100_500sims.R
