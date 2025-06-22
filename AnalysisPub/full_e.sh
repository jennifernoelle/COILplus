#!/bin/bash
#SBATCH --output=full_e.out
#SBATCH -J full_e
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 3e_analysis_COIL_exp_500sims.R
