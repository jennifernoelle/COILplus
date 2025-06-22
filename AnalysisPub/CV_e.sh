#!/bin/bash
#SBATCH --output=CVe.out
#SBATCH -J CVe
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4e_cv_COIL_exp_500sims.R
