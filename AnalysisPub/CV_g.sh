#!/bin/bash
#SBATCH --output=CVg.out
#SBATCH -J CVg
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4g_cv_COILp_nb_exp_500sims.R
