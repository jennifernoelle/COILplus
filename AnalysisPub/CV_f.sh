#!/bin/bash
#SBATCH --output=CVf.out
#SBATCH -J CVf
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4f_cv_COILp_b_exp_500sims.R
