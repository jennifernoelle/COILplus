#!/bin/bash
#SBATCH --output=full_f_expcoilp_nb.out
#SBATCH -J f_f_cexp
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3f_analysis_COILp_b_exp_500sims.R
