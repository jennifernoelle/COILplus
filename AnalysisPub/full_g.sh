#!/bin/bash
#SBATCH --output=full_g_p75coilp_nb.out
#SBATCH -J full_g_cpexp_nb
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3g_analysis_COILp_nb_exp_500sims.R
