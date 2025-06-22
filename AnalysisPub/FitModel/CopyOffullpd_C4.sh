#!/bin/bash
#SBATCH --output=fullpd_C4_alt.out
#SBATCH -J pd_C4_alt
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH Analysis_pd_C4_b_400sims_alt.R
