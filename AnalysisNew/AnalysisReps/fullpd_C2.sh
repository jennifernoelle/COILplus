#!/bin/bash
#SBATCH --output=fullpd_C2.out
#SBATCH -J pd_C2
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH Analysis_pd_C2_b_100sims.R
