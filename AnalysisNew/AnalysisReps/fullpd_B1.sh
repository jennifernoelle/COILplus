#!/bin/bash
#SBATCH --output=fullpd_B1.out
#SBATCH -J pd_B1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH Analysis_pd_B1_nb_20sims.R
