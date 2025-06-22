#!/bin/bash
#SBATCH --output=full_star.out
#SBATCH -J full_star
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH 3star_Analysis_pd_B4_nb_400sims.R
