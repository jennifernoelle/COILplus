#!/bin/bash
#SBATCH --output=CVpdB4.out
#SBATCH -J CVpdB4
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_B4_nb_400sims.R
