#!/bin/bash
#SBATCH --output=CVpdB3.out
#SBATCH -J CVpdB3
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_B3_nb_200sims.R
