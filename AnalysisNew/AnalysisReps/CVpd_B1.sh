#!/bin/bash
#SBATCH --output=CVpdB1.out
#SBATCH -J CVpdB1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH CV_pd_B1_nb_20sims.R
