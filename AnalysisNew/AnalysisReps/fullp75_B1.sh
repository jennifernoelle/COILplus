#!/bin/bash
#SBATCH --output=fullp75_B1.out
#SBATCH -J p75_B1
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH Analysis_p75_B1_nb_20sims.R
