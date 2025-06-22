#!/bin/bash
#SBATCH --output=fullp0_A4.out
#SBATCH -J p0_A4
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

module load R/4.1.1-rhel8 
R CMD BATCH Analysis_p0_A4_old_400sims.R
