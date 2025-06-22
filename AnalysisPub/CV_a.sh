#!/bin/bash
#SBATCH --output=CVa.out
#SBATCH -J CVa
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 4a_cv_COIL_0_100_500sims.R
