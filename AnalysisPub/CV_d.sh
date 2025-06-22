#!/bin/bash
#SBATCH --output=CVd.out
#SBATCH -J CVd
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=15GB

module load R/4.1.1-rhel8 
R CMD BATCH 4d_cv_COILp_nb_75_100_500sims.R
