#!/bin/bash
#SBATCH --output=full_d.out
#SBATCH -J full_d
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 3d_analysis_COILp_nb_75_100_500sims.R
