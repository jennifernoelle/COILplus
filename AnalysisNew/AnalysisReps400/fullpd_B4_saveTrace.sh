#!/bin/bash
#SBATCH --output=fullpd_B4.out
#SBATCH -J pd_B4_test
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 4
#SBATCH --mem-per-cpu=8GB

export OMP_NUM_THREADS=1 
export MKL_NUM_THREADS=1 

module load R/4.1.1-rhel8 

R CMD BATCH Analysis_pd_B4_nb_400sims_saveTrace.R
