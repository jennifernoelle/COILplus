#!/bin/bash
#SBATCH --output=traits_save100.out
#SBATCH -J traits100
#SBATCH --partition=scavenger --account=dunsonlab
#SBATCH -c 10
#SBATCH --mem-per-cpu=10GB

module load R/4.1.1-rhel8 
R CMD BATCH 6a_trait_matching_save_100lighter.R 
