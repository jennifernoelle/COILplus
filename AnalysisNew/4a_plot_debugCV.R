# File to plot results for cv and full run

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
# wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
# setwd(wd_path)

# Save results using convention: res_date_i.rda
sampleP <- TRUE

# Where the processed data are saved:
data_path <- 'ProcessedDataNew/'
# Where you want to save MCMC results:
results_path <- paste0('ResultsNew/', date, '/')
# Where the functions are available:
source_path <- 'HelperScriptsJKNew/'


# ------ STEP 0: Some functions. --------- #


# We want the github version of superheat
if(!("superheat" %in% rownames(installed.packages()))){
  install.packages("devtools")
  remotes::install_github("rlbarter/superheat")
}

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(superheat)
library(abind)
library(gridExtra)
library(grid)
library(dplyr)
library(caret)
library(pROC)
library(bayesplot)
library(tidyverse)

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'F_obs.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'OP_full.dat')) # site level obs plants
load(paste0(data_path, 'OV_full.dat')) # site level obs verts

v.taxa <- read.csv(paste0(data_path, 'v_taxa.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa.csv'))

sum(v.taxa$Animal_Species_Corrected == rownames(obs_A))
sum(p.taxa$Plant_Species_Corrected == colnames(obs_A))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_W <- Obs_W
obs_X <- Obs_X
obs_A <- A_obs
obs_F <- F_obs

## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

## Improved guess: Expert 1 modified
# Same study: 1, same site: 0.75 -> 0.85
# Same country and habitat: 0.5 -> 0.65, same region and habitat: 0.45 -> 0.35, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05 

## Improved default guess: 0.75/1
O_P <- ifelse(O_P == 1, 1, 0.75)
O_V <- ifelse(O_V == 1, 1, 0.75)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]
s.names <- dimnames(obs_A)[[3]]


# -------------- STEP 1: Specifications. ------------ #

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 2

# Number of cv samples per repetition
n.cv <- 100


# --------------- STEP 4: Cross validation results ----------------- #
date <- 'pd_B1_nb_20sims_debug'
results_path <- paste0('ResultsNew/', date, '/')


# Getting the results together (held out indicies and predictions)
all_indices_nb <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds_nb <- array(NA, dim = c(repetitions, nV, nP))
for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
  pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
  all_indices_nb[rr, , ] <- cv_indices
  our_preds_nb[rr, , ] <- pred
}


# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred_nb <- array(NA, dim = c(repetitions, n.cv))
for (rr in 1 : repetitions) {
  for (ii in 1 : n.cv) {
    pred_nb[rr, ii] <- our_preds_nb[rr, all_indices_nb[rr, ii, 1], all_indices_nb[rr, ii, 2]]
  }
}

# Average and median in the held out data.
pred_mean <- apply(pred_nb, 1, mean)
pred_median <- apply(pred_nb, 1, median)


date <- 'pd_C1_b_20sims_debug'
results_path <- paste0('ResultsNew/', date, '/')

# Getting the results together (held out indicies and predictions)
all_indices_b <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds_b <- array(NA, dim = c(repetitions, nV, nP))
for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
  pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
  all_indices_b[rr, , ] <- cv_indices
  our_preds_b[rr, , ] <- pred
}


# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred_b <- array(NA, dim = c(repetitions, n.cv))
for (rr in 1 : repetitions) {
  for (ii in 1 : n.cv) {
    pred_b[rr, ii] <- our_preds_b[rr, all_indices_b[rr, ii, 1], all_indices_b[rr, ii, 2]]
  }
}

### Compare
# Average and median in the held out data.
apply(pred_b, 1, mean)
apply(pred_nb, 1, mean)

# The arrays with all post samples all equal between blocked and unblocked
all.equal(our_preds_b, our_preds_nb)

# Are the chains equal?? NO, thank goodness, just coincidence we got the same means?
c(all_indices_b[1,,])
c(all_indices_b[2,,])

c(pred_b[1,])
c(pred_b[2,])

# So, the CV blocked and unblocked are the same when i run in parallel
# Now go see if the results are the same if I just run the simple cv code - yes
# Are the results the same if I run the full sample, no cv code - no
# Found an issue: block_sampleOccP = block_sampleOccP needed in an argument
# The function now works as expected in a loop, but not in mclapply





