# Estimating variable importance based on interaction model
# This code is slow, run overnight and consider reducing the number of permutations for exploratory analysis
# Plot and analyze in a separate file

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
setwd(wd_path)

data_path <- 'ProcessedDataPub/'
save_plots_path <- 'ResultsPub/'
source_path <- 'HelperScriptsPub/'

# Set models to loop over for full results
model_names <- c('COIL_0_100_500sims', # this one i named cv differently
                 'COIL_75_100_500sims', 
                 'COILp_b_75_100_500sims', 
                 'COILp_nb_75_100_500sims', 
                 'COIL_exp_500sims', 
                 'COILp_b_exp_500sims', 
                 'COILp_nb_exp_500sims')


# Save results using convention: res_date_i.rda
date <- model_names[6]
results_path <- paste0('ResultsPub/', date , '/')
save_path <- paste0('ResultsPub/', date, '/')

# ------ STEP 0: Some functions. --------- #

library(tidyverse)
library(grid)
library(gridExtra)

source(paste0(source_path, 'TraitMatching3_function.R'))
source(paste0(source_path, 'useful_functions.R'))

# Define a useful function
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Loading the data:
obs_X <- loadRData(paste0(data_path, 'Obs_X.dat')) # vert traits
obs_W <- loadRData(paste0(data_path, 'Obs_W.dat')) # plant traits
obs_A <- loadRData(paste0(data_path, 'A_obs.dat'))

# Check for missing values in the covariates:
any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)

# Getting the sample sizes:
nV <- nrow(obs_X)
nP <- nrow(obs_W)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# --------------- STEP 1: Getting the results together ----------------- #

# MCMC chains saved from full model fit:
nchains <- 4

# Putting together the predictions from the chains:

### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
}

all_pred <- all_X <- all_W <- NULL
for (ii in 1 : nchains) {
  load(paste0(results_path, 'res_', date, '_', ii, '.dat'))
  all_pred[[ii]] <- res$all_pred
  all_X[[ii]] <-res$Xs
  all_W[[ii]] <- res$Ws
}

# Number of posterior samples by chain:
Nsims <- dim(all_pred[[1]])[1]
pV <- length(all_X[[1]])
pP <- length(all_W[[1]])

# Using the linear predictor of the interaction model (not corrected for bias due to detection, number of studies):
mod_pL1s <- array(NA, dim = c(nchains * Nsims, nV, nP))
for (cc in 1 : nchains) {
  wh_entries <- Nsims * (cc - 1) + 1 : Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}

# Extracting the imputed missing values
Xs <- list()
Ws <- list()

if(any_X_miss){
for(p in 1:pV){
  nmiss.p <- length(all_X[[1]][[p]])
  Xs[[p]] <- matrix(NA, nchains*Nsims, nmiss.p)
  for (cc in 1 : nchains) {
    wh_entries <- Nsims * (cc - 1) + 1 : Nsims
    Xs[[p]][wh_entries,  ] <- (all_X[[cc]])[[p]]
  }
}
}

if(any_W_miss){
for(p in 1:pP){
  nmiss.p <- length(all_W[[1]][[p]])
  Ws[[p]] <- matrix(NA, nchains*Nsims, nmiss.p)
  for (cc in 1 : nchains) {
    wh_entries <- Nsims * (cc - 1) + 1 : Nsims
    Ws[[p]][wh_entries,  ] <- (all_W[[cc]])[[p]]
  }
}
}


# --------------- STEP 2: Performing trait matching ----------------- #

# Dealing with extreme values for which logit(x) is infinite.
mod_pL1s[mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}

# Line by line running of the function
B <- 500 # permutations
# obs_X = Obs_X
# obs_W = Obs_W
# obs_only <- TRUE
# 

t1 <- Sys.time()
trait_match <- TraitMatching3(B = B, mod_pL1s = mod_pL1s,
                              Xs = NULL, Ws = NULL,  # Imputed values not used.
                              obs_X = obs_X, obs_W = obs_W, obs_only = TRUE)

print(Sys.time() - t1)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W
corr_obs_X <- trait_match$corr_obs_X
corr_obs_W <- trait_match$corr_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X_', date, '.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W_', date, '.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W_', date, '.dat'))
save(corr_obs_X, file = paste0(save_path, 'corr_obs_X_', date, '.dat'))
save(corr_obs_W, file = paste0(save_path, 'corr_obs_W_', date, '.dat'))
