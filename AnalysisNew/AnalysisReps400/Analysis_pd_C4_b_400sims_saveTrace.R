# This version doesn't save full probabilities but rather running mean
# Doesn't save marginal probabilities or latent factors at all

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'pd_C4_b_400sims_saveTrace'

# Where the processed data are saved:
data_path <- 'ProcessedDataNew/'
# Where you want to save MCMC results:
save_path_base <- 'ResultsNew/'
# Where the functions are available:
source_path <- 'HelperScriptsJKNew/'

# Create the results folder
ifelse(!dir.exists(file.path(save_path_base, date)), dir.create(file.path(save_path_base, date)), FALSE)
save_path <- paste0('ResultsNew/', date, '/')

# ------ STEP 0: Some functions. --------- #

source(paste0(source_path, 'UpdExtraVar_function.R'))
source(paste0(source_path, 'UpdTraitCoef_function.R'))
source(paste0(source_path, 'UpdLatFac_function.R'))
source(paste0(source_path, 'UpdProbObs_function.R'))
source(paste0(source_path, 'UpdOccur_function.R'))
source(paste0(source_path, 'UpdOccurP_function.R'))
source(paste0(source_path, 'UpdOccurP_function_blocked.R'))
source(paste0(source_path, 'UpdRho_function.R'))
source(paste0(source_path, 'OmegaFromV_function.R'))
source(paste0(source_path, 'useful_functions.R'))
source(paste0(source_path, 'CorrMat_function.R'))
source(paste0(source_path, 'PredictInteractions_function.R'))
source(paste0(source_path, 'GetPredLatFac_function.R'))
source(paste0(source_path, 'GetPredWeights_function.R'))
source(paste0(source_path, 'Utils_OccurP.R'))
source(paste0(source_path, 'MCMC_function_working_saveRhoSomePis.R'))

library(parallel)
library(foreach)
library(abind)
library(magrittr)
library(truncnorm)

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'F_obs.dat'))
load(paste0(data_path, 'Obs_X.dat')) # mammal traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'OP_full.dat')) # site level obs plants
load(paste0(data_path, 'OV_full.dat')) # site level obs verts

## Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_W <- Obs_W
obs_X <- Obs_X
obs_A <- A_obs
obs_F <- F_obs

## The following assignments were used in creating obs_OP and obs_OV
# Species occurs in the same study -> 1
# Species occurs at the same site, but a different study -> 0.75
# Species occurs in the same country, but a different site -> 0.5
# Species occurs in the same zone, but different ... -> 0.25

## Old Expert guess: 
# Same study: 1, same site: 0.75 -> 0.85
# Same country and habitat: 0.5 -> 0.65, same region and habitat: 0.45 -> 0.35, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05 

## New Expert guess: 
# Same study and site: 1, same site diff study: 0.75 -> 0.85
# Same country different site: 0.5 -> 0.75, 
# Same zone: 0.25 -> 0.5

## Improved default guess: 0.75/1
# O_P <- ifelse(O_P == 1, 1, 0)
# O_V <- ifelse(O_V ==1, 1, 0)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]


# -------------- STEP 1: Specifications. ------------ #

Nsims <- 400 #10000 #500 #1000 #5000 # original 10000, reasonable 2500
thin <-  10 #5 # original 40
burn <-  floor(Nsims*thin*5) #5 #40000 #2000 #5000 # 2000 # 22000 # original 40000, reasonable 2500
use_H <- 10 # original 10
theta_inf <- 0.01


# Prior distributions:
stick_alpha <- 5
prior_theta <- c(1, 1)
prior_tau <- c(5, 5)
prior_rho <- c(5, 5)  # I do not update this for now.
prior_mu0 <- 0
prior_sigmasq0 <- 10
prior_sigmasq <- c(1, 1)

# sampling <- NULL
sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                 gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                 delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                 z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                 rV = TRUE, miss_X = TRUE, miss_W = TRUE, O_B = TRUE,
                 O_P = TRUE, p_OB = TRUE, p_OP = TRUE)

start_values <- NULL
block_sampleOccP <- TRUE
bias_cor <- TRUE # Performing bias correction.

# Setup some random pis to save all posterior samples of for traceplots
# Matrix that chooses 9 non-recorded interactions to save for trace plots
set_out <- matrix(0, nrow = nV, ncol = nP)   
set_out[sample(which(comb_A == 0), 9)] <- 1

# Getting the indices
set.seed(1234)
pi_trace_indices <- matrix(NA, nrow = 9, ncol = 2)
wh <- which(set_out == 1)
for (ii in 1 : 9) {
  row_ii <- wh[ii] %% nV
  row_ii <- ifelse(row_ii == 0, nV, row_ii)
  col_ii <- ceiling(wh[ii] / nV)
  pi_trace_indices[ii, ] <- c(row_ii, col_ii)
}

# # Line by line assignments
# use_shrinkage <- TRUE
# cut_feed <- FALSE
# p_occur_B <- O_V
# p_occur_P <- O_P
# focus <- obs_F
# 
# 
# # --------------- TEST WITHOUT PARALLELIZATION -------------------#
# 
# t1 <- Sys.time()
# mcmc <- MCMC(obs_A, focus = obs_F, p_occur_B = O_V, p_occur_P = O_P, obs_X, obs_W, Cu, Cv,
#                           Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
#                           bias_cor = TRUE, theta_inf = 0.01,
#                           mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
#                           mh_pprior_sd = 0.1, mh_p_step = 0.1,
#                           stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
#                           prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
#                           prior_sigmasq = c(1, 1), start_values = NULL,
#                           sampling = sampling,
#                           cut_feed = FALSE,
#                           block_sampleOccP = block_sampleOccP)
# Sys.time() - t1
# 
# 
# # Binding different predictions of interest: Posterior samples of the
# # interaction indicators, the linear predictor of the interaction model,
# all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
# 
# # Phylogenetic correlation parameter for bird and plant correlation matrices.
# correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
# 
# # Running mean of detection probabilities
# p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
# 
# # Running mean of latent factors
# factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
# 
# # Occurrence probabilities for plants
# occ_plants <- list(p_OPs = mcmc$p_OP_mean, p_accept = mcmc$p_OP_accepted)
# 
# # Combining the results we are interested in to a list and saving:
# res <- list(all_pred = all_pred, correlations = correlations,
#             p_detect = p_detect, factors = factors, occ_plants = occ_plants)
# 
# ## Compare prior to posterior mean occurrence probs
# post.pi <- occ_plants$p_OPs
# plot(x = c(p_occur_P), y = c(post.pi) , main = "Posterior vs Prior",
#      xlab = "Prior Pi", ylab = "Post Mean Pi")
# 
# 

# --------------- STEP 2: SETUP PARALLEL ----------------- #

# # Set up function to execute in parallel
n.chains <- 4
mcmc.parallel <- function(cc, obs_A, focus = obs_F, p_occur_B = O_V, p_occur_P = O_P, obs_X, obs_W, Cu, Cv,
                          Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
                          bias_cor = TRUE, theta_inf = 0.01,
                          mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                          mh_pprior_sd = 0.1, mh_p_step = 0.1,
                          stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                          prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                          prior_sigmasq = c(1, 1), start_values = NULL,
                          sampling = sampling,
                          cut_feed = FALSE,
                          block_sampleOccP = block_sampleOccP, 
                          pi_trace_indices = pi_trace_indices){
  

set.seed(cc)

  ## Create a unique filename for each interation of the parallel loop
  each_filename <- paste0('res_', date, '_', as.character(cc), '.dat')
  each_filepath <- file.path(save_path, each_filename)
  
  mcmc <- MCMC(obs_A, focus, p_occur_B, p_occur_P, obs_X, obs_W, Cu, Cv,
           Nsims, burn, thin, use_H = 10, use_shrinkage = TRUE,
           bias_cor = TRUE, theta_inf = 0.01,
           mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
           mh_pprior_sd = 0.1, mh_p_step = 0.1,
           stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
           prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
           prior_sigmasq = c(1, 1), start_values = NULL,
           sampling = sampling, cut_feed = FALSE, 
           block_sampleOccP = block_sampleOccP, 
           pi_trace_indices = pi_trace_indices) 
  
  
  # Binding different predictions of interest: Posterior samples of the
  # interaction indicators, the linear predictor of the interaction model,
  #all_pred <- abind::abind(pred_L = mcmc$Ls, probL = mcmc$mod_pL1s, along = 4)
  
  # Phylogenetic correlation parameter for bird and plant correlation matrices.
  correlations <- cbind(U = mcmc$rU, V = mcmc$rV)
  correlations_accept <- list(ru_accepted = mcmc$ru_accepted, rv_accepted = mcmc$rv_accepted)
  
  random_pis <- list(save_pis = mcmc$save_pis, save_pi_indices = mcmc$pi_trace_indices)
  
  # # Running mean of detection probabilities
  # p_detect <- list(pis = mcmc$pi_mean, pjs = mcmc$pj_mean)
  # 
  # # Running mean of latent factors
  # factors <- list(U = mcmc$U_mean, V = mcmc$V_mean)
  # 
  # # Imputed values of missing covariates
  # Xs <- mcmc$Xs
  # Ws <- mcmc$Ws
  # 
  # # Occurrence indicators, probabilities, acceptance rates (last two not relevant if we don't sample P)
  # occ_plants <- list(OP_mean = mcmc$OP_mean, p_OPs = mcmc$p_OP_mean, p_accept = mcmc$p_OP_accepted) 
  # occ_verts <- list(OB_mean = mcmc$OB_mean, p_BPs = mcmc$p_OB_mean, p_accept = mcmc$p_OB_accepted) 
  # 
  # Combining the results we are interested in to a list and saving:
  res <- list(#all_pred = all_pred, 
              correlations = correlations, accept = correlations_accept, 
              random_pis = random_pis
              # p_detect = p_detect, factors = factors, occ_plants = occ_plants, 
              # occ_verts = occ_verts,  Xs = Xs, Ws = Ws
              )
  save(res, file = each_filepath)
  
  rm(res)
}

#---------------------- STEP 3: RUN THE SAMPLER -------------------------------------------#
t1 <- Sys.time()
mclapply(1:n.chains, function(i) mcmc.parallel(cc=i,
                                                obs_A = obs_A, focus = obs_F, p_occur_B = O_V, p_occur_P = O_P,
                                                obs_X = obs_X, obs_W = obs_W, Cu = Cu, Cv = Cv,
                                                Nsims = Nsims, burn = burn, thin = thin,
                                                use_H = use_H, bias_cor = bias_cor,use_shrinkage = TRUE,
                                                theta_inf = theta_inf, mh_n_pis = mh_n_pis,
                                                mh_n_pjs = mh_n_pjs, mh_n_rho = mh_n_rho,
                                                mh_pprior_sd = 0.1, mh_p_step = 0.1,
                                                stick_alpha = stick_alpha, prior_theta = prior_theta,
                                                prior_tau = prior_tau, prior_rho = prior_rho,
                                                prior_mu0 = prior_mu0, prior_sigmasq0 = prior_sigmasq0,
                                                prior_sigmasq = prior_sigmasq, start_values = start_values,
                                                sampling = sampling, cut_feed = FALSE, 
                                                block_sampleOccP = block_sampleOccP, 
                                               pi_trace_indices = pi_trace_indices),
                           mc.cores = n.chains)



Sys.time() - t1
detectCores()
