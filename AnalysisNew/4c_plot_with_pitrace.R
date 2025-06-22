# File to plot results for cv and full run
# This results file only saves info for trace to be lighter

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
# wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
# setwd(wd_path)

# Save results using convention: res_date_i.rda

#date <- 'p0_A4_old_400sims_savemore' 
#date <- 'p75_A4_old_400sims'
#date <- 'p75_B4_nb_400sims'
#date <- 'p75_C4_b_400sims' # CV is right for this one, but accidentallly saved re-ran S200 to save more, so don't have occurrence indicators
date <- 'pd_A4_old_400sims'
date <- 'pd_B4_nb_400sims'
date <- 'pd_C4_b_400sims'

date <- 'pd_B4_nb_400sims_saveRho_test'
date <- 'pd_B4_nb_400sims_saveRho_MT'
date <- 'pd_C4_b_400sims_saveTrace'


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
library(rstan)

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

# Read in taxonomical data
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


## The following assignments were used in creating obs_OP
# Same study: 1, same site: 0.75
# Same country and habitat: 0.5, same region and habitat: 0.45, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05

## Improved guess: Expert 1 modified
# Same study: 1, same site: 0.75 -> 0.85
# Same country and habitat: 0.5 -> 0.65, same region and habitat: 0.45 -> 0.35, same habitat only: 0.25,
# Same country not habitat: 0.1, same region not habitat: 0.05 
# 
# ## Improved default guess: 0.75/1
# O_P <- ifelse(O_P == 1, 1, 0.75)
# O_V <- ifelse(O_V == 1, 1, 0.75)

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]
s.names <- dimnames(obs_A)[[3]]

# Set sampleP and occurrence prior probs according to the version of the model fit
sampleP <- 1- grepl('A4', date) # A4 indicates the old sampler and so occurrence probs are not sampled

if(grepl('p0', date)){ # Unobserved species are assumed not present
  O_P <- ifelse(O_P == 1, 1, 0)
  O_V <- ifelse(O_V ==1, 1, 0)
}else if(grepl('p75', date)){ # Unobserved species are assumed present with p = 0.75
  O_P <- ifelse(O_P == 1, 1, 0.75)
  O_V <- ifelse(O_V ==1, 1, 0.75)
} # Otherwise keep expert defined probabilities loaded above


# -------------- STEP 1: Specifications. ------------ #

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 10

# Number of cv samples per repetition
n.cv <- 100

# --------------- STEP 1: Getting the results together - full sample  ----------------- #

### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
}

### Bind together predicted interactions and detection probabilities 

# Number of posterior samples used:
use_Nsims <- dim(all_res[[1]]$random_pis[[1]])[2]
use_Npis <- dim(all_res[[1]]$random_pis[[1]])[1]

# Creating an array to bind results across chains:
random_pis <- array(NA, dim = c(use_Npis, use_Nsims, nchains))
random_pi_indices <- matrix(NA, nrow = use_Npis, 2)

for (ii in 1 : nchains) {
 random_pis[,,ii] <- all_res[[ii]]$random_pis[[1]]
 random_pi_indices <- all_res[[ii]]$random_pis[[2]] # these are actually the same in all chains
}

# --------------- STEP 3: Taxonomic correlation of latent factors ----------------- #
burn <- 5000 # I saved all samples in this run so I still have to discard
if(nchains > 1){
  all_cor <- abind::abind(all_res[[1]]$correlations, all_res[[2]]$correlations, along = 3)
  if(nchains>2){for (cc in 3 : nchains) {
    all_cor <- abind::abind(all_cor, all_res[[cc]]$correlations, along = 3)
  }
  }
}else{
  all_cor <- all_res[[1]]$correlations
}

# Posterior means and 95% credible intervals for the rho parameters in the
# latent factors for bird and plant species:
apply(all_cor, 2, mean)
apply(all_cor[,,1], 2, quantile, probs = c(0.025, 0.975)) #U
quantile(c(all_cor[,,1]), probs = c(0.025, 0.975)) #U
apply(all_cor[,,2], 2, quantile, probs = c(0.025, 0.975)) #V
quantile(c(all_cor[,,2]), probs = c(0.025, 0.975)) #V


# Mixing of the rho parameters
#all_cor <- array(all_cor, dim = dim(all_cor), dimnames =  list("Iterations", "Parameter", "Chain"))
all_cor <- provideDimnames(all_cor, base = list("Iterations", "Parameter", "Chain"))
all_cor <- aperm(all_cor, c(1, 3, 2))

u_df <- as.data.frame(all_cor[,,1]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
        mutate(Chain = gsub("Chain", "", Chain)) %>% 
        mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
        rename("U" = "value")

u_df_save <- u_df[-c(1:burn), ] # DEBUG, we actually need to drop the first burn*chains obs doing it this way

mcmc_trace(u_df_save)

v_df <- as.data.frame(all_cor[,,2]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("V" = "value")
v_df_save <- v_df[-c(1:burn),]
mcmc_trace(v_df_save)# + theme(legend.position = "bottom")

# Look at Gelman Rubin Rhat
Rhat(all_cor[-c(1:burn),,1])
ess_bulk(all_cor[,,1])
ess_tail(all_cor[,,1])

Rhat(all_cor[-c(1:burn),,2])
ess_bulk(all_cor[,,2])
ess_tail(all_cor[,,2])



#plot(all_cor[,1], type = "l")

# Trace plots for random probabilities: now the same indices are heldout across chains
if(nchains > 1){
  trace_pis <- abind::abind(t((all_res[[1]]$random_pis)$save_pis), t((all_res[[2]]$random_pis)$save_pis), along = 3)
  #pi_indices <- rbind((all_res[[1]]$random_pis)$save_pi_indices, (all_res[[2]]$random_pis)$save_pi_indices)
  if(nchains>2){for (cc in 3 : nchains) {
    trace_pis <- abind::abind(trace_pis, t(all_res[[cc]]$random_pis$save_pis), along = 3)
    #pi_indices <- rbind(pi_indices, all_res[[cc]]$random_pis$save_pi_indices) # don't need, same across chains
  }
  }
}else{
  trace_pis <- (all_res[[1]]$random_pis)$save_pis
  #pi_indices <- (all_res[[1]]$random_pis)$save_pi_indices
}

pi_indices <- (all_res[[1]]$random_pis)$save_pi_indices
colnames(trace_pis) <- paste0("V", pi_indices[,1], "-P", pi_indices[,2] )
trace_pis <- provideDimnames(trace_pis, base = list("Iterations", "Parameter", "Chain"))
trace_pis <- aperm(trace_pis, c(1, 3, 2))

mcmc_trace(trace_pis)

# Traceplots for running mean of pis
compute_ma <- function(x){
  cx <- cumsum(x)
  rmean <- cx/1:length(cx)
  return(rmean)
}

compute_ma(test)

ma_pis <- apply(trace_pis, c(2,3), compute_ma)
mcmc_trace(ma_pis[5000:nrow(ma_pis), ,])



# --------------- STEP 4: Cross validation results ----------------- #

# Getting the results together (held out indicies and predictions)
all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
our_preds <- array(NA, dim = c(repetitions, nV, nP))
for (rr in 1 : repetitions) {
  cv_indices <- loadRData(paste0(results_path, 'cv_indices_', date, '_', rr, '.dat'))
  pred <- loadRData(paste0(results_path, 'pred_', date, '_', rr, '.dat'))
  all_indices[rr, , ] <- cv_indices
  our_preds[rr, , ] <- pred
}


# Predictions of the held out data from model: 100 indices held out each time - always interactions
pred <- array(NA, dim = c(repetitions, n.cv))
for (rr in 1 : repetitions) {
  for (ii in 1 : n.cv) {
    pred[rr, ii] <- our_preds[rr, all_indices[rr, ii, 1], all_indices[rr, ii, 2]]
  }
}

# Average and median probability of interaction based on the overall data
overall_mean <- cbind(apply(our_preds, 1, mean))
overall_median <- cbind(apply(our_preds, 1, median))

# Average and median in the held out data.
pred_mean <- apply(pred, 1, mean)
pred_median <- apply(pred, 1, median)

# Pseudo precision: 
mean(pred_mean/overall_mean)
mean(pred_median/overall_median)

# Pseudo accuracy: 
#mean(pred_mean)
sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75
#sum(pred>0.85)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75
#sum(pred>0.95)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75
#pred_mean


# Creating the data frame we will plot:
plot_dta <- data.frame(value = rbind(pred_mean / overall_mean, pred_median / overall_median), 
                       stat = rep(c('Pred:Overall (mean)', 'Pred:Overall (median)'), repetitions))

# Plotting cross validation results:
png(filename = paste0(results_path, "cv_res_", date, ".png"))
ggplot(data = plot_dta) +
  geom_boxplot(aes(x = stat, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  ggtitle('Out of sample performance', subtitle = date) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) + 
  theme(text = element_text(size = 20))
dev.off()



# --------------- STEP 5: Performing trait matching ----------------- #
# Using the linear predictor of the interaction model (not corrected for bias due to detection, number of studies):
mod_pL1s <- array(NA, dim = c(nchains * Nsims, nV, nP)) 
for (cc in 1 : nchains) {
  wh_entries <- Nsims * (cc - 1) + 1 : Nsims
  mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2]
}

# Dealing with extreme values for which logit(x) is infinite.
mod_pL1s[mod_pL1s > 1 - 10^{-10}] <- 1 - 10^{-10}

t1 <- Sys.time()
trait_match <- TraitMatching2(B = 100, mod_pL1s = mod_pL1s,
                              Xs = NULL, Ws = NULL,  # Imputed values not used.
                              obs_X = Obs_X, obs_W = Obs_W, obs_only = TRUE)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X_', date, '.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W_', date, '.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W_', date, '.dat'))

print(Sys.time() - t1)
#--------------------------- MORE LATENT FACTORS --------------------------------#

## Overview of latent factors for animals
if(nchains > 1){
  all_u <- abind::abind(all_res[[1]]$factors[[1]], all_res[[2]]$factors[[1]], along = 3)
  all_v <- abind::abind(all_res[[1]]$factors[[2]], all_res[[2]]$factors[[2]], along = 3)  
  if(nchains>2){for (cc in 3 : nchains) {
    all_u <- abind::abind(all_u, all_res[[cc]]$factors[[1]], along = 3)
    all_v <- abind::abind(all_v, all_res[[cc]]$factors[[2]], along = 3)
  }
  }
}else{
  all_u <- all_res[[1]]$factors[[1]]
  all_v <- all_res[[1]]$factors[[2]]
}

# Look at animal latent factors
mean_u <- apply(all_u, c(1,2), mean)
u_df <- data.frame(factor = paste0("Factor ", 1:10), MeanAbsValue = colMeans(abs(mean_u)), MeanValue = colMeans(mean_u))
u_df$factor <- factor(u_df$factor, levels = u_df$factor)
ggplot(u_df, aes(x = factor, y= MeanAbsValue)) + 
  geom_bar(stat = "identity") +
  theme_minimal() + 
  ggtitle("Animal Latent Factors") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 


u_df_animals <- data.frame(mean_u)
colnames(u_df_animals) <- paste0("Factor_", 1:10)
u_df_animals$Mammal <- m.names
ggplot(u_df_animals, aes(x = Factor_1, y = Factor_2, label = Mammal)) + 
  geom_point() + 
  geom_label(data = subset(u_df_animals, Mammal %in% c("Papio_anubis", "Papio_cynocephalus"))) + 
  theme_minimal() + 
  ggtitle("Animal Latent Factors") + 
  theme(text = element_text(family = "serif", size = 20)) 


# Look at plant latent factors
mean_v <- apply(all_v, c(1,2), mean)
v_df <- data.frame(factor = paste0("Factor ", 1:10), MeanAbsValue = colMeans(abs(mean_v)), MeanValue = colMeans(mean_v))
v_df$factor <- factor(v_df$factor, levels = v_df$factor)
ggplot(v_df, aes(x = factor, y= MeanAbsValue)) + 
  geom_bar(stat = "identity") +
  theme_minimal() + 
  ggtitle("Plant Latent Factors") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90), text = element_text(family = "serif", size = 20)) 


v_df_plants <- data.frame(mean_v)
colnames(v_df_plants) <- paste0("Factor_", 1:10)
v_df_plants$Species <- plant.names 
v_df_plants <- left_join(v_df_plants, traits.plants)
ggplot(v_df_plants, aes(x = Factor_1, y = Factor_2, color = Order)) + 
  geom_point() + 
  #geom_label(data = subset(u_df_animals, Mammal %in% c("Papio_anubis", "Papio_cynocephalus"))) + 
  theme_minimal() + 
  ggtitle("Plant Latent Factors") + 
  theme(text = element_text(family = "serif", size = 20)) 

