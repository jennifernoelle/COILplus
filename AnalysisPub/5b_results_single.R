# File to plot results for cv and full run
# Includes new code to plot trace quantities and investigate occurrences
# This version is to examine the primary results, which sample occurrence probs
# Use the compare_results file to look at all results

# --------- TO DO: set  your directories and name the current results files using the date ---------#


# Save results using convention: res_date_i.rda
date <- 'COILp_b_exp_500sims'


# Where the processed data are saved:
data_path <- 'ProcessedDataNew/'
# Where you want to save MCMC results:
results_path <- paste0('ResultsPub/', date, '/')
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


# ------------------------- STEP 1: Specifications. -------------------------- #

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 10

# Number of cv samples per repetition
n.cv <- 100

# ------------- STEP 1: Getting the results together - full sample ---------- #

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
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nV, nP))
occ_p <- array(NA, dim = c(nchains, nP, nS))
occ_v <- array(NA, dim = c(nchains, nV, nS))
p_occ_p <- array(NA, dim = c(nchains, nP, nS))
p_occ_v <- array(NA, dim = c(nchains, nV, nS))
p_accept_p <- array(NA, dim = c(nchains, nP, nS)) 
p_accept_v <- array(NA, dim = c(nchains, nV, nS))

for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
 # Save posterior mean of occurrence indicators for each chain
  occ_p[ii, , ] <- (all_res[[ii]]$occ_plants)$OP_mean
  occ_v[ii, , ] <- (all_res[[ii]]$occ_verts)$OB_mean
   if(sampleP){
   p_occ_p[ii, , ] <- all_res[[ii]]$occ_plants$p_OPs
   p_accept_p[ii, , ] <- all_res[[ii]]$occ_plants$p_accept
   p_occ_v[ii, , ] <- all_res[[ii]]$occ_verts$p_BPs
   p_accept_v[ii, , ] <- all_res[[ii]]$occ_verts$p_accept
 }
}

dimnames(pred_ours)[2 : 3] <- list(vertebrate = rownames(obs_A), plant = colnames(obs_A))
dimnames(p_occ_p)[2:3] <- dimnames(occ_p)[2:3] <- list(plant = colnames(obs_A), dimnames(obs_A)[[3]])
dimnames(p_occ_v)[2:3] <- dimnames(occ_v)[2:3] <- list(vert = rownames(obs_A), dimnames(obs_A)[[3]])

# Calculating the posterior means across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
mean_occ_v <- apply(occ_v, c(2,3), mean) # occ indicators
mean_occ_p <- apply(occ_p, c(2,3), mean)
mean_p_occ_v <- apply(p_occ_v, c(2,3), mean)
mean_p_occ_p <- apply(p_occ_p, c(2,3), mean)
mean_accept_v <- apply(p_accept_v, c(2,3), mean)
mean_accept_p <- apply(p_accept_p, c(2,3), mean)

# Sanity checking occurrence probs 
mean(mean_occ_p[O_P ==1]) # should be 1 regardless of sampleP
mean(mean_occ_v[O_V ==1]) # should be 1 regardless of sampleP
mean(mean_p_occ_p[O_P ==1]) # should be 1 if sampleP, otherwise we don't fill in the array 
mean(mean_p_occ_v[O_V ==1]) # should be 1 even if sampleP == FALSE because occurrence probs are still 0/1
c1 <- p_occ_p[1,,]
mean(c1[O_P==1])
c2 <- p_occ_p[2,,]
mean(c2[O_P==1])


# Setting the recorded interactions to NA (so that they don't overpower the colors)
mean_pred_na[comb_A == 1] <- NA

# EDA with occurrence indicators and probabilities
# We want to see if these update more when we sampleP
plot(c(mean_occ_p), c(O_P)) # Should update except in p=0/1 scenario
plot(c(mean_occ_v), c(O_V)) # Should update except in p=0/1 scenario
plot(c(mean_p_occ_p), c(O_P)) # Should only update if sampleP
plot(c(mean_p_occ_v), c(O_V))  # Should only update if sampleP


# Save the posterior network 
write.csv(mean_pred, paste0(results_path, "post_network", date, ".csv"))

#-------------------- STEP 2: SANITY CHECK: interactions -----------------------------#

# Note that because post prob is taken to be the mean across posterior samples of L
# Zero values are possible when the number of posterior samples is relativley low
# And unique values are limited to multiples of 1/n_samples
mean(mean_pred)
summary(c(mean_pred))
sort(c(mean_pred))[1000] # Might be actually zero: plausible given small number of post samples


# Check: model should be outputting 1 whenever there is a known interaction 
sum(comb_A)
sum(mean_pred == 1)
mean(mean_pred[comb_A == 1])

#write.csv(mean_pred, paste0(results_path, "post_network", date, ".csv"))

sum(comb_A) # observed interactions
sum(comb_A) /length(comb_A)
sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.75, na.rm = TRUE) # v likely non-observed interactions




# ----------------- STEP 3: PLOTTING THE HEATMAP ----------------------------#

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
v_group <- v.taxa$Animal_Family
p_group <- p.taxa$Plant_Family

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set plot_pred to pred_ours for results based on our method 
plot_pred <- mean_pred_na

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_v_size <- 10
min_p_size <- 25

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_v_groups <- names(which(v_size_cluster >= min_v_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(results_path, "heatmap_", date, ".png"), width = 1000, height =1000)
superheat(X = plot_pred[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          bottom.label.size = 0.2, left.label.size = 0.12,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          title = date, 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 20)
#dev.off()




# Plot observed network no plant legend
comb_A_na <- ifelse(comb_A ==1, NA, 0)
#png(paste0(results_path, "Figures/heatmap_observed.png"), width = 1000, height = 510)
superheat(X = comb_A[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 9,
          #bottom.label = "none",
          bottom.label.text.size =10,
          bottom.label.size = 0.6, left.label.size = 0.3,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"),
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
          #title.alignment = "center",
          legend = TRUE
)



# Plot imputed network without title
#png(paste0(results_path, "Figures/heatmap_observed.png"), width = 1000, height = 510)
superheat(X = plot_pred[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 8,
          #bottom.label = "none",
          bottom.label.text.size = 8,
          bottom.label.size = 0.40, left.label.size = 0.18,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"),
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 15, 
          #title.alignment = "center",
          legend = TRUE
)

#--------------------- MCMC DIAGNOSTICS ---------------------------------------#
# Look at mixing etc of log likelihood
burn <- 10000
use_Nsims_unthinned <- length(all_res[[1]]$logL) - burn
logL_full <- matrix(data = NA, nrow = use_Nsims_unthinned + burn, ncol = nchains)
logL <- matrix(data = NA, nrow = use_Nsims_unthinned, ncol = nchains)
for (cc in 1 : nchains) {
  logL[, cc]  <- all_res[[cc]]$logL[(burn + 1): (use_Nsims_unthinned + burn)]
  logL_full[, cc]  <- all_res[[cc]]$logL
}


L_df <- as.data.frame(logL) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("V", "", Chain)) %>% 
  mutate(Chain = as.numeric(Chain)) %>%
  rename("logL" = "value")

L_df_full <- as.data.frame(logL_full) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("V", "", Chain)) %>% 
  mutate(Chain = as.numeric(Chain)) %>%
  rename("logL" = "value")

# Save traceplot
png(paste0(save_plots_path, "trace_full", m, ".png"), width = 1500, height =1000)
mcmc_trace(L_df_full)
dev.off()


# Save traceplot
png(paste0(save_plots_path, "trace", m, ".png"), width = 1500, height =1000)
mcmc_trace(L_df)
dev.off()


# Look at Gelman Rubin Rhat and save
mcmc_diag <- c(Rhat(logL[,]), ess_bulk(logL), ess_tail(logL))
mcmc_diag_full <- c(Rhat(logL_full[,]), ess_bulk(logL_full), ess_tail(logL_full))

mcmc_diag
mcmc_diag_full

# --------------- STEP 4: Taxonomic correlation of latent factors ----------------- #

burn <- 10000 # I saved all samples in this run so I still have to discard

if(nchains > 1){
  all_cor_full <- abind::abind(all_res[[1]]$correlations, all_res[[2]]$correlations, along = 3)
  if(nchains>2){for (cc in 3 : nchains) {
    all_cor_full <- abind::abind(all_cor_full, all_res[[cc]]$correlations, along = 3)
  }
  }
}else{
  all_cor_full <- all_res[[1]]$correlations
}

all_cor <- all_cor_full[-c(1:burn),,]
# Posterior means and 95% credible intervals for the rho parameters in the
# latent factors for bird and plant species:
apply(all_cor, 2, mean)
quantile(c(all_cor[,1,]), probs = c(0.025, 0.975)) #U
quantile(c(all_cor[,2,]), probs = c(0.025, 0.975)) #V

# 
# apply(all_cor[,,1], 2, quantile, probs = c(0.025, 0.975)) #U
# quantile(c(all_cor[,,1]), probs = c(0.025, 0.975)) #U
# apply(all_cor[,,2], 2, quantile, probs = c(0.025, 0.975)) #V
# quantile(c(all_cor[,,2]), probs = c(0.025, 0.975)) #V


# Mixing of the rho parameters
#all_cor <- array(all_cor, dim = dim(all_cor), dimnames =  list("Iterations", "Parameter", "Chain"))
all_cor_full <- provideDimnames(all_cor_full, base = list("Iterations", "Parameter", "Chain"))
all_cor_full <- aperm(all_cor_full, c(1, 3, 2))
all_cor <- all_cor_full[(burn + 1): (use_Nsims_unthinned + burn),,]

u_df <- as.data.frame(all_cor[,,1]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("U" = "value")
u_df_full <- as.data.frame(all_cor_full[,,1]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("U" = "value")

mcmc_trace(u_df)
mcmc_trace(u_df_full)

v_df <- as.data.frame(all_cor[,,2]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("V" = "value")
v_df_full <- as.data.frame(all_cor_full[,,2]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("V" = "value")

mcmc_trace(v_df)# + theme(legend.position = "bottom")
mcmc_trace(v_df_full)# + theme(legend.position = "bottom")

# Look at Gelman Rubin Rhat
Rhat(all_cor_full[-c(1:burn),,1])
ess_bulk(all_cor_full[,,1])
ess_tail(all_cor_full[,,1])

Rhat(all_cor_full[-c(1:burn),,2])
ess_bulk(all_cor_full[,,2])
ess_tail(all_cor_full[,,2])



#---------------- STEP 5: MCMC Diagnostics ------------------------------------#

# Look at mixing etc of log likelihood
use_Nsims_unthinned <- length(all_res[[1]]$logL)
logL <- matrix(data = NA, nrow = use_Nsims_unthinned, ncol = nchains)
for (cc in 1 : nchains) {
  logL[, cc] <- all_res[[cc]]$logL
}


L_df <- as.data.frame(logL) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("V", "", Chain)) %>% 
  mutate(Chain = as.numeric(Chain)) %>%
  rename("logL" = "value")

mcmc_trace(L_df) # remember to remove burn-in


# --------------- STEP 6: Cross validation results ----------------- #

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

