# File to plot results for cv and full run
# We have strange results where the CV performance is identical but full run performance is not

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
# wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
# setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'p75_B4_nb_400sims'
date <- 'p75_C4_b_400sims'


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
use_Nsims <- dim(all_res[[1]]$all_pred)[1]

# Creating an array to bind results across chains:
pred_ours <- array(NA, dim = c(nchains * use_Nsims, nV, nP))
#raw_probs <- array(NA, dim = c(nchains * use_Nsims, nV, nP))
p_occ_p <- array(NA, dim = c(nchains, nP, nS))
p_occ_v <- array(NA, dim = c(nchains, nV, nS))
p_accept_p <- array(NA, dim = c(nchains, nP, nS)) 
p_accept_v <- array(NA, dim = c(nchains, nV, nS))

for (ii in 1 : nchains) {
  # Using the posterior samples of the L matrix:
  pred_ours[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 1] # this is the imputed L's from the model output (included bias corr)
  #raw_probs[1 : use_Nsims + use_Nsims * (ii - 1), , ] <- all_res[[ii]]$all_pred[, , , 2] # these are the raw probs, we don't use, I'm just debugging

   if(sampleP){
   p_occ_p[ii, , ] <- all_res[[ii]]$occ_plants[[1]] 
   p_accept_p[ii, , ] <- all_res[[ii]]$occ_plants[[2]]
   p_occ_v[ii, , ] <- all_res[[ii]]$occ_verts[[1]] 
   p_accept_v[ii, , ] <- all_res[[ii]]$occ_verts[[2]]
 }
}

dimnames(pred_ours)[2 : 3] <- list(vertebrate = rownames(obs_A), plant = colnames(obs_A))
dimnames(p_occ_p)[2:3] <- list(plant = colnames(obs_A), dimnames(obs_A)[[3]])
dimnames(p_occ_v)[2:3] <- list(vert = rownames(obs_A), dimnames(obs_A)[[3]])

# Calculating the posterior means across mcmc iterations:
mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
#mean_prob_raw <- apply(raw_probs, c(2, 3), mean) # just for debugging
mean_p_occ_v <- apply(p_occ_v, c(2,3), mean)
mean_p_occ_p <- apply(p_occ_p, c(2,3), mean)
mean_accept_v <- apply(p_accept_v, c(2,3), mean)
mean_accept_p <- apply(p_accept_p, c(2,3), mean)

# Debugging occurrence probs 
mean(mean_p_occ_p[O_P ==1])
mean(mean_p_occ_v[O_V ==1])
c1 <- p_occ_p[1,,]
mean(c1[O_P==1])
c2 <- p_occ_p[2,,]
mean(c2[O_P==1])


# Setting the recorded interactions to NA (so that they don't overpower the colors)
mean_pred_na[comb_A == 1] <- NA

#-------------------- STEP 2: SANITY CHECK: interactions -----------------------------#

# Note that because post prob is taken to be the mean across posterior samples of L
# Zero values are possible when the number of posterior samples is relativley low
# And unique values are limited to multiples of 1/n_samples
mean(mean_pred)
summary(c(mean_pred))
#summary(c(mean_prob_raw))
#sort(c(mean_prob_raw))[1000] # Small but nonzero 
sort(c(mean_pred))[1000] # Might be actually zero: plausible given small number of post samples


# Check: model should be outputting 1 whenever there is a known interaction 
sum(comb_A)
sum(mean_pred == 1)
mean(mean_pred[comb_A == 1])

write.csv(mean_pred, paste0(results_path, "post_network", date, ".csv"))

sum(comb_A) # observed interactions
sum(comb_A) /length(comb_A)
sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.75, na.rm = TRUE) # v likely non-observed interactions

sum(mean_pred > 0.5) # observed and unobserved likely interactions
sum(mean_pred > 0.5)/length(mean_pred)
sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
sum(mean_pred_na > 0.5, na.rm = TRUE)/sum(!(is.na(mean_pred_na))) # predicted prevalence among non-observed interactions

sum(mean_pred > 0.75) # observed and unobserved of v likely interactions
sum(mean_pred > 0.75)/length(mean_pred)
sum(mean_pred_na > 0.75, na.rm = TRUE) # v likely non-observed interactions
sum(mean_pred_na > 0.75, na.rm = TRUE)/sum(!(is.na(mean_pred_na))) # predicted prevalence among non-observed interactions

#png(file = paste0(results_path, "Histogram", date, ".png"))
hist(c(mean_pred), main = paste0("African Frugivory fit ", date), 
     xlab = "Posterior Interaction Probability")
abline(v = mean(mean_pred))
#dev.off()


#---------------- STEP 3: SANITY CHECK - OCCURRENCE PROBABILITIES ------------------#
#Look at occurrence probabilities
# DEBUG why are we sampling non-1 values when OP == 1?

if(sampleP){
# Occurrence probabilty should be 1 if species was actually observed
mean(mean_p_occ_p[O_P == 1])
mean(mean_p_occ_v[O_V == 1])
par(mfrow = c(1,2))
plot(c(mean_p_occ_p), c(O_P))
plot(c(mean_p_occ_v), c(O_V))
}

par(mfrow = c(1,1))

# I'm no longer saving all samples of occurrence probabilities
# if(sampleP){
# # Look at mixing for a few select probabilities for plants
# png(paste0(results_path, date, "_mixing.png"), width = 2000, height =1000)
# par(mfrow = c(3,4))
# for(p in 1:12){
#   plant <- sample(1:nP,1)
#   study <- sample(1:nS, 1)
#   plot(p_occ_p[,plant, study], type = "l", main = paste0(colnames(obs_A)[plant], " x ", s.names[study]))
#   abline(h = O_P[plant, study], col = "blue", lwd = 3, lty = 2)
#   abline(h = mean_p_occ_p[plant, study], col = "red", lwd = 3, lty = 3)
# }
# dev.off()
# }



# --------------- STEP 4: Cross validation results ----------------- #

# Check to see if these are identical between two different models
# Save results using convention: res_date_i.rda

# Unblocked sampler
date <- 'p75_B4_nb_400sims'
results_path <- paste0('ResultsNew/', date, '/')

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

predB <- pred

# Blocked sampler
date <- 'p75_C4_b_400sims'
results_path <- paste0('ResultsNew/', date, '/')

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

predC <- pred

### Compare: pred is 10x100 matrix storing post probs for each
sum(predC == predB)



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

