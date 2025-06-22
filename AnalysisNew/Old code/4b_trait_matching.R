# Estimating variable importance based on interaction model
# This code is slow, so best to just run in terminal and plot using another interactive file

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
#wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
#setwd(wd_path)


# Paths that are the same for all models 
data_path <- 'ProcessedDataNew/'
save_path <- 'ResultsNew/'
source_path <- 'HelperScriptsJKNew/'

model_names_full <- c('p0_A4_old_400sims_savemore', #for results with occurrence probs
                      'p75_A4_old_400sims_savemore', # for results with occurrence probs
                      'p75_B4_nb_400sims',
                      'p75_C4_b_400sims', # CV is right for this one, but accidentally saved re-ran S200 to save more, so don't have occurrence indicators
                      'pd_A4_old_400sims',
                      'pd_B4_nb_400sims',
                      'pd_C4_b_400sims'
)


# Save results using convention: res_date_i.rda
date <- model_names_full[7]
results_path <- paste0('ResultsNew/', date , '/')



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
B <- 5 #100 permutations
# obs_X = Obs_X
# obs_W = Obs_W
# obs_only <- TRUE
# 

t1 <- Sys.time()
trait_match <- TraitMatching3(B = B, mod_pL1s = mod_pL1s,
                              Xs = NULL, Ws = NULL,  # Imputed values not used.
                              obs_X = obs_X, obs_W = obs_W, obs_only = TRUE)

(Sys.time() - t1)

rsq_resampling_X <- trait_match$rsq_resampling_X
rsq_resampling_W <- trait_match$rsq_resampling_W
rsq_obs_X <- trait_match$rsq_obs_X
rsq_obs_W <- trait_match$rsq_obs_W
corr_obs_X <- trait_match$corr_obs_X
corr_obs_W <- trait_match$corr_obs_W

save(rsq_obs_X, file = paste0(save_path,  'rsq_obs_X_', date, '.dat'))
save(rsq_obs_W, file = paste0(save_path, 'rsq_obs_W_', date, '.dat'))
save(corr_obs_X, file = paste0(save_path,  'corr_obs_X_', date, '.dat'))
save(corr_obs_W, file = paste0(save_path, 'corr_obs_W_', date, '.dat'))
save(rsq_resampling_X, file = paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
save(rsq_resampling_W, file = paste0(save_path, 'rsq_resampling_W_', date, '.dat'))


# Calculating the number of permuted standard deviations away from the mean.

# Starting from the vert covariates:
wh_obs <- rsq_obs_X
wh_resampling <- rsq_resampling_X
sd_awayX <- rep(NA,  length(wh_obs))
names(sd_awayX) <- colnames(Obs_X) #good_namesX
for  (cc in 1 : length(wh_obs)) {
  sd_awayX[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}

# And for the plant covariates:
wh_obs <- rsq_obs_W
wh_resampling <- rsq_resampling_W
sd_awayW <- rep(NA,  length(wh_obs))
names(sd_awayW) <- colnames(Obs_W) #good_namesW
for  (cc in 1 : length(wh_obs)) {
  sd_awayW[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}


# Plotting the tiles of variable importance ordering the variables in
# decreasing importance:

# For the bird species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))

# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#D4DEF7', high = '#425075') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 8))


# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction across the three chains,
# and setting the recorded interactions to NA:
use_Ls <- do.call(abind, c(lapply(all_res, function(x) x$all_pred[, , , 1]), along = 1))
use_mean_Ls <- apply(use_Ls, c(2, 3), mean)
use_mean_Ls[comb_A == 1] <- NA

# Which covariate is to be plotted. The ones we want are listed first.
wh_X <- 1
wh_W <- 1

# Showing only the species that have the covariate measured.
keep_birds <- which(!is.na(obs_X[, wh_X]))
keep_plants <- which(!is.na(obs_W[, wh_W]))
use_out <- use_mean_Ls[keep_birds, keep_plants]

# Because some species have identical values for the covariate, in order for
# plot to show all of them, we need to slightly pertube their values. That way,
# the increasing or decreasing order is not altered, but there is no overlap in
# the covariate values:
bird_cov <- obs_X[keep_birds, wh_X]
bird_cov <- bird_cov + rnorm(length(keep_birds), sd = sd(bird_cov) * 0.0001)
plant_cov <- obs_W[keep_plants, wh_W]
plant_cov <- plant_cov + rnorm(length(plant_cov), sd = sd(plant_cov) * 0.0001)

# Creating a data frame in which the species are ordered by their covariate
# values. This will allow us to plot the probability of interaction across the
# covariates in an interpretable way. We also note that we need to turn the
# covariates to factors in order for them to be plotted in the correct order.
plot_dta <- data.frame(cov_bird = rep(bird_cov, length(keep_plants)),
                       cov_plant = rep(plant_cov, each = length(keep_birds)),
                       probability = as.numeric(use_out))
plot_dta$use_cov_bird <- factor(as.numeric(factor(plot_dta$cov_bird)))
plot_dta$use_cov_plant <- factor(as.numeric(factor(plot_dta$cov_plant)))

g <- ggplot() +
  geom_raster(aes(x = use_cov_plant, y = use_cov_bird, fill = probability), data = plot_dta) +
  scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1)) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Bird information: Increasing Body Mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12)),
                        bottom = textGrob("Plant information: Increasing Fruit Diameter", 
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12)),
                        vp=viewport(width=0.5, height=0.6))
