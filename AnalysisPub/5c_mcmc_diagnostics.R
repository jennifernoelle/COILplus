# This file just looks at the extra logL results saved for 3star

# --------- TO DO: set  your directories and name the current results files using the date ---------#


# Save results using convention: res_date_i.rda
date <- 'pd_B4_nb_400sims_v2' # only saved logL for fastest run
#date <- 'COILp_b_exp_500sims'


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



# ------------- STEP 1: Getting the results together - full sample ---------- #

# Number of chains for the full sample runs: old numbers 20000 burn, 4000 post thinned to 400
nchains <- 4
burn <- 20000
thin <- 10

### Put all chains together in a list
all_res <- NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  all_res[[ii]] <- res
}

# Look at mixing etc of log likelihood
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