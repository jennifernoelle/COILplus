# File to get CV data from the longest CV runs (400 posterior samples)
# And plot the results
# This version adds nicer heatmaps and results for occurrence with the one where they were missing before

# ---------- TO DO: set  your directories and list of models to analyze ---------#

# Paths that are the same for all models 
data_path <- 'ProcessedDataPub/'
save_plots_path <- 'ResultsPub/'
source_path <- 'HelperScriptsPub/'

# Set models to loop over for full results
model_names <- c(#'COIL_0_100_500sims', 
                 'COIL_75_100_500sims', 
                 'COILp_b_75_100_500sims',
                 #'COILp_nb_75_100_500sims', 
                 'COIL_exp_500sims', 
                 'COILp_b_exp_500sims'
                 #'COILp_nb_exp_500sims'
                 )


# Some runs failed so we need to temporarily skip:
model_names_skip <- NULL

# Set models to loop over for cross validation
model_names_cv <- c(#'COIL_p0_500sims', 
                 'COIL_75_100_500sims', 
                 'COILp_b_75_100_500sims',
                 #'COILp_nb_75_100_500sims', 
                 'COIL_exp_500sims', 
                 'COILp_b_exp_500sims'
                 #'COILp_nb_exp_500sims'
                 )

# Some CV runs failed so we temporarily skip over
model_names_cv_skip <- NULL #c('COIL_75_100_500sims') # Rerunning due to typo in .sh, should be done 5/29

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 10

# Number of cv samples per repetition
n.cv <- 100

# Number of burn-in iterations: use this because we save all logL
burn <- 10000


# -------------------- STEP 0: Load functions, packages, data ---------------- #

# We want the github version of superheat
if(!("superheat" %in% rownames(installed.packages()))){
  install.packages("devtools")
  remotes::install_github("rlbarter/superheat")
}

library(ggplot2)
library(ggridges)
library(ggplotify)
library(RColorBrewer)
library(ggnewscale)
library(gplots)
library(xtable)
library(superheat)
library(abind)
library(gridExtra)
library(grid)
library(dplyr)
library(caret)
library(pROC)
library(rstan)
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
O_P_expert <- O_P
O_V_expert <- O_V

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

# Species occurs in the same study -> 1
# Species occurs at the same site, but a different study -> 0.75
# Species occurs in the same country, but a different site -> 0.5
# Species occurs in the same zone, but different ... -> 0.25
O_labels_expert <- c("Different zone", "Same zone only", 
                       "Same country only","Same site", 
                       "Same study")

#---------------------- STEP 1. SET UP USEFUL THINGS ---------------------------#

# Getting the combined network for the interactions recorded in any study
comb_A <- apply(obs_A, c(1, 2), sum)
comb_A <- (comb_A > 0) * 1

# Useful values
nP <- ncol(obs_A)
nV <- nrow(obs_A)
nS <- dim(obs_A)[3]
s.names <- dimnames(obs_A)[[3]]

# Storage for performance metrics
precision_df <- data.frame(matrix(ncol = 4,  nrow = 0))
recall_df <- data.frame(matrix(ncol = 4,  nrow = 0))
postmean_df <- data.frame(matrix(ncol = 4,  nrow = 0))

# Storage for new predictions
new_preds_df <- data.frame(matrix(ncol = 4, nrow = 0))

# Storage for density plotting data
dplot_list <- histplot_list <- list()
preds_combined <- mean_density_combined <- data.frame(matrix(ncol = 4, nrow = 0))
Blues10 <- colorRampPalette(c("lightblue", "blue"))(10)

# Storage for occurrence indicators summaries
p_occ_df <- v_occ_df <- data.frame(matrix(ncol = 4, nrow = 0))

# Storage for mcmc diagnostics
mcmc_df <- data.frame(matrix(ncol = 4, nrow = 0))
mcmc_df_full <- data.frame(matrix(ncol = 4, nrow = 0))

## Heatmap setup

# Start by sorting the taxa by family if they aren't already (animals need some cleaning, plants are fine)
v.taxa.ordered <- v.taxa %>% mutate(original_row = row_number()) %>% 
  arrange(Animal_Class, Animal_Family, Animal_Genus, Animal_Species_Corrected)
use.codes <- sort(paste(LETTERS, rep(1:3, 26), sep = "-")) # Create codes that will alphabetize nicely
v.family.codebook <- data.frame(Animal_Family = unique(v.taxa.ordered$Animal_Family)) %>% 
  mutate(v.show.order = paste0("Family ", use.codes[row_number()])) %>% 
  right_join(., v.taxa.ordered)


# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
v_group <- v.family.codebook$v.show.order # use code names so they'll be in the right order when sorted alph by superheat
p_group <- p.taxa$Plant_Family
v.taxa.order <- v.taxa.ordered$original_row

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size
v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
v_size_df <- data.frame(v_size_cluster) %>% rownames_to_column(., var = "family.code")
v.family.codebook <- left_join(v.family.codebook, v_size_df, join_by(v.show.order == family.code))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_v_size to 10, and min_plant_size to 25. Additionally, 
# we retain the smaller clusters Hominidae and Elephantidae because they are of
# scientific interest. Setting both to 0 will produce the full results.
min_v_size <- 12
min_p_size <- 30

which_v_groups_extra <- unique(v.family.codebook[which(v.family.codebook$Animal_Family == "Elephantidae" 
                                                       | v.family.codebook$Animal_Family == "Hominidae"), "v.show.order"])
keep_v_groups <- sort(c(names(which(v_size_cluster >= min_v_size)), which_v_groups_extra))
keep_p_groups <- names(which(p_size_cluster >= min_p_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


#--------------------- STEP 2. LOAD EACH MODEL AND ANALYZE---------------------#
counter <- 1
for(ii in seq_along(model_names_cv)){
  m <- model_names[ii]
  m_cv <- model_names_cv[ii]
  skip_full <- m %in% model_names_skip
  skip_cv <- m_cv %in% model_names_cv_skip
  
  cat("\n Analyzing results for ", m)
  
  #### Basic setup for analyzing this model 

  # Where are the model results stored
  results_path  <- paste0('ResultsPub/', m , '/')
  
  # Set sampleP and occurrence prior probs according to the version of the model fit
  sampleP <- 1 - grepl('COIL_', m) # COIL_ indicates the old sampler and so occurrence probs are not sampled
  
  if(grepl('0_100', m)){ # Unobserved species are assumed not present
    O_P <- ifelse(O_P_expert == 1, 1, 0)
    O_V <- ifelse(O_V_expert ==1, 1, 0)
    O_labels <- c("Different study","Same study")
  }else if(grepl('75_100', m)){ # Unobserved species are assumed present with p = 0.75
    O_P <- ifelse(O_P_expert == 1, 1, 0.75)
    O_V <- ifelse(O_V_expert ==1, 1, 0.75)
    O_labels <- c("Different study","Same study")
  }else{
    O_P <- O_P_expert
    O_V <- O_V_expert
    O_labels <- O_labels_expert
  }
  
  # Set which occurrence probs were used
  scenario <- ifelse(grepl('0_100', m), "0/100%", 
                     ifelse(grepl('75_100', m), "75/100%", "Expert"))
  
  # Set which sampler was used
  sampler <- ifelse(grepl("COIL_", m), "Old", ifelse(grepl("COILp_b_", m), "New, blocked", "New, unblocked"))
  
  
  #### Load the full-sample results to analyze occurrence indicators and heatmamps
  if(!skip_full){
  
  all_res <- NULL
  for (ii in 1 : nchains) {
    cat("\n Chain = ", ii)
    res <- loadRData(paste0(results_path, 'res_',  m, '_', ii, '.dat')) 
    cat(" Dim is ", dim(res[[1]]))
    all_res[[ii]] <- res
  }
  
  # Bind together predicted interactions and detection probabilities 
  use_Nsims <- dim(all_res[[1]]$all_pred)[1]  # Number of posterior samples used
  
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
    #if(m_full != 'p75_C4_b_400sims'){ # start big if statement to exclude this bad model from occurrence plots
    occ_p[ii, , ] <- (all_res[[ii]]$occ_plants)$OP_mean
    occ_v[ii, , ] <- (all_res[[ii]]$occ_verts)$OB_mean
    if(sampleP){
      p_occ_p[ii, , ] <- all_res[[ii]]$occ_plants$p_OPs
      p_accept_p[ii, , ] <- all_res[[ii]]$occ_plants$p_accept
      p_occ_v[ii, , ] <- all_res[[ii]]$occ_verts$p_BPs
      p_accept_v[ii, , ] <- all_res[[ii]]$occ_verts$p_accept
    }
   # }# End if statement for occurrences to omit bad model
  }
  
  dimnames(pred_ours)[2 : 3] <- list(vertebrate = rownames(obs_A), plant = colnames(obs_A))
  dimnames(p_occ_p)[2:3] <- dimnames(occ_p)[2:3] <- list(plant = colnames(obs_A), dimnames(obs_A)[[3]])
  dimnames(p_occ_v)[2:3] <- dimnames(occ_v)[2:3] <- list(vert = rownames(obs_A), dimnames(obs_A)[[3]])
  
  # Calculating the posterior means across mcmc iterations:
  mean_pred <- mean_pred_na <- apply(pred_ours, c(2, 3), mean)
  mean_pred_na[comb_A == 1] <- NA
  mean_occ_v <- apply(occ_v, c(2,3), mean) # occ indicators
  mean_occ_p <- apply(occ_p, c(2,3), mean)
  mean_p_occ_v <- apply(p_occ_v, c(2,3), mean)
  mean_p_occ_p <- apply(p_occ_p, c(2,3), mean)
  mean_accept_v <- apply(p_accept_v, c(2,3), mean)
  mean_accept_p <- apply(p_accept_p, c(2,3), mean)
  #rm(all_res)
  
  ## Summarize new predictions
  pred_likely <- sum(mean_pred_na > 0.5, na.rm = TRUE) # likely non-observed interactions
  pred_vlikely <- sum(mean_pred_na > 0.75, na.rm = TRUE) # v likely non-observed interactions
  
  dta_new_preds <- data.frame(scenario = scenario, sampler = sampler, stat = c('Phat > 50%', 'Phat > 75%'), 
                              value = c(pred_likely, pred_vlikely))
  new_preds_df <- rbind(new_preds_df, dta_new_preds)
  
  
  
  ## Plant occurrences
  
  # Posterior mean occurrence indicators by prior 
  df1a <- data.frame(obs_op = c(O_P), post_op = c(mean_occ_p), 
                    color = factor(O_P, labels = O_labels), 
                    type = "Occurrence baseline probabilities") %>% 
    filter(obs_op != 1) # Don't care about known present species
  # Posterior mean occurrence indicators by actual scenario
  df1b <- data.frame(obs_op = c(O_P_expert), post_op = c(mean_occ_p), 
                     color = factor(O_P_expert, labels = O_labels_expert), 
                     type = "Occurrence scenario") %>% 
    filter(obs_op != 1)
 
  # Save results to the summary data frame
  occ_means_plants <- data.frame(mean_post_occ = c(mean_occ_p), where_observed = c(O_P_expert)) %>%
                      group_by(where_observed) %>%
                      summarize_all(mean) %>% 
                      mutate(sampler = sampler, prior = scenario)
  
  p_occ_df <- rbind(p_occ_df, occ_means_plants)
  
  ## Vertebrate occurrences
   # Posterior mean occurrence indicators by prior 
  df1a <- data.frame(obs_op = c(O_V), post_op = c(mean_occ_v), 
                     color = factor(O_V, labels = O_labels), 
                     type = "Occurrence baseline probabilities") %>% 
    filter(obs_op != 1) # Don't care about known present species
  # Posterior mean occurrence indicators by actual scenario
  df1b <- data.frame(obs_op = c(O_V_expert), post_op = c(mean_occ_v), 
                     color = factor(O_V_expert, labels = O_labels_expert), 
                     type = "Occurrence scenario") %>% 
    filter(obs_op != 1)
  
  # Save results to the summary data frame
  occ_means_verts <- data.frame(mean_post_occ = c(mean_occ_v), where_observed = c(O_V_expert)) %>%
    group_by(where_observed) %>%
    summarize_all(mean) %>% 
    mutate(sampler = sampler, prior = scenario)
  
  v_occ_df <- rbind(v_occ_df, occ_means_verts)
  
  
  ### HEATMAPS
  
  #  New heatmaps
    # Keep hominids and elephants too
    # Sort animal families intuitively
  
  
  # Start by sorting the taxa by family if they aren't already (animals need some cleaning, plants are fine)
  v.taxa.ordered <- v.taxa %>% mutate(original_row = row_number()) %>% 
    arrange(Animal_Class, Animal_Family, Animal_Genus, Animal_Species_Corrected)
  use.codes <- sort(paste(LETTERS, rep(1:3, 26), sep = "-")) # Create codes that will alphabetize nicely
  v.family.codebook <- data.frame(Animal_Family = unique(v.taxa.ordered$Animal_Family)) %>% 
    mutate(v.show.order = paste0("Family ", use.codes[row_number()])) %>% 
    right_join(., v.taxa.ordered)
  
  
  # Creating the clusters that will be used
  # The following two lines specify that horizontal and vertical lines in our
  # plot will separate species by taxonomic families:
  v_group <- v.family.codebook$v.show.order # use code names so they'll be in the right order when sorted alph by superheat
  p_group <- p.taxa$Plant_Family
  v.taxa.order <- v.taxa.ordered$original_row

  plot_pred <- mean_pred_na[v.taxa.order,]
  
  # Calculating the size of each cluster, will be used when plotting results
  # for families of certain size
  
  v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
  v_size_df <- data.frame(v_size_cluster) %>% rownames_to_column(., var = "family.code")
  v.family.codebook <- left_join(v.family.codebook, v_size_df, join_by(v.show.order == family.code))
  p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))
  
  # Set the minimum cluster size that should be plotted. For the results of the
  # manuscript, we set min_v_size to 10, and min_plant_size to 25. Additionally, 
  # we retain the smaller clusters Hominidae and Elephantidae because they are of
  # scientific interest. Setting both to 0 will produce the full results.
  min_v_size <- 12
  min_p_size <- 30
  
  which_v_groups_extra <- unique(v.family.codebook[which(v.family.codebook$Animal_Family == "Elephantidae" 
                                                         | v.family.codebook$Animal_Family == "Hominidae"), "v.show.order"])
  keep_v_groups <- sort(c(names(which(v_size_cluster >= min_v_size)), which_v_groups_extra))
  keep_p_groups <- names(which(p_size_cluster >= min_p_size))
  
  keep_v_index <- which(v_group %in% keep_v_groups)
  keep_p_index <- which(p_group %in% keep_p_groups)
  
  # Plot with legends etc
  # Plotting those with minimum size as specified:
  # We replaced observed interactions with black NA
  png(paste0(save_plots_path, "heatmap_withlegend_", m, ".png"), width = 1500, height =1000)
  superheat(X = plot_pred[keep_v_index, keep_p_index],
            pretty.order.cols = FALSE, 
            pretty.order.rows = FALSE,
            membership.rows = v_group[keep_v_index],
            membership.cols = p_group[keep_p_index],
            grid.hline.col = "#00257D", grid.vline.col = '#00257D',
            grid.hline.size = 0.75, grid.vline.size = 0.75,
            left.label.text.size =15,
            bottom.label.text.size = 15,
            bottom.label.size = 0.55, left.label.size = 0.3,
            bottom.label.text.angle = 90,
            #legend.breaks = seq(0, 1, by = 0.2),
            #legend.vspace = 0.05,
            legend.width = 2.5,
            legend.height = 0.15,
            legend.text.size = 35,
            heat.col.scheme = "grey", heat.na.col = 'black',
            heat.pal.values = seq(0, 1, by = 0.05)
            #title = date, 
            #left.label.text.size = 20, 
            #bottom.label.text.size = 20
  )
  dev.off()
  
  # Plot without legends etc
  # Plotting those with minimum size as specified:
  # We replaced observed interactions with black NA
  png(paste0(save_plots_path, "heatmap_nolegend_", m, ".png"), width = 1500, height =1000)
  superheat(X = plot_pred[keep_v_index, keep_p_index],
            pretty.order.cols = FALSE, 
            pretty.order.rows = FALSE,
            membership.rows = v_group[keep_v_index],
            membership.cols = p_group[keep_p_index],
            grid.hline.col = "#00257D", grid.vline.col = '#00257D',
            grid.hline.size = 0.75, grid.vline.size = 0.75,
            bottom.label = "none", left.label = "none", legend = FALSE,
            # left.label.text.size =15,
            # bottom.label.text.size = 15,
            # bottom.label.size = 0.55, left.label.size = 0.3,
            # bottom.label.text.angle = 90,
            #legend.breaks = seq(0, 1, by = 0.2),
            #legend.vspace = 0.05,
            # legend.width = 2.5,
            # legend.height = 0.15,
            # legend.text.size = 35,
            heat.col.scheme = "grey", heat.na.col = 'black',
            heat.pal.values = seq(0, 1, by = 0.05)
            #title = date, 
            #left.label.text.size = 20, 
            #bottom.label.text.size = 20
  )
  dev.off()
  
  ### MCMC Diagnostics
  
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
  
  # Save traceplot: only saved these quantities for main models
  png(paste0(save_plots_path, "trace_full", m, ".png"), width = 1500, height =1000)
  if(sum(is.na(L_df_full$logL))==0) {mcmc_trace(L_df_full)}
  dev.off()
  
  
  # Save traceplot
  png(paste0(save_plots_path, "trace", m, ".png"), width = 1500, height =1000)
  if(sum(is.na(L_df$logL))==0) mcmc_trace(L_df)
  dev.off()
  
  
  # Look at Gelman Rubin Rhat and save
  mcmc_diag <- c(Rhat(logL[,]), ess_bulk(logL), ess_tail(logL))
  mcmc_diag_full <- c(Rhat(logL_full[,]), ess_bulk(logL_full), ess_tail(logL_full))
  
  mcmc_df <- rbind(mcmc_df, c(m, mcmc_diag))
  mcmc_df_full <- rbind(mcmc_df_full, c(m, mcmc_diag_full))
  }

  #### Load the CV results 
  
  if(!skip_cv){
  
  # Getting the results together (held out indices and predictions)
  all_indices <- array(NA, dim = c(repetitions, n.cv, 2))
  our_preds <- array(NA, dim = c(repetitions, nV, nP))
  for (rr in 1 : repetitions) {
    cv_indices <- loadRData(paste0(results_path, 'cv_indices_', m_cv, '_', rr, '.dat'))
    pred <- loadRData(paste0(results_path, 'pred_', m_cv , '_', rr, '.dat'))
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
  
  # Average and median in the held out data for each cv replicate.
  pred_mean <- apply(pred, 1, mean)
  pred_median <- apply(pred, 1, median)
  
  
  ## Plot density of predictions
  
  # Compute mean density curve across CV replicates
  den_allreps <- apply(pred, 1, density)
  z <- seq(0, 1, length.out=100)
  zavg <- rowMeans(sapply(den_allreps, function(d) approx(d$x, d$y, z, yleft=0, yright=0)$y))
  df_dens_avg <- data.frame(z = z, zavg = zavg, scenario = scenario, sampler = sampler)
  
  # Prepare data frame with predictions for each replicate 
  preds_df <- as.data.frame(t(pred)) %>% 
    pivot_longer(., cols = everything(), names_to = "replicate") %>%
    mutate(scenario = scenario, sampler = sampler)
  
  # this.plot <-  ggplot() + 
  #   stat_density(aes(x = value, color = replicate), geom = 'line', alpha = 0.25, data = preds_df, position = 'identity') + 
  #   geom_line(aes(x = z, y = zavg), data = df_dens_avg) +
  #   ggtitle(paste0(scenario, ", Sampler: ", sampler)) + 
  #   theme_minimal() + 
  #   theme(legend.position = "none")
  # 
  # hist.plot <- ggplot() + 
  #   geom_histogram(aes(x = value, fill = as.factor(replicate)), alpha = 0.25, 
  #                  data = preds_df)  +
  #   ggtitle(paste0(scenario, ", Sampler: ", sampler)) + 
  #   theme_minimal() + 
  #   theme(legend.position = "none")
  # 
  # dplot_list[[counter]] <- this.plot
  # histplot_list[[counter]] <- hist.plot
  # 
  # Save the post preds and mean density curve for this model and scenario
  preds_combined <- rbind(preds_combined, preds_df)
  mean_density_combined <- rbind(mean_density_combined, df_dens_avg)
  
    ## True precision is: TP / (TP + FP) but we have presence-only so all of our FP's are uncertain
  # Pseudo precision: 
  mean(pred_mean/overall_mean)
  mean(pred_median/overall_median)
  
  ## True recall is TP/(TP + FN) - we can actually compute this for our heldout true interactions
  # Recall combing all repetitions:
  sum(pred>0.5)/length(pred) # what proportion of true interactions are predicted as "likely" >0.5
  sum(pred>0.75)/length(pred) # what proportion of true interactions are predicted as "v likely" >0.75

  # Recall for each repetition:
  pred_recall_50 <- apply(pred, 1, function(x) sum(x > 0.5)/100) # just out of curiosity, looking at a few thresholds
  pred_recall_75 <- apply(pred, 1, function(x) sum(x > 0.75)/100) # just out of curiosity, looking at a few thresholds
  pred_recall_binom <- apply(pred, 1, function(x) sum(rbinom(n = 1, size = 100, prob = x))/100)
  
  # Creating the data frame we will plot:
  plot_dta_precision <- data.frame(value = rbind(pred_mean / overall_mean), 
                         stat = rep(c('Pred:Overall (mean)'), repetitions), 
                         scenario = scenario, sampler = sampler)
  overallmean_dta <- data.frame(value = c(overall_mean, pred_mean), 
                                   stat = rep(c('Overall mean', 'Heldout Mean'), c(repetitions,repetitions)), 
                                   scenario = scenario, sampler = sampler)
  
  plot_dta_recall <- data.frame(value = c(pred_recall_50, pred_recall_75, pred_recall_binom), 
                                stat = rep(c('Recall at 50%', 'Recall at 75%', 'Recall sampled'), times = c(repetitions, repetitions, repetitions)), 
                                scenario = scenario, sampler = sampler)
  
  precision_df <- rbind(precision_df, plot_dta_precision)
  postmean_df <- rbind(postmean_df, overallmean_dta)
  recall_df <- rbind(recall_df, plot_dta_recall)
  }
  
  counter <- counter + 1
}


#------------------------- OOS Performance Metrics ----------------------------#


### Main Text Plots

# Main text: give big picture overview by comparing new-blocked sampler expert scenario, old sampler expert scenario, 0/1 old sampler

# We see big advantage for the new sampler w.r.t pseudo-precision
png(filename = paste0(save_plots_path, 'precision_compared.png'), width = 500, height = 500)
precision_df %>% 
  filter((scenario == "Expert" | scenario == "0/100%") &
           (sampler == "Old" | sampler == "New, blocked")) %>% 
  mutate(model = paste0(scenario, ", ", sampler)) %>% 
  mutate(model = ifelse(model == "0/100%, Old", "COIL 0/100", 
                        ifelse(model == "Expert, Old", "COIL Exp", "COIL+ Exp"))) %>%
  ggplot() +
  geom_boxplot(aes(x = model, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = c(0,15), n.breaks = 6) +
  theme(text = element_text(size = 24, family = 'serif'))
dev.off()

# Repeat with unblocked sampler to see if it does better now - NOPE
precision_df %>% 
  filter((scenario == "Expert" | scenario == "0/100%") &
           (sampler == "Old" | sampler == "New, blocked" | sampler == "New, unblocked")) %>% 
  mutate(model = paste0(scenario, ", ", sampler)) %>% 
  mutate(model = ifelse(model == "0/100%, Old", "COIL 0/100", 
                        ifelse(model == "Expert, Old", "COIL Exp", 
                        ifelse(model == "Expert, New, unblocked", "COIL+ ExpUB", "COIL+ Exp")))) %>%
  ggplot() +
  geom_boxplot(aes(x = model, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = c(0,15), n.breaks = 6) +
  theme(text = element_text(size = 24, family = 'serif'))

# Note that old sampler appears good via recall, we see it has terrible precision, old and new sampler are very similar
png(filename = paste0(save_plots_path, 'recall_compared.png'), width = 500, height = 500)
recall_df %>% 
  filter((scenario == "Expert-defined occurrence probabilities" | scenario == "0/100% occurrence probabilities") &
           (sampler == "Old" | sampler == "New, blocked") & 
           stat == "Recall at 75%") %>% 
  mutate(model = paste(scenario, ", ", sampler)) %>% 
  mutate(model = ifelse(model == "0/100% occurrence probabilities ,  Old", "A", 
                        ifelse(model == "Expert-defined occurrence probabilities ,  Old", "B", "C"))) %>%
  ggplot() +
  geom_boxplot(aes(x = model, y = value)) +
  theme_bw() +
  ylab('Recall at 75%') +
  xlab('') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = c(0,15), n.breaks = 6) +
  theme_bw() + 
  theme(text = element_text(size = 24, family = 'serif'))
dev.off()


### SM Plots: add performance for 75/100 prior
# Put it next to the Expert prior so we can see that the recall is really inferior with default prior vs expert

## Recall: proportion of true heldout observations recovered

# Combined plot: similar performance across models w.r.t recall but expert scenarios give much better recall 
png(filename = paste0(save_plots_path, 'recall_compared_SM.png'), width = 1000, height = 1000)
recall_df %>% filter(scenario != "0/100%" & # these are too far off they make the plot harder to read
                      stat != "Recall sampled" & # The models are too close so this is too much noise
                      sampler != "New, unblocked") %>% # We're going with the blocked sampler
  mutate(sampler = ifelse(sampler == "Old", "COIL", "COIL+"), 
         scenario = ifelse(scenario == "75/100%", "75/100% occurrence prior", "Expert occurrence prior"), 
         stat = paste0(stat, " threshold")) %>% 
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ stat + scenario) + 
  #ggtitle('Out of sample performance: recall') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 24, family = 'serif'))
dev.off()


# Combined plot: similar peformance with the 75/100 scenario, but strong advantage with expert scenarios; 75/100 gives better precision, 
# but after a certain point, does it really matter
png(filename = paste0(save_plots_path, 'precision_compared_SM.png'), width = 1000, height = 500)
precision_df %>% filter(scenario != "0/100%" & # these are too far off they make the plot harder to read
                           stat != "Recall sampled" & # The models are too close so this is too much noise
                           sampler != "New, unblocked") %>% # We're going with the blocked sampler
                mutate(sampler = ifelse(sampler == "Old", "COIL", "COIL+"), 
                       scenario = ifelse(scenario == "75/100%", "75/100% occurrence prior", "Expert occurrence prior"), 
                       stat = paste0(stat, " threshold")) %>% 
ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  facet_wrap(~ scenario) + 
  #ggtitle('Out of sample performance: precision') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'))
dev.off()


# OOS Performance Metrics table: summarized in SM Table 4.
recall_summary <- recall_df %>% group_by(stat, scenario, sampler) %>% 
              summarize(mean_value = mean(value)) %>%
              filter(stat != "Recall sampled")
precision_summary <- precision_df%>% 
  group_by(stat, scenario, sampler) %>% 
  summarize(mean_value = round(mean(value), 3), .groups = "drop")

precision_summary <- precision_df %>% # convert to character to stop rounding
  group_by(stat, scenario, sampler) %>%
  summarize(mean_value = sprintf("%.3f", mean(value)), .groups = "drop")

            
postmean_df %>% group_by(scenario, sampler, stat) %>% 
  summarize(mean_pi = mean(value))


# Number of new predictions for SM Table 4
new_preds_df

# Posterior density plots for CV interactions
grid.arrange(grobs = dplot_list[2:7], nrow = 2)
grid.arrange(grobs = histplot_list[2:7], nrow = 2)



#------------------------ OCCURRENCE PLOTS ------------------------------------#

p_occs_summary <- p_occ_df %>% pivot_wider(names_from = c(sampler, prior), values_from = mean_post_occ) %>%
             mutate(where_observed = ifelse(where_observed == 0, 'Different zone', 
                                     ifelse(where_observed == 0.25, 'Same zone only', 
                                     ifelse(where_observed == 0.5, 'Same country only', 
                                     ifelse(where_observed == 0.75, 'Same site', 'Same study'))))) %>% 
            mutate(taxa = "plants")
v_occs_summary <- v_occ_df %>% pivot_wider(names_from = c(sampler, prior), values_from = mean_post_occ) %>%
  mutate(where_observed = ifelse(where_observed == 0, 'Different zone', 
                                 ifelse(where_observed == 0.25, 'Same zone only', 
                                        ifelse(where_observed == 0.5, 'Same country only', 
                                               ifelse(where_observed == 0.75, 'Same site', 'Same study'))))) %>% 
  mutate(taxa = "frugivores")

# Make nice summary tables for only the key models
print(xtable(p_occs_summary))
print(xtable(v_occs_summary))
print(xtable(rbind(p_occs_summary, v_occs_summary)))

# Compare new and old sampler under expert occurrence scenario
grid.arrange(grobs = occ_p_plot_box_list[c(5,7)], nrow = 1)
grid.arrange(grobs = occ_v_plot_box_list[c(5,7)], nrow = 1)

grid.arrange(grobs = occ_p_plot_ridges_list[c(5,7)], nrow = 1)
occ_p_plot_ridges_list[7]







 # OTHER STUFF BELOW


# --------------- STEP 3: Taxonomic correlation of latent factors ----------------- #
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
apply(all_cor, 2, quantile, probs = c(0.025, 0.975))


# Mixing of the rho parameters
#all_cor <- array(all_cor, dim = dim(all_cor), dimnames =  list("Iterations", "Parameter", "Chain"))
all_cor <- provideDimnames(all_cor, base = list("Iterations", "Parameter", "Chain"))
all_cor <- aperm(all_cor, c(1, 3, 2))

u_df <- as.data.frame(all_cor[,,1]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("U" = "value")

mcmc_trace(u_df)

v_df <- as.data.frame(all_cor[,,2]) %>% pivot_longer(., everything(), names_to = "Chain") %>% 
  mutate(Chain = gsub("Chain", "", Chain)) %>% 
  mutate(Chain = as.numeric(ifelse(Chain == "", 4, Chain))) %>%
  rename("V" = "value")
mcmc_trace(v_df)

#plot(all_cor[,1], type = "l")




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

#-------------------------- MISCELLANEOUS EXTRA PLOTS ----------------------------#

#### Density plots for CV preds


M1Colours <- colorRampPalette(c("lightblue", "blue"))(10)
M2Colours <- colorRampPalette(c("tomato", "red"))(10)
M3Colours <- colorRampPalette(c("lightgreen", "green"))(10)
names(M1Colours) <- names(M2Colours) <-  names(M3Colours) <- levels(preds_combined$replicate)

preds_combined$replicate <- as.factor(preds_combined$replicate)

# First plot for Expert defined scenario
scenario <- "Expert-defined occurrence probabilities"
ggplot() + 
  # New Blocked
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "New, blocked" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M1", values = M1Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "blue", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "New, blocked" & 
                                                                                                  mean_density_combined$scenario == scenario,]) + 
  # New Unblocked
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "New, unblocked" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M2", values = M2Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "red", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "New, unblocked" & 
                                                                                                 mean_density_combined$scenario == scenario,]) + 
  # Old
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "Old" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M3", values = M3Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "green", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "Old" & 
                                                                                                   mean_density_combined$scenario == scenario,]) + 
  theme_minimal() +
  theme(legend.position = "none")


# Now plot for default scenario
scenario <- "75/100% occurrence probabilities"
ggplot() + 
  # New Blocked
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "New, blocked" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M1", values = M1Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "blue", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "New, blocked" & 
                                                                                                  mean_density_combined$scenario == scenario,]) + 
  # New Unblocked
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "New, unblocked" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M2", values = M2Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "red", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "New, unblocked" & 
                                                                                                 mean_density_combined$scenario == scenario,]) + 
  # Old
  stat_density(aes(x = value, color = as.factor(replicate)), geom = 'line', alpha = 0.25, 
               data = preds_combined[preds_combined$sampler == "Old" & 
                                       preds_combined$scenario == scenario,], position = 'identity') +
  scale_color_manual(name = "M3", values = M3Colours) + 
  new_scale_color() +
  geom_line(aes(x = z, y = zavg), color = "green", linewidth = 1.5, data = mean_density_combined[mean_density_combined$sampler == "Old" & 
                                                                                                   mean_density_combined$scenario == scenario,]) + 
  theme_minimal() +
  theme(legend.position = "none")


# Combined plot: similar performance across models w.r.t recall but expert scenarios give much better recall 
recall_df %>% filter(scenario != "0/100% occurrence probabilities") %>% # these are too far off they make the plot harder to read
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ scenario + stat) + 
  ggtitle('Out of sample performance: recall') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))




# Combined plot, only using proper, sampled recall
recall_df %>% filter(stat == 'Recall sampled' & scenario!="0/100% occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ scenario) + 
  ggtitle('Out of sample performance: recall') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))

# Expert scenario: slight advantage to old at 50%, 75%
recall_df %>%   filter(scenario == "Expert-defined occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ stat) + 
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))


# Main text: overview comparing 0/1 to new and old models under expert occurrence scenario: 
# Note there is a slight advantage to old at 50%, 75%
recall_df %>%   filter(scenario == "Expert-defined occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ stat) + 
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))



# New default scenario: slight advantage to old sampler
recall_df %>%   filter(scenario == "75/100% occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('') +
  xlab('') +
  facet_wrap(~ stat) + 
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))


### Pseudo-precision: the ratio of means

# Combined plot: similar peforamance with the 75/100 scenario, but strong advantage with expert scenarios; 75/100 gives better precision, 
# but after a certain point, does it really matter
ggplot(data = precision_df) +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  facet_wrap(~ scenario) + 
  ggtitle('Out of sample performance: precision') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))

# Expert-defined occurrence scenarios - only
precision_df %>% 
  filter(scenario == "Expert-defined occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))

# New default occurrence scenarios 75/100 - for supplement
precision_df %>% 
  filter(scenario == "75/100% occurrence probabilities") %>%
  ggplot() +
  geom_boxplot(aes(x = sampler, y = value)) +
  theme_bw() +
  ylab('Prediction Mean / Overall Mean') +
  xlab('') +
  theme(legend.position = 'none') +
  #scale_y_continuous(limits = function(x) c(0.9, x[2]), n.breaks = 6) 
  theme(text = element_text(size = 20, family = 'serif'), axis.text.x = element_text(angle = 90, hjust = 1))


