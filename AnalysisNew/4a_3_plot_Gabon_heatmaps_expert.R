# File to plot results for cv and full run

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
# wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
# setwd(wd_path)


# Where the raw data are saved:
data_path_raw <- 'RawDataNew/'
# Where the processed data are saved:
data_path <- 'ProcessedDataNew/'
# Where the functions are available:
source_path <- 'HelperScriptsJKNew/'
save_path <- 'ResultsNew/Gabon/'


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

# Supplementary data on species:
v.taxa <- read.csv(paste0(data_path, 'v_taxa.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa.csv'))

# Species list for Gabon
v.gabon <- read.csv(paste0(data_path_raw, 'frugivores_gabon.csv'))[,2]
p.gabon <- read.csv(paste0(data_path_raw, 'plants_gabon.csv'))[,2]
p.gabon <- paste(str_split_i(p.gabon, pattern = " ", i = 1), 
                 str_split_i(p.gabon, pattern = " ", i = 2), sep = "_")

# Rename for convenience
Cu <- Cu_phylo
Cv <- Cv_phylo
obs_W <- Obs_W
obs_X <- Obs_X
obs_A <- A_obs
obs_F <- F_obs

# Check for consistency
sum(v.taxa$Animal_Species_Corrected != rownames(obs_A))
sum(p.taxa$Plant_Species_Corrected != colnames(obs_A))

sum(!(v.gabon %in% rownames(obs_A)))
sum(!(p.gabon %in% colnames(obs_A)))

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

## Improved default guess: 0.75/1
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


# -------------- STEP 1: Specifications. ------------ #

# Number of chains for the full sample runs
nchains <- 4

# Number of cross validation repetitions:
repetitions <- 10

# Number of cv samples per repetition
n.cv <- 100

# --------------- STEP 1: Getting the results together - full sample  ----------------- #

# Read the posterior networks for each of the dates we're interested in

# Where you want to save MCMC results:
# Save results using convention: res_date_i.rda

date <-  'pd_B4_nb_400sims' # Expert prior using unblocked sampling of occurrence indicators and probabilities

# date_list <- list('p0_A4_old_400sims','p75_A4_old_400sims', 'p75_B4_nb_400sims',
#                   'p75_C4_b_400sims','pd_A4_old_400sims','pd_B4_nb_400sims',
#                   'pd_C4_b_400sims')
# 
# date_desc_list <- list('1/0 Prior, Old Model', '1/0.75 Prior, Old Model', 
#                        '1/0.75 Prior, New Model Blocked', '1/0.75 Prior, New Model Unblocked', 
#                        'Expert Prior, Old Model', 'Expert Prior, New Model Blocked', 
#                        'Expert Prior, New Model Unblocked')
# 
# 
# ### Read all posterior networks into a list
# net_list <- list()
# 
# for (ii in 1 : length(date_list)) {
#   this.version <- date_list[[ii]]
#   cat("\n Version = ", this.version)
#   results_path <- paste0('ResultsNew/', this.version , '/')
#   net_list[[ii]] <- read.csv(paste0(results_path, 'post_network',  this.version, '.csv')) 
# }

results_path <- paste0('ResultsNew/', date , '/')
this_net <- read.csv(paste0(results_path, 'post_network',  date, '.csv'))

#---------------------- Subset to the Gabon species --------------------------#

# Subset observed data for Gabon, but keep ordering
keep_vs <- which(rownames(obs_A) %in% v.gabon)
keep_ps <- which(colnames(obs_A) %in% p.gabon)

v.taxa.g <- v.taxa[keep_vs,]
p.taxa.g <- p.taxa[keep_ps,]

v.g.names <- rownames(obs_A)[keep_vs]
p.g.names <- colnames(obs_A)[keep_ps]
nV.g <- length(v.g.names)
nP.g <- length(p.g.names)

this_net <- this_net_na <- this_net[keep_vs, keep_ps]
# net_list_g <- lapply(net_list, function(x) x[keep_vs, keep_ps])
comb_A_subset <- comb_A[keep_vs, keep_ps]

this_net_na[comb_A_subset ==1] <- NA # set observed interactions to NA



# ----------- STEP 2: PLOTTING THE HEATMAP --------------------#

# New code to reorder the vertebrates and keep only the families of interest


## Heatmap setup

# First, we want the following categories: ungulate, monkey, elephant, carnivore, bird, bat, ape
v.taxa.g <- v.taxa.g %>%
            mutate(animal_type = ifelse(Animal_Taxa_Type == "Large_birds" | Animal_Taxa_Type == "Small_birds", "Bird", Animal_Taxa_Type)) %>%
            mutate(animal_type = gsub("s$", "", animal_type))

v_group <- v.taxa.g$animal_type
p_group <- p.taxa.g$Plant_Family
v.taxa.order <- v.taxa.ordered$original_row

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_v_size to 10, and min_plant_size to 25. Additionally, 
# we retain the smaller clusters Hominidae and Elephantidae because they are of
# scientific interest. Setting both to 0 will produce the full results.
min_v_size <- 1
min_p_size <- 5

keep_v_groups <- c("Ape", "Bat", "Bird", "Carnivore", "Elephant", "Monkey", "Ungulate")
keep_p_groups <- names(which(p_size_cluster >= min_p_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


superheat(X = this_net_na[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 4,
          bottom.label.text.size = 4,
          bottom.label.size = 0.2, left.label.size = 0.12,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05), 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 10)


png(paste0(save_path, "heatmap_gabon_goodtaxa_", date, ".png"), 
    width = 1500, height =1000)
superheat(X = this_net_na[keep_v_index, keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.size = 0.75, grid.vline.size = 0.75,
          left.label.text.size =10,
          bottom.label.text.size = 10,
          bottom.label.size = 0.35, left.label.size = 0.2,
          bottom.label.text.angle = 90,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          legend.width = 2.5,
          legend.height = 0.15,
          legend.text.size = 35,
          heat.col.scheme = "grey", heat.na.col = 'black',
          heat.pal.values = seq(0, 1, by = 0.05)
          #title = date
)
dev.off()

# 
# # Creating the clusters that will be used
# # The following two lines specify that horizontal and vertical lines in our
# # plot will separate species by taxonomic families:
# v_group <- v.taxa.g$Animal_Family
# p_group <- p.taxa.g$Plant_Family
# 
# # Calculating the size of each cluster, will be used when plotting results
# # for families of certain size:
# v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
# p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))
# 
# # Set the minimum cluster size that should be plotted. For the results of the
# # manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# # to 0 will produce the full results.
# min_v_size <- 5
# min_p_size <- 5
# 
# keep_p_groups <- names(which(p_size_cluster >= min_p_size))
# keep_v_groups <- names(which(v_size_cluster >= min_v_size))
# 
# keep_v_index <- which(v_group %in% keep_v_groups)
# keep_p_index <- which(p_group %in% keep_p_groups)


for(p in 1:length(net_list_g_na)){
  this.version <- date_list[[p]]
  this.description <- date_desc_list[[p]]
  plot_pred <- net_list_g_na[[p]]
  
  # Plotting those with minimum size as specified:
  # We replaced observed interactions with black NA
  png(paste0(save_path, "heatmap_gabon", this.version, ".png"), 
      width = 1000, height =1000)
  superheat(X = plot_pred[keep_v_index, keep_p_index],
            membership.rows = v_group[keep_v_index],
            membership.cols = p_group[keep_p_index],
            grid.hline.col = "#00257D", grid.vline.col = '#00257D',
            grid.hline.size = 0.3, grid.vline.size = 0.3,
            bottom.label.text.angle = 90,
            left.label.text.size = 4,
            bottom.label.text.size = 4,
            bottom.label.size = 0.2, left.label.size = 0.12,
            #legend.breaks = seq(0, 1, by = 0.2),
            #legend.vspace = 0.05,
            heat.col.scheme = "grey", heat.na.col = 'black',
            heat.pal.values = seq(0, 1, by = 0.05), 
            title = paste0("Post Probs: ", this.description), 
            #left.label.text.size = 20, 
            #bottom.label.text.size = 20, 
            title.size = 10)
  dev.off()
}


# Plot observed network no plant legend
#comb_A_na <- ifelse(comb_A ==1, NA, 0)
png(paste0(save_path, "heatmap_gabon_observed.png"), width = 1000, height = 510)
superheat(X = comb_A[keep_v_index, keep_p_index],
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
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"),
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          title.size = 10, 
          title = "Observed",
          #title.alignment = "center",
          legend = TRUE
)
dev.off()


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



# ----------- STEP 2b: PLOTTING THE HEATMAP for printing --------------------#

# Creating the clusters that will be used
# The following two lines specify that horizontal and vertical lines in our
# plot will separate species by taxonomic families:
v_group <- v.taxa.g$Animal_Family
p_group <- p.taxa.g$Plant_Family

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size:
v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_bird_size to 10, and min_plant_size to 20. Setting both
# to 0 will produce the full results.
min_v_size <- 5
min_p_size <- 5

keep_p_groups <- names(which(p_size_cluster >= min_p_size))
keep_v_groups <- names(which(v_size_cluster >= min_v_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)


p <- 6
  this.version <- date_list[[p]]
  this.description <- date_desc_list[[p]]
  plot_pred <- net_list_g_na[[p]]
  
  # Plotting those with minimum size as specified:
  # We replaced observed interactions with black NA
  png(paste0(save_path, "heatmap_gabon", this.version, ".png"), 
      width = 1000, height =1000)
  superheat(X = plot_pred[keep_v_index, keep_p_index],
            membership.rows = v_group[keep_v_index],
            membership.cols = p_group[keep_p_index],
            grid.hline.col = "#00257D", grid.vline.col = '#00257D',
            grid.hline.size = 0.3, grid.vline.size = 0.3,
            bottom.label.text.angle = 90,
            left.label.text.size = 6,
            bottom.label.text.size = 6,
            bottom.label.size = 0.2, left.label.size = 0.15,
            #legend.breaks = seq(0, 1, by = 0.2),
            #legend.vspace = 0.05,
            heat.col.scheme = "grey", heat.na.col = 'black',
            heat.pal.values = seq(0, 1, by = 0.05), 
            title = "",
            legend.text.size = 15,
            #title = paste0("Post Probs: ", this.description), 
            #left.label.text.size = 20, 
            #bottom.label.text.size = 20, 
            title.size = 10)
  dev.off()



# Plot observed network no plant legend
#comb_A_na <- ifelse(comb_A ==1, NA, 0)
png(paste0(save_path, "heatmap_gabon_observed.png"), width = 1000, height = 510)
superheat(X = comb_A[keep_v_index, keep_p_index],
          membership.rows = v_group[keep_v_index],
          membership.cols = p_group[keep_p_index],
          grid.hline.col = "#00257D", grid.vline.col = '#00257D',
          grid.hline.size = 0.3, grid.vline.size = 0.3,
          bottom.label.text.angle = 90,
          left.label.text.size = 4,
          bottom.label.text.size = 4,
          bottom.label.size = 0.2, left.label.size = 0.12,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          #heat.col.scheme = "grey", 
          #heat.na.col = 'black',
          heat.pal = c("white", "black"),
          heat.pal.values = seq(0, 1, by = 0.05), 
          #title = "(a) Observed frugivory network", 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20, 
          #title.size = 10, 
          #title = "Observed",
          title = "",
          #title.alignment = "center",
          legend = FALSE
)
dev.off()


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

