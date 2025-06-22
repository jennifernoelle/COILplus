# Estimating variable importance based on interaction model
# This code is just for plotting

# --------- TO DO: set  your directories and name the current results files using the date ---------#

# The directory where the analysis is performed: should already be your WD if you cloned the repo
#wd_path <- '/Users/camilledesisto/Documents/GitHub/African-Frugivory'
#wd_path <- "/home/grad/jnk21/projects/African-Frugivory-V2"
# wd_path<- "/hpc/group/dunsonlab/jkampe/African-Frugivory-V3"
# setwd(wd_path)

# Save results using convention: res_date_i.rda
date <- 'pd_C4_b_400sims'

# Where the processed data are saved:
data_path <- 'ProcessedDataNew/'
# Where you put the MCMC results:
results_path <- paste0('ResultsNew/', date, '/')
# Where you want to save the trait matching
save_path <- paste0('ResultsNew/', date, '/')
# Where the functions are available:
source_path <- 'HelperScriptsJKNew/'

# ------ STEP 0: Some functions. --------- #

library(tidyverse)
library(grid)
library(gridExtra)
library(abind)
library(superheat)

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

# Read in taxonomical data
v.taxa <- read.csv(paste0(data_path, 'v_taxa.csv'))
p.taxa <- read.csv(paste0(data_path, 'p_taxa.csv'))


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

# all_res <- <-NULL
all_Ls <- all_X <- all_W <-NULL
for (ii in 1 : nchains) {
  cat("\n Chain = ", ii)
  res <- loadRData(paste0(results_path, 'res_',  date, '_', ii, '.dat')) 
  cat(" Dim is ", dim(res[[1]]))
  #all_res[[ii]] <- res
  #all_pred[[ii]] <- res$all_pred
  all_Ls[[ii]] <- res$all_pred[,,,1] # selects posterior samples of L
  all_X[[ii]] <-res$Xs
  all_W[[ii]] <- res$Ws
}

# all_pred <- all_X <- all_W <- NULL
# for (ii in 1 : nchains) {
#   load(paste0(results_path, 'res_', date, '_', ii, '.dat'))
#   all_pred[[ii]] <- res$all_pred
#   all_X[[ii]] <-res$Xs
#   all_W[[ii]] <- res$Ws
# }

# Number of posterior samples by chain:
Nsims <- dim(all_Ls[[1]])[1]
pV <- length(all_X[[1]])
pP <- length(all_W[[1]])

# # Using the linear predictor of the interaction model (not corrected for bias due to detection, number of studies):
# # Used for performing the trait matching analysis but not plotting
# mod_pL1s <- array(NA, dim = c(nchains * Nsims, nV, nP))
# for (cc in 1 : nchains) {
#   wh_entries <- Nsims * (cc - 1) + 1 : Nsims
#   mod_pL1s[wh_entries, , ] <- all_pred[[cc]][, , , 2] # Selects prob_L
# }

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


# --------------- STEP 2: Load and plot the trait matching ----------------- #


# Load the trait matching results
load(paste0(save_path,  'rsq_obs_X_', date, '.dat'))
load(paste0(save_path, 'rsq_obs_W_', date, '.dat'))
load(paste0(save_path, 'rsq_resampling_X_', date, '.dat'))
load(paste0(save_path, 'rsq_resampling_W_', date, '.dat'))
load(paste0(save_path, 'corr_obs_X_', date, '.dat'))
load(paste0(save_path, 'corr_obs_W_', date, '.dat'))

good_namesX <- c("Generation Length", "Log Body Mass", "IUCN Status", "Habitat")
good_namesW <- c("Fruit Length", "Fruit Width", "Wood Density")

# Calculating the number of permuted standard deviations away from the mean.

# Starting from the vert covariates:
wh_obs <- rsq_obs_X
wh_resampling <- rsq_resampling_X
sd_awayX <- rep(NA,  length(wh_obs))
names(sd_awayX) <- good_namesX #colnames(obs_X) #good_namesX
for  (cc in 1 : length(wh_obs)) {
  sd_awayX[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc])
}

# And for the plant covariates:
wh_obs <- rsq_obs_W # this gives the squared correlation coefficient between covars and interactions
wh_resampling <- rsq_resampling_W # this is the squared correlation coeff between shuffled covars and interactions
sd_awayW <- rep(NA,  length(wh_obs))
names(sd_awayW) <- good_namesW
for  (cc in 1 : length(wh_obs)) { # for each covar, compute how many sds away from resampled mean the obs rsq is
  sd_awayW[cc] <- (wh_obs[cc] - mean(wh_resampling[, cc])) / sd(wh_resampling[, cc]) # note these sds are tiny
}


# Plotting the tiles of variable importance ordering the variables in
# decreasing importance:

# For the vertebrate species:
xx <- data.frame(value = sd_awayX, covariate = names(sd_awayX), y = 1)
xx <- xx[order(- xx$value), ]
xx$covariate <- factor(xx$covariate, levels = xx$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif"))

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = xx) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'bottom', axis.title = element_blank()) +
  guides(fill = guide_legend("Frugivore Variable Importance")) +
  scale_fill_gradient(low = '#BFF0B6', high = '#3B6E32') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif", color = "black"), 
        legend.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 15, family = "serif", color = "black"), 
        legend.title = element_text(size = 20, family = "serif"))


# For the plant species:
ww <- data.frame(value = sd_awayW, covariate = names(sd_awayW), y = 1)
ww <- ww[order(- ww$value), ]
ww$covariate <- factor(ww$covariate, levels = ww$covariate)

ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'none', axis.title = element_blank()) +
  scale_fill_gradient(low = '#D4DEF7', high = '#425075') +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif"))


ggplot() + geom_tile(aes(x = covariate, y = y, fill = value), color = 'white', data = ww) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank()) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(legend.position = 'bottom', axis.title = element_blank()) +
  scale_fill_gradientn(colors = c('#D4DEF7','#425075'), breaks = c(600, 1000, 1650)) +
  guides(fill = guide_legend("Plant Variable Importance")) +
  theme(axis.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 20, family = "serif", color = "black"), 
        legend.text = element_text(angle = 0, hjust = 0.5, vjust = 0, size = 15, family = "serif", color = "black"), 
        legend.title = element_text(size = 20, family = "serif"))


# ------- PART B: Posterior probabilities based on the important covariates.

# Taking the posterior probabilities of interaction by averaging R across the three chains,
# and setting the recorded interactions to NA:
# use_Ls <- do.call(abind, c(lapply(all_res, function(x) x$all_pred[, , , 1]), along = 1)) # all_pred[,,,1] is pred_L
use_Ls <- abind(all_Ls, along = 1) # all_pred[,,,1] is pred_L
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
  #scale_fill_gradient(low = "#F5D4C7", high = "#02A65F", na.value = '#016B3B',
  #                    name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  scale_fill_gradient(low = "#FDE725FF", high = "#02A65F", na.value = '#016B3B',
                      name = 'Posterior\ninteraction\nprobability\n', limits = c(0, 1)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 30),
        axis.title.y = element_text(vjust = - 1, family = 'serif'), text = element_text(family = 'serif')) +
  ylab(expression(symbol('\256'))) + xlab(expression(symbol('\256')))

gridExtra::grid.arrange(g, left = textGrob("Frugivores: increasing log body mass", rot = 90,
                                           x = 1.3, y = 0.57, gp = gpar(fontsize = 12, fontfamily = 'serif')),
                        bottom = textGrob("Plants: increasing wood density",
                                          x = 0.435, y = 1.3, gp = gpar(fontsize = 12, fontfamily = 'serif')),
                        vp=viewport(width=0.5, height=0.6))


#------ PART C. PLOTTING IMPORTANT COVARIATES WITH SIGNED CORRELATION ---------#

# TO DO HERE: I came up with a workaround to get the families in the right order
# Next step is to save the images, and edit the labels in BioRendr, adding graphics for animal families
# Restrict to fewer families
# Save images, open in photos, paste into paint

# Now group taxa by group for these plots and create pretty names
colnames(corr_obs_W) <- c("Fruit length", "Fruit width", "Wood density") #colnames(obs_W) # plant covars

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
plot_pred_W <- corr_obs_W[v.taxa.order,]

# Calculating the size of each cluster, will be used when plotting results
# for families of certain size


v_size_cluster <- sapply(unique(v_group), function(x) sum(v_group == x))
v_size_df <- data.frame(v_size_cluster) %>% rownames_to_column(., var = "family.code")
v.family.codebook <- left_join(v.family.codebook, v_size_df, join_by(v.show.order == family.code))
# Look at which families might be interesting
v.family.codebook %>% select(Animal_Family, Animal_Taxa_Type_Courser, v_size_cluster) %>% unique()

p_size_cluster <- sapply(unique(p_group), function(x) sum(p_group == x))

# Set the minimum cluster size that should be plotted. For the results of the
# manuscript, we set min_v_size to 10, and min_plant_size to 25. Additionally, 
# we retain the smaller clusters Hominidae and Elephantidae because they are of
# scientific interest. Setting both to 0 will produce the full results.
table(v_size_cluster)
table(p_size_cluster)
min_v_size <- 12
min_p_size <- 30

which_v_groups_extra <- unique(v.family.codebook[which(v.family.codebook$Animal_Family == "Elephantidae" 
                                | v.family.codebook$Animal_Family == "Hominidae"), "v.show.order"])
keep_v_groups <- sort(c(names(which(v_size_cluster >= min_v_size)), which_v_groups_extra))
keep_p_groups <- names(which(p_size_cluster >= min_p_size))

keep_v_index <- which(v_group %in% keep_v_groups)
keep_p_index <- which(p_group %in% keep_p_groups)

### Plant covariates versus animal groups

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(save_path, "corr_w_heatmap_", date, ".png"), width = 1500, height =1000)

# Plot without labels for adding pretty silhouettes
superheat(X = plot_pred_W[keep_v_index, ],
          pretty.order.cols = FALSE, 
          pretty.order.rows = FALSE,
          membership.rows = v_group[keep_v_index],
          #membership.cols = p_group[keep_p_index],
          grid.hline.col = "black", grid.vline.col = 'black', # #00257D
          grid.hline.size = 0.5, grid.vline.size = 0.5,
          #left.label.text.size = 10,
          #left.label = "none", bottom.label = "none", legend =  FALSE
          #bottom.label.text.size = 10,
          #bottom.label.size = 0.075, left.label.size = 0.2,
          #bottom.label.text.angle = 90,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          # legend.width = 2.5,
          # legend.height = 0.15,
          # legend.text.size = 25,
          # padding = 2
          #heat.col.scheme = "grey", heat.na.col = 'black',
          #heat.pal.values = seq(0, 1, by = 0.05)
          #title = date, 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20
)
#dev.off()

# Get family names for pretty plots
unique(v.family.codebook[v.family.codebook$v_size_cluster >= min_v_size | 
                         v.family.codebook$v.show.order %in% which_v_groups_extra, 1:2])



### Frugivore covariates versus plant groups
plot_pred_X <- corr_obs_X[, -c(3,4)]
colnames(plot_pred_X) <- c("Generation length", "Log body mass") #colnames(obs_X) # frugivore covars

# Plotting those with minimum size as specified:
# We replaced observed interactions with black NA
#png(paste0(save_path, "corr_x_heatmap_", date, ".png"), width = 1500, height =1000)
superheat(X = plot_pred_X[keep_p_index, ],
          pretty.order.cols = FALSE, 
          pretty.order.rows = FALSE,
          membership.rows = p_group[keep_p_index],
          left.label = "none", bottom.label = "none", legend =  FALSE,
          grid.hline.col = "black", grid.vline.col = 'black', # #00257D
          grid.hline.size = 0.5, grid.vline.size = 0.5,
          # left.label.text.size = 5,
          # bottom.label.text.size = 10,
          # bottom.label.size = 0.075, left.label.size = 0.2,
          #bottom.label.text.angle = 90,
          #legend.breaks = seq(0, 1, by = 0.2),
          #legend.vspace = 0.05,
          legend.width = 2.5,
          legend.height = 0.15,
          legend.text.size = 25,
          #heat.col.scheme = "grey", heat.na.col = 'black',
          #heat.pal.values = seq(0, 1, by = 0.05)
          #title = date, 
          #left.label.text.size = 20, 
          #bottom.label.text.size = 20
)
#dev.off()

