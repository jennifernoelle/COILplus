# Builds plant phylogenies and also adds additional taxonomical info to traits files
# Uses taxonomical info to impute missing fruit lengths

# -------- TO DO --------- #

# The directory where the original data are:
data_path <- 'RawDataNew/'
# The directory where the processed data should be saved:
save_path <- 'ProcessedDataNew/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# --------- BEGINNING -------- #

# Setting the working directory.
#setwd(wd_path)
if(!("V.PhyloMaker" %in% rownames(installed.packages()))){
  install.packages("remotes")
  remotes::install_github("jinyizju/V.PhyloMaker")
}

# Loading libraries.
library(ape)
library(V.PhyloMaker)
library(gplots)
library(magrittr)
library(dplyr)


# Read taxonomical data
p.taxa <- read.csv(paste0(save_path, "p_taxa.csv")) %>%
          rename(Species = Plant_Species_Corrected, Genus = Plant_Genus, Family = Plant_Family)


n_occur <- data.frame(table(p.taxa$Species))
n_occur[n_occur$Freq > 1, ]
p.taxa[p.taxa$Species %in% n_occur[n_occur$Freq > 1, 1],]

nP <- nrow(p.taxa)

#------------------- GET PHYLOGENY ----------------------#

y <- phylo.maker(p.taxa, scenarios = "S3") # default scenario

#create a random tree to test later
y_rand <- rtree(nP)

# Check to see which species are in the tree
tree.species <- y$species.list
filter(tree.species, status != "bind")
p.taxa$Species[!(p.taxa$Species %in% tree.species$species)]

# Using the phylogenetic tree in the ape R package to get a phylogenetic
# correlation matrix.
Cv_phylo <- round(ape::vcv(y[[1]], corr = TRUE), 10)
p.taxa$Species[!(p.taxa$Species %in% rownames(Cv_phylo))]
rownames(Cv_phylo)[!(rownames(Cv_phylo) %in% p.taxa$Species)]


Cv_phylo1 <- round(ape::vcv(y_rand, corr = TRUE), 10)


use_order <- unlist(sapply(p.taxa$Species,function(r) which(rownames(Cv_phylo) == r)))
Cv_phylo <- Cv_phylo[use_order, use_order]


# Most plants are not very closely related
gplots::heatmap.2(Cv_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE,
                  labRow = FALSE, labCol = FALSE)
hist(c(Cv_phylo))

if (save_files) {
  save(Cv_phylo, file = paste0(save_path, 'Cv_phylo.dat'))
}



