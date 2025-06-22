# Builds plant phylogenies and also adds additional taxonomical info to traits files
# Uses taxonomic info to impute missing fruit lengths

# ----------------------------------- TO DO ---------------------------------- #

# The directory where the original data are:
data_path <- 'RawDataPub/'
# The directory where the processed data should be saved:
save_path <- 'ProcessedDataPub/'
save_path_generic <- 'RawDataGeneric/'
# Whether the processed data should be saved or not:
save_files <- TRUE


# ------------------------------- STEP O. LOAD PACKAGES ---------------------- #

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


# Read taxonomic data
p.taxa <- read.csv(paste0(save_path, "p_taxa.csv")) %>%
          rename(Species = Plant_Species_Corrected, Genus = Plant_Genus, Family = Plant_Family)


n_occur <- data.frame(table(p.taxa$Species))
n_occur[n_occur$Freq > 1, ]
p.taxa[p.taxa$Species %in% n_occur[n_occur$Freq > 1, 1],]

nP <- nrow(p.taxa)

#------------------------------- STEP 1. GET PHYLOGENY ------------------------#

y <- phylo.maker(p.taxa, scenarios = "S3") # default scenario

# Create a random tree to test later
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

gplots::heatmap.2(Cv_phylo, dendrogram = 'none', trace = 'none', Rowv = FALSE, Colv = FALSE,
                  labRow = FALSE, labCol = FALSE)
hist(c(Cv_phylo))



# Create generic version of the matrix to use with anonymized data
codebook <- read.csv(paste0(save_path, 'p_taxa_codebook.csv')) %>% 
  select(Plant_Species_Corrected, Plant_Species_Generic)
sum(rownames(Cv_phylo) != codebook$Plant_Species_Corrected)
Cv_phylo_generic <- Cv_phylo
rownames(Cv_phylo_generic) <- colnames(Cv_phylo_generic) <- codebook$Plant_Species_Generic

if (save_files) {
  save(Cv_phylo, file = paste0(save_path, 'Cv_phylo.dat'))
  save(Cv_phylo_generic, file = paste0(save_path_generic, 'Cv_phylo.dat')) 
}



