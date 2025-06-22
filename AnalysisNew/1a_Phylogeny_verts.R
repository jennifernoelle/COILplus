
# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawDataNew/'
save_path <- 'ProcessedDataNew/' # Create this directory

library(tidyverse)
save_files <- TRUE

# --------- BEGINNING -------- #

# Setting the working directory.
#setwd(wd_path)

# Loading libraries.
library(ape)
library(tidyverse)
library(data.table)
library(superheat)

if(!("U.PhyloMaker" %in% rownames(installed.packages()))){
  install.packages("remotes")
  devtools::install_github("jinyizju/U.PhyloMaker") 
}
library(U.PhyloMaker)

# ---------------------------  PART 1: Load the data --------------------------#

# Read the cleaned taxonomical data corresponding to our data
v.taxa <- read.csv(paste0(save_path, "v_taxa.csv")) 
v.names <- v.taxa$Animal_Species_Corrected

# Read the downloaded tree data from VertLife for birds, mammals and create mega-tree
megatree_mammal <- read.tree(paste0(data_path, "mammal_megatree.tre")) # mammal megatree from VertLife
megatree_bird <- read.tree(paste0(data_path, "bird_megatree.tre")) # bird megatree from VertLife
megatree_frug <- bind.tree(megatree_mammal, megatree_bird, where = "root", #bind mammal and bird consensus megatrees into one combined megatree
                           position = 0, interactive = FALSE)

# Replace species names for eight species which are not included at the species level and tree placement seems suspect
# We're replacing them with species in the same genus for phylo, then will rename the correlation matrix
v.taxa <- v.taxa %>% 
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Myonycteris_leptodon","Myonycteris_angolensis" ,Animal_Species_Corrected )) %>% 
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Epomophorus_minor" ,"Epomops_franqueti", Animal_Species_Corrected2 )) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Scotonycteris_ophiodon" ,"Scotonycteris_bergmansi", Animal_Species_Corrected2)) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Galagoides_demidoff" ,"Galagoides_granti",  Animal_Species_Corrected2)) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Piliocolobus_tephrosceles", "Procolobus_verus", Animal_Species_Corrected2)) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Piliocolobus_badius","Procolobus_verus" , Animal_Species_Corrected2)) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected =="Piliocolobus_tholloni", "Procolobus_verus", Animal_Species_Corrected2)) %>%
  mutate(Animal_Species_Corrected2 = ifelse(Animal_Species_Corrected == "Grammomys_poensis",  "Grammomys_caniceps", Animal_Species_Corrected2)) 


# Prepare the taxonomical data
frugs_sp <- v.taxa %>% rename(species = Animal_Species_Corrected2, genus = Animal_Genus, family = Animal_Family) %>% 
                       mutate(species = gsub("_", " ", species)) %>% 
                       select(species, genus, family) #get U.PhyloMaker "sp.list"
frugs_gn <- v.taxa %>% rename(genus = Animal_Genus, family = Animal_Family) %>% 
                       select(genus, family) #get U.PhyloMaker "gen.list"

#---------------------------- PART 2: Create phylogenetic hypothesis----------#

# Get phylogenetic tree of species in frugivore database
frug_tree <- U.PhyloMaker::phylo.maker(sp.list = frugs_sp$species, 
                                       tree = megatree_frug, gen.list = frugs_gn, 
                                       nodes.type = 1, scenario = 3)
# Extract the tree
frug_tree_phylo <- frug_tree$phylo

# Checking how many species from our interaction data exist in the phylogenetic trees: all!
v.names2 <- gsub(" ", "_", frugs_sp$species)
x_names <- frug_tree_phylo$tip.label
sort(v.names2[which(!(v.names2 %in% x_names))])
sort(x_names[which(!(x_names %in% v.names2))])

# Using the phylogenetic tree and the ape R package to get a phylogenetic
# correlation matrix.
Cu_phylo <- round(ape::vcv(frug_tree_phylo, corr = TRUE), 10)
sum(Cu_phylo >1)/ length(Cu_phylo)
heatmap(Cu_phylo)

#### Rename the species in the phylogenetic correlation matrix
# Note that in some cases we actually have the place holder species in our collection, so we need both

## Non-duplicate species: we can just rename these
x_names <- ifelse(x_names == "Scotonycteris_bergmansi", "Scotonycteris_ophiodon", x_names)
x_names <- ifelse(x_names == "Galagoides_granti","Galagoides_demidoff", x_names) 
x_names <- ifelse(x_names == "Grammomys_caniceps", "Grammomys_poensis", x_names)
rownames(Cu_phylo) <- colnames(Cu_phylo) <- x_names

# Duplicate species: Epomops_franqueti, Myonycteris_angolensis,Procolobus_verus
# We have to create a new row
w.Ef <- which(row.names(Cu_phylo) == "Epomops_franqueti")
w.Ma <- which(row.names(Cu_phylo) == "Myonycteris_angolensis")
w.Pv <- which(row.names(Cu_phylo) == "Procolobus_verus")

Cu.Ef <- Cu_phylo[, w.Ef]
Cu_phylo <- rbind(Cu.Ef, Cu_phylo)
Cu_phylo <- cbind(c(1, Cu.Ef), Cu_phylo) # Self-corr = 1
rownames(Cu_phylo)[1] <- colnames(Cu_phylo)[1] <- "Epomophorus_minor"

Cu.Ma <- Cu_phylo[, w.Ma]
Cu_phylo <- rbind(Cu.Ma, Cu_phylo)
Cu_phylo <- cbind(c(1, Cu.Ma), Cu_phylo) # Self-corr = 1
rownames(Cu_phylo)[1] <- colnames(Cu_phylo)[1] <- "Myonycteris_leptodon"

Cu.Pv <- Cu_phylo[, w.Pv]
Cu_phylo <- rbind(Cu.Pv, rbind(Cu.Pv, rbind(Cu.Pv, Cu_phylo))) # Three species replaced with this
Cu_phylo <- cbind(c(1,1,1, Cu.Pv), cbind(c(1,1,1, Cu.Pv), cbind(c(1,1,1, Cu.Pv), Cu_phylo))) # Self-corr = 1
rownames(Cu_phylo)[1:3] <- colnames(Cu_phylo)[1:3] <- c("Piliocolobus_tephrosceles","Piliocolobus_badius","Piliocolobus_tholloni" )

# Re-order the rows and columns to match the order of species in our data:
use_order <- sapply(v.names, function(r) which(rownames(Cu_phylo) == r))
table(sapply(use_order, length))
use_order <- use_order[sapply(use_order, length) > 0]
use_order <- as.numeric(use_order)
Cu_phylo <- Cu_phylo[use_order, use_order]


heatmap(Cu_phylo)
# # Abbreviate names for plotting
# plot.names <- gsub("_", " ", colnames(Cu_phylo))
# g <- substr(gsub(" .*", "", plot.names), 1,3)
# sp <- gsub(".* " , "", plot.names)
# plot.names <- paste0(g, ". ", sp)
# plot.cu <- Cu_phylo
# colnames(plot.cu) <- rownames(plot.cu) <- plot.names
# 
#png(filename = "EDA/Phylo_mammals.png", width = 1000, height =1000)
superheat(X = Cu_phylo,
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
          title = "Phylogenetic correlation of vertebrates",
          #left.label.text.size = 20,
          #bottom.label.text.size = 20,
          title.size = 15)
# dev.off()


if(save_files){
  save(Cu_phylo, file = paste0(save_path, 'Cu_phylo.dat'))
}
