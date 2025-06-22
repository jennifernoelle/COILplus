

library(ape);library(tidyverse);library(U.PhyloMaker)

megatree_mammal <- read.tree("mammal_megatree.tre")#mammal megatree from VertLife
megatree_bird <- read.tree("bird_megatree.tre")#bird megatree from VertLife
frug_dat_total2 <- read.csv("frug_dat_total2.csv")#read combined frugivory database

megatree_frug <- bind.tree(megatree_mammal, megatree_bird, where = "root", position = 0, interactive = FALSE)#bind mammal and bird consensus megatrees into one combined megatree

frugs_taxo <- frug_dat_total2 %>% select(Animal_Species_Corrected, Animal_Family) %>% mutate(genus = Animal_Species_Corrected) %>% separate(genus, into = c("genus", "sp")) %>%  select(Animal_Species_Corrected, genus, Animal_Family) %>% mutate(species = str_replace(Animal_Species_Corrected, "_", " ")) %>%distinct() %>% rename(family= Animal_Family)%>% select(species, genus, family) #get U.PhyloMaker "sp.list"

frugs_gn <- frugs_taxo %>% select(genus, family)#get U.PhyloMaker "gen.list"

frug_tree <- U.PhyloMaker::phylo.maker(frugs_taxo, megatree_frug, frugs_gn, nodes.type = 1, scenario = 3)#get phylogenetic tree of species in frugivore database
#four species fail to bind to the tree

frug_tree_phylo2 <- frug_tree$phylo #get "phylo" class tree from the U.PhyloMaker output
