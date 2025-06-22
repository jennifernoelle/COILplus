# This file performs basic data cleaning on the new combined adjacency matrices
# This file processes the raw data to create the following key data sources: 
# Adj_matrix_new.rds, O_A, O_P
# Saves useful info about studies

# MISC to do don't forget to do before running: 
  # X Reorder adjacency matrices phylogenetically maybe just by family
  # Make sure to fix taxa in phylogeny and traits files too if needed

# MISC to do while running
  # Create a new EDA file to put all of this into
  # Save some of the processed data with study info from before as raw data to do the occurrence processing
  # Include info that will allow us to produce map with animal focus
  # Replace missing lat-long with country centroids for plotting?
  

# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawDataNew/'
save_path <- 'ProcessedDataNew/' # Create this directory

save_files <- FALSE

library(tidyverse)
library(foreach)
library(doParallel)
library(abind)

# Load data
frug  <- read.csv(paste0(data_path, "Frugivory_combined_clean.csv"))[, -1]
zones <- read.csv(paste0(data_path, "StudyZone_clean2.csv"))[, -c(8,9)] 


#----------------------------- PART A: EDA ------------------------------------#

### 1. Interactions data
apply(frug, 2, function(x) sum(is.na(x)))
table(frug$Zone) # 354 missing
table(frug$Country) # 668 missing
table(frug$Focus) # 293 missing & 268 both limited: treat these as v. narrow studies (i.e. all 0s noninformative)
table(frug$Study_Site) # 962 missing
table(frug$Study_ID)
nrow(frug) # 10635 interactions

frug_unique <- frug %>% group_by(Animal_Species_Corrected, Plant_Species_Corrected) %>% 
         dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected) %>% 
         unique() 
nrow(frug_unique) # 5938 unique interactions

### 2. Site metadata
z.sites <- sort(zones$StudyZone) # has lat long
f.sites <- sort(unique(frug$Study_Site))
f.sites[!(f.sites %in% z.sites)]

### 3. check for multi-site studies and multicountry studies5938
table(length(frug$Study_Site)) # Each observation has only one site associated to it
# 
# frug_studies <- frug %>% group_by(Study_ID) %>% # But many studies cover multiple sites # <--------- DEBUG: try restoring this line if anything breaks
#   summarize(n_SitesThisStudy = length(unique(Study_Site)), 
#             n_CountriesThisStudy = length(unique(Country)), 
#             n_ZonesThisStudy = length(unique(Zone)))

frug_studies <- frug %>% group_by(Study_ID) %>% # But many studies cover multiple sites
  summarize(n_SitesThisStudy = length(unique(Study_Site)), 
            n_CountriesThisStudy = length(unique(Country)), 
            n_ZonesThisStudy = length(unique(Zone)), 
            focus = unique(Focus), 
            n_AnimalsThisStudy = length(unique(Animal_Species_Corrected)), 
            n_PlantsThisStudy =  length(unique(Plant_Species_Corrected)))
filter(frug_studies, n_SitesThisStudy > 1 | n_CountriesThisStudy > 1 | n_ZonesThisStudy > 1)
filter(frug_studies, is.na(n_SitesThisStudy))
table(frug_studies$focus)
table(frug_studies[frug_studies$focus == "Animal", 'n_AnimalsThisStudy'])
table(frug_studies[frug_studies$focus == "Plant", 'n_PlantsThisStudy'])

# Look at multi-sites studies:
sum(frug_studies$n_SitesThisStudy>1) # 17 multi-site studies
sum(frug_studies[frug_studies$n_SitesThisStudy>1, ]$n_SitesThisStudy) # 39 study-sites from multi-site studies
sum(frug_studies$n_CountriesThisStudy>1 & frug_studies$n_SitesThisStudy == 1) # 4 multi-country studies with 1 site due to NA site
sum(frug_studies$n_ZonesThisStudy>1 & frug_studies$n_SitesThisStudy == 1) # 11 multi-country studies with 1 site due to NA site


#-------------------------- PART B: Cleaning and extracting useful info -------------------#

#### Recode Study ID
# Note that if one study contains NA and non-NA study sites, this will create new ID for all NA, which will just be the Study_ID
# Also clean some plant families in here
frug <- frug %>%        
        left_join(., frug_studies) %>% 
        # Differentiate meta analysis sources
        mutate(Study_ID = ifelse(grepl("D", Study_ID) == FALSE, paste0("CD", Study_ID), Study_ID))%>% 
         # Create new study identifier for multi-site (or zone or country) studies 
        mutate(Study_ID_Plus = ifelse(n_SitesThisStudy>1 & !is.na(Study_Site), paste0(Study_ID, Study_Site), 
                                      ifelse(n_CountriesThisStudy > 1 & !is.na(Country), paste0(Study_ID, Country), 
                                      ifelse(n_ZonesThisStudy >1 & !is.na(Zone), paste0(Study_ID, Zone),
                                      ifelse((n_SitesThisStudy>1 | n_CountriesThisStudy > 1 | n_ZonesThisStudy >1) 
                                             & (is.na(Study_Site) &is.na(Country) & is.na(Zone)), paste0(Study_ID, "NA"), 
                                      Study_ID))))) %>%  
        # Some species are known by different scientific names
        mutate(Animal_Species_Corrected = ifelse(Animal_Species_Corrected == "Piliocolobus_temmincki",
                                           "Piliocolobus_badius", Animal_Species_Corrected)) %>%
        mutate(Animal_Species_Corrected = ifelse(Animal_Species_Corrected == "Dyaphorophyia_castanea",
                                                 "Platysteira_hormophora", Animal_Species_Corrected)) %>%
        mutate(Animal_Species_Corrected = ifelse(Animal_Species_Corrected == "Iduna_pallida", 
                                           "Hippolais_pallida", Animal_Species_Corrected)) %>% 
        mutate(Animal_Species_Corrected = ifelse(Animal_Species_Corrected == "Melaniparus_albiventris", 
                                                 "Parus_albiventris", Animal_Species_Corrected)) %>%
        mutate(Animal_Species_Corrected = ifelse(Animal_Species_Corrected == "Melaniparus_funereus", 
                                                 "Parus_funereus", Animal_Species_Corrected)) %>%
        mutate(Plant_Family = ifelse(Plant_Species_Corrected == "Poga_oleosa", "Anisophylleaceae", Plant_Family)) %>%
        mutate(Plant_Family = ifelse(Plant_Species_Corrected == "Harungana_madagascariensis", "Hypericaceae", Plant_Family)) %>%
        mutate(Plant_Family = ifelse(Plant_Species_Corrected == "Mammea_africana", "Calophyllaceae", Plant_Family)) %>%
        mutate(Plant_Family = ifelse(Plant_Species_Corrected == "Olea_capensis", "Oleaceae", Plant_Family)) %>%
        mutate(Plant_Family = ifelse(Plant_Species_Corrected == "Olea_welwitschii", "Oleaceae", Plant_Family)) %>% 
        # Misc typos: Plants
        mutate(Plant_Species_Corrected = ifelse(Plant_Species_Corrected == "Vachellia xanthophloea","Vachellia_xanthophloea", Plant_Species_Corrected)) %>% 
        mutate(Plant_Species_Corrected = ifelse(Plant_Species_Corrected == "Senegalia_ nigrescens","Senegalia_nigrescens", Plant_Species_Corrected)) %>% 
        mutate(Plant_Species_Corrected = trimws(Plant_Species_Corrected)) %>%   # Remove leading or trailing whitespace
        mutate(Animal_Genus = gsub("_.*", "", Animal_Species_Corrected), 
               Plant_Genus = gsub("_.*", "", Plant_Species_Corrected)) %>% 

        # Include the original genus as an alternative for reclassified species
        mutate(Animal_Genus2 = ifelse(Animal_Species_Corrected == "Platysteira_hormophora", "Dyaphorophyia", Animal_Genus)) %>%
        mutate(Animal_Genus2 = ifelse(Animal_Species_Corrected == "Hippolais_pallida", "Iduna", Animal_Genus2)) %>% 
        mutate(Animal_Genus2 = ifelse(Animal_Species_Corrected == "Parus_albiventris","Melaniparus",  Animal_Genus2)) %>%
        mutate(Animal_Genus2 = ifelse(Animal_Species_Corrected =="Parus_funereus","Melaniparus",  Animal_Genus2)) 
        
        
### Clean taxonomic information 
# Extract plant and animal list for frugivory only data: alphabetize now, sort phylo later

v.taxa <- frug %>%
  dplyr::select(Animal_Species_Corrected, Animal_Genus, Animal_Genus2, Animal_Family, Animal_Taxa_Type) %>%
  unique() %>%
  mutate(Animal_Class = ifelse(Animal_Taxa_Type == "Large_birds" | Animal_Taxa_Type == "Small_birds", "Aves", "Mammalia")) %>% 
  mutate(Animal_Taxa_Type_Courser = ifelse(Animal_Taxa_Type == "Apes" | Animal_Taxa_Type == "Monkeys" |Animal_Taxa_Type == "BushBaby", 
                               "Primates", Animal_Taxa_Type)) %>%
  arrange(Animal_Class, Animal_Taxa_Type_Courser, Animal_Taxa_Type, Animal_Family, Animal_Genus, Animal_Species_Corrected) 
v.names <- v.taxa$Animal_Species_Corrected
  
p.taxa <- frug %>% 
    dplyr::select(Plant_Species_Corrected, Plant_Genus, Plant_Family) %>%
    unique() %>%
    mutate(Plant_Family = ifelse(Plant_Genus == "Bersama", "Francoaceae", Plant_Family)) %>% 
    mutate(Plant_Family = ifelse(Plant_Genus == "Carduus" | Plant_Genus == "Sphaeranthus", "Asteraceae", Plant_Family)) %>% 
    mutate(Plant_Family = ifelse(Plant_Genus == "Cleistanthus" | Plant_Genus == "Spondianthus", "Phyllanthaceae", Plant_Family)) %>% 
    mutate(Plant_Family = ifelse(Plant_Genus == "Memecylon" | Plant_Genus == "Warneckea", "Melastomataceae", Plant_Family)) %>% 
    mutate(Plant_Family = ifelse(Plant_Genus == "Tricliceras", "Passifloraceae", Plant_Family)) %>% 
    arrange(Plant_Family, Plant_Genus, Plant_Species_Corrected) 


# Look for duplicate plants
n_occur <- data.frame(table(p.taxa$Plant_Species_Corrected))
n_occur[n_occur$Freq > 1, ]
p.taxa[p.taxa$Plant_Species_Corrected %in% n_occur[n_occur$Freq > 1, 1],]

p.names <- p.taxa$Plant_Species_Corrected

# v.names <- frug %>%
#   dplyr::select(Animal_Species_Corrected)%>%
#   #arrange(Animal_Species_Corrected) %>%
#   unique() %>%
#   pull(Animal_Species_Corrected) 

# p.names <- frug %>%
#   dplyr::select(Plant_Species_Corrected) %>%
#   arrange(Plant_Species_Corrected) %>%
#   unique() %>%
#   pull(Plant_Species_Corrected)

s.names <- frug %>%
  dplyr::select(Study_ID)%>%
  unique() %>%
  pull(Study_ID)

sp.names <- frug %>%
  dplyr::select(Study_ID_Plus)%>%
  unique() %>%
  pull(Study_ID_Plus)

nV <- length(v.names)
nP <- length(p.names)
nS <- length(s.names)
nSP <- length(sp.names) 


# Look at vertebrate and plant families
animal_families <- frug %>%
  dplyr::select(Animal_Species_Corrected, Animal_Family) %>%
  unique()

# Note that some studies are multi-site: my counts don't match up with the ones in the Excel file, weird 
site_info <- frug %>% group_by(Study_Site) %>% 
  summarize(nV = length(unique(Animal_Species_Corrected)), 
            nP = length(unique(Plant_Species_Corrected)), 
            nS = length(unique(Study_ID))) %>% 
  left_join(., zones, join_by(Study_Site == StudyZone)) %>% 
  dplyr::select(Study_Site, nV, nP, nS, Lat, Long)

# Summarize study and site info
study_site_info <- frug %>% group_by(Study_ID, Study_Site, Country, Zone) %>% 
                   summarize(n_animal_species = length(unique(Animal_Species_Corrected)), 
                             n_plant_species = length(unique(Plant_Species_Corrected)), 
                             n_interactions = n()) %>%
                   left_join(., site_info) %>% 
                   dplyr::select(-c(nV, nP, nS))

sum(is.na(study_site_info$Lat))

## Create data frame summarizing study information
# How many times is each study site repeated? 
# animal_families <- frug %>% 
#   select(Animal_Species_Corrected, Animal_Family) %>% 
#   unique() 
# study_focus1 <- data.frame(Study_ID = double(), nfocus = integer(), Species = character())
# study_focus2 <- data.frame(Study_ID = double(), nfocus = integer(), Family1 = character(), Family2 = character(),
#                            Family3 = character(), 
#                            Species1 = character(), 
#                            Species2 = character(), Species3 = character(), Species4 = character(), 
#                            Species5 = character(), Species6 = character(), Species7 = character())
# 
# counter1 <- 1
# counter2 <- 1
# for(ss in 1:nSP){ # maybe redo with Study_ID_Plus
#   this_study <- sp.names[ss]
#   these_obs <- frug[frug$Study_ID_Plus == this_study, ]
#   these_animals <- unique(these_obs$Animal_Species_Corrected)
#   these_families <- unique(these_obs$Animal_Family)
#   cat("\n animals", these_animals)
#   cat("\n families", these_families)
#   # for(a in 1:length(these_animals)){
#   #   study_focus1[counter1,] <- c(this_study, length(these_animals), these_animals[a])
#   #   counter1 <- counter1 + 1
#   # }
#   # study_focus2[counter2, ] <- c(this_study, length(these_animals), these_families, rep(NA, 3 - length(these_families)), 
#   #                               these_animals, rep(NA, 7-length(these_animals)))
#   # counter2 <- counter2 + 1
# }
# study_focus1$Study_ID <- as.numeric(study_focus1$Study_ID)
# study_focus2$Study_ID <- as.numeric(study_focus2$Study_ID)
# study_focus1 <- left_join(study_focus1, animal_families, join_by("Species" == "Animal_Species_Corrected"))


#---------------------- PART C: CREATE THE ADJACENCY MATRICES --------------------------#
# Convert raw data into a list of adjacency matrices
# Adj_matrix <- array( 0, dim= c(nrow(animal_list), nrow(plant_list), nrow(study_list)))
# interactions <- matrix(0, nrow= nrow(animal_list), ncol= nrow(plant_list))
# Adj_matrix_new <- lapply(seq_len(length(studies)), function(X) interactions)
# Note that we use Study_ID_Plus rather than Study_ID to distinguish between sites at multi-site studies

# First: create a list of adjacency matrices of varying dimension
frug_new <- frug %>% dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected, Study_ID, Study_ID_Plus)
frug_new$Interaction <- 1

Adj_matrix_new <- list()
for(ss in 1:nSP){
  Adj <- frug_new %>% filter(Study_ID_Plus==sp.names[ss]) %>% 
    dplyr::select(Animal_Species_Corrected, Plant_Species_Corrected, Interaction) %>% 
    distinct() %>% pivot_wider(names_from = Plant_Species_Corrected, values_from=Interaction)
  # Case 1: one plant species only (species names will be lost)
  #if(sum(dim(as.matrix(as.data.frame(Adj)[, -1]))) == 2){
  if(dim(Adj)[2] == 2){
    this.v <- Adj$Animal_Species_Corrected
    this.p <- colnames(Adj)[2]
    Adj_df  <- as.data.frame(Adj)[,-1]
    Adj_matrix <- as.matrix(Adj_df)
    dimnames(Adj_matrix) <- list(this.v, this.p)
  # Case 2: multiple species of at least one type
  }else{
    Adj_df  <- as.data.frame(Adj)
    rownames(Adj_df) <- Adj_df$Animal_Species_Corrected
    Adj_df <- Adj_df[,-1]
    Adj_df[is.na(Adj_df)] <- 0
    Adj_matrix <- as.matrix(Adj_df)
  }
  s.id <- sp.names[ss]
  Adj_matrix_new[[s.id]] <- Adj_matrix # add back rownames and colnames for 1,1 studies
}

# Now create adjacency array by augmenting with zeros
acomb <- function(...) abind(..., along=3)
registerDoParallel(10)

A_obs <- foreach(i=1:nSP, .combine = 'acomb') %dopar% {
  A <- Adj_matrix_new[[i]]
  wh.v <- sapply(dimnames(A)[[1]], function(x) which(x==v.names))
  wh.p <- sapply(dimnames(A)[[2]], function(x) which(x==p.names))
  
  A.full <- matrix(0, nrow = nV, ncol = nP)
  A.full[wh.v, wh.p] <- A
  A.full
}

dimnames(A_obs) <- list(v.names, p.names, sp.names)

# Double check: can delete later
this.A <- A_obs[,,1]
rownames(this.A)[which(rowSums(this.A) != 0)]
colnames(this.A)[which(colSums(this.A)!= 0)]
frug_new %>% filter(Study_ID_Plus==sp.names[1])

this.A <- A_obs[,,5]
rownames(this.A)[which(rowSums(this.A) != 0)]
colnames(this.A)[which(colSums(this.A)!= 0)]
frug_new %>% filter(Study_ID_Plus==sp.names[5])

#----------------- PART D: CREATING THE OCCURRENCE MATRICES -------------------#

## Extract occurrence information using a variety of methods, determining occurrence probability when
# Species occurs in the same study -> 1
# Species occurs at the same site, but a different study -> 0.75
# Species occurs in the same country, but a different site -> 0.5
# Species occurs in the same zone, but different ... -> 0.25

## Create occurrence matrices
O_P <- matrix(0, nrow = nP, ncol = nSP) 
rownames(O_P) <- p.names
O_V <- matrix(0, nrow = nV, ncol = nSP) 
rownames(O_V) <- v.names
colnames(O_P) <- colnames(O_V) <- sp.names


## EDA: check for multi-site studies in the combined data
frug %>% group_by(Study_ID) %>% 
  summarize(n_Sites = length(unique(Study_Site)), 
            n_Countries = length(unique(Country)), 
            n_Zones = length(unique(Zone))) %>% 
  filter(n_Sites > 1 | n_Countries > 1 | n_Zones > 1)

for(ss in 1:nSP){
  IDplus <- sp.names[ss]
  study_site <- unique(frug[frug$Study_ID_Plus == IDplus, "Study_Site"])
  study_zone <- unique(frug[frug$Study_ID_Plus == IDplus, "Zone"])
  study_country <- unique(frug[frug$Study_ID_Plus == IDplus, "Country"]) # site, country, zone can all be NA, think about how to handle that
  
  study_site_dta <- frug[which(frug$Study_ID_Plus == IDplus ),] # same study and site: this is always nonempty
  site_dta <- frug[which(frug$Study_Site == study_site & !(is.na(frug$Study_Site))),] # same non-NA site: may be empty
  country_dta <- frug[which(frug$Country == study_country & !(is.na(frug$Country))), ] # same non-NA country: may be empty
  zone_dta <- frug[which(frug$Zone == study_zone & !(is.na(frug$Zone))), ] # same non-NA region: may be empty
  
  for(p in 1:nP){
    this.plant <- p.names[p]
    if(length(study_site_dta$Plant_Species_Corrected[study_site_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 1
    } else if(length(site_dta$Plant_Species_Corrected[site_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.75
    }  else if(length(country_dta$Plant_Species_Corrected[country_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.5
    } else if(length(zone_dta$Plant_Species_Corrected[zone_dta$Plant_Species_Corrected==this.plant])>0){
      O_P[p,ss] <- 0.25
    } 
  }  
  for(v in 1:nV){
    this.animal <- v.names[v]
    if(length(study_site_dta$Animal_Species_Corrected[study_site_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 1
    } else if(length(site_dta$Animal_Species_Corrected[site_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.75
    } else if(length(country_dta$Animal_Species_Corrected[country_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.5
    } else if(length(zone_dta$Animal_Species_Corrected[zone_dta$Animal_Species_Corrected==this.animal])>0){
      O_V[v,ss] <- 0.25
    }
  }
  
}


# ------ PART E: Preparing and subsetting focus -------- #

### USING DEFAULT FOCUS

# Creating two arrays that correspond to whether the studies would have
# recorded an observed interaction or not, and whether species co-occur in the
# study area. We will assume that these arrays are 0/1. For the focus array,
# this makes sense as we have access to whether the study is animal or plant
# focused or a network study. For the species occurrence, we know for a fact
# that a species with a recorded interaction occurs in the area. The only
# simplifying assumption is therefore when we set probability of occurrence to
# zero for species without any recorded interaction.


# We have to categorize studies first: all animal-focused except for maybe one (Study 13)
F_obs <- array(0, dim = c(nV, nP, nSP))
dimnames(F_obs) <- list(v.names, p.names, sp.names)
methods <- frug %>% group_by(Study_ID_Plus) %>% 
                    summarize(method = unique(Focus)) %>% 
                    mutate(method = ifelse(is.na(method), "Unknown", method))

# EDA code (to remove possibly)
# investigate bothlimited a little more, maybe it's actually obvious which is which
bl <- methods %>% filter(method =="BothLimited") %>% 
  pull("Study_ID_Plus")
num <- 1
Adj_matrix_new[[which(bl[num] == sp.names)]] # run this multiple times for an idea
num <- num + 1
         
for (ss in 1 : nSP) {
  this_study <- sp.names[ss]
  this_A <- A_obs[,,ss]
  this_method <- methods[ss, "method"]
  these_v <- which(rowSums(this_A)>0)
  these_p <- which(colSums(this_A)>0)
  # Observability depends on Focus
  if(this_method == "Animal"){
    F_obs[these_v, , ss] <- 1 # Animal focused: all interactions between target vertebrates and any plants recorded
  } else if(this_method == "Plant"){
    F_obs[, these_p , ss] <- 1 # Plant focused: all interactions between target plants and any vertebrates recorded    
  } else if(this_method == "BothLimited" | this_method == "Unknown"){
    F_obs[,,ss] <- this_A # BothLimited, only observed interactions are known to be observable b/c we can't tell which plants vs which animals were focal
  } else { # Both or network studies: all interactions are observable
    F_obs[,,ss] <- 1 # All observations recorded in a network study
  }
}




if(save_files){
  # The cleaned core data files
  save(A_obs, file = paste0(save_path, "A_obs.dat"))
  save(O_P, file = paste0(save_path, 'OP_full.dat'))
  save(O_V, file = paste0(save_path, 'OV_full.dat'))
  save(F_obs, file = paste0(save_path, 'F_obs.dat'))
  # Cleaned site info and study methods
  write.csv(site_info, file = paste0(save_path, 'site_info.csv'), row.names = FALSE)
  write.csv(methods, file = paste0(save_path,'study_methods.csv'), row.names = FALSE)
  write.csv(study_site_info, file = paste0(save_path, 'study_site_info.csv'), row.names = FALSE)
  # Cleaned plant and vertebrate taxa information
  write.csv(v.taxa, file = paste0(save_path, "v_taxa.csv"), row.names = FALSE)
  write.csv(p.taxa, file = paste0(save_path, "p_taxa.csv"), row.names = FALSE)
  # save(study_info2, file = paste0(save_path, 'study_info2.dat'))
  # save(unique_locs, file = paste0(save_path, 'loc_info.dat'))
}

