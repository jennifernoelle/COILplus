# This file performs basic data cleaning on the interaction data
# This file processes the raw data to create the following key data sources: 
  # A_obs, F_obs, O_A, O_P
# Saves useful info about studies: 
  # site_info, study_methods, study_site_info
# Creates taxa databases:
  # v.taxa, p.taxa
# Creates anonymized data for distribution: generic species names and sites  (country and zone ok)

# Uses these files
  # Frugivory_combined_cleaned.csv
  # StudyZone_clean2.csv


# ---------------------------------- TO DO ----------------------------------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawDataPub/'
save_path <- 'ProcessedDataPub/' # Create this directory
save_path_generic <- 'RawDataGeneric/'

save_files <- TRUE

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



#--------------------- PART B2. REPLACE TAXA NAMES AND SITE NAMES WITH GENERICS ---------------#

### Vertebrates/frugivores

# Set up useful data structures
nV_fams <- length(unique(v.taxa$Animal_Family))
V_sp <- paste0("Frugivore_", 1:nV)

# Create generic-real codebook
v_family_generic <- data.frame(Animal_Family_Generic = paste0("Family_", 1:nV_fams), 
                               Animal_Family = unique(v.taxa$Animal_Family))
v.taxa$Animal_Species_Generic = V_sp
v_taxa_generic <- left_join(v.taxa, v_family_generic)


### Plants

# Set up useful data structures
nP_fams <- length(unique(p.taxa$Plant_Family))
P_sp <- paste0("Plant_", 1:nP)

# Create generic-real codebook
p_family_generic <- data.frame(Plant_Family_Generic = paste0("Family_", 1:nP_fams), 
                               Plant_Family = unique(p.taxa$Plant_Family))
p.taxa$Plant_Species_Generic = P_sp

# Exclude family from the data used for merging because frug contains 
# taxonomical errors, corrected above in the p.taxa
p_taxa_generic <- left_join(p.taxa, p_family_generic) %>% 
  select(-Plant_Family) 


### Site names
site_names <- unique(frug$Study_Site)
nSites <- length(site_names)

sites_generic <- data.frame(Study_Site_generic = paste0("Site_", 1:nSites), 
                            Study_Site = site_names) %>% 
                 mutate(Study_Site_generic = ifelse(is.na(Study_Site), NA, Study_Site_generic)) 
sites_codebook <- frug %>% select(Study_ID, Study_ID_Plus, Study_Site, Country, Zone, 
                                  n_SitesThisStudy, n_CountriesThisStudy, n_ZonesThisStudy) %>% 
                  unique() %>% 
                  left_join(., sites_generic) %>%
                  mutate(Study_ID_Plus_generic = ifelse(n_SitesThisStudy>1 & !is.na(Study_Site), paste0(Study_ID, Study_Site_generic), 
                                ifelse(n_CountriesThisStudy > 1 & !is.na(Country), paste0(Study_ID, Country), 
                                       ifelse(n_ZonesThisStudy >1 & !is.na(Zone), paste0(Study_ID, Zone),
                                              ifelse((n_SitesThisStudy>1 | n_CountriesThisStudy > 1 | n_ZonesThisStudy >1) 
                                                     & (is.na(Study_Site) &is.na(Country) & is.na(Zone)), paste0(Study_ID, "NA"), 
                                                     Study_ID)))))
                


### Create the generic data to save

frug_generic <- frug %>% 
                left_join(., sites_codebook) %>%
                left_join(., p_taxa_generic) %>%
                left_join(., v_taxa_generic) %>%
                select(Study_ID, Study_ID_Plus_generic, Country, Zone, Focus, Study_Site_generic, 
                      Animal_Species_Generic, Plant_Species_Generic, 
                       Animal_Family_Generic, Plant_Family_Generic) %>%
                rename(Animal_Species_Corrected = Animal_Species_Generic, 
                       Animal_Family = Animal_Family_Generic, Plant_Family = Plant_Family_Generic, 
                       Plant_Species_Corrected = Plant_Species_Generic, 
                       Study_Site = Study_Site_generic, Study_ID_Plus = Study_ID_Plus_generic)

p_taxa_generic <- p_taxa_generic %>% 
  select(Plant_Species_Generic, Plant_Family_Generic)

v_taxa_generic <- v_taxa_generic %>% 
  select(Animal_Species_Generic, Animal_Family_Generic)

# Save codebook versions for going between anonymous and labelled data
p.taxa.codebook <- p.taxa
v.taxa.codebook <- v.taxa

# Remove generics from labelled data
p.taxa <- p.taxa %>% select(., -Plant_Species_Generic)
v.taxa <- v.taxa %>% select(.,-Animal_Species_Generic)


### Save generics
if(save_files){
  write.csv(frug_generic, file = paste0(save_path_generic, 'frug_generic.csv'), row.names = FALSE)
  write.csv(v_taxa_generic, file = paste0(save_path_generic, "v_taxa_generic.csv"), row.names = FALSE)
  write.csv(p_taxa_generic, file = paste0(save_path_generic, "p_taxa_generic.csv"), row.names = FALSE)
}

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
# study area. 

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
  write.csv(v.taxa.codebook, file = paste0(save_path, "v_taxa_codebook.csv"), row.names = FALSE)
  write.csv(p.taxa.codebook, file = paste0(save_path, "p_taxa_codebook.csv"), row.names = FALSE)
  
}

