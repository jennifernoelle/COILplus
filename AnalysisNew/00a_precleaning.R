# This file performs pre-cleaning
# We may not want to include it in the final repo


# -------- TO DO --------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawDataNew/'

library(tidyverse)
save_files <- TRUE

# Load data
frug  <- read.csv(paste0(data_path, "Frugivory_combined.csv"))[, -1]
zones <- read.csv(paste0(data_path, "StudyZone_clean.csv"))[, -c(8,9)] 
zones <- zones %>% mutate(StudyZone = ifelse(StudyZone == "", NA, StudyZone))
zones <- zones[which(rowSums(is.na(zones)) != ncol(zones)), ]

apply(frug, 2, function(x) sum(is.na(x)))
apply(zones, 2, function(x) sum(is.na(x)))

z.sites <- sort(zones$StudyZone)
f.sites <- sort(unique(frug$Study_Site))
f.sites[!(f.sites %in% z.sites)]

frug <- frug %>%
           mutate(Study_Site = ifelse(Study_Site == "Adiopodoumé","Adiopodoume", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Bossematié", "Bossematie", Study_Site)) %>% 
           mutate(Study_Site = ifelse(Study_Site == "Banyang-bo", "Banyang-Mbo", Study_Site)) %>% 
           mutate(Study_Site = ifelse(Study_Site == "Bioko_Island", "Bioko", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Lwamundo", "Lwamunda", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Manovo Gounda St Floris National Park- Gounda ", "Manovo Gounda St Floris National Park", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Mawambi_Hills", "Mawambi-Hills", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Ndoki_Odzala", "Ndoki", Study_Site)) %>%
           mutate(Study_Site = ifelse(Study_Site == "Odobullu Forest", "Bale Mountains", Study_Site)) %>% 
           mutate(Study_Site = gsub("_", "-", Study_Site))
  
# Can't find: Kala Malou National Park, Logging concession

# Make sure country and region are consistent for each site
site_country <- frug %>% group_by(Study_Site) %>%
         summarize(n_countries = length(unique(Country)), 
                   n_zones = length(unique(Zone))) %>%
         filter(n_countries > 1 | n_zones > 1)
site_country
  
frug <- frug %>% mutate(Country = ifelse(Study_Site == "Budongo" & !is.na(Study_Site), "Uganda", Country)) %>% 
                 mutate(Zone = ifelse(Study_Site == "Budongo" & !is.na(Study_Site), "AE", Zone))

if(save_files){
  write.csv(frug, file = paste0(data_path, 'Frugivory_combined_clean.csv'))
  write.csv(zones, file = paste0(data_path, 'StudyZone_clean2.csv'))
}
  
  