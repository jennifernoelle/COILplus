# Clean the supplementary traits data: Obs_X, Obs_W

# ----------------------------------- TO DO ---------------------------------- #

# Set the directories below to correspond to paths on your machine:

# The directory within the working directory where the data are:
data_path <- 'RawDataNew/'
save_path <- 'ProcessedDataNew/' # Create this directory
save_files <- TRUE

library(tidyverse)

#------------------------------ PART 1: EDA -----------------------------------#

# Read the data
traits.v.raw <- read.csv(paste0(data_path, "frugivore_trait.csv"))[, -1]
traits.p.raw <- read.csv(paste0(data_path, "plant_trait.csv"))[, -1]
v.taxa <- read.csv(paste0(save_path, "v_taxa.csv"))
p.taxa <- read.csv(paste0(save_path, "p_taxa.csv"))

# Use cleaned species info: problem - this leads to mismatches due to cleaning of the taxa df
v.taxa$Animal_Species_Corrected[!(v.taxa$Animal_Species_Corrected %in% traits.v.raw$Species)]
traits.v.raw <- traits.v.raw %>% mutate(Species= ifelse(Species == "Piliocolobus_temmincki",
                                         "Piliocolobus_badius", Species)) %>%
  mutate(Species = ifelse(Species == "Dyaphorophyia_castanea",
                                           "Platysteira_hormophora", Species)) %>%
  mutate(Species = ifelse(Species == "Iduna_pallida", 
                                           "Hippolais_pallida", Species)) %>% 
  mutate(Species = ifelse(Species == "Melaniparus_albiventris", 
                                           "Parus_albiventris", Species)) %>%
  mutate(Species = ifelse(Species == "Melaniparus_funereus", 
                                           "Parus_funereus", Species)) 

p.taxa$Plant_Species_Corrected[!(v.taxa$Plant_Species_Corrected %in% traits.p.raw$Species)]
traits.v <- left_join(v.taxa, traits.v.raw, join_by(Animal_Species_Corrected == Species)) %>% 
             select(-c(Genus, Family, Taxa)) %>% 
             rename(Species = Animal_Species_Corrected, Genus = Animal_Genus, 
                    Family = Animal_Family, 
                    Taxa = Animal_Taxa_Type)

traits.p <- left_join(p.taxa, traits.p.raw[, -4]) %>% 
  select(-c(Genus)) %>% 
  rename(Species = Plant_Species_Corrected, Genus = Plant_Genus, Family =Plant_Family)


apply(traits.v.raw, 2, function(x) sum(is.na(x)))
apply(traits.v, 2, function(x) sum(is.na(x)))
apply(traits.p.raw, 2, function(x) sum(is.na(x)))
apply(traits.p, 2, function(x) sum(is.na(x)))

ggplot(traits.v[traits.v$Taxa != "Elephant",], aes(x = Taxa, y = Body_Mass)) + 
  geom_boxplot()

ggplot(traits.v, aes(x = Taxa, y = Generation_Length)) + 
  geom_boxplot()

# Check for inconsistent information: vertebrates
n_occur <- data.frame(table(traits.v$Species))
n_occur[n_occur$Freq > 1, ]
traits.v[traits.v$Species %in% n_occur[n_occur$Freq > 1, 1],]
traits.v <- traits.v %>% 
            mutate(Body_Mass = ifelse(Species == "Piliocolobus_badius", (8059 + 7750)/2, Body_Mass)) %>% 
            mutate(IUCN_Status = ifelse(Species == "Piliocolobus_badius", "EN", IUCN_Status)) %>%   
            mutate(Habitat = ifelse(Species == "Piliocolobus_badius", "Mixed", Habitat)) %>%      
            mutate(Body_Mass = ifelse(Species == "Cercopithecus_erythrogaster", (4000 + 8059)/2, Body_Mass)) %>% 
            mutate(Generation_Length = ifelse(Species == "Cercopithecus_erythrogaster", (9+10)/2, Generation_Length)) %>% 
            mutate(Habitat = ifelse(Species == "Cercopithecus_erythrogaster", "Mixed", Habitat)) %>% 
            mutate(Body_Mass = ifelse(Species == "Lophocebus_albigena", (7096 + 6000)/2, Body_Mass)) %>% 
            mutate(IUCN_Status = ifelse(Species == "Lophocebus_albigena", "VU", IUCN_Status)) %>% 
            mutate(IUCN_Status = ifelse(Species == "Mandrillus_leucophaeus", "EN", IUCN_Status)) %>% 
            unique()
  
            
# Inconsistent information: plants
n_occur <- data.frame(table(traits.p$Species))
n_occur[n_occur$Freq > 1, ]
traits.p[traits.p$Species %in% n_occur[n_occur$Freq > 1, 1],]

# ------ PART 2: Preparing covariates for use in the model -------- #

# Remove taxonomic info from traits data
# Put continuous covariates first
# Only binary and continuous covariates allowed
Obs_X <- traits.v %>%
  mutate(IUCN_Endangered = ifelse(IUCN_Status == "CR" | IUCN_Status=="EN", 1, 0), 
         Habitat_Forest = ifelse(Habitat == "Forest", 1, 0), 
         Log_Body_Mass = log(Body_Mass)) %>%
  select(c(Generation_Length, Log_Body_Mass, IUCN_Endangered, Habitat_Forest))
rownames(Obs_X) <- traits.v$Species

Obs_W <- traits.p %>%
  select(FruitLength, FruitWidth, WoodDensity) # %>%
  #mutate(FruitWidth = as.numeric(ifelse(FruitWidth == ".1.1", 1.1, FruitWidth)))
rownames(Obs_W) <- traits.p$Species


# Impute missing fruit lengths, widths, wood density by family or genus
apply(traits.p, 2, function(x) sum(is.na(x)))
LtoW <- mean(Obs_W$FruitLength/Obs_W$FruitWidth, na.rm = TRUE)

fruit.df <- traits.p %>% 
  # If only one of length or width is missing
  mutate(FruitWidth = ifelse(is.na(FruitWidth) & !(is.na(FruitLength)), 
                             FruitLength/LtoW, FruitWidth), 
         FruitLength = ifelse(is.na(FruitLength) & !(is.na(FruitWidth)), 
                              FruitWidth*LtoW,  FruitLength)) %>%
  group_by(Genus) %>% # First get genus mean fruit widths
  mutate(GenusMean.FL = mean(FruitLength, na.rm = TRUE)) %>%
  mutate(GenusMean.FW = mean(FruitWidth, na.rm = TRUE)) %>%
  mutate(GenusMean.WD = mean(WoodDensity, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Family) %>% # Next get family mean fruit widths
  mutate(FamilyMean.FL = mean(FruitLength, na.rm = TRUE)) %>%
  mutate(FruitLength = ifelse(is.na(FruitLength), 
                              ifelse(is.na(GenusMean.FL), FamilyMean.FL, 
                                     GenusMean.FL), FruitLength)) %>%
  mutate(FamilyMean.FW = mean(FruitWidth, na.rm = TRUE)) %>%
  mutate(FruitWidth = ifelse(is.na(FruitWidth), 
                              ifelse(is.na(GenusMean.FW), FamilyMean.FW, 
                                     GenusMean.FW), FruitWidth)) %>%
  mutate(FamilyMean.WD = mean(WoodDensity, na.rm = TRUE)) %>% 
  mutate(WoodDensity = ifelse(is.na(WoodDensity), 
                             ifelse(is.na(GenusMean.WD), FamilyMean.WD, 
                                    GenusMean.WD), WoodDensity)) %>%
  select(FruitLength, FruitWidth, WoodDensity, Species) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

fruits.ordered <- left_join(traits.p, fruit.df, by="Species")
Obs_W$FruitLength <- as.numeric(fruits.ordered$FruitLength.y)
Obs_W$FruitWidth <- as.numeric(fruits.ordered$FruitWidth.y)
Obs_W$WoodDensity <- as.numeric(fruits.ordered$WoodDensity.y)


# Impute missing generation lengths by family
apply(traits.v, 2, function(x) sum(is.na(x)))
gen.df <- traits.v %>% 
  group_by(Genus) %>% # First get genus mean fruit widths
  mutate(GenusMean.g = mean(Generation_Length, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Family) %>% # Next get family mean fruit widths
  mutate(FamilyMean.g = mean(Generation_Length, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(Generation_Length = ifelse(is.na(Generation_Length), 
                              ifelse(is.na(GenusMean.g), FamilyMean.g, 
                                     GenusMean.g), Generation_Length)) %>%
  mutate_all(~ifelse(is.nan(.), NA, .))

gl.ordered <- left_join(traits.v, gen.df, by="Species")
Obs_X$Generation_Length <- as.numeric(gl.ordered$Generation_Length.y)

# How many missing values are left
apply(Obs_X, 2, function(x) sum(is.na(x)))
apply(Obs_W, 2, function(x) sum(is.na(x)))

if (save_files) {
  save(Obs_X, file = paste0(save_path, 'Obs_X.dat'))
  save(Obs_W, file = paste0(save_path, 'Obs_W.dat'))
}




