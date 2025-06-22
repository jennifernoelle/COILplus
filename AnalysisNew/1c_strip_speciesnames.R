# This file takes the processed data and replaces species names with generics
# We should probably actually include the code to create O and F from the data
# Think about how to clean sites info <---------------- TO DO 
# What we really  need is to strip names as soon as frug is cleaned enough to process
# into adjacency matrices...do that later


#---------------------------- TO DO -------------------------------------------#

# The directory where the processed data should be saved:
data_path <- 'ProcessedDataPub/'
save_path <- 'GenericData/'
# Whether the processed data should be saved or not:
save_files <- TRUE

#---------------------- STEP 1: LOAD DATA AND PACKAGES ------------------------#

# Loading the data:
load(paste0(data_path, 'Cu_phylo.dat'))
load(paste0(data_path, 'Cv_phylo.dat'))
load(paste0(data_path, 'A_obs.dat'))
load(paste0(data_path, 'F_obs.dat'))
load(paste0(data_path, 'Obs_X.dat')) # vertebrate traits
load(paste0(data_path, 'Obs_W.dat')) # plant traits
load(paste0(data_path, 'OP_full.dat')) # site level obs plants
load(paste0(data_path, 'OV_full.dat')) # site level obs verts

# Taxa database
v_taxa <- read.csv(paste0(data_path, "v_taxa.csv"))
p_taxa <- read.csv(paste0(data_path, "p_taxa.csv"))

# Study info



#--------------------- STEP 2. REPLACE TAXA NAMES WITH GENERICS ---------------#

### Vertebrates/frugivores

# Set up useful data structures
nV <- nrow(v_taxa)
nV_fams <- length(unique(v_taxa$Animal_Family))
V_sp <- paste0("Frugivore_", 1:nV)

# Create generic-real codebook
v_family_generic <- data.frame(Animal_Family_Generic = paste0("Family_", 1:nV_fams), 
                               Animal_Family = unique(v_taxa$Animal_Family))
v_taxa$Animal_Species_Generic = V_sp
v_taxa_generic <- left_join(v_taxa, v_family_generic) %>% 
                  select(Animal_Species_Generic, Animal_Family_Generic, Animal_Family)

# Replace species names with generics in all the key files
rownames(Cu_phylo) <- colnames(Cu_phylo) <- V_sp
rownames(A_obs) <- rownames(F_obs) <- rownames(Obs_X) <- rownames(O_V) <- V_sp


### Plants

# Set up useful data structures
nP <- nrow(p_taxa)
nP_fams <- length(unique(p_taxa$Plant_Family))
P_sp <- paste0("Plant_", 1:nP)

# Create generic-real codebook
p_family_generic <- data.frame(Plant_Family_Generic = paste0("Family_", 1:nP_fams), 
                               Plant_Family = unique(p_taxa$Plant_Family))
p_taxa$Plant_Species_Generic = P_sp
p_taxa_generic <- left_join(p_taxa, p_family_generic) %>% 
  select(Plant_Species_Generic, Plant_Family_Generic, Plant_Family)

# Replace species names with generics in all the key files
rownames(Cv_phylo) <- colnames(Cv_phylo) <- P_sp
colnames(A_obs) <- colnames(F_obs) <- rownames(Obs_W) <- rownames(O_P) <- P_sp

#-------------------- REPLACE STUDY LOCATION WITH GENERICS --------------------#


#--------------------------- SAVE FILES ---------------------------------------#

if(save_files){
  write.csv(v_taxa_generic, paste0(save_path,"v_taxa_generic.csv"), row.names = FALSE)
  write.csv(p_taxa_generic, paste0(save_path,"p_taxa_generic.csv"), row.names = FALSE)
  save(A_obs, paste0(save_path, "A_obs.dat"))
  save(Cu_phlyo, paste0(save_path, "Cu_phylo.dat"))
  save(Cu_phlyo, paste0(save_path, "Cv_phylo.dat"))
  save(Obs_X, paste0(save_path, "Obs_X.dat"))
  save(Obs_W, paste0(save_path, "Obs_W.dat"))
}




                  





