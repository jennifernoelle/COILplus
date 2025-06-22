# ---------- TO DO: set  your directories and list of models to analyze ---------#

# Paths that are the same for all models 
data_path_raw <- 'RawDataPub/'
data_path_processed <- 'ProcessedDataPub/'
save_plots_path <- 'ResultsPub'

# Load packages 
library(ggplot2)
library(sf)
library(raster)
library(rnaturalearth)


#### Load data: the data is very large so it'll take a few minutes
# Specify the URL of the file you want to download
url <- "https://zenodo.org/records/7764460/files/ps_africa_treecover_2019_100m_v0.1.tif?download=1"
# Specify the file name and location where you want to save the file on your computer
file_name <- "tree_cover.csv"
# Call the download.file() function, passing in the URL and file name/location as arguments
download.file(url, paste(data_path_raw, file_name, sep = ""), mode = "wb")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sites <- read.csv(paste0(data_path_processed,"site_info.csv"))
tree_cover <- raster(paste0(data_path_raw, file_name))

#### Setup data: this part  takes about ten minutes
africa <- world[world$continent == "Africa", ]
crop_extent <- extent(xmin(tree_cover), xmax(tree_cover), -16, 18)
tree_cover_clipped <- crop(tree_cover, crop_extent)
tree_cover_downsampled <- aggregate(tree_cover_clipped, fact = 50, fun = mean)
tree_cover_df <- as.data.frame(tree_cover_downsampled, xy = TRUE)
colnames(tree_cover_df) <- c("x", "y", "tree_cover")

#### Plots

# Number of plant species
ggplot(data = africa) +
  geom_raster(data = tree_cover_df, aes(x = x, y = y, fill = tree_cover)) +
  geom_sf(fill = NA, color = "black") +
  scale_fill_gradient(name = "Tree Cover", low = "white", high = "darkgreen", na.value = "transparent") +
  geom_point(data = sites, aes(x = Long, y = Lat, size = nP), 
             color = "black", fill = "goldenrod", shape = 21, stroke = 0.8) +
  theme_void() +
  coord_sf(xlim = c(crop_extent[1], crop_extent[2]), ylim = c(-14, 14)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_size_continuous(name = "# Plant Species")

# Number of studies
png(filename=(paste0(save_plots_path, "studies_map.png")), width = 1500, height = 700)

ggplot(data = africa) +
  geom_raster(data = tree_cover_df, aes(x = x, y = y, fill = tree_cover)) +
  geom_sf(fill = NA, color = "grey35") +
  scale_fill_gradient(name = "Tree Cover", low = "white", high = "darkgreen", 
                      na.value = "transparent", guide = "none") +
  geom_point(data = sites, aes(x = Long, y = Lat, size = nS), 
             color = "darkgoldenrod", fill = "goldenrod", shape = 21, stroke = 0.8, alpha = 0.9) +
  theme_void() +
  coord_sf(xlim = c(crop_extent[1], crop_extent[2]), ylim = c(-14, 14)) +
  theme(panel.border = element_rect(color = "grey35", fill = NA, size = 1)) +
  scale_size_continuous(name = "Number of studies", range = c(2,12)) + 
  theme(legend.position = "bottom", text = element_text(family = "serif", size = 16))

dev.off()


ggplot(data = africa) +
  geom_raster(data = tree_cover_df, aes(x = x, y = y, fill = tree_cover)) +
  geom_sf(fill = NA, color = "grey35") +
  scale_fill_gradient(name = "Tree Cover", low = "white", high = "#117733", 
                      na.value = "transparent", guide = "none") +
  geom_point(data = sites, aes(x = Long, y = Lat, size = nS), 
             color = "grey35", fill ="#882255", shape = 21, stroke = 0.8, alpha = 0.9) +
  theme_void() +
  coord_sf(xlim = c(crop_extent[1], crop_extent[2]), ylim = c(-14, 14)) +
  theme(panel.border = element_rect(color = "grey35", fill = NA, size = 1)) +
  scale_size_continuous(name = "Number of Studies") + 
  theme(legend.position = "bottom", text = element_text(family = "serif", size = 24))


