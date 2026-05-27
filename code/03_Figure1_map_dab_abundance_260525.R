# ==============================================================================
# DAB Transcript Analysis in Arctic Samples
# Author: Patrick White
# Date: August 11, 2025
# Description: Analysis of domoic acid biosynthesis (DAB) transcripts and 
#              Pseudo-nitzschia abundance across Arctic locations
# ==============================================================================

# Load required libraries
library(dplyr)        # Data manipulation
library(stringr)      # String manipulation
library(tidyr)        # Data reshaping
library(ggplot2)      # Plotting
library(sf)           # Spatial data handling
library(rnaturalearth) # World map data
library(ggspatial)     # scale bar

# ================================================================
# DATA LOADING AND PREPARATION
# ==============================================================================

# Read data (same file used previously)
comprehensive_data <- read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')


# Optional: ensure Loocation is a factor (keeps consistent coloring / ordering if needed)
comprehensive_data$Loocation <- factor(comprehensive_data$Loocation, 
                                          levels = c("BELCHER", "TALBOT", "SVERDRUP", "STATION_OG",
                                                     'GRISE_FIORD', 'JAKEMAN'))




# Calculate mean and standard deviation by location, year, and DAB gene
long_dab_count_summ3_mean = comprehensive_data %>% 
    group_by(Loocation, year, dab_name) %>% 
    summarise(MEAN_DABs = mean(count_dab, na.rm = TRUE),
              SD_DABs = sd(count_dab, na.rm = TRUE),
              .groups = 'drop') %>%
    mutate(MEAN_DABs = ifelse(is.na(MEAN_DABs), 0, MEAN_DABs))

# Define coordinate variables for each location
belcher_lon <- -80.5
belcher_lat <- 75.788031
talbot_lon <- -77.901711
talbot_lat <- 77.911486
sverdrup_lon <- -83.230252
sverdrup_lat <- 75.764418
grise_lon <- -82.997032
grise_lat <- 76.733404
jakeman_lon <- -81.020115
jakeman_lat <- 76.409441
station_og_lon <- -82.989778
station_og_lat <- 76.188338

# Define coordinates for Arctic locations in Jones Sound and Nares Strait
arctic_coords <- data.frame(
    Loocation = c("BELCHER", "TALBOT", "SVERDRUP", "GRISE_FIORD", "JAKEMAN", 'STATION_OG'),
    longitude = c(belcher_lon, talbot_lon, sverdrup_lon, grise_lon, jakeman_lon, station_og_lon),
    latitude = c(belcher_lat, talbot_lat, sverdrup_lat, grise_lat, jakeman_lat, station_og_lat)
)

# Filter data for Arctic locations only
arctic_dab_data <- long_dab_count_summ3_mean %>%
    filter(Loocation %in% c("BELCHER", "TALBOT", "SVERDRUP", "GRISE_FIORD", "JAKEMAN","STATION_OG")) %>%
    left_join(arctic_coords, by = "Loocation") %>%
    # Offset dabC points slightly to the right
    mutate(longitude_adj = ifelse(dab_name == "dabC", longitude + .75, longitude)) %>% 
    mutate(longitude_adj = ifelse(dab_name == "dabA", longitude - .75, longitude)) %>% 
    mutate(Loocation = ifelse(Loocation == "STATION_OG", "STATION OG", Loocation)) %>% 
    mutate(Loocation = ifelse(Loocation == "GRISE_FIORD", "GRISE FIORD", Loocation)) 

# Get higher resolution world map data
world <- ne_countries(scale = "large", returnclass = "sf")

# Define map bounds for Canadian Arctic (Jones Sound and Nares Strait region)
map_bounds <- list(
    xmin = -90, xmax = -75,
    ymin = 75.5, ymax = 78.2
)

# Add Auyuittuq coordinates
auyuittuq_coords <- data.frame(
    Loocation = "AUSUITTUQ",
    longitude = -82.8945,
    latitude = 76.4185
)

# Add Qikiqtaaluk coordinates
qikiqtaaluk_coords <- data.frame(
    Loocation = "QIKIQTAALUK",
    longitude = -87,
    latitude = 76.25
)


# Create grid positions for each location
grid_offset <- 0.3
arctic_dab_data_grid <- arctic_dab_data %>%
    mutate(
        # Create grid positions: DabA/DabC as rows, 2019/2021 as columns
        x_offset = case_when(
            year == 2019 ~ -grid_offset/2,
            year == 2021 ~ grid_offset/2
        ),
        y_offset = case_when(
            dab_name == "DabA" ~ grid_offset/8,
            dab_name == "DabB" ~ -grid_offset/8,
            dab_name == "DabC" ~ -grid_offset/8*3,
            dab_name == "DabD" ~ -grid_offset/8*5
        ),
        longitude_grid = longitude + x_offset,
        latitude_grid = latitude + y_offset
    )

# Add dummy rows for missing data
dummy_rows <- data.frame(
    # GRISE FIORD 2021 - DabA, DabB, DabC, DabD
    Loocation = rep("GRISE FIORD", 4),
    year = rep(2021, 4),
    dab_name = c("DabA", "DabB", "DabC", "DabD"),
    MEAN_DABs = rep(NA, 4),
    longitude = rep(grise_lon, 4),
    latitude = rep(grise_lat, 4),
    longitude_adj = rep(grise_lon, 4),
    longitude_grid = rep(grise_lon + grid_offset/2, 4),
    latitude_grid = c(grise_lat + grid_offset/8, grise_lat - grid_offset/8, grise_lat - grid_offset/8*3, grise_lat - grid_offset/8*5),
    x_offset = rep(grid_offset/2, 4),
    y_offset = c(grid_offset/8, -grid_offset/8, -grid_offset/8*3, -grid_offset/8*5)
) %>%
rbind(data.frame(
    # JAKEMAN 2021 - DabA, DabB, DabC, DabD  
    Loocation = rep("JAKEMAN", 4),
    year = rep(2021, 4),
    dab_name = c("DabA", "DabB", "DabC", "DabD"),
    MEAN_DABs = rep(NA, 4),
    longitude = rep(jakeman_lon, 4),
    latitude = rep(jakeman_lat, 4),
    longitude_adj = rep(jakeman_lon, 4),
    longitude_grid = rep(jakeman_lon + grid_offset/2, 4),
    latitude_grid = c(jakeman_lat + grid_offset/8, jakeman_lat - grid_offset/8, jakeman_lat - grid_offset/8*3, jakeman_lat - grid_offset/8*5),
    x_offset = rep(grid_offset/2, 4),
    y_offset = c(grid_offset/8, -grid_offset/8, -grid_offset/8*3, -grid_offset/8*5)
)) %>%
rbind(data.frame(
    # TALBOT 2019 - DabA, DabB, DabC, DabD
    Loocation = rep("TALBOT", 4),
    year = rep(2019, 4),
    dab_name = c("DabA", "DabB", "DabC", "DabD"),
    MEAN_DABs = rep(NA, 4),
    longitude = rep(talbot_lon, 4),
    latitude = rep(talbot_lat, 4),
    longitude_adj = rep(talbot_lon, 4),
    longitude_grid = rep(talbot_lon - grid_offset/2, 4),
    latitude_grid = c(talbot_lat + grid_offset/8, talbot_lat - grid_offset/8, talbot_lat - grid_offset/8*3, talbot_lat - grid_offset/8*5),
    x_offset = rep(-grid_offset/2, 4),
    y_offset = c(grid_offset/8, -grid_offset/8, -grid_offset/8*3, -grid_offset/8*5)
)) %>%
rbind(data.frame(
    # STATION OG (SSW) 2019 - DabA, DabB, DabC, DabD
    Loocation = rep("STATION OG", 4),
    year = rep(2019, 4),
    dab_name = c("DabA", "DabB", "DabC", "DabD"),
    MEAN_DABs = rep(NA, 4),
    longitude = rep(station_og_lon, 4),
    latitude = rep(station_og_lat, 4),
    longitude_adj = rep(station_og_lon, 4),
    longitude_grid = rep(station_og_lon - grid_offset/2, 4),
    latitude_grid = c(station_og_lat + grid_offset/8, station_og_lat - grid_offset/8, station_og_lat - grid_offset/8*3, station_og_lat - grid_offset/8*5),
    x_offset = rep(-grid_offset/2, 4),
    y_offset = c(grid_offset/8, -grid_offset/8, -grid_offset/8*3, -grid_offset/8*5)
))

# Combine original data with dummy rows
arctic_dab_data_grid <- rbind(arctic_dab_data_grid %>% select(-SD_DABs), dummy_rows)

names(arctic_dab_data_grid)

# ==============================================================================
# CUSTOM POLYGON SYMBOLS FOR MAP VISUALIZATION
# ============================================================================


# We'll keep other locations as square points; compute color range from point data
value_min <- min(arctic_dab_data_grid$MEAN_DABs, na.rm = TRUE)
value_max <- max(arctic_dab_data_grid$MEAN_DABs, na.rm = TRUE)



# Create map with averaged DAB values (one square per location-DAB combination)
arctic_map_averaged <- ggplot() +
    geom_sf(data = world, fill = "lightgray", color = "white", size = 0.2) +
    coord_sf(xlim = c(map_bounds$xmin, map_bounds$xmax), 
             ylim = c(map_bounds$ymin, map_bounds$ymax),
             expand = FALSE, crs = st_crs(4326)) +
    # Add squares colored by averaged DAB abundance
    geom_point(data = arctic_dab_data_grid, 
               aes(x = longitude_grid, y = latitude_grid, fill = MEAN_DABs), 
               shape = 22, alpha = 1, size = 8, color = "black", stroke = 0.5) +
    # Add location labels
    geom_text(data = distinct(arctic_dab_data_grid, Loocation, longitude, latitude), 
              aes(x = longitude, y = latitude + 0.1, label = Loocation),
              size = 4, fontface = "bold", color = "black") +
    # Add Auyuittuq star
    geom_point(data = auyuittuq_coords,
               aes(x = longitude, y = latitude),
               shape = 21, size = 4, color = "black", fill = "black") +
    # Add Auyuittuq label
    geom_label(data = auyuittuq_coords,
               aes(x = longitude-.99, y = latitude+.00, label = Loocation),
               size = 4, fontface = "bold", color = "black", 
               fill = "white", label.size = 0.25) +
    geom_label(data = qikiqtaaluk_coords,
               aes(x = longitude+.4, y = latitude - 0.2, label = Loocation),
               size = 6, fontface = "bold", color = "black", 
               fill = "white", label.size = 1) +
    ggspatial::annotation_scale(
        location = "br",
        width_hint = 0.25,
        text_cex = 1,
        pad_x = grid::unit(0.25, "cm"),
        pad_y = grid::unit(0.25, "cm")
    ) +
    scale_fill_gradient(low = "white", high = "purple", name = "Domoic acid production quantity\n(Averaged across years)",
                       labels = scales::trans_format("identity", scales::scientific_format()),
                       na.value = "grey") +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14),
          legend.title = element_text(vjust = 0.8, size = 10),
          legend.text = element_text(size = 7),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = 'bottom') +
    labs(title = "",
         x = "Longitude",
         y = "Latitude")

print(arctic_map_averaged)

ggsave(arctic_map_averaged, 
       filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/Dab_fig1_map_dab_averaged_260526.pdf", 
       width = 2300, height = 1490, units = "px", dpi = 150, device = "pdf")
