# ===============================================================================
# DAB Sampling Locations - Supplementary Map
# Author: Patrick White
# Description: Builds Jones Sound and Nares Strait station maps from the
#              environmental metadata table, labels stations, and outputs a
#              combined figure with a shared year legend.
# ===============================================================================

# Load required libraries
library(dplyr)        # Data manipulation
library(stringr)      # String manipulation
library(tidyr)        # Data reshaping
library(ggplot2)      # Plotting
library(sf)           # Spatial data handling
library(rnaturalearth) # World map data
library(cowplot)      # Plot arrangements
library(ggrepel)      # Improved label placement to avoid overlaps
library(ggspatial)    # Scale bars and north arrows

c25 = c(
  "dodgerblue2", "#E31A1C")

# Load the metadata table used for station locations and years
comprehensive_data_meta = read.csv('/Users/patrickwhite/Bertrand Lab Dropbox/Bertrand Lab shared workspace/User/Patrick_White_DB/corsair_drive/Patrick_W/Transcriptomics/dab_paper_writing/review1/code/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')

# Get higher resolution world map data
world <- ne_countries(scale = "large", country = "Canada", returnclass = "sf")

# Define map bounds for the Jones Sound plot
map_bounds <- list(
    xmin = -84, xmax = -80,
    ymin = 75.6, ymax = 76.7
)

# Prepare station table used for plotting (unique station-year locations)
station_tbl <- comprehensive_data_meta %>%
  mutate(station_id = ifelse(grepl('VIO', sample), str_sub(sample, 1, 6), str_sub(sample, 1, 3))) %>%
  select(Longitude, Latitude, Loocation, year, station_id) %>%
  distinct() %>% 
  mutate(Longitude = ifelse(Loocation == 'STATION_OG', -82.922133, Longitude)) %>% 
  mutate(Latitude = ifelse(Loocation == 'STATION_OG', 76.379117, Latitude))

# Helper to build a station map for a given bounding box
make_map <- function(bounds, title = NULL, df) {
  # compute year-based scales so legend can combine shape + fill
  yrs <- as.character(sort(unique(df$year)))
  n_yrs <- length(yrs)
  # shapes that support fill
  shape_palette <- rep(c(21,22,23,24,25), length.out = n_yrs)
  # fill colors from existing palette (c25)
  fill_palette <- c25[seq_len(n_yrs)]
  shape_vals <- setNames(shape_palette, yrs)
  fill_vals  <- setNames(fill_palette, yrs)

  ggplot() +
    geom_sf(data = world, fill = "lightgray", color = "white", size = 0.2) +
    coord_sf(xlim = c(bounds$xmin, bounds$xmax),
             ylim = c(bounds$ymin, bounds$ymax),
             expand = FALSE, crs = st_crs(4326)) +
    geom_point(data = df,
               aes(x = Longitude, y = Latitude, shape = as.character(year), fill = as.character(year)),
               size = 4, color = "black", stroke = 0.5, alpha = 0.75) +
    # Use ggrepel to offset labels and draw connecting segments to the points.
    geom_label_repel(
      data = df,
      aes(x = Longitude, y = Latitude, label = station_id),
      size = 4,
      color = "black",
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = "gray40",
      segment.size = 0.3,
      force = 1,
      max.iter = 2000,
      seed = 42,
      show.legend = FALSE
    ) +
    # combined legend: use shape scale to produce the legend and override fill in legend keys,
    # hide the separate fill legend so only one "year" legend remains
    scale_shape_manual(name = "year", values = shape_vals,
                       guide = guide_legend(override.aes = list(fill = unname(fill_vals), colour = "black"))) +
    scale_fill_manual(name = "year", values = fill_vals, guide = "none") +
    annotation_scale(
      location = "tr",
      width_hint = 0.25,
      text_cex = 1.2,
      style = "ticks"
    ) +
    theme_minimal() +
    # Add a black border around the plotting panel so each map has a clear outline.
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.6),
      # keep the panel background neutral (no fill) so geom_sf fill is visible
      panel.background = element_rect(fill = NA, colour = NA),
      # increased font sizes for axes, legend, and title
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 13),
      legend.text  = element_text(size = 12),
      legend.key.size = unit(1.2, "lines"),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    labs(title = title)
}

# Compute Talbot bounds from data (fallback to small buffer around dataset median)
talbot_df <- comprehensive_data_meta %>% filter(str_detect(Loocation, regex("TALBOT", ignore_case = TRUE)))
if (nrow(talbot_df) > 0) {
  buffer <- 1
  talbot_bounds <- list(
    xmin = min(talbot_df$Longitude, na.rm = TRUE) - buffer*2,
    xmax = max(talbot_df$Longitude, na.rm = TRUE) + buffer*2,
    ymin = min(talbot_df$Latitude, na.rm = TRUE) - buffer/2,
    ymax = max(talbot_df$Latitude, na.rm = TRUE) + buffer/2
  )
} else {
  # fallback: small box around dataset median
  lon_m <- median(comprehensive_data_meta$Longitude, na.rm = TRUE)
  lat_m <- median(comprehensive_data_meta$Latitude, na.rm = TRUE)
  delta <- 0.5
  talbot_bounds <- list(xmin = lon_m - delta*2, xmax = lon_m + delta*2, ymin = lat_m - delta, ymax = lat_m + delta)
}

# Build maps: create one plot that contains the legend (from p_main), then extract that legend
p_main_legend <- make_map(
  map_bounds,
  title = "Jones Sound",
  df = station_tbl %>% filter(Loocation != 'TALBOT') %>% mutate(station_id = ifelse(Loocation == 'JAKEMAN', 'VIO_6', station_id))
) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 315, hjust = 0, size = 12)
  )

# Create a legend-free version of the main plot
p_main_nolegend <- p_main_legend + theme(legend.position = "none")

# Nares Strait plot without its own legend (we'll use the legend from p_main_legend)
p_talbot_nolegend <- make_map(
  talbot_bounds,
  title = "Nares Strait",
  df = station_tbl %>% filter(Loocation == 'TALBOT')
) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 315, hjust = 0, size = 12)
  )

# Extract legend from the p_main_legend object
legend_grob <- cowplot::get_legend(p_main_legend)

# Arrange the two plots and the shared legend (legend derived from p_main)
gridded <- plot_grid(
  p_main_nolegend, p_talbot_nolegend, legend_grob,
  ncol = 3, rel_widths = c(1, 0.94, 0.28)
)
gridded

# Save combined figure
ggsave(
  gridded,
  filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/dab_sampling_map_supp_260113.pdf",
  device = 'pdf', width = 3000, height = 3000, units = 'px'
)

