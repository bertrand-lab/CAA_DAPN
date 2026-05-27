# ==============================================================================
# DAB Transcript Analysis in Arctic Samples
# Author: Patrick White
# Date: August 11, 2025
# Description: Analysis of domoic acid biosynthesis (DAB) transcripts and 
#              Pseudo-nitzschia abundance across Arctic locations
# ==============================================================================

# Load required libraries
library(dplyr)        # Data manipulation
library(data.table)   # Fast data reading
library(stringr)      # String manipulation
library(tidyr)        # Data reshaping
library(ggplot2)      # Plotting
library(sf)           # Spatial data handling
library(rnaturalearth) # World map data
library(rnaturalearthdata) # Additional map data
library(cowplot)      # Plot arrangements

# ==============================================================================
# DATA TRANSFORMATION TO LONG FORMAT
# ==============================================================================

dabs_counts = fread('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/arctic_dabs_transcripts_only_raw_counts_without_multiple_annotations_260113.csv')

names(dabs_counts)

dabs_counts$dab_name 

# Convert wide format count data to long format for analysis
# Filter for only transcripts with DAB annotations
long_dab_count <- dabs_counts %>%
    filter(!is.na(dab_name)) %>%
    pivot_longer(cols = `VSM21-SV7-1_A`:`VSM21-SV1-2_A`, 
                 names_to = "sample", 
                 values_to = "count")

# Summarize DAB transcript counts by sample and gene type
long_dab_count_summ <- long_dab_count %>%
    group_by(sample, dab_name) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = 'drop')# %>%
    #mutate(sample = str_replace(sample, "SSW_B", "SSW")) %>%
    #group_by(sample, dab_name) %>%
    #summarise(count = mean(count, na.rm = TRUE), .groups = 'drop')



# ==============================================================================
# RELATIVE ABUNDANCE CALCULATIONS
# ==============================================================================

# Calculate total transcript counts per sample for normalization
df <- read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/arctic_dabs_total_counts_per_sample_260113.csv')

# Calculate relative abundance (proportion of total transcripts)
long_dab_count_summ2 <- left_join(long_dab_count_summ, df, by = "sample") %>%
    mutate(rel_abund = count / count_all) 


# ==============================================================================
# SAMPLE GROUPING AND METADATA INTEGRATION
# ==============================================================================

# Add sample grouping based on naming patterns
long_dab_count_summ2 <- long_dab_count_summ2 %>%
    mutate(sample_group = case_when(
        str_detect(sample, "VIO") ~ "2019",               # 2019 samples
        str_detect(sample, "SSW") ~ "2021",               # 2021 seawater samples
        TRUE ~ "2021"                                     # Default to 2021
    ))

# Extract Pseudo-nitzschia transcript data from the main dataset
long_pseudo_count <- fread('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/arctic_dabs_PN_only_raw_counts_without_multiple_annotations_260113.csv')

# Summarize total Pseudo-nitzschia counts per sample
long_pseudo_count_summ <- long_pseudo_count %>%
    group_by(sample) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
    mutate(classification_type = "Pseudo-nitzschia")
# ==============================================================================
# PSEUDO-NITZSCHIA RELATIVE ABUNDANCE CALCULATIONS
# ==============================================================================

# Calculate Pseudo-nitzschia percentage of total transcripts
pseudo_nitzschia_comparison <- long_pseudo_count_summ %>%
    left_join(df, by = "sample") %>%
    mutate(pn_percentage = (count / count_all) * 100,
           rel_abund = count / count_all) %>%
    select(sample, count, count_all, pn_percentage, rel_abund, classification_type) %>% 
    rename('count_pn' = 'count',
           'count_pn_total' = 'count_all',
           'pn_rel_abund' = 'rel_abund')

# Create comprehensive combined dataset with all data points
comprehensive_data <- long_dab_count_summ2 %>%
    # Add Pseudo-nitzschia data
    left_join(pseudo_nitzschia_comparison %>% select(sample, count_pn )
    , by = 'sample') %>%
    rename('count_dab' = 'count') %>%
    # Calculate DAB percentages and ratios
    mutate(
        dab_to_pn_ratio = count_dab / count_pn) %>% 
        select(-rel_abund, -count_all)

comprehensive_data

location_information = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/transcriptomic_location_information.csv')

location_information$SampleName = str_remove_all(location_information$SampleName, " ")
comprehensive_data_meta = left_join(comprehensive_data, location_information, by = c("sample" = "SampleName"))

comprehensive_data_meta$Loocation = ifelse(comprehensive_data_meta$sample == "SSW", "STATION_OG", comprehensive_data_meta$Loocation)
comprehensive_data_meta$Loocation = ifelse(comprehensive_data_meta$sample == "SSW_B", "STATION_OG", comprehensive_data_meta$Loocation)


# Write comprehensive dataset
write.csv(comprehensive_data_meta, 
          '/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/arctic_comprehensive_abundanceTPM_summary_dabs_260113_filtered.csv', 
          row.names = FALSE )
