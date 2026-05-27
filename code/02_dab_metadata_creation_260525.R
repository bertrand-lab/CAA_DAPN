# ===============================================================================
# this code is for taking the different environmental datasets from 2019 and 2021, combining them, and then linking them to the comprehensive dataset of Dab transcripts.
# ===============================================================================

# Load all required libraries
library(dplyr)        # Data manipulation
library(stringr)      # String manipulation
library(tidyr)        # Data reshaping


# Template environmental data structure
env_data_2021 = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/env_metadata_2021_white2025.csv') ## this dataset is from the data file 'devon_metadata_microbial_2021_5m_meaned_complete_06_05_2024.2.csv' https://github.com/Patrick-L-White/CAA_PPCC/blob/main/data/devon_metadata_microbial_2021_5m_meaned_complete_06_05_2024.2.csv Publisdhed in the data repository for the paper White et al. (2025) in ISME communications.

env_data_2019 = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/env_metadata_2019_bhatia2021.csv') ## taken from the supplemental data from the paper Bhatia et al. (2021) in Global Biogeochemical Cycles, which is available here: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021GB006976 and specifically: https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2021GB006976&file=2021GB006976-sup-0002-Table+SI-S03.pdf


env_data_2021_filt = env_data_2021 %>% select(SAMPLE_ID, Depth, Distance_from_glacier, Latitude, Longitude, 
             SiO4, PO4, NO2, NO3, NH3, Temperature, Salinity, 
             Dissolved.O2, Turbidity, Chlorophyl.a)%>% 
    mutate(year = 2021)

env_data_2019_filt = env_data_2019 %>% select(CTD_ID, Depth_m, Distance, Latitude, Longitude, 
  SiO4_uM, PO4_uM, NO2_uM, NO3_uM, NH3_uM, Temp_C, Salinity_PSU, 
  Dis_O2_Conc, Turbidity_NTU, Chl_A_ug_L) %>%
  rename(
    SAMPLE_ID = CTD_ID,
    Depth = Depth_m,
    Distance_from_glacier = Distance,
    SiO4 = SiO4_uM,
    PO4 = PO4_uM,
    NO2 = NO2_uM,
    NO3 = NO3_uM,
    NH3 = NH3_uM,
    Temperature = Temp_C,
    Salinity = Salinity_PSU,
    Dissolved.O2 = Dis_O2_Conc,
    Turbidity = Turbidity_NTU,
    Chlorophyl.a = Chl_A_ug_L) %>% 
    mutate(year = 2019)

## file that links the Dab trancript IDs with the sample IDs from the environmental samples. These were created from communications with authors on the Bhatia et al. (2021) paper and the White et al. (2025) paper.
link2019_ids = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/env_metadata_sample_IDs_link_2019_and_2021.csv')

link2019_ids$link_ID[!link2019_ids$link_ID %in% env_data_2019_filt$SAMPLE_ID ]
env_data_2019_filt = env_data_2019_filt[env_data_2019_filt$SAMPLE_ID %in% link2019_ids$link_ID, ]
link2019_ids = link2019_ids[link2019_ids$link_ID %in% env_data_2019_filt$SAMPLE_ID, ]
link2019_ids= link2019_ids %>% distinct()

env_data_2019_filt = left_join(env_data_2019_filt, link2019_ids %>% select(link_ID, sample_id), by = c('SAMPLE_ID' = 'link_ID'))

env_data_2019_filt = env_data_2019_filt %>% mutate(SAMPLE_ID = sample_id) %>% select(-sample_id)

env_meta_combined = rbind(env_data_2019_filt, env_data_2021_filt)

# Manually add environmental data for SSW sample
ssw_env <- data.frame(
  SAMPLE_ID = "SSW",
  Depth = 5,
  Distance_from_glacier = NA,  # Not provided
  Latitude = NA,  # Not provided
  Longitude = NA,  # Not provided
  SiO4 = 1.57,
  PO4 = 0.42,
  NO2 = 0.04,
  NO3 = 0.22,
  NH3 = 0.23,
  Temperature = 2.58977318,
  Salinity = 31.0710169,
  Dissolved.O2 = 358.269372,
  Turbidity = 0.28344727,
  Chlorophyl.a = 0.13052188,
  year = 2021
)


# Add SSW_B data to the combined environmental data
env_meta_combined = rbind(env_meta_combined, ssw_env)

# Load the comprehensive dataset created in the previous script
comprehensive_data <- read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/arctic_comprehensive_abundanceTPM_summary_dabs_260113_filtered.csv')

# id processing to match across datasets
env_meta_combined$SAMPLE_ID = str_replace_all(env_meta_combined$SAMPLE_ID, '-', '_')
env_meta_combined$SAMPLE_ID = str_replace_all(env_meta_combined$SAMPLE_ID, ' ', '')

comprehensive_data$sample = str_remove_all(comprehensive_data$sample, 'VSM21-')
comprehensive_data$sample = str_replace_all(comprehensive_data$sample, '-', '_')
comprehensive_data$sample = str_replace_all(comprehensive_data$sample, ' ', '')
comprehensive_data$sample = str_replace_all(comprehensive_data$sample, '1_A', '1')
comprehensive_data$sample = str_replace_all(comprehensive_data$sample, '2_A', '2')
comprehensive_data$sample = str_replace_all(comprehensive_data$sample, '3_A', '3')

comprehensive_data_meta = left_join(comprehensive_data, env_meta_combined, by = c("sample" = "SAMPLE_ID"))

tms = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')

comprehensive_data_meta_tms = left_join(comprehensive_data_meta, tms, by = 'sample') 

comprehensive_data_meta_tms = comprehensive_data_meta_tms %>% 
  mutate(sample_group = ifelse(sample == "SSW", "2021", sample_group)) 

comprehensive_data_meta_tms = comprehensive_data_meta_tms %>% select(-Distance_from_glacier)

comprehensive_data_meta = comprehensive_data_meta %>% select(-Distance_from_glacier)

## mean across the duplicated SSW filter
comprehensive_data_meta_tms_ssw  = comprehensive_data_meta_tms %>%
  filter(grepl('SSW', sample)) %>%
  group_by(dab_name) %>%
  summarise(across(count_dab, mean), across(-count_dab, first), .groups = "drop") 

comprehensive_data_meta_tms_not_sww = comprehensive_data_meta_tms  %>% filter(!grepl('SSW', sample))

comprehensive_data_meta_tms = rbind(comprehensive_data_meta_tms_not_sww, comprehensive_data_meta_tms_ssw)

comprehensive_data_meta_tms$Loocation = ifelse(comprehensive_data_meta_tms$Loocation == "SSW", "STATION_OG", comprehensive_data_meta_tms$Loocation)

comprehensive_data_meta_tms = comprehensive_data_meta_tms %>% filter(Depth < 50)

#write.csv(comprehensive_data_meta, file = '/Volumes/CORSAIR/Patrick_W/Transcriptomics/env_metadata/env_metadata_and_dabs_TPM_260113.csv', row.names = FALSE, quote = TRUE)

write.csv(comprehensive_data_meta_tms, file = '/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv', row.names = FALSE)
