# Load only used libraries
library(dplyr)        # Data manipulation
library(ggplot2)      # Plotting
library(cowplot)      # Plot arrangements
library(tidyr)        # Data reshaping
library(stringr)      # String manipulation (for str_split)
library(data.table)   # Efficient data reading
library(vegan)        # PERMANOVA analysis

common_font_theme <- theme_bw(base_size = 14) + theme(
  # axis text and titles (shared)
  axis.text.x = element_text(size = 14),   # match p2
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  
  # legend and plot title text sizing (shared defaults)
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  plot.title = element_text(size = 18, face = "bold")
)

# Download Brunson et al. 2024 supplementary data from PNAS. Cant be automatically downloaded with curl. 
brunson_url   = "https://www.pnas.org/doi/suppl/10.1073/pnas.2319177121/suppl_file/pnas.2319177121.sd01.csv"
brunson_local = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/brunson/pnas.2319177121.sd01.csv"

counts_pn_only = read.csv(brunson_local) ## this is from supplementary data of brunson et al. 2024 follow this link to downlaod it: https://www.pnas.org/doi/suppl/10.1073/pnas.2319177121/suppl_file/pnas.2319177121.sd01.csv


## this is from supplementary data of brunson et al. 2024 Table S3. FOllow this link to downlaod it and create yur own table. https://www.pnas.org/doi/suppl/10.1073/pnas.2319177121/suppl_file/pnas.2319177121.sapp.pdf. DabB was additionally identified by blasting against the brunson data with the sequence ID 'contig_393632_37_744_-'

dabs <- data.frame(
  orf_id = c(
    "contig_393632_37_744_-",
    "contig_392124_1_1473_+",
    "contig_469576_120_1547_-",
    "contig_519263_1_1377_-",
    "contig_390219_292_1467_-",
    "contig_392175_3_581_-",
    "contig_698932_2_493_-",
    "contig_519989_125_1249_-",
    "contig_383882_249_1943_-"
  ),
  dab = c(
    "dabB",
    "dabA",
    "dabA",
    "dabC",
    "dabC",
    "dabC",
    "dabC",
    "dabC",
    "dabD"
  ),
  stringsAsFactors = FALSE
)


# remove X and _MWII_R1 from column names
names(counts_pn_only)[34:85] = gsub("X", "", names(counts_pn_only)[34:85])
names(counts_pn_only)[34:85] = gsub("_MWII_R1", "", names(counts_pn_only)[34:85])

# Extract column names and convert to Date format to identify out-of-range samples
sample_dates <- as.Date(names(counts_pn_only)[34:85], format="%m%d%y")
out_of_range_samples <- names(counts_pn_only)[34:85][
  sample_dates < as.Date("040115", format="%m%d%y") | 
    sample_dates > as.Date("100115", format="%m%d%y")
]
counts_pn_only_filtered = counts_pn_only

counts = counts_pn_only_filtered
# drop metadata columns, keep sample count columns
counts = counts %>% select(-(orf_contam_type:transmembrane_domains))

# map ORFs to DAB IDs and compute per-sample totals
dabs_counts = left_join(counts, dabs, by = c('orf_id' = 'orf_id'))

dabs_counts_colsums = colSums(dabs_counts  %>% select(-orf_id, -dab))
dabs_counts_colsums = as.data.frame(dabs_counts_colsums)
dabs_counts_colsums$sample = rownames(dabs_counts_colsums)

# aggregate DAB counts by sample and DAB
dabs_counts_summed = dabs_counts %>% filter(!is.na(dab))  %>% 
  pivot_longer(cols = `011415`:`123014`,
               names_to = "sample",
               values_to = "count") %>% 
  group_by(sample, dab) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = 'drop')

# aggregate total counts per sample for normalization
counts_summed = dabs_counts %>% 
  pivot_longer(cols = `011415`:`123014`,
               names_to = "sample",
               values_to = "count") %>% 
  group_by(sample) %>%
  summarise(count = sum(count, na.rm = TRUE), .groups = 'drop')

# compute percent of reads per DAB
dab_percentages_brunson = left_join(dabs_counts_summed, counts_summed, by = 'sample') %>% 
  mutate(percent_dab = (count.x / count.y))  %>% 
  select(sample, dab, percent_dab)

dabs_wide_perc_brunson = dab_percentages_brunson %>% pivot_wider(names_from = dab, values_from = percent_dab) 

#==============================================================================
# calculate ratios among DABs in Brunson dataset
#==============================================================================
dabs_wide_perc_brunson$dab_a_to_d <- ifelse(dabs_wide_perc_brunson$dabD == 0, NA, dabs_wide_perc_brunson$dabA / dabs_wide_perc_brunson$dabD)
dabs_wide_perc_brunson$dab_b_to_d <- ifelse(dabs_wide_perc_brunson$dabD == 0, NA, dabs_wide_perc_brunson$dabB / dabs_wide_perc_brunson$dabD)
dabs_wide_perc_brunson$dab_c_to_d <- ifelse(dabs_wide_perc_brunson$dabD == 0, NA, dabs_wide_perc_brunson$dabC / dabs_wide_perc_brunson$dabD)
dabs_wide_perc_brunson$dab_a_to_b <- ifelse(dabs_wide_perc_brunson$dabB == 0, NA, dabs_wide_perc_brunson$dabA / dabs_wide_perc_brunson$dabB)
dabs_wide_perc_brunson$dab_a_to_c <- ifelse(dabs_wide_perc_brunson$dabC == 0, NA, dabs_wide_perc_brunson$dabA / dabs_wide_perc_brunson$dabC)
dabs_wide_perc_brunson$dab_b_to_c <- ifelse(dabs_wide_perc_brunson$dabC == 0, NA, dabs_wide_perc_brunson$dabB / dabs_wide_perc_brunson$dabC)

# reshape to long format for boxplots
ratio_data <- dabs_wide_perc_brunson %>%
  select(sample, dab_a_to_d, dab_b_to_d, dab_c_to_d, dab_a_to_b, dab_a_to_c, dab_b_to_c) %>%
  pivot_longer(cols = -sample, names_to = "ratio", values_to = "value")

comprehensive_data_meta = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')

# focus on stations used in this comparison
comprehensive_data_meta = comprehensive_data_meta %>% filter(Loocation %in% c('BELCHER','SVERDRUP','JAKEMAN', 'STATION_OG'))
# Convert comprehensive_data_meta to wide format like dabs_wide_perc_brunson
comprehensive_data_meta_wide <- comprehensive_data_meta %>%
  select(sample, dab_name, dab_to_pn_ratio) %>%
  pivot_wider(names_from = dab_name, values_from = dab_to_pn_ratio)

## perform dab ratios
comprehensive_data_meta_wide$dab_a_to_d <- ifelse(comprehensive_data_meta_wide$DabD == 0, NA, comprehensive_data_meta_wide$DabA / comprehensive_data_meta_wide$DabD)
comprehensive_data_meta_wide$dab_b_to_d <- ifelse(comprehensive_data_meta_wide$DabD == 0, NA, comprehensive_data_meta_wide$DabB / comprehensive_data_meta_wide$DabD)
comprehensive_data_meta_wide$dab_c_to_d <- ifelse(comprehensive_data_meta_wide$DabD == 0, NA, comprehensive_data_meta_wide$DabC / comprehensive_data_meta_wide$DabD)
comprehensive_data_meta_wide$dab_a_to_b <- ifelse(comprehensive_data_meta_wide$DabB == 0, NA, comprehensive_data_meta_wide$DabA / comprehensive_data_meta_wide$DabB)
comprehensive_data_meta_wide$dab_a_to_c <- ifelse(comprehensive_data_meta_wide$DabC == 0, NA, comprehensive_data_meta_wide$DabA / comprehensive_data_meta_wide$DabC)
comprehensive_data_meta_wide$dab_b_to_c <- ifelse(comprehensive_data_meta_wide$DabC == 0, NA, comprehensive_data_meta_wide$DabB / comprehensive_data_meta_wide$DabC)

# reshape to long format for boxplots
ratio_data_meta <- comprehensive_data_meta_wide %>%
  select(sample, dab_a_to_d, dab_b_to_d, dab_c_to_d, dab_a_to_b, dab_a_to_c, dab_b_to_c) %>%
  pivot_longer(cols = -sample, names_to = "ratio", values_to = "value") 

ratio_data$dataset = 'Bloom Brunson et al. (2024)'
ratio_data_meta$dataset = 'This study' 

combined_ratio_data = rbind(ratio_data, ratio_data_meta)

combined_ratio_data = combined_ratio_data %>% 
  mutate(ratio = gsub("_to_", " : ", ratio)) %>% 
  mutate(ratio = gsub("^dab_", "", ratio)) %>% 
  mutate(ratio = toupper(ratio)) %>%
  mutate(ratio = gsub("\\b([A-Za-z])\\b", "dab\\1", ratio, perl = TRUE))

# Convert to wide format for PCA (one row per sample)
combined_ratio_data_wide <- combined_ratio_data %>%
  pivot_wider(names_from = ratio, values_from = value, values_fill = NA) 

# boxplot of DAB ratios, colored by dataset and bloom status

combined_ratio_data$dataset = ifelse(combined_ratio_data$sample %in% out_of_range_samples & combined_ratio_data$dataset != 'This study', "Non bloom Brunson et al. (2024)", ifelse(combined_ratio_data$sample %in% out_of_range_samples & combined_ratio_data$dataset == 'This study', "Non bloom This study", combined_ratio_data$dataset))

boxplot_data <- combined_ratio_data %>%
  filter(!is.na(value))


p_ratios <- ggplot(boxplot_data, aes(x = ratio, y = value, fill = dataset)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), 
              alpha = 0.5, size = 2, shape = 21, color = 'black') +
  scale_fill_manual(values = c("Non bloom Brunson et al. (2024)" = "#CBA6FF","Bloom Brunson et al. (2024)" = "purple", "This study" = "#FFB6D9")) +
  labs(x = "", y = "Ratio Value", fill = "Dataset") +
  common_font_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap(~ratio, scales = "free", strip.position = "top") +
  theme(strip.background = element_rect(fill = "white", color = "black"))

print(p_ratios)


## gg save this 
ggsave(
  p_ratios,
  filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/Dab_figS2_dab_brunson_ratios_boxplot_260114.pdf",
  device = 'pdf', width = 12, height = 8
)


# filter out non-bloom Brunson samples for PCA comparison
combined_ratio_data_wide = combined_ratio_data_wide %>% filter(sample %in% out_of_range_samples == FALSE | dataset == 'This study')   %>% na.omit()

combined_ratio_data_wide$dataset = ifelse(combined_ratio_data_wide$dataset == 'Bloom Brunson et al. (2024)', "Brunson et al. (2024)", combined_ratio_data_wide$dataset)

# filter out non bloom brunson samples

# PCA without dabD (avoid dominance by a single DAB)
pca_result_no_dabD <- prcomp(combined_ratio_data_wide %>% select(-sample, -dataset, -contains("dabD")), scale. = TRUE)


pca_data_no_dabD <- as.data.frame(pca_result_no_dabD$x) %>%
  mutate(dataset = combined_ratio_data_wide$dataset)

pca_plot_no_dabD <- ggplot(pca_data_no_dabD, aes(x = PC1, y = PC2, fill = dataset)) +
  geom_point(size = 4, shape = 21, color = 'black', alpha = 0.7) +
  geom_segment(data = as.data.frame(pca_result_no_dabD$rotation[, 1:2] * 3),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", alpha = 0.7, inherit.aes = FALSE) +
  geom_text(data = as.data.frame(pca_result_no_dabD$rotation[, 1:2] * 3.2),
            aes(x = PC1, y = PC2, label = rownames(pca_result_no_dabD$rotation)),
            color = "black", fontface = "bold", size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Brunson et al. (2024)" = "#CBA6FF", "This study" = "#FFB6D9")) +
  labs(x = paste0("PC1 (", round(summary(pca_result_no_dabD)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result_no_dabD)$importance[2,2]*100, 1), "%)"),
       fill = "Dataset") +
  common_font_theme +
  theme(legend.position = "none")

# PCA with dabD included
pca_result_with_dabD <- prcomp(combined_ratio_data_wide %>% select(-sample, -dataset), scale. = TRUE)

pca_data_with_dabD <- as.data.frame(pca_result_with_dabD$x) %>%
  mutate(dataset = combined_ratio_data_wide$dataset)

pca_plot_with_dabD <- ggplot(pca_data_with_dabD, aes(x = PC1, y = PC2, fill = dataset)) +
  geom_point(size = 4, shape = 21, color = 'black', alpha = 0.7) +
  geom_segment(data = as.data.frame(pca_result_with_dabD$rotation[, 1:2] * 3),
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black", alpha = 0.7, inherit.aes = FALSE) +
  geom_text(data = as.data.frame(pca_result_with_dabD$rotation[, 1:2] * 3.2),
            aes(x = PC1, y = PC2, label = rownames(pca_result_with_dabD$rotation)),
            color = "black", fontface = "bold", size = 3, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Brunson et al. (2024)" = "#CBA6FF", "This study" = "#FFB6D9")) +
  labs(x = paste0("PC1 (", round(summary(pca_result_with_dabD)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result_with_dabD)$importance[2,2]*100, 1), "%)"),
       fill = "Dataset") +
  common_font_theme +
  theme(legend.position = "none")

# Combine PCA plots side by side
combined_pca_plot <- plot_grid(pca_plot_no_dabD, pca_plot_with_dabD,
                               ncol = 2, align = "hv")

print(combined_pca_plot)

# PERMANOVA without dabD (Euclidean on scaled ratios)
combined_ratio_data_no_dabD <- combined_ratio_data_wide %>% select(-sample, -contains("dabD"))
mcombined_ratio_data_no_dabD = as.matrix(combined_ratio_data_no_dabD %>% select(-dataset))
permanova_result_no_dabD <- adonis2(mcombined_ratio_data_no_dabD ~ dataset, data = combined_ratio_data_no_dabD, permutations = 999, method = "euclidean")
print("PERMANOVA results WITHOUT dabD:")
print(permanova_result_no_dabD)



# PERMANOVA with dabD
combined_ratio_datawithdabD <- combined_ratio_data_wide %>% select(-sample)
mcombined_ratio_datawithdabD = as.matrix(combined_ratio_datawithdabD %>% select(-dataset))
permanova_result_with_dabD <- adonis2(mcombined_ratio_datawithdabD ~ dataset, data = combined_ratio_datawithdabD, permutations = 999, method = "euclidean")
print("PERMANOVA results WITH dabD:")
print(permanova_result_with_dabD)

#==============================================================================
# Proportion of Pseudo-nitzschia reads comparison
#==============================================================================

euk_reads =  read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/eukaryotic_read_counts_summary.csv')

# Download Brunson et al. 2024 supplementary data from Zenodo
brunson_url <- "https://zenodo.org/records/10728894/files/annotation_all.filtered.tab.gz?download=1"
brunson_local <- "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/public_datasets/brunson/annotation_all.filtered.tab"


# use this to download brunson data if not already downloaded, or go to path specified and download manually if you have trouble with curl method in R

 if (!file.exists(brunson_local)) {
   download.file(brunson_url, destfile = brunson_local, method = "curl")
}

counts_pn_only = read.csv(brunson_local, sep = '\t') ## this is from supplementary data of brunson et al. 2024 https://www.pnas.org/doi/suppl/10.1073/pnas.2319177121/suppl_file/pnas.2319177121.sd01.csv

brunson_counts = read.csv(brunson_local, sep = '\t', header = TRUE) # this data is from brunson et al. 2024, https://zenodo.org/records/10728894/files/annotation_all.filtered.tab.gz?download=1

# tag PN vs non-PN for proportion calculation
brunson_counts$is_pn = grepl('Pseudo-nitzschia', brunson_counts$best_hit_species, ignore.case = TRUE)

counts_long = brunson_counts %>% 
pivot_longer(cols = contains('_MWII_R1'), names_to = 'sample', values_to = 'counts')

counts_summary = counts_long %>%
  group_by(sample, is_pn) %>%
  dplyr::summarise(total_counts = sum(counts)) %>%
  ungroup() %>%
  pivot_wider(names_from = is_pn, values_from = total_counts, values_fill = 0) %>%
  dplyr::rename(pn_counts = `TRUE`, non_pn_counts = `FALSE`) %>%
  mutate(total_counts = pn_counts + non_pn_counts,
         pn_proportion = pn_counts / total_counts)

# remove X and _MWII_R1 from column names
counts_summary$sample = gsub("X", "", counts_summary$sample)
counts_summary$sample = gsub("_MWII_R1", "", counts_summary$sample)

# keep only samples not in the out-of-range sample list
counts_summary_filtered = counts_summary %>% 
  filter(!sample %in% out_of_range_samples)

violin_plot <- ggplot() +
  geom_violin(data = counts_summary_filtered, aes(x = "Brunson", y = pn_proportion), fill = '#CBA6FF', alpha = 0.5) +
  geom_jitter(data = counts_summary_filtered, aes(x = "Brunson", y = pn_proportion), 
              width = 0.15, alpha = 0.7, size = 3, shape = 21, fill = '#CBA6FF', color = 'black') +
  geom_violin(data = euk_reads, aes(x = "Ours", y = pn_proportion), fill = '#FFB6D9', alpha = 0.5) +
  geom_jitter(data = euk_reads, aes(x = "Ours", y = pn_proportion), 
              width = 0.15, alpha = 0.7, size = 3, shape = 21, fill = '#FFB6D9', color = 'black') +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(y = 'Proportion of eukaryotic\nreads mapping to\nPseudo-nitzschia (%)', x = 'Dataset') +
  common_font_theme +
  theme(
    plot.margin = grid::unit(c(0, 0.2, 0, 0.7), "cm")
  )

# Rename to p2 for consistency with grid code
p2 <- violin_plot 
p2
# NOTE: This script assumes comprehensive_data_meta has been created
# by running dab_metadata_creation_09032025.R first

# shared legend for combined figure
shared_legend <- cowplot::get_legend(ggplot() +
  geom_point(aes(x = c("Brunson et al. (2024)", "This study"), 
                 y = c(1, 1),
                 fill = c("Brunson et al. (2024)", "This study")),
             shape = 21, size = 5, alpha = .7) +
             theme_minimal() +
  scale_fill_manual(values = c("Brunson et al. (2024)" = "#CBA6FF", "This study" = "#FFB6D9"), 
                    name = "Dataset") +
  theme(legend.position = "right")+ common_font_theme)  


p2_nolegend <- p2 + theme(legend.position = "none") + 
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12)
  )
p1_nolegend <- combined_pca_plot + theme(legend.position = "none") + theme(plot.margin = grid::unit(c(.2, 0.2, .2, 1.6), "cm")) # top, right, bottom, left

# arrange plots vertically and place shared legend to the right
plots <- cowplot::plot_grid(p2_nolegend, p1_nolegend, ncol = 1, labels = c('A', 'B'), label_size = 20,rel_heights = c(.7, 1))
gridded = cowplot::plot_grid(plots, shared_legend, rel_widths = c(.4, 0.15))

gridded

ggsave(
  gridded,
  filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/Dab_fig2_dab_brunson_ratios_260113.pdf",
  device = 'pdf', width = 3000, height = 2000, units = 'px'
)
