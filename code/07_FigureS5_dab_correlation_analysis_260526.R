# DAB Transcripts and Environmental Data Integration - PCA Analysis
# Author: Patrick White
# Date: August 29, 2025
# Description: Combines comprehensive DAB and Pseudo-nitzschia transcript data
#              with environmental variables for correlation analysis
# ===============================================================================

# Load all required libraries
library(dplyr)        # Data manipulation
library(tidyr)        # Data reshaping
library(cowplot)      # Plot arrangements
library(ggplot2)      # Plotting
library(stringr)      # String manipulation
library(corrplot)      # Correlation plots

# NOTE: This script assumes comprehensive_data_meta has been created
# by running dab_metadata_creation_09032025.R first

# Read combined metadata + DAB/PN TPMs (assumes this file exists)
comprehensive_data_meta = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')

comprehensive_data_meta$dab_name <- gsub("Dab", "dab", comprehensive_data_meta$dab_name)

comprehensive_data_meta = comprehensive_data_meta %>% select(sample, Loocation, dab_name, dab_to_pn_ratio, count_dab, Depth, SiO4, PO4, NO2, NO3, NH3, Temperature, Salinity, Dissolved.O2, Turbidity, Chlorophyl.a, Mn55.MR., Fe56.MR., Cd111.LR.)

# convert each numeric column to z score normailzed
comprehensive_data_meta <- comprehensive_data_meta %>%
  mutate(across(where(is.numeric), ~ (.-mean(., na.rm=TRUE))/sd(., na.rm=TRUE)))

# Pivot comprehensive_data_meta so each DAB becomes its own set of columns:
# - dab_count_<DAB>        : raw count
# - dab_count_all_<DAB>    : count_all (if present)
# - dab_to_pn_ratio_<DAB>  : dab_to_pn_ratio (if present)
# Clean DAB names to safe column names and fill missing counts with 0.


dab_ratio_wide <- comprehensive_data_meta %>%
  select(sample, dab_name, dab_to_pn_ratio) %>%
  pivot_wider(
    names_from = dab_name,
    values_from = dab_to_pn_ratio,
    values_fill = NA,
    names_glue = "{dab_name} normalized to PN"  # suffix instead of prefix
  )

names(comprehensive_data_meta)
names(dab_ratio_wide)

# Environmental metadata (one-row-per-sample)
env_data <- comprehensive_data_meta %>%
  select(-Loocation, -count_dab,-dab_name, -dab_to_pn_ratio,  -Depth) %>%
  distinct(sample, .keep_all = TRUE)

# Combine environmental and DAB wide data
total_data <- env_data %>%
  inner_join(dab_ratio_wide, by = "sample")

names(total_data)

numeric_vars = total_data %>%
  select(where(is.numeric)) %>%
  na.omit()

# Select numeric variables only; safely drop known numeric metadata if present
tm_data <- total_data %>%
  select(where(is.numeric)) %>%
  select(Cd111.LR., Mn55.MR., Fe56.MR., contains('to PN')) %>%
  na.omit()

names(tm_data)

plot(tm_data)

  # Hard-coded desired ordering for metals + PN-normalized DAB ratios only
  desired_order <- c(
    "Mn55.MR.", "Fe56.MR.", "Cd111.LR.",
    "dabA normalized to PN", "dabB normalized to PN", "dabC normalized to PN", "dabD normalized to PN",
    "dab_PN_ratio_DabA", "dab_PN_ratio_DabB", "dab_PN_ratio_DabC", "dab_PN_ratio_DabD"
  )

  # Keep only variables that actually exist in tm_data (preserve requested order)
  desired_order <- intersect(desired_order, names(tm_data))

  # validate and reorder
  numeric_vars_ord <- tm_data[, desired_order, drop = FALSE]

  # correlation matrix & p-values (pairwise complete)
  cor_mat <- cor(numeric_vars_ord, method = "spearman")
  res1 <- cor.mtest(numeric_vars_ord, conf.level = .95)

  # palette: blue for negative, red for positive
  col_palette <- colorRampPalette(c("blue", "white", "red"))(200)

# reusable function to create ggplot correlation tile plot
make_corr_ggplot <- function(cor_mat, p_mat = NULL, order = colnames(cor_mat),
                             upper = TRUE, show_all_labels = FALSE,
                             show_signif = FALSE, sig_level = 0.05,
                             palette = col_palette, label_size =  2) {
   vars <- order
   corr_df_full <- expand.grid(Var1 = vars, Var2 = vars, stringsAsFactors = FALSE)
   corr_df_full$cor <- as.vector(cor_mat)
   if (!is.null(p_mat)) corr_df_full$p <- as.vector(p_mat) else corr_df_full$p <- NA
   corr_df_full$Var1f <- factor(corr_df_full$Var1, levels = vars)
   corr_df_full$Var2f <- factor(corr_df_full$Var2, levels = vars)
   if (upper) {
     corr_df <- corr_df_full[as.numeric(corr_df_full$Var1f) <= as.numeric(corr_df_full$Var2f), ]
   } else {
     corr_df <- corr_df_full
   }

   p <- ggplot2::ggplot(corr_df, ggplot2::aes(x = Var2f, y = Var1f, fill = cor)) +
     ggplot2::geom_tile(color = "white") +
     ggplot2::scale_fill_gradientn(colors = palette, limits = c(-1, 1), na.value = "grey90", name = "Spearman rho") +
     ggplot2::coord_fixed() +
     ggplot2::theme_minimal() +
     ggplot2::theme(
       axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
       axis.title = ggplot2::element_blank(),
       panel.grid = ggplot2::element_blank()
     )

   if (show_all_labels) {
     p <- p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", cor)), color = "black", size = label_size)
   } else if (show_signif && !is.null(p_mat)) {
     p <- p + ggplot2::geom_text(ggplot2::aes(label = ifelse(p <= sig_level, sprintf("%.2f", cor), "")),
                                 color = "black", size = label_size)
   }

   return(p)
 }

# create cor1: upper triangle with all coefficients shown
cor1 <- make_corr_ggplot(cor_mat, p_mat = res1$p, order = colnames(cor_mat),
                         upper = TRUE, show_all_labels = FALSE, show_signif = TRUE,
                         sig_level = 0.05, palette = col_palette, label_size = 2)

## without metals

not_tm_vars <- total_data %>%
  select(-Fe56.MR., -Cd111.LR., -Mn55.MR.) %>%
  select(where(is.numeric)) %>%
  na.omit()

names(not_tm_vars)


  # Replace these names with the exact column names you want in that order
  desired_order <- c(
    # environmental (Depth removed as metals/Depth excluded in this section)
    "SiO4", "PO4", "NO2", "NO3", "NH3",
    "Temperature", "Salinity", "Dissolved.O2", "Turbidity", "Chlorophyl.a",
    # DAB TPM column name variants (created with names_glue = "{dab_name} (TPM)")
    "dabA (TPM)", "dabB (TPM)", "dabC (TPM)", "dabD (TPM)",
    # raw DAB count variants (if present)
    "dab_count_DabA", "dab_count_DabB", "dab_count_DabC", "dab_count_DabD",
    # pn count and PN-normalized DAB ratios (names from earlier pivots)
    "PN abundance (TPM)",
    "dabA normalized to PN", "dabB normalized to PN", "dabC normalized to PN", "dabD normalized to PN"
  )

  # Keep only variables that actually exist in numeric_vars (preserve requested order)
  desired_order <- intersect(desired_order, names(not_tm_vars))

  # validate and reorder
  numeric_vars_ord <- not_tm_vars[, desired_order, drop = FALSE]
  numeric_vars_ord = numeric_vars_ord %>% select(-contains("TPM")) 

  # correlation matrix & p-values (pairwise complete)
  cor_mat <- cor(numeric_vars_ord, method = "spearman")
  res1 <- cor.mtest(numeric_vars_ord, conf.level = .95)

  # palette: blue for negative, red for positive
  col_palette <- colorRampPalette(c("blue", "white", "red"))(200)

# create cor2: upper triangle with labels only for significant correlations
cor2 <- make_corr_ggplot(cor_mat, p_mat = res1$p, order = colnames(cor_mat),
                         upper = TRUE, show_all_labels = FALSE, show_signif = TRUE,
                         sig_level = 0.05, palette = col_palette, label_size = 2)


gridded <- cowplot::plot_grid(cor1, cor2, ncol = 2, labels = c("A", "B"))

gridded

ggsave(
  gridded,
  filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/corellation_supp_PN_dab_macros_CTD_260113.png",
  device = 'png', width = 3000, height = 3000, units = 'px'
)
