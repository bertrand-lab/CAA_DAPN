# ==============================================================================
# DAB Transcript Analysis in Arctic Samples
# Author: Patrick White
# Date: August 11, 2025
# Description: Analysis of domoic acid biosynthesis (DAB) transcripts and 
#              Pseudo-nitzschia abundance across Arctic locations
# ==============================================================================

# Load minimal required libraries
library(dplyr)
library(ggplot2)
library(stringr)

# Color palette used in the original correlation plots
c25 =  c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#CC79A7",  # reddish purple
  "#D55E00", # vermillion
  "#0072B2" # blue
)


# Read data (same file used previously)
comprehensive_data = read.csv('/Users/patrickwhite/Documents/GitHub/CAA_DAPN/formatting_metadata/env_metadata_and_dabs_TPM_TMs_260113_filtered.csv')


# Optional: ensure Loocation is a factor (keeps consistent coloring / ordering if needed)
comprehensive_data$Loocation <- factor(comprehensive_data$Loocation, 
                                          levels = c("BELCHER", "TALBOT", "SVERDRUP", "STATION_OG",
                                                     'GRISE_FIORD', 'JAKEMAN'))

 comprehensive_data$dab_name <- gsub("Dab", "dab", comprehensive_data$dab_name)

# ==============================================================================
# CORRELATION ANALYSIS
# ==============================================================================


# split the data into a list of data.frames, one per dab_name and run lm on each sequentially
df_by_dab <- split(comprehensive_data, comprehensive_data$dab_name)

# convert the lapply version to an explicit for-loop
names_dabs <- names(df_by_dab)
lm_list <- vector("list", length(names_dabs))

for (i in seq_along(names_dabs)) {
  nm <- names_dabs[i]
  sub <- df_by_dab[[nm]]
  
  m <- lm(count_dab ~ count_pn, data = sub)
  s <- summary(m)
  coefm <- coef(m)
  pval <- s$coefficients["count_pn", "Pr(>|t|)"]
  
  lm_list[[i]] <- data.frame(
    dab_name = nm,
    intercept = coefm[1],
    slope = coefm[2],
    r.squared = s$r.squared,
    p.value = pval,
    stringsAsFactors = FALSE
  )
}

lm_stats <- do.call(rbind, lm_list)

lm_stats$eq_label <- paste0(
  "y = ", signif(lm_stats$slope, 3),
  " x + ", signif(lm_stats$intercept, 3),
  "\nR² = ", signif(lm_stats$r.squared, 3)
)

comprehensive_data$Loocation = str_replace_all(comprehensive_data$Loocation, '_', ' ')

# df_by_dab is available if you want to "send" each filtered dataframe elsewhere:
# names(df_by_dab) gives the dab_names, df_by_dab[[nm]] is the dataframe for that dab.
# plot with regression line and per-facet equation / R^2 annotation
p_corr2 = ggplot(comprehensive_data, aes(x = count_pn, y = count_dab, fill = Loocation, group = dab_name)) +
  geom_point(size = 5, alpha = .7 , shape = 21, color = 'black') +
  labs(
    x = "Total pseudo-nitzschia transcripts (TPM)", 
    y = "Dab transcripts (TPM)",
    fill = "Site"
  ) +
  theme_bw() +
  # don't add a legend entry for the smooth (removes line from legend)
  geom_smooth(method = "lm", se = FALSE, color = "black", show.legend = FALSE) +
  theme(
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = 'black'),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) +
  scale_fill_manual(values = c25) +
  # ensure legend keys show points only (no line)
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 21))) +
  facet_wrap(~dab_name, scales = 'free') +
  geom_text(
    data = lm_stats,
    aes(x = -Inf, y = Inf, label = eq_label),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.1,
    size = 4
  )
p_corr2

# lm_stats contains the slope, intercept, R^2 and p-value for each dab_name
print(lm_stats)


ggsave(
  p_corr2,
  filename = "/Users/patrickwhite/Documents/GitHub/CAA_DAPN/plots/dab_vs_pn_correlation_filtered_260113.pdf",
  device = 'pdf', width = 3000, height = 2000, units = 'px'
)

