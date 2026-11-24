# plot_sample_qc.R
# Supplementary figure: quality control of the samples, every region.
#   a  Library size (millions of reads counted on genes) of every sample
#   b  Detected genes (CPM > 1) of every sample
#   c  Standardized connectivity of every sample; below the dashed line
#      (max_connectivity_z) the sample is an outlier and is removed
#   d  Distribution of the genes of every sample (kept samples, ordered by
#      median), before and after removing the technical covariates: line =
#      5-95% of the genes, thick line = 25-75%, point = median. After the
#      correction the samples should be aligned.
#   e  GC and length biases (NOISeq) of every sample, before and after the
#      EDASeq correction: R^2 of the expression explained by GC / length.
# Outlier samples (sample_qc) in red in a-c, with their ID in c.
# Run by Snakemake (rule plot_sample_qc).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  source(snakemake@input[["theme"]])
})

max_z <- snakemake@params[["max_connectivity_z"]]

read_tables <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(specimenID = col_character(),
                            individualID = col_character(),
                            batch = col_character(),
                            model_batch = col_character(),
                            apoe4_allele = col_character())) %>%
    bind_rows()
}

qc <- read_tables(snakemake@input[["sample_qc"]]) %>%
  set_pub_levels() %>%
  mutate(status = factor(if_else(outlier, "outlier", "kept"), levels = c("kept", "outlier")))
distribution <- read_tables(snakemake@input[["distribution"]]) %>%
  set_pub_levels() %>%
  mutate(stage = factor(stage, levels = c("before", "after"),
                        labels = c("Before correction", "After correction")))

cat("Outliers per region:\n")
print(as.data.frame(qc %>% group_by(region) %>%
                      summarise(samples = n(), outliers = sum(outlier))))

# Outliers in red on top of the rest; same look in a, b and c
status_layers <- function() {
  list(
    geom_point(aes(colour = status), size = 0.5, alpha = 0.7, stroke = 0,
               position = position_jitter(width = 0.25, height = 0, seed = 1)),
    scale_colour_manual(values = c(kept = COL_NEUTRAL, outlier = COL_UP),
                        labels = c(kept = "Kept", outlier = "Outlier (removed)"),
                        name = NULL, drop = FALSE),
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))),
    labs(x = NULL),
    theme_pub()
  )
}

ordered_status <- arrange(qc, status) # outliers drawn last

# ---- a. Library size --------------------------------------------------------

p_library <- ggplot(ordered_status, aes(region, library_size / 1e6)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = PUB_LINE_WIDTH,
               colour = INK_SECONDARY) +
  status_layers() +
  labs(y = "Library size (M reads)")

# ---- b. Detected genes ------------------------------------------------------

p_detected <- ggplot(ordered_status, aes(region, detected_genes)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = PUB_LINE_WIDTH,
               colour = INK_SECONDARY) +
  status_layers() +
  labs(y = "Detected genes (CPM > 1)")

# ---- c. Standardized connectivity -------------------------------------------

p_connectivity <- ggplot(ordered_status, aes(region, connectivity_z)) +
  geom_hline(yintercept = max_z, linetype = "dashed", linewidth = PUB_LINE_WIDTH) +
  status_layers() +
  geom_text_repel(data = filter(qc, outlier), aes(label = specimenID),
                  size = PUB_LABEL_SIZE * 0.8, segment.size = 0.2,
                  min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  labs(y = "Standardized connectivity (z)")

# ---- d. Distribution of every sample before and after the correction -------

distribution <- distribution %>%
  group_by(region, stage) %>%
  mutate(rank = rank(q50, ties.method = "first")) %>%
  ungroup()

p_distribution <- ggplot(distribution, aes(rank)) +
  geom_linerange(aes(ymin = q05, ymax = q95), colour = "grey80", linewidth = 0.15) +
  geom_linerange(aes(ymin = q25, ymax = q75), colour = INK_SECONDARY, linewidth = 0.3) +
  geom_point(aes(y = q50), size = 0.15, colour = INK) +
  facet_grid(stage ~ region, scales = "free_x") +
  labs(x = "Samples (ordered by median)", y = "Expression (VST, log2 scale)") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ---- e. GC and length biases before and after the correction ----------------

bias <- read_tables(snakemake@input[["bias"]]) %>%
  set_pub_levels() %>%
  mutate(stage = factor(stage, levels = c("before", "after"),
                        labels = c("Before correction", "After correction")),
         bias = factor(bias, levels = c("GC", "length"),
                       labels = c("GC content bias", "Length bias")))

cat("\nMedian R^2 of the biases:\n")
print(as.data.frame(bias %>% group_by(region, bias, stage) %>%
                      summarise(median_r2 = median(r2), .groups = "drop")))

p_bias <- ggplot(bias, aes(region, r2, fill = stage)) +
  geom_boxplot(outlier.size = 0.2, linewidth = PUB_LINE_WIDTH, width = 0.7,
               position = position_dodge(width = 0.8)) +
  facet_grid(. ~ bias) +
  scale_fill_manual(values = c(`Before correction` = "grey80",
                               `After correction` = NPG[["blue"]]),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = expression(R^2 ~ "of the bias (per sample)")) +
  theme_pub()

# ---- Figure -----------------------------------------------------------------

figure <- ((p_library | p_detected) / p_connectivity / p_distribution / p_bias) +
  plot_layout(heights = c(1, 1, 1.4, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)

cat("\nDone\n")
