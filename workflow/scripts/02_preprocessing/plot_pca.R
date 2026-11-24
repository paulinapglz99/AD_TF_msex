# plot_pca.R
# One panel with the PCA of every brain region, one row per region:
#   left   PC1 vs PC2: color = diagnosis, shape = sex, one ellipse per sex
#   right  Variance of each PC explained by each covariate (R^2 of PC ~ covariate),
#          written where the BH-adjusted p < 0.05
# Letters: first the PCAs (a-e), then the heatmaps (f-j), in region order.
# Generic: draws the tables written by pca_pre.R (rule plot_pca_pre) and,
# later, by the PCA of the network QC.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  source(snakemake@input[["theme"]])
})

covariates <- unlist(snakemake@params[["covariates"]])

read_all <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(specimenID = col_character(),
                            individualID = col_character(), # numeric in Mayo
                            batch = col_character(),
                            model_batch = col_character(),
                            apoe4_allele = col_character())) %>%
    bind_rows()
}

scores <- read_all(snakemake@input[["scores"]]) %>% set_pub_levels()
variance <- read_all(snakemake@input[["variance"]])
association <- read_all(snakemake@input[["association"]])

pc_label <- function(region_name, i) {
  v <- variance %>% filter(region == region_name, pc == paste0("PC", i)) %>% pull(variance)
  sprintf("PC%d (%.1f%%)", i, v)
}

# ---- One row per region -----------------------------------------------------

plot_region <- function(region_name) {
  p_pca <- scores %>%
    filter(region == region_name) %>%
    ggplot(aes(PC1, PC2)) +
    stat_ellipse(aes(linetype = sex, group = sex), colour = INK,
                 linewidth = PUB_LINE_WIDTH) +
    geom_point(aes(colour = diagnosis, shape = sex), size = 0.8, alpha = 0.8) +
    scale_colour_dx() +
    scale_shape_sex() +
    scale_linetype_sex() +
    labs(x = pc_label(region_name, 1), y = pc_label(region_name, 2), title = region_name) +
    theme_pub()

  region_association <- association %>% filter(region == region_name)
  n_pcs <- n_distinct(region_association$pc)
  pc_axis <- region_association %>%
    distinct(pc, variance) %>%
    mutate(label = sprintf("%s\n%.0f%%", pc, variance))

  # R^2 is written only where BH-adjusted p < 0.05; same 0-1 scale in every region
  p_association <- region_association %>%
    mutate(pc = factor(pc, levels = paste0("PC", seq_len(n_pcs))),
           covariate = factor(covariate, levels = rev(covariates)),
           label = if_else(padj < 0.05, sprintf("%.2f", r2), "")) %>%
    ggplot(aes(pc, covariate, fill = r2)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(aes(label = label), size = PUB_LABEL_SIZE * 0.8, colour = INK) +
    # NA (white): covariate with a single value in the region, e.g. protocol in Mayo
    scale_fill_gradient(low = "white", high = NPG[["blue"]], limits = c(0, 1),
                        breaks = c(0, 0.5, 1), name = expression(R^2), na.value = "white",
                        guide = guide_colourbar(barwidth = unit(15, "mm"))) +
    scale_x_discrete(labels = setNames(pc_axis$label, pc_axis$pc)) +
    scale_y_discrete(labels = label_variable) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())

  list(pca = p_pca, association = p_association)
}

rows <- lapply(REGION_ORDER, plot_region)

# One letter per plot: first all the PCAs (a, b, ...), then all the
# covariate heatmaps, in the order of the regions
n_regions <- length(rows)
for (i in seq_len(n_regions)) {
  rows[[i]]$pca <- rows[[i]]$pca + labs(tag = letters[i])
  rows[[i]]$association <- rows[[i]]$association + labs(tag = letters[n_regions + i])
}

panel <- rows %>%
  lapply(\(r) list(r$pca, r$association)) %>%
  unlist(recursive = FALSE) %>%
  wrap_plots(ncol = 2, widths = c(1, 1.6), guides = "collect") &
  # Legends at the foot of the panel (collected once by patchwork)
  theme(legend.position = "bottom", legend.justification = "center")

save_pub(panel, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)

cat("Done\n")
