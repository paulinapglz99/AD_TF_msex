# plot_network_topology.R
# Supplementary figure: topology of the 20 networks (region x diagnosis x sex).
# Encoding of every panel: color = sex; open circle = control, filled = AD.
#   a-f  Edges, density, mean degree, clustering, modularity and number of
#        Leiden communities. In each region women on the left and men on the
#        right; the segment joins control and AD, so its length and direction
#        are the change with AD in each sex.
#   g    Degree distribution P(k) of every network, log-log, with logarithmic
#        bins (density per unit of degree). Rows: control and AD.
#   h    Edges of every network vs its number of samples: with more samples
#        ARACNe finds more edges, so comparisons between networks mix biology
#        with sample size.
# All panels share one legend.
# Run by Snakemake (rule plot_network_topology).

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

read_all <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(gene = col_character(), specimenID = col_character(),
                            individualID = col_character(), batch = col_character(),
                            model_batch = col_character(), apoe4_allele = col_character())) %>%
    bind_rows()
}

# region, group -> region, group, diagnosis, sex (with the plot order)
split_group <- function(df) {
  df %>%
    separate_wider_delim(group, "_", names = c("diagnosis", "sex"), cols_remove = FALSE) %>%
    set_pub_levels()
}

metrics <- read_all(snakemake@input[["metrics"]]) %>% split_group()
nodes <- read_all(snakemake@input[["nodes"]]) %>% dplyr::select(region, group, degree)
n_samples <- read_all(snakemake@input[["samples"]]) %>%
  dplyr::count(region, group, name = "n_samples")

# The same scales and point size in every panel, so patchwork keeps one legend
POINT_SIZE <- 1.6
SHAPE_DX <- c(CT = 1, AD = 16)
scale_colour_sex <- function(...) {
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex", ...)
}
scale_shape_dx <- function(...) {
  scale_shape_manual(values = SHAPE_DX, labels = DX_LABELS, name = "Diagnosis", ...)
}
legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

# ---- a-f. Topology of the 20 networks ---------------------------------------

topology_metrics <- c("n_edges", "density", "degree_mean", "clustering_global",
                      "modularity", "n_communities")
regions <- levels(metrics$region)
offset <- c(female = -0.17, male = 0.17) # women left, men right of each region

metric_long <- metrics %>%
  pivot_longer(all_of(topology_metrics), names_to = "metric", values_to = "value") %>%
  mutate(x = as.numeric(region) + offset[as.character(sex)])

# One row per region, sex and metric: control and AD side by side
metric_change <- metric_long %>%
  dplyr::select(region, sex, metric, x, diagnosis, value) %>%
  pivot_wider(names_from = diagnosis, values_from = value)

plot_metric <- function(m) {
  ggplot(filter(metric_long, metric == m), aes(x, value, colour = sex)) +
    geom_segment(data = filter(metric_change, metric == m),
                 aes(x = x, xend = x, y = CT, yend = AD),
                 linewidth = PUB_LINE_WIDTH * 2, show.legend = FALSE) +
    geom_point(aes(shape = diagnosis), size = POINT_SIZE, stroke = 0.6) +
    scale_colour_sex() +
    scale_shape_dx() +
    scale_x_continuous(breaks = seq_along(regions), labels = regions,
                       expand = expansion(add = 0.4)) +
    # 800K instead of 8e+05; small values are left as they are
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale(),
                                                     drop0trailing = TRUE)) +
    labs(x = NULL, y = NULL, title = label_variable(m)) +
    theme_pub() +
    legend_right +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

metric_plots <- lapply(topology_metrics, plot_metric)

# ---- g. Degree distributions ------------------------------------------------

# Logarithmic bins of degree; P(k) = genes in the bin / (all genes x bin width),
# so the curve does not depend on the width of the bins
n_bins <- 25
breaks <- 10^seq(0, log10(max(nodes$degree) + 1), length.out = n_bins + 1)

degree_distribution <- nodes %>%
  filter(degree > 0) %>%
  mutate(bin = cut(degree, breaks = breaks, include.lowest = TRUE, labels = FALSE)) %>%
  dplyr::count(region, group, bin) %>%
  group_by(region, group) %>%
  mutate(degree = sqrt(breaks[bin] * breaks[bin + 1]), # geometric middle of the bin
         p = n / (sum(n) * (breaks[bin + 1] - breaks[bin]))) %>%
  ungroup() %>%
  split_group() %>%
  mutate(diagnosis = factor(DX_LABELS[as.character(diagnosis)], levels = DX_LABELS))

p_degree <- ggplot(degree_distribution, aes(degree, p, colour = sex)) +
  geom_line(linewidth = PUB_LINE_WIDTH * 1.5) +
  facet_grid(diagnosis ~ region) +
  # Ticks on powers of 10 only, so the labels are readable
  scale_x_log10(breaks = 10^(0:4), labels = scales::label_log()) +
  scale_y_log10(breaks = 10^seq(-8, 0, by = 2), labels = scales::label_log()) +
  # Same colors as the other panels: their legend already explains them
  scale_colour_sex(guide = "none") +
  labs(x = "Degree (k)", y = "P(k)") +
  theme_pub()

# ---- h. Edges vs number of samples ------------------------------------------

size_vs_n <- metrics %>%
  left_join(n_samples, by = c("region", "group"))

cat("Spearman correlation of edges and number of samples:",
    round(cor(size_vs_n$n_samples, size_vs_n$n_edges, method = "spearman"), 2), "\n")

# Control and AD of the same region and sex joined by a segment
size_change <- size_vs_n %>%
  dplyr::select(region, sex, diagnosis, n_samples, n_edges) %>%
  pivot_wider(names_from = diagnosis, values_from = c(n_samples, n_edges))

p_size <- ggplot(size_vs_n, aes(n_samples, n_edges, colour = sex)) +
  geom_segment(data = size_change,
               aes(x = n_samples_CT, xend = n_samples_AD, y = n_edges_CT, yend = n_edges_AD),
               linewidth = PUB_LINE_WIDTH * 1.5, alpha = 0.5, show.legend = FALSE) +
  geom_point(aes(shape = diagnosis), size = POINT_SIZE, stroke = 0.6) +
  # One label per pair (on the AD network); repelled from every point and from
  # each other, so none overlap
  geom_text_repel(data = filter(size_vs_n, diagnosis == "AD"), aes(label = region),
                  size = PUB_LABEL_SIZE * 0.75, show.legend = FALSE,
                  box.padding = 0.35, point.padding = 0.25, max.overlaps = Inf,
                  segment.size = 0.2, min.segment.length = 0, seed = 1) +
  scale_colour_sex() +
  scale_shape_dx() +
  # Log scale: most networks have 25-60 samples and would crowd on a linear axis
  scale_x_log10(breaks = c(25, 50, 100, 200)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = "Samples of the network", y = "Edges") +
  theme_pub() +
  legend_right

# ---- Figure -----------------------------------------------------------------

# One letter per plot: a-f the metrics, g the degree distributions, h the edges
# vs samples. Every legend is the same, so guides = "collect" keeps only one.
top <- wrap_plots(metric_plots, ncol = 3)
bottom <- (p_degree | p_size) + plot_layout(widths = c(2.2, 1))
figure <- (top / bottom) +
  plot_layout(heights = c(1.3, 1), guides = "collect") +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 170)

cat("Done\n")
