# plot_networks_supplementary.R
# Supplementary figure of the networks:
#   a  Topology of the 20 networks: edges, density, mean degree, clustering,
#      modularity and number of Leiden communities. Color = diagnosis,
#      shape = sex.
#   b  Every pair of networks of each region with every measure: edge and node
#      Jaccard, Spearman rho of degree and NMI of communities.
#   c  Sensitivity of the rewired x DEG overlap to the top fraction of rewired
#      genes (log2 odds ratio with 5%, 10% and 20%).
# WARNING: each network has its own number of samples; with more samples
# ARACNe finds more edges.
# Run by Snakemake (rule plot_networks_supplementary).

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

read_all <- function(files) {
  files %>% lapply(read_tsv, show_col_types = FALSE) %>% bind_rows()
}

metrics <- read_all(snakemake@input[["metrics"]])
pairs <- read_all(snakemake@input[["pairs"]])
tests <- read_all(snakemake@input[["tests"]])

legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

# ---- a. Topology of the 20 networks -----------------------------------------

topology_metrics <- c("n_edges", "density", "degree_mean", "clustering_global",
                      "modularity", "n_communities")

topology <- metrics %>%
  separate_wider_delim(group, "_", names = c("diagnosis", "sex"), cols_remove = FALSE) %>%
  pivot_longer(all_of(topology_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = topology_metrics,
                         labels = label_variable(topology_metrics))) %>%
  set_pub_levels()

p_topology <- ggplot(topology, aes(region, value, colour = diagnosis, shape = sex)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.3) +
  facet_wrap(~metric, nrow = 2, scales = "free_y") +
  scale_colour_dx() +
  scale_shape_sex() +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  legend_right +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- b. Every pair, every measure -------------------------------------------

measures <- c("edge_jaccard", "node_jaccard", "degree_spearman", "community_nmi")

# Pairs in the same order in every region: x before y in GROUP_ORDER
ordered_pairs <- pairs %>%
  mutate(swap = match(x, GROUP_ORDER) > match(y, GROUP_ORDER),
         first = if_else(swap, y, x),
         second = if_else(swap, x, y),
         pair = label_pair(first, second))
pair_levels <- ordered_pairs %>%
  distinct(first, second, pair) %>%
  arrange(match(first, GROUP_ORDER), match(second, GROUP_ORDER)) %>%
  pull(pair)

comparison <- ordered_pairs %>%
  pivot_longer(all_of(measures), names_to = "measure", values_to = "value") %>%
  mutate(measure = factor(measure, levels = measures, labels = label_variable(measures)),
         pair = factor(pair, levels = rev(pair_levels))) %>%
  set_pub_levels()

p_pairs <- ggplot(comparison, aes(region, pair, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", value)), size = PUB_LABEL_SIZE * 0.7) +
  facet_wrap(~measure, nrow = 1) +
  scale_fill_gradient(low = "white", high = NPG[["blue"]], limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "Similarity") +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  legend_right +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank())

# ---- c. Sensitivity of the rewired x DEG overlap ----------------------------

sensitivity <- tests %>%
  mutate(log2_or = log2(odds_ratio),
         top = paste0(100 * top_fraction, "%")) %>%
  set_pub_levels()

p_sensitivity <- ggplot(sensitivity, aes(top_fraction, log2_or, colour = sex)) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH) +
  geom_line(linewidth = PUB_LINE_WIDTH) +
  geom_point(aes(shape = main), size = 1.3) +
  facet_grid(. ~ region) +
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "Main", `FALSE` = "Sensitivity"), name = NULL) +
  scale_x_continuous(breaks = sort(unique(sensitivity$top_fraction)),
                     labels = \(x) paste0(100 * x, "%")) +
  labs(x = "Top fraction of rewired genes",
       y = expression(log[2] ~ "odds ratio, rewired" %*% "DEG")) +
  theme_pub() +
  legend_right

# ---- Figure -----------------------------------------------------------------

figure <- (p_topology / p_pairs / p_sensitivity) +
  plot_layout(heights = c(1.3, 1.1, 0.8)) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 210)

cat("Done\n")
