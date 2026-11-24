# plot_networks_supplementary.R
# Supplementary figure of the network comparisons (the topology of each network
# is in plot_network_topology.R):
#   a-d  Every pair of networks of each region with one measure each: edge
#        Jaccard, node Jaccard, Spearman rho of degree and NMI of communities.
#        Each heatmap has its own color scale (the measures have very
#        different ranges); the value is written in every tile.
#   e    Sensitivity of the rewired x DEG overlap to the top fraction of
#        rewired genes (log2 odds ratio with 5%, 10% and 20%).
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

pairs <- read_all(snakemake@input[["pairs"]])
tests <- read_all(snakemake@input[["tests"]])

legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

# ---- a-d. Every pair, one measure per heatmap -------------------------------

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
  mutate(pair = factor(pair, levels = rev(pair_levels))) %>%
  set_pub_levels()

# Grey scale: blue, red, green and orange already mean diagnosis and sex
plot_pairs <- function(m) {
  df <- comparison %>%
    filter(measure == m) %>%
    # Scale of each heatmap: from its own minimum to its own maximum
    mutate(shade = (value - min(value, na.rm = TRUE)) /
             max(diff(range(value, na.rm = TRUE)), 1e-9))
  ggplot(df, aes(region, pair, fill = shade)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", value), colour = shade > 0.6),
              size = PUB_LABEL_SIZE * 0.8) +
    scale_fill_gradient(low = "grey97", high = "grey35", guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
    labs(x = NULL, y = NULL, title = label_variable(m)) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(), axis.ticks = element_blank())
}

# Pair names on every heatmap: without them the letter touches the title
pair_plots <- lapply(measures, plot_pairs)

# ---- e. Sensitivity of the rewired x DEG overlap ----------------------------

sensitivity <- tests %>%
  mutate(log2_or = log2(odds_ratio),
         top = paste0(100 * top_fraction, "%")) %>%
  set_pub_levels()

p_sensitivity <- ggplot(sensitivity, aes(top_fraction, log2_or, colour = sex)) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH) +
  geom_line(linewidth = PUB_LINE_WIDTH) +
  geom_point(aes(shape = main), size = 1.3) +
  facet_grid(. ~ region) +
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex",
                      guide = guide_legend(order = 1)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "Main", `FALSE` = "Sensitivity"), name = NULL,
                     guide = guide_legend(order = 2)) +
  # Numbers only (the % is in the axis title), so neighbouring facets do not touch
  scale_x_continuous(breaks = sort(unique(sensitivity$top_fraction)),
                     labels = \(x) 100 * x, expand = expansion(mult = 0.08)) +
  labs(x = "Top % of rewired genes",
       y = expression(log[2] ~ "odds ratio, rewired" %*% "DEG")) +
  theme_pub() +
  legend_right +
  theme(panel.spacing.x = unit(4, "mm"))

# ---- Figure -----------------------------------------------------------------

# One letter per measure (a-d) and e for the sensitivity
top <- wrap_plots(pair_plots, ncol = 2)
figure <- (top / p_sensitivity) +
  plot_layout(heights = c(2, 0.8)) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 175)

cat("Done\n")
