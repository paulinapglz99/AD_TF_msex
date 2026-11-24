# plot_networks_by_sex.R
# Main figure: effect of AD on the co-expression networks of women and men.
#   a  Convergence/divergence index of every region and measure: similarity of
#      the networks of women and men with AD minus in controls. > 0 (red):
#      with AD the networks of the two sexes are more alike; < 0 (blue): less.
#   b  Edge Jaccard of AD vs CT in women and in men, every region: the lower,
#      the more AD rewires the network. Under the bars, the samples of each
#      network (AD / CT): with more samples ARACNe finds more edges.
#   c  Are the rewired genes the DEGs of the same sex? log2 odds ratio of the
#      overlap (main top fraction); * = BH-adjusted permutation p < 0.05.
#   d  Adjusted rewiring of every gene by AD in women vs in men, colored by the
#      sex in which it is among the top rewired genes; rho = Spearman.
# Letters follow the reading order: a on top, b and c side by side, d below.
# Run by Snakemake (rule plot_networks_by_sex).

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
  files %>%
    lapply(read_tsv, show_col_types = FALSE, col_types = cols(gene = col_character())) %>%
    bind_rows()
}

convergence <- read_all(snakemake@input[["convergence"]]) %>% set_pub_levels()
pairs <- read_all(snakemake@input[["pairs"]]) %>% set_pub_levels()
genes <- read_all(snakemake@input[["genes"]]) %>% set_pub_levels()
tests <- read_all(snakemake@input[["tests"]]) %>% set_pub_levels()

legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

measures <- c("edge_jaccard", "node_jaccard", "degree_spearman", "community_nmi")

# ---- a. Convergence/divergence index ----------------------------------------

p_convergence <- convergence %>%
  mutate(measure = factor(measure, levels = measures, labels = label_variable(measures))) %>%
  ggplot(aes(region, index, fill = direction)) +
  geom_hline(yintercept = 0, linewidth = PUB_LINE_WIDTH) +
  geom_col(width = 0.7) +
  facet_wrap(~measure, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(convergence = COL_UP, divergence = COL_DOWN),
                    labels = c(convergence = "Convergence", divergence = "Divergence"),
                    name = "With AD") +
  labs(x = NULL, y = "Similarity between sexes\n(AD − CT)") +
  theme_pub() +
  legend_right +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- b. Edge Jaccard of AD vs CT in each sex --------------------------------

within_sex <- pairs %>%
  mutate(sex = case_when(
    x %in% c("AD_female", "CT_female") & y %in% c("AD_female", "CT_female") ~ "female",
    x %in% c("AD_male", "CT_male") & y %in% c("AD_male", "CT_male") ~ "male"
  )) %>%
  filter(!is.na(sex)) %>%
  # n of the AD network and of the CT network, whatever the order of the pair
  mutate(n_ad = if_else(startsWith(x, "AD"), n_x, n_y),
         n_ct = if_else(startsWith(x, "AD"), n_y, n_x),
         label = paste0(n_ad, "/", n_ct)) %>%
  set_pub_levels()

p_jaccard <- ggplot(within_sex, aes(region, edge_jaccard, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = label, y = 0), position = position_dodge(width = 0.8),
            vjust = 1.4, size = PUB_LABEL_SIZE * 0.8) +
  scale_fill_sex() +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  labs(x = "Region (samples AD/CT under each bar)", y = "Edge Jaccard, AD vs CT") +
  theme_pub() +
  legend_right

# ---- d. Rewiring in women vs in men -----------------------------------------

rewiring <- genes %>%
  filter(!is.na(rewired_in)) %>%
  mutate(rewired_in = factor(rewired_in, levels = SHARED_ORDER)) %>%
  arrange(rewired_in)

correlations <- rewiring %>%
  group_by(region) %>%
  summarise(rho = cor(rewiring_adj_female, rewiring_adj_male, method = "spearman"),
            .groups = "drop") %>%
  mutate(label = sprintf("ρ = %.2f", rho))
cat("Spearman correlation of the rewiring in women and in men:\n")
print(as.data.frame(correlations))

limit <- unname(quantile(abs(c(rewiring$rewiring_adj_female, rewiring$rewiring_adj_male)),
                         0.999, na.rm = TRUE))

p_rewiring <- ggplot(rewiring, aes(rewiring_adj_female, rewiring_adj_male)) +
  geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH / 2) +
  geom_vline(xintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH / 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = PUB_LINE_WIDTH) +
  geom_point(aes(colour = rewired_in), size = 0.3, alpha = 0.6, stroke = 0) +
  geom_text(data = correlations, aes(label = label), x = -limit, y = limit,
            hjust = 0, vjust = 1, size = PUB_LABEL_SIZE) +
  facet_grid(. ~ region) +
  scale_x_continuous(limits = c(-limit, limit), oob = scales::squish) +
  scale_y_continuous(limits = c(-limit, limit), oob = scales::squish) +
  scale_colour_shared(name = "Rewired in") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  coord_equal() +
  labs(x = "Adjusted rewiring by AD in women", y = "Adjusted rewiring\nby AD in men") +
  theme_pub() +
  legend_right

# ---- c. Rewired genes x DEGs ------------------------------------------------

overlap <- tests %>%
  filter(main) %>%
  mutate(padj = p.adjust(p_perm, method = "BH"),
         log2_or = log2(odds_ratio),
         star = if_else(!is.na(padj) & padj < 0.05, "*", "")) %>%
  set_pub_levels()

cat("\nRewired x DEG overlap (main fraction):\n")
print(as.data.frame(overlap %>% dplyr::select(region, sex, n_both, expected_both, odds_ratio,
                                              p_perm, padj)))

p_overlap <- ggplot(overlap, aes(region, log2_or, colour = sex)) +
  geom_hline(yintercept = 0, linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  geom_text(aes(label = star), position = position_dodge(width = 0.5),
            vjust = -0.6, size = PUB_LABEL_SIZE * 1.5, show.legend = FALSE) +
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex") +
  labs(x = NULL, y = expression(log[2] ~ "odds ratio, rewired" %*% "DEG")) +
  theme_pub() +
  legend_right

# ---- Figure -----------------------------------------------------------------

figure <- (p_convergence / (p_jaccard | p_overlap) / p_rewiring) +
  plot_layout(heights = c(1, 1, 1.1)) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 170)

cat("\nDone\n")
