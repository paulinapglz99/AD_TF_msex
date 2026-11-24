# plot_networks_by_sex.R
# Main figure: effect of AD on the co-expression networks of women and men,
# from the whole network to the genes and their function.
#   a  How much does AD rewire each sex? Edge Jaccard of AD vs CT in women and
#      in men, every region: the lower, the more AD rewires the network. Under
#      each lollipop, the samples of the AD / CT networks (with more samples
#      ARACNe finds more edges).
#   b  Are the networks of the two sexes more or less alike with AD?
#      Convergence index of every region and measure: similarity of the
#      networks of women and men with AD minus in controls. > 0: they converge;
#      < 0: they diverge.
#   c  Adjusted rewiring of every gene by AD in women vs in men, colored by the
#      sex in which it is among the top rewired genes; rho = Spearman.
#   d  Are the rewired genes the DEGs of the same sex? log2 odds ratio of the
#      overlap (main top fraction); * = BH-adjusted permutation p < 0.05.
#   e  Function of the genes rewired in one sex only: top terms (one
#      collection) of the genes rewired in women only and in men only, per
#      region. Size = -log10(padj).
# Color = sex in every panel (women green, men orange). Facets by region share
# one letter; each measure or question has its own. Every legend appears once.
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

collection_name <- snakemake@params[["collection"]]
top_terms <- snakemake@params[["top_terms"]]

read_all <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(gene = col_character(), rewired_in = col_character())) %>%
    bind_rows()
}

convergence <- read_all(snakemake@input[["convergence"]]) %>% set_pub_levels()
pairs <- read_all(snakemake@input[["pairs"]]) %>% set_pub_levels()
genes <- read_all(snakemake@input[["genes"]]) %>% set_pub_levels()
tests <- read_all(snakemake@input[["tests"]]) %>% set_pub_levels()
enrichment <- read_all(snakemake@input[["enrichment"]]) %>%
  filter(collection == collection_name, rewired_in %in% c("female", "male")) %>%
  set_pub_levels()

legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

regions <- levels(pairs$region)
measures <- c("edge_jaccard", "node_jaccard", "degree_spearman", "community_nmi")

# Women on the left and men on the right of each region (panels a and d)
offset <- c(female = -0.17, male = 0.17)
x_by_sex <- function(region, sex) as.numeric(region) + offset[as.character(sex)]
scale_x_regions <- function() {
  scale_x_continuous(breaks = seq_along(regions), labels = regions,
                     expand = expansion(add = 0.45))
}
scale_colour_sex <- function(...) {
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex", ...)
}

# Heights (mm) of the three rows: a and b; c; d and e
# (c has square facets, so its row only needs their height)
heights <- c(60, 42, 85)

# ---- a. Edge Jaccard of AD vs CT in each sex --------------------------------

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
  set_pub_levels() %>%
  mutate(x = x_by_sex(region, sex))

p_jaccard <- ggplot(within_sex, aes(x, edge_jaccard, colour = sex)) +
  geom_segment(aes(xend = x, y = 0, yend = edge_jaccard), linewidth = PUB_LINE_WIDTH * 2) +
  geom_point(size = 1.8) +
  # n under each lollipop, vertical so the two sexes of a region do not overlap
  geom_text(aes(label = label, y = 0), angle = 90, hjust = 1.1,
            size = PUB_LABEL_SIZE * 0.8, show.legend = FALSE) +
  scale_colour_sex() +
  scale_x_regions() +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.05))) +
  labs(x = NULL, y = "Edge Jaccard, AD vs CT", caption = "Samples AD/CT under each network") +
  theme_pub() +
  legend_right +
  theme(plot.caption = element_text(size = PUB_BASE_SIZE - 1, hjust = 0))

# ---- b. Convergence/divergence index ----------------------------------------

limit_index <- max(abs(convergence$index), na.rm = TRUE)

p_convergence <- convergence %>%
  mutate(measure = factor(measure, levels = rev(measures),
                          labels = label_variable(rev(measures)))) %>%
  ggplot(aes(region, measure, fill = index)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  # "0.00" instead of "-0.00" or "+0.00"
  geom_text(aes(label = if_else(round(index, 2) == 0, "0.00", sprintf("%+.2f", index))),
            size = PUB_LABEL_SIZE * 0.8) +
  # Purple/brown, not the red/blue of the diagnosis
  scale_fill_gradient2(low = "#8C510A", mid = "white", high = "#762A83", midpoint = 0,
                       limits = c(-limit_index, limit_index),
                       breaks = c(-limit_index, 0, limit_index),
                       labels = c("Diverge", "0", "Converge"),
                       name = "Sexes with AD\n(AD − CT)") +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  legend_right +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        # Room between the title and the "Converge" label at the top of the bar
        legend.title = element_text(margin = margin(b = 8)))

# ---- c. Rewiring in women vs in men -----------------------------------------

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
  scale_colour_shared(name = "Rewired by AD in") +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  coord_equal() +
  labs(x = "Adjusted rewiring by AD in women", y = "Adjusted rewiring\nby AD in men") +
  theme_pub() +
  legend_right

# ---- d. Rewired genes x DEGs ------------------------------------------------

overlap <- tests %>%
  filter(main) %>%
  mutate(padj = p.adjust(p_perm, method = "BH"),
         log2_or = log2(odds_ratio),
         star = if_else(!is.na(padj) & padj < 0.05, "*", "")) %>%
  set_pub_levels() %>%
  mutate(x = x_by_sex(region, sex))

cat("\nRewired x DEG overlap (main fraction):\n")
print(as.data.frame(overlap %>% dplyr::select(region, sex, n_both, expected_both, odds_ratio,
                                              p_perm, padj)))

p_overlap <- ggplot(overlap, aes(x, log2_or, colour = sex)) +
  geom_hline(yintercept = 0, linewidth = PUB_LINE_WIDTH, colour = COL_NEUTRAL) +
  geom_segment(aes(xend = x, y = 0, yend = log2_or), linewidth = PUB_LINE_WIDTH * 2) +
  geom_point(size = 1.8) +
  # The star goes on the far side of the point from 0
  geom_text(aes(label = star, vjust = if_else(log2_or >= 0, -0.3, 1.1)),
            size = PUB_LABEL_SIZE * 1.5, show.legend = FALSE) +
  # Same colors as panel a, whose legend explains them
  scale_colour_sex(guide = "none") +
  scale_x_regions() +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = NULL, y = expression(log[2] ~ "odds ratio, rewired" %*% "DEG")) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- e. Function of the genes rewired in one sex only -----------------------

# About 3 mm per term in the bottom row; lower top_terms until they fit
rows_that_fit <- floor((heights[3] - 25) / 3)
top_ids <- function(n) {
  enrichment %>%
    group_by(region, rewired_in) %>%
    slice_min(padj, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(id) %>%
    pull(id)
}
n_shown <- top_terms
while (n_shown > 1 && length(top_ids(n_shown)) > rows_that_fit) {
  n_shown <- n_shown - 1
}
cat("\nPanel e:", collection_name, "| top", n_shown, "terms per region and sex ->",
    length(top_ids(n_shown)), "terms; fit:", rows_that_fit, "\n")

terms <- enrichment %>%
  filter(id %in% top_ids(n_shown)) %>%
  mutate(rewired_in = factor(rewired_in, levels = c("female", "male"),
                             labels = c("Rewired in women only", "Rewired in men only")))

# Rows grouped by the sex where the term is most significant
term_order <- terms %>%
  group_by(description) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(rewired_in, desc(padj)) %>%
  pull(description)

p_function <- if (nrow(terms) == 0) {
  ggplot() +
    annotate("text", x = 0, y = 0, size = PUB_LABEL_SIZE,
             label = "No significant terms for the genes rewired in one sex only") +
    theme_void()
} else {
  terms %>%
    mutate(description = factor(description, levels = term_order)) %>%
    ggplot(aes(region, description)) +
    geom_point(aes(size = pmin(-log10(padj), 10), colour = rewired_in)) +
    facet_grid(. ~ rewired_in) +
    scale_colour_manual(values = c(`Rewired in women only` = COL_SEX[["female"]],
                                   `Rewired in men only` = COL_SEX[["male"]]),
                        guide = "none") +
    scale_size_continuous(range = c(0.5, 3), name = expression(-log[10](padj))) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(labels = \(x) shorten(x, 45)) +
    labs(x = NULL, y = NULL, title = gsub("_", " ", collection_name)) +
    theme_pub() +
    legend_right +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = PUB_BASE_SIZE - 1),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2))
}

# ---- Figure -----------------------------------------------------------------

# Each row on its own: in R, `/` binds tighter than `+`, so a plot_layout()
# written inside the stack would apply to the wrong part
row_ab <- (p_jaccard | p_convergence) + plot_layout(widths = c(1, 1.1))
row_de <- (p_overlap | p_function) + plot_layout(widths = c(1, 2.4))

figure <- (row_ab / p_rewiring / row_de) +
  plot_layout(heights = heights) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = sum(heights) + 15)

cat("\nDone\n")
