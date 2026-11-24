# plot_de_supplementary.R
# Supplementary figure of the differential expression:
#   a  Volcano plots of AD vs CT, women (top row) and men (bottom row), every
#      region. The 3 most significant DEGs of each panel are named.
#   b  Volcano plots of the sex x diagnosis interaction, every region. Up = the
#      AD effect is larger (more positive) in women than in men. Genes expressed
#      in one sex only (chrY, XIST...) were not tested and are shown in black;
#      the 2 most extreme of them and every significant gene are named.
#   c  Spearman correlation of the log2FC of AD between every region and sex.
#   d  Top pathways of the interaction (GSEA): for each region, the top_n
#      significant pathways (lowest p). Color = NES: positive = the pathway goes
#      more up (or less down) with AD in women than in men.
# Run by Snakemake (rule plot_de_supplementary).

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

collection_name <- snakemake@params[["collection"]]
top_n <- snakemake@params[["top_n"]]
padj_max <- snakemake@params[["padj"]]
log2fc_min <- snakemake@params[["log2fc"]]

read_tables <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(chromosome = col_character())) %>%
    bind_rows()
}

de_sex <- bind_rows(
  read_tables(snakemake@input[["de_female"]]) %>% mutate(sex = "female"),
  read_tables(snakemake@input[["de_male"]]) %>% mutate(sex = "male")
)
de_interaction <- read_tables(snakemake@input[["de_interaction"]])
gsea_interaction <- read_tables(snakemake@input[["gsea_interaction"]]) %>%
  filter(collection == collection_name)

# One panel per sex and region: "Female · DLPFC" ... "Male · CRB"
panel_levels <- paste(rep(SEX_LABELS[SEX_ORDER], each = length(REGION_ORDER)),
                      rep(REGION_ORDER, times = length(SEX_ORDER)), sep = " · ")
panel_of <- function(sex, region) {
  factor(paste(SEX_LABELS[as.character(sex)], region, sep = " · "), levels = panel_levels)
}

# Volcano data: -log10(p), directions in drawing order (significant on top)
volcano_data <- function(table) {
  table %>%
    filter(!is.na(pvalue)) %>%
    mutate(log10p = -log10(pvalue),
           name = coalesce(symbol, gene_id),
           direction = factor(direction, levels = DIRECTION_ORDER)) %>%
    arrange(direction)
}

# Legends to the right, one item per line, in every panel
legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

# Panels c and d share the bottom row: their legends go below them, so the
# heatmap and the pathway names keep the whole width
legend_bottom <- theme(legend.position = "bottom", legend.direction = "horizontal",
                       legend.justification = "center")

fc_limit <- unname(quantile(abs(c(de_sex$log2FoldChange, de_interaction$log2FoldChange)),
                            0.999, na.rm = TRUE))

volcano_layers <- function() {
  list(
    geom_vline(xintercept = c(-log2fc_min, log2fc_min), linetype = "dashed",
               colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH / 2),
    geom_point(aes(colour = direction), size = 0.3, alpha = 0.7, stroke = 0),
    # Few breaks: the small panels sit side by side and their labels would overlap
    scale_x_continuous(limits = c(-fc_limit, fc_limit), breaks = c(-1, 0, 1),
                       oob = scales::squish),
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))),
    labs(y = expression(-log[10](p))),
    theme_pub(),
    legend_right
  )
}

# ---- a. Volcano plots of AD vs CT -------------------------------------------

volcano_sex <- de_sex %>%
  volcano_data() %>%
  mutate(panel = panel_of(sex, region))

named_sex <- volcano_sex %>%
  filter(direction %in% c("up", "down")) %>%
  group_by(panel) %>%
  slice_min(pvalue, n = 3, with_ties = FALSE) %>%
  ungroup()

p_volcano <- ggplot(volcano_sex, aes(log2FoldChange, log10p)) +
  volcano_layers() +
  geom_text_repel(data = named_sex, aes(label = name), size = PUB_LABEL_SIZE * 0.8,
                  segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  facet_wrap(~panel, nrow = 2, scales = "free_y") +
  scale_colour_direction() +
  labs(x = expression(log[2] * "FC of AD (vs CT)"))

# ---- b. Volcano plots of the interaction ------------------------------------

volcano_interaction <- de_interaction %>%
  volcano_data() %>%
  set_pub_levels()

named_interaction <- bind_rows(
  volcano_interaction %>% filter(direction %in% c("up", "down")),
  volcano_interaction %>%
    filter(direction == "not_tested") %>%
    group_by(region) %>%
    slice_min(pvalue, n = 2, with_ties = FALSE) %>%
    ungroup()
)

cat("Interaction: significant and not tested genes per region:\n")
print(as.data.frame(dplyr::count(volcano_interaction, region, direction)))

p_interaction <- ggplot(volcano_interaction, aes(log2FoldChange, log10p)) +
  volcano_layers() +
  geom_text_repel(data = named_interaction, aes(label = name), size = PUB_LABEL_SIZE * 0.8,
                  segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  facet_grid(. ~ region) +
  scale_colour_direction(labels = c(ns = "Not significant",
                                    not_tested = "Not tested (one sex only)",
                                    down = "AD effect larger in men",
                                    up = "AD effect larger in women")) +
  labs(x = expression(log[2] * "FC of the interaction"))

# ---- c. Correlation of the AD effect between regions and sexes --------------

effects_wide <- de_sex %>%
  mutate(panel = panel_of(sex, region)) %>%
  dplyr::select(gene_id, panel, log2FoldChange) %>%
  pivot_wider(names_from = panel, values_from = log2FoldChange)

correlation <- effects_wide %>%
  dplyr::select(all_of(panel_levels)) %>%
  as.matrix() %>%
  cor(method = "spearman", use = "pairwise.complete.obs")

correlation_long <- correlation %>%
  as_tibble(rownames = "row") %>%
  pivot_longer(-row, names_to = "column", values_to = "rho") %>%
  mutate(row = factor(row, levels = rev(panel_levels)),
         column = factor(column, levels = panel_levels))

p_correlation <- ggplot(correlation_long, aes(column, row, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = PUB_LABEL_SIZE * 0.8) +
  scale_fill_effect(name = "Spearman ρ", limit = 1, breaks = c(-1, 0, 1),
                    guide = guide_colourbar(barwidth = unit(25, "mm"))) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  legend_bottom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank())

# ---- d. Top pathways of the interaction -------------------------------------

# Heights (mm): a and b fixed; c and d share the rest of the page.
# About 3 mm per pathway in d.
heights <- c(a = 75, b = 45)
height_bottom <- FIG_MAX_HEIGHT - sum(heights)
rows_that_fit <- floor((height_bottom - 25) / 3)

# Ids of the n most significant interaction pathways of each region
top_pathways <- function(n) {
  gsea_interaction %>%
    filter(padj < padj_max) %>%
    group_by(region) %>%
    slice_min(pvalue, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(id) %>%
    pull(id)
}

n_shown <- top_n
while (n_shown > 1 && length(top_pathways(n_shown)) > rows_that_fit) {
  n_shown <- n_shown - 1
}
top_ids <- top_pathways(n_shown)
cat("\nPanel d: top", n_shown, "interaction pathways per region ->", length(top_ids),
    "pathways; fit in the page:", rows_that_fit, "\n")

top <- gsea_interaction %>%
  filter(id %in% top_ids) %>%
  set_pub_levels() %>%
  mutate(significant = padj < padj_max)

# Rows ordered by similarity of their NES across regions
nes_matrix <- top %>%
  dplyr::select(description, region, NES) %>%
  pivot_wider(names_from = region, values_from = NES, values_fill = 0) %>%
  tibble::column_to_rownames("description") %>%
  as.matrix()
row_order <- if (nrow(nes_matrix) > 1) {
  rownames(nes_matrix)[hclust(dist(nes_matrix))$order]
} else {
  rownames(nes_matrix)
}
top <- mutate(top, description = factor(description, levels = row_order))

p_pathways <- ggplot(top, aes(region, description)) +
  geom_point(aes(colour = NES, size = pmin(-log10(padj), 10), alpha = significant)) +
  scale_colour_effect(name = "NES\n(interaction)", limit = max(abs(top$NES))) +
  scale_size_continuous(range = c(0.3, 3), name = expression(-log[10](padj))) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_y_discrete(labels = \(x) shorten(x, 42)) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  legend_bottom +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = PUB_BASE_SIZE - 1),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2))

# ---- Figure -----------------------------------------------------------------

# No free() inside this nested row: patchwork fails with it ("subscript out of bounds")
bottom <- (p_correlation | p_pathways) + plot_layout(widths = c(1.25, 1))

figure <- free(p_volcano) / free(p_interaction) / bottom +
  plot_layout(heights = c(heights, height_bottom)) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)

cat("\nDone\n")
