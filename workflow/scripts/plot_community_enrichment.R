# plot_community_enrichment.R
# Supplementary figure: functions of the main Leiden communities of the 20
# networks. For each network, the top_communities largest communities and the
# most significant term of each (one collection). Each row is one of those
# terms; each column one network (F-CT, F-AD, M-CT, M-AD within each region).
# A point means the term is significant in some community of that network;
# its size is -log10(padj) in the community where it is most enriched.
# Same encoding as the other figures: color = sex; open circle = control,
# filled = AD. So the figure shows which functions form a module in which
# region, diagnosis and sex.
# Run by Snakemake (rule plot_community_enrichment).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  source(snakemake@input[["theme"]])
})

collection_name <- snakemake@params[["collection"]]
top_communities <- snakemake@params[["top_communities"]]

enrichment <- snakemake@input[["enrichment"]] %>%
  lapply(read_tsv, show_col_types = FALSE,
         col_types = cols(community = col_integer(), padj = col_double())) %>%
  bind_rows() %>%
  filter(collection == collection_name)

# Communities are numbered by size (1 = largest): best term of the largest ones
top_terms <- enrichment %>%
  filter(community <= top_communities) %>%
  group_by(region, group, community) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(id, description)

cat("Collection:", collection_name, "| terms shown:", nrow(top_terms), "\n")

# Each term in each network: its best community (lowest padj)
dots <- enrichment %>%
  filter(id %in% top_terms$id) %>%
  group_by(region, group, id, description) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(log10_padj = pmin(-log10(padj), 20)) %>%
  separate_wider_delim(group, "_", names = c("diagnosis", "sex"), cols_remove = FALSE) %>%
  mutate(group = factor(group, levels = GROUP_ORDER)) %>%
  set_pub_levels()

# Rows ordered by similarity of their pattern across networks
pattern <- dots %>%
  mutate(column = paste(region, group)) %>%
  dplyr::select(description, column, log10_padj) %>%
  pivot_wider(names_from = column, values_from = log10_padj, values_fill = 0) %>%
  tibble::column_to_rownames("description") %>%
  as.matrix()
row_order <- if (nrow(pattern) > 1) {
  rownames(pattern)[hclust(dist(pattern))$order]
} else {
  rownames(pattern)
}

p_dots <- dots %>%
  mutate(description = factor(description, levels = row_order)) %>%
  ggplot(aes(group, description)) +
  geom_point(aes(size = log10_padj, colour = sex, shape = diagnosis), stroke = 0.5) +
  facet_grid(. ~ region) +
  scale_x_discrete(labels = GROUP_LABELS, drop = FALSE) +
  scale_y_discrete(labels = \(x) shorten(x, 55)) +
  scale_colour_manual(values = COL_SEX, labels = SEX_LABELS, name = "Sex",
                      guide = guide_legend(order = 1, override.aes = list(size = 2))) +
  scale_shape_manual(values = c(CT = 1, AD = 16), labels = DX_LABELS, name = "Diagnosis",
                     guide = guide_legend(order = 2, override.aes = list(size = 2))) +
  scale_size_continuous(range = c(0.4, 2.8), name = expression(-log[10](padj)),
                        guide = guide_legend(order = 3)) +
  labs(x = NULL, y = NULL, title = gsub("_", " ", collection_name)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = PUB_BASE_SIZE - 1),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.2),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.3),
        axis.line = element_blank(),
        legend.position = "right", legend.direction = "vertical")

# Height grows with the number of terms, up to the page
height <- min(40 + 2.8 * length(row_order), FIG_MAX_HEIGHT)
save_pub(p_dots, snakemake@output[["pdf"]], width = "double", height = height)

cat("Done\n")
