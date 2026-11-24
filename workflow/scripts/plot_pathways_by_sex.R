# plot_pathways_by_sex.R
# Supplementary figure: top pathways (GSEA) of AD in women and in men, every
# region, one panel per gene set collection (a = first collection, b = second...).
# In each panel: for each region and sex, the top_n significant pathways with
# the largest |NES|, all shown in every region and sex.
# Color = NES, size = -log10(padj), faded = not significant.
# Run by Snakemake (rule plot_pathways_by_sex).

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

collections <- unlist(snakemake@params[["collections"]])
top_n <- snakemake@params[["top_n"]]
padj_max <- snakemake@params[["padj"]]

read_tables <- function(files) {
  files %>% lapply(read_tsv, show_col_types = FALSE) %>% bind_rows()
}

gsea <- bind_rows(
  read_tables(snakemake@input[["gsea_female"]]) %>% mutate(sex = "female"),
  read_tables(snakemake@input[["gsea_male"]]) %>% mutate(sex = "male")
)

# Each panel gets the same share of the page; about 3 mm per pathway
panel_height <- FIG_MAX_HEIGHT / length(collections)
rows_that_fit <- floor((panel_height - 25) / 3)

dot_plot <- function(collection_name) {
  pathways <- filter(gsea, collection == collection_name)

  # Ids of the n significant pathways with the largest |NES| of each region and sex
  top_pathways <- function(n) {
    pathways %>%
      filter(padj < padj_max) %>%
      group_by(region, sex) %>%
      slice_max(abs(NES), n = n, with_ties = FALSE) %>%
      ungroup() %>%
      distinct(id) %>%
      pull(id)
  }

  # Lower top_n until the pathways fit in the panel
  n_shown <- top_n
  while (n_shown > 1 && length(top_pathways(n_shown)) > rows_that_fit) {
    n_shown <- n_shown - 1
  }
  top_ids <- top_pathways(n_shown)
  cat(collection_name, ": top", n_shown, "pathways per region and sex (config top_n =",
      top_n, ") ->", length(top_ids), "pathways; fit in the panel:", rows_that_fit, "\n")

  top <- pathways %>%
    filter(id %in% top_ids) %>%
    set_pub_levels() %>%
    mutate(significant = padj < padj_max)

  # Rows ordered by similarity of their NES across regions and sexes
  nes_matrix <- top %>%
    mutate(column = paste(region, sex)) %>%
    dplyr::select(description, column, NES) %>%
    pivot_wider(names_from = column, values_from = NES, values_fill = 0) %>%
    tibble::column_to_rownames("description") %>%
    as.matrix()
  row_order <- if (nrow(nes_matrix) > 1) {
    rownames(nes_matrix)[hclust(dist(nes_matrix))$order]
  } else {
    rownames(nes_matrix)
  }
  top <- mutate(top, description = factor(description, levels = row_order))

  ggplot(top, aes(sex, description)) +
    geom_point(aes(colour = NES, size = pmin(-log10(padj), 10), alpha = significant)) +
    facet_grid(. ~ region) +
    scale_colour_effect(name = "NES", limit = max(abs(top$NES))) +
    scale_size_continuous(range = c(0.3, 3), name = expression(-log[10](padj))) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_x_discrete(labels = SEX_LABELS) +
    # One line per pathway: long names are cut
    scale_y_discrete(labels = \(x) shorten(x, 60)) +
    labs(x = NULL, y = NULL, title = gsub("_", " ", collection_name)) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = PUB_BASE_SIZE - 1),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.2),
          legend.position = "right", legend.direction = "vertical")
}

figure <- collections %>%
  lapply(dot_plot) %>%
  wrap_plots(ncol = 1) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)

cat("\nDone\n")
