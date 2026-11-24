# plot_degs_ad_by_sex.R
# Main figure: effect of AD in women and in men, in every region.
#   a  Number of DEGs (AD vs CT) in women and in men: up above 0, down below 0.
#      Under each bar, the samples of that sex (AD / CT): women have 2-3 times
#      more AD samples, which alone gives them more DEGs.
#   b  log2FC of AD in women vs in men, one point per gene, colored by where the
#      gene is a DEG. Uses every gene, so it does not depend on the DEG threshold
#      or on the sample size. Diagonal = same effect in both sexes;
#      rho = Spearman correlation of the two effects.
#   c  Same at the pathway level (GSEA): NES of AD in women vs in men, one point
#      per pathway. Colored: pathways significant in the GSEA of the sex x
#      diagnosis interaction, red if the AD effect is larger in women, blue if
#      larger in men. Their names are in the supplementary pathway figure.
# Pathways: one gene set collection (config: figures/pathways/collection).
# Run by Snakemake (rule plot_degs_ad_by_sex).

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
padj_max <- snakemake@params[["padj"]]

read_tables <- function(files) {
  files %>%
    lapply(read_tsv, show_col_types = FALSE,
           col_types = cols(specimenID = col_character(),
                            individualID = col_character(),
                            batch = col_character(),
                            model_batch = col_character(),
                            apoe4_allele = col_character(),
                            chromosome = col_character())) %>%
    bind_rows()
}

read_gsea <- function(files) {
  read_tables(files) %>% filter(collection == collection_name)
}

de <- bind_rows(
  read_tables(snakemake@input[["de_female"]]) %>% mutate(sex = "female"),
  read_tables(snakemake@input[["de_male"]]) %>% mutate(sex = "male")
)
samples <- read_tables(snakemake@input[["samples"]])

gsea <- bind_rows(
  read_gsea(snakemake@input[["gsea_female"]]) %>% mutate(sex = "female"),
  read_gsea(snakemake@input[["gsea_male"]]) %>% mutate(sex = "male")
)
gsea_interaction <- read_gsea(snakemake@input[["gsea_interaction"]])

# Legends to the right, one item per line, in every panel
legend_right <- theme(legend.position = "right", legend.direction = "vertical",
                      legend.justification = "center")

# Same look for the two scatter plots (b, c)
scatter_layers <- function(limit) {
  list(
    geom_hline(yintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH / 2),
    geom_vline(xintercept = 0, colour = COL_NEUTRAL, linewidth = PUB_LINE_WIDTH / 2),
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = PUB_LINE_WIDTH),
    facet_grid(. ~ region),
    scale_x_continuous(limits = c(-limit, limit), oob = scales::squish),
    scale_y_continuous(limits = c(-limit, limit), oob = scales::squish),
    guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))),
    coord_equal(),
    theme_pub(),
    legend_right
  )
}

# rho of each region, written in the top-left corner of its panel
correlation_labels <- function(data, x, y) {
  data %>%
    group_by(region) %>%
    summarise(rho = cor({{ x }}, {{ y }}, method = "spearman"), .groups = "drop") %>%
    mutate(label = sprintf("ρ = %.2f", rho))
}

# ---- a. Number of DEGs ------------------------------------------------------

# x axis label of each bar: sex and its samples, AD / CT
n_samples <- samples %>%
  dplyr::count(region, sex, diagnosis) %>%
  pivot_wider(names_from = diagnosis, values_from = n) %>%
  mutate(label = paste0(SEX_LABELS[sex], "\n", AD, " / ", CT))

degs <- de %>%
  filter(direction %in% c("up", "down")) %>%
  dplyr::count(region, sex, direction) %>%
  complete(region, sex, direction = c("up", "down"), fill = list(n = 0)) %>%
  mutate(signed = if_else(direction == "up", n, -n),
         vjust = if_else(direction == "up", -0.3, 1.3)) %>%
  left_join(n_samples, by = c("region", "sex")) %>%
  set_pub_levels() %>%
  arrange(sex) %>%
  mutate(label = factor(label, levels = unique(label)))

cat("DEGs per region and sex:\n")
print(as.data.frame(degs %>% dplyr::select(region, sex, direction, n)))

p_counts <- ggplot(degs, aes(label, signed, fill = sex)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = PUB_LINE_WIDTH) +
  geom_text(aes(label = n, vjust = vjust), size = PUB_LABEL_SIZE) +
  facet_grid(. ~ region, scales = "free_x") +
  scale_fill_sex() +
  # Room above and below the bars for their numbers
  scale_y_continuous(labels = abs, expand = expansion(mult = 0.25)) +
  labs(x = "Samples (AD / CT)", y = "DEGs (down | up)") +
  theme_pub() +
  legend_right

# ---- b. log2FC in women vs log2FC in men ------------------------------------

effects <- de %>%
  filter(tested) %>%
  dplyr::select(region, gene_id, sex, log2FoldChange, direction) %>%
  pivot_wider(names_from = sex, values_from = c(log2FoldChange, direction)) %>%
  filter(!is.na(log2FoldChange_female), !is.na(log2FoldChange_male)) %>%
  mutate(shared = factor(case_when(
    direction_female != "ns" & direction_male != "ns" ~ "both",
    direction_female != "ns" ~ "female",
    direction_male != "ns" ~ "male",
    .default = "none"
  ), levels = SHARED_ORDER)) %>%
  set_pub_levels() %>%
  arrange(shared) # DEGs drawn on top of the rest

fc_correlations <- correlation_labels(effects, log2FoldChange_female, log2FoldChange_male)
cat("\nSpearman correlation of the AD effect (log2FC) in women and in men:\n")
print(as.data.frame(fc_correlations))

# Same symmetric limits in every region; the few genes beyond them are drawn on the edge
fc_limit <- unname(quantile(abs(c(effects$log2FoldChange_female,
                                  effects$log2FoldChange_male)), 0.999))

p_effects <- ggplot(effects, aes(log2FoldChange_female, log2FoldChange_male)) +
  scatter_layers(fc_limit) +
  geom_point(aes(colour = shared), size = 0.3, alpha = 0.6, stroke = 0) +
  geom_text(data = fc_correlations, aes(label = label), x = -fc_limit, y = fc_limit,
            hjust = 0, vjust = 1, size = PUB_LABEL_SIZE) +
  scale_colour_shared(name = "DEG in") +
  labs(x = expression(log[2] * "FC of AD in women"),
       y = expression(log[2] * "FC of AD in men"))

# ---- c. NES in women vs NES in men ------------------------------------------

cat("\nCollection:", collection_name, "\nPathways significant in the interaction (padj <",
    padj_max, "):\n")
print(as.data.frame(gsea_interaction %>%
                      group_by(region) %>%
                      summarise(significant = sum(padj < padj_max))))

# Interaction NES > 0: the pathway goes more up (or less down) with AD in women
interaction <- gsea_interaction %>%
  transmute(region, id,
            interaction = case_when(padj >= padj_max ~ "ns",
                                    NES > 0 ~ "women",
                                    .default = "men"))

nes <- gsea %>%
  dplyr::select(region, id, sex, NES) %>%
  pivot_wider(names_from = sex, values_from = NES, names_prefix = "NES_") %>%
  filter(!is.na(NES_female), !is.na(NES_male)) %>%
  left_join(interaction, by = c("region", "id")) %>%
  mutate(interaction = factor(coalesce(interaction, "ns"), levels = c("ns", "men", "women"))) %>%
  set_pub_levels() %>%
  arrange(interaction) # significant pathways drawn on top

nes_correlations <- correlation_labels(nes, NES_female, NES_male)
nes_limit <- max(abs(c(nes$NES_female, nes$NES_male)))

p_nes <- ggplot(nes, aes(NES_female, NES_male)) +
  scatter_layers(nes_limit) +
  geom_point(aes(colour = interaction), size = 0.6, alpha = 0.8, stroke = 0) +
  geom_text(data = nes_correlations, aes(label = label), x = -nes_limit, y = nes_limit,
            hjust = 0, vjust = 1, size = PUB_LABEL_SIZE) +
  scale_colour_manual(values = c(ns = "grey80", men = COL_DOWN, women = COL_UP),
                      labels = c(ns = "Not significant",
                                 men = "Larger in men",
                                 women = "Larger in women"),
                      name = "AD effect\n(interaction)") +
  labs(x = "NES of AD in women", y = "NES of AD in men")

# ---- Figure -----------------------------------------------------------------

figure <- (p_counts / p_effects / p_nes) +
  plot_layout(heights = c(1, 1.1, 1.1)) +
  plot_annotation(tag_levels = "a")

save_pub(figure, snakemake@output[["pdf"]], width = "double", height = 150)

cat("\nDone\n")
