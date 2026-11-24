# pca_pre.R
# PCA of one brain region before any correction, to choose the covariates of
# the DESeq2 model. Genes: Sage's filter (CPM > min_cpm in >= min_fraction of
# the samples of one group), variance stabilizing transformation (DESeq2::vst)
# and the most variable genes.
# Saves tables only (plots: plot_pca.R), so they can be reused:
#   scores       PC1..n_pcs of every sample, with its covariates
#   variance     % of variance of every PC
#   association  R^2 and p of PC ~ covariate, for each PC and covariate
# Run by Snakemake (rule pca_pre).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

# Bioconductor first and the tidyverse last: Bioconductor packages mask dplyr
# functions such as count() and select() if they are loaded after it
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  source(snakemake@input[["gene_filter"]])
})

region_name <- snakemake@wildcards[["region"]]
min_cpm <- snakemake@params[["min_cpm"]]
min_fraction <- snakemake@params[["min_fraction"]]
top_genes <- snakemake@params[["top_genes"]]
n_pcs <- snakemake@params[["n_pcs"]]
covariates <- unlist(snakemake@params[["covariates"]])

samples <- read_tsv(snakemake@input[["samples"]],
                    col_types = cols(batch = col_character(),
                                     model_batch = col_character(),
                                     apoe4_allele = col_character()))
counts <- readRDS(snakemake@input[["counts"]])

count_matrix <- counts %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

# ---- Gene filter ------------------------------------------------------------

keep <- filter_genes(count_matrix, samples, min_cpm, min_fraction)

cat("Genes passing the filter:", sum(keep), "of", length(keep), "\n")

# ---- PCA --------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(round(count_matrix[keep, ]),
                              colData = tibble::column_to_rownames(samples, "specimenID"),
                              design = ~1)
expression <- assay(vst(dds, blind = TRUE))

most_variable <- head(order(apply(expression, 1, var), decreasing = TRUE), top_genes)
pca <- prcomp(t(expression[most_variable, ]))

variance <- tibble(region = region_name,
                   pc = paste0("PC", seq_along(pca$sdev)),
                   variance = 100 * pca$sdev^2 / sum(pca$sdev^2))

scores <- pca$x[, seq_len(n_pcs)] %>%
  as_tibble(rownames = "specimenID") %>%
  left_join(samples, by = "specimenID")

# ---- Variance of each PC explained by each covariate ------------------------

r_squared <- function(pc, covariate) {
  data <- tibble(score = scores[[pc]], value = scores[[covariate]]) %>%
    filter(!is.na(value))
  if (n_distinct(data$value) < 2) return(tibble(r2 = NA_real_, p = NA_real_))
  fit <- lm(score ~ value, data = data)
  tibble(r2 = summary(fit)$r.squared, p = anova(fit)[["Pr(>F)"]][1])
}

association <- expand_grid(pc = paste0("PC", seq_len(n_pcs)), covariate = covariates) %>%
  rowwise() %>%
  mutate(r_squared(pc, covariate)) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  left_join(variance, by = "pc") %>%
  relocate(region, pc, variance)

write_tsv(scores, snakemake@output[["scores"]])
write_tsv(variance, snakemake@output[["variance"]])
write_tsv(association, snakemake@output[["association"]])

cat("\nDone\n")
