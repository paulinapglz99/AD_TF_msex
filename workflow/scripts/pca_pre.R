# pca_pre.R
# PCA of one brain region before any correction (all samples, outliers
# included), to choose the covariates of the DESeq2 model. Genes: Sage's filter
# (CPM > min_cpm in >= min_fraction of the samples of one group), variance
# stabilizing transformation (DESeq2::vst, blind) and the most variable genes.
# Saves tables only (plots: plot_pca.R), so they can be reused:
#   scores       PC1..n_pcs of every sample, with its covariates
#   variance     % of variance of every PC
#   association  R^2 and p of PC ~ covariate, for each PC and covariate
# Run by Snakemake (rule pca_pre).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  source(snakemake@input[["gene_filter"]])
  source(snakemake@input[["expression_utils"]])
})

region_name <- snakemake@wildcards[["region"]]

samples <- read_tsv(snakemake@input[["samples"]],
                    col_types = cols(individualID = col_character(),
                                     batch = col_character(),
                                     model_batch = col_character(),
                                     apoe4_allele = col_character()))
count_matrix <- readRDS(snakemake@input[["counts"]]) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

keep <- filter_genes(count_matrix, samples,
                     snakemake@params[["min_cpm"]], snakemake@params[["min_fraction"]])
cat("Genes passing the filter:", sum(keep), "of", length(keep), "\n")

expression <- vst_counts(count_matrix[keep, samples$specimenID], samples, blind = TRUE)

pca <- pca_tables(expression, samples, region_name,
                  top_genes = snakemake@params[["top_genes"]],
                  n_pcs = snakemake@params[["n_pcs"]],
                  covariates = unlist(snakemake@params[["covariates"]]))

write_tsv(pca$scores, snakemake@output[["scores"]])
write_tsv(pca$variance, snakemake@output[["variance"]])
write_tsv(pca$association, snakemake@output[["association"]])

cat("\nDone\n")
