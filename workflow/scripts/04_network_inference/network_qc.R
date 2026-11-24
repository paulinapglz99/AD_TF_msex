# network_qc.R
# Expression matrix of one brain region for the co-expression networks (ARACNe).
#   1. Samples: those kept by sample_qc (no outliers).
#   2. Genes: Sage's filter (CPM > min_cpm in >= min_fraction of the samples)
#      in ALL the groups (diagnosis x sex), so the networks of a region share
#      the same genes and every gene is really expressed in every group
#      (chrY genes, XIST... are left out).
#   3. VST (design ~ group) with the GC/length normalization factors of
#      gc_normalize, and removal of the technical covariates of the DESeq2
#      model with limma::removeBatchEffect, protecting the groups.
#      The same correction as the differential expression.
#   4. PCA after the correction (tables, plotted by plot_pca_post), and the
#      distribution of every sample before and after (for the QC panel).
# Run by Snakemake (rule network_qc).

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
covariates <- unlist(snakemake@params[["covariates"]])

samples <- read_tsv(snakemake@input[["samples"]],
                    col_types = cols(individualID = col_character(),
                                     batch = col_character(),
                                     model_batch = col_character(),
                                     apoe4_allele = col_character()))
count_matrix <- readRDS(snakemake@input[["counts"]]) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()
count_matrix <- count_matrix[, samples$specimenID]

# ---- Genes expressed in every group -----------------------------------------

passes <- passes_by_group(count_matrix, samples,
                          snakemake@params[["min_cpm"]], snakemake@params[["min_fraction"]])
in_all_groups <- Reduce(`&`, passes)

# GC/length normalization factors: only genes with GC and length annotation have them
factors <- readRDS(snakemake@input[["factors"]])
keep <- rownames(count_matrix)[in_all_groups] %>% intersect(rownames(factors))

cat("Region:", region_name, "| samples:", ncol(count_matrix), "\n")
cat("Genes passing the filter in each group:\n")
print(vapply(passes, sum, integer(1)))
cat("Genes passing in ALL groups:", sum(in_all_groups), "of", length(in_all_groups),
    "| with GC/length correction (kept):", length(keep), "\n")

# ---- Corrected expression ---------------------------------------------------

before <- vst_counts(count_matrix[keep, ], samples, blind = FALSE,
                     normalization_factors = factors)
expression <- remove_covariates(before, samples, covariates)
saveRDS(expression, snakemake@output[["expression"]])

# Distribution of every sample before and after removing the covariates
# (quantiles of its genes), for the QC panel
sample_quantiles <- function(matrix, stage) {
  quantiles <- apply(matrix, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  as_tibble(t(quantiles), rownames = "specimenID") %>%
    setNames(c("specimenID", "q05", "q25", "q50", "q75", "q95")) %>%
    mutate(region = region_name, stage = stage, .before = 1)
}
bind_rows(sample_quantiles(before, "before"), sample_quantiles(expression, "after")) %>%
  write_tsv(snakemake@output[["distribution"]])

# ---- PCA after the correction -----------------------------------------------

pca <- pca_tables(expression, samples, region_name,
                  top_genes = snakemake@params[["top_genes"]],
                  n_pcs = snakemake@params[["n_pcs"]],
                  covariates = unlist(snakemake@params[["pca_covariates"]]))

write_tsv(pca$scores, snakemake@output[["scores"]])
write_tsv(pca$variance, snakemake@output[["variance"]])
write_tsv(pca$association, snakemake@output[["association"]])

cat("\nDone\n")
