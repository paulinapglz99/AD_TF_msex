# Expression helpers shared by pca_pre, sample_qc and network_qc, so the three
# use exactly the same normalization, correction and PCA.
# Scripts load it with:
#   source(snakemake@input[["expression_utils"]])
# Bioconductor functions are called with pkg::fun, so this file does not load
# (and mask) anything.

library(dplyr)
library(tidyr)

# Sample table ready to be a model's colData: group as factor, categorical
# covariates as factors, numeric covariates centered and scaled.
# The same transformation as in deseq2.R.
model_coldata <- function(samples, covariates) {
  samples %>%
    mutate(group = factor(group),
           across(all_of(covariates) & where(is.character), factor),
           across(all_of(covariates) & where(is.numeric), \(x) as.numeric(scale(x)))) %>%
    tibble::column_to_rownames("specimenID")
}

# GC/length normalization factors (rule gc_normalize) for some genes and
# samples, centered again per gene (geometric mean 1) after the subset
normalization_factors_for <- function(factors, genes, sample_ids) {
  subset <- factors[genes, sample_ids, drop = FALSE]
  subset / exp(rowMeans(log(subset)))
}

# Variance stabilizing transformation (DESeq2::vst) of a raw count matrix
# (genes x samples, columns in the order of `samples`).
# blind = TRUE: ignores the groups (exploration, pca_pre);
# blind = FALSE: design ~ group, recommended before removing covariates.
# normalization_factors: GC/length factors (gene x sample) from gc_normalize;
# NULL = DESeq2's own size factors (library size only).
vst_counts <- function(count_matrix, samples, blind = TRUE, normalization_factors = NULL) {
  col_data <- samples %>%
    mutate(group = factor(group)) %>%
    tibble::column_to_rownames("specimenID")
  dds <- DESeq2::DESeqDataSetFromMatrix(round(count_matrix), colData = col_data,
                                        design = ~group)
  if (!is.null(normalization_factors)) {
    DESeq2::normalizationFactors(dds) <- normalization_factors_for(
      normalization_factors, rownames(count_matrix), colnames(count_matrix)
    )
  }
  SummarizedExperiment::assay(DESeq2::vst(dds, blind = blind))
}

# VST (with the GC/length factors, if given) + removal of the technical
# covariates (see remove_covariates). Returns log-scale expression, genes x samples.
corrected_expression <- function(count_matrix, samples, covariates,
                                 normalization_factors = NULL) {
  vst_counts(count_matrix, samples, blind = FALSE,
             normalization_factors = normalization_factors) %>%
    remove_covariates(samples, covariates)
}

# Removal of the technical covariates from log-scale expression with
# limma::removeBatchEffect, protecting the biological groups (diagnosis x sex).
# Same covariates as the DESeq2 model: one categorical (model_batch, if
# present) and the numeric ones.
remove_covariates <- function(expression, samples, covariates) {
  col_data <- model_coldata(samples, covariates)

  categorical <- covariates[vapply(col_data[covariates], is.factor, logical(1))]
  numeric <- setdiff(covariates, categorical)
  stopifnot(length(categorical) <= 2) # removeBatchEffect takes up to 2 batch factors

  limma::removeBatchEffect(
    expression,
    batch = if (length(categorical) >= 1) col_data[[categorical[1]]] else NULL,
    batch2 = if (length(categorical) == 2) col_data[[categorical[2]]] else NULL,
    covariates = if (length(numeric) > 0) as.matrix(col_data[numeric]) else NULL,
    design = model.matrix(~group, data = col_data)
  )
}

# PCA of an expression matrix (genes x samples, log scale) on its most variable
# genes, and variance of each PC explained by each covariate (R^2 of PC ~ covariate).
# Returns three tables: scores (PCs + sample columns), variance, association.
pca_tables <- function(expression, samples, region_name, top_genes, n_pcs, covariates) {
  most_variable <- head(order(apply(expression, 1, var), decreasing = TRUE), top_genes)
  pca <- prcomp(t(expression[most_variable, ]))

  variance <- tibble(region = region_name,
                     pc = paste0("PC", seq_along(pca$sdev)),
                     variance = 100 * pca$sdev^2 / sum(pca$sdev^2))

  scores <- pca$x[, seq_len(n_pcs)] %>%
    as_tibble(rownames = "specimenID") %>%
    left_join(samples, by = "specimenID")

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

  list(scores = scores, variance = variance, association = association)
}
