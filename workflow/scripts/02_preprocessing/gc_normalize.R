# gc_normalize.R
# GC content and length correction of one brain region, as in the old
# 2.pre-promRNA script: EDASeq corrects, NOISeq diagnoses.
#   1. Genes: Sage's filter (at least one group) with GC and length annotation.
#   2. EDASeq, on the raw counts (all AD/CT samples of the region):
#        within-lane length correction (full quantile) ->
#        within-lane GC correction (full quantile) ->
#        between-lane normalization (upper quartile, library size)
#      Order: GC and length are correlated across genes (r ~ -0.47), so the
#      correction done last brings back part of the other bias. Tested on TC
#      (median NOISeq R^2, GC / length): uncorrected 0.36 / 0.66;
#      GC -> length 0.63 / 0.04; length -> GC 0.17 / 0.25, the only order that
#      lowers both.
#   3. Normalization factors for DESeq2 (EDASeq's recommendation): the counts
#      stay raw integers and the correction goes in normalizationFactors(), so
#      deseq2, sample_qc and network_qc all use the same corrected scale.
#      factor = (count + 0.1) / (normalized + 0.1), centered per gene.
#   4. NOISeq diagnostics before and after: R^2 and p of the GC bias and of the
#      length bias of every sample (cubic spline regression of the expression
#      on GC / length bins). High R^2 with low p = the expression depends on it.
# Run by Snakemake (rule gc_normalize).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

# Bioconductor first and the tidyverse last: Bioconductor packages mask dplyr
# functions such as count() and select() if they are loaded after it
suppressPackageStartupMessages({
  library(EDASeq)
  library(NOISeq)
  library(dplyr)
  library(tidyr)
  library(readr)
  source(snakemake@input[["gene_filter"]])
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
count_matrix <- round(count_matrix[, samples$specimenID])

annotation <- read_tsv(snakemake@input[["annotation"]],
                       col_types = cols(chromosome = col_character())) %>%
  filter(!is.na(gc_content), !is.na(length))

# ---- Genes ------------------------------------------------------------------

passes <- filter_genes(count_matrix, samples,
                       snakemake@params[["min_cpm"]], snakemake@params[["min_fraction"]])
genes <- intersect(rownames(count_matrix)[passes], annotation$gene_id)

cat("Region:", region_name, "| samples:", ncol(count_matrix), "\n")
cat("Genes passing the filter:", sum(passes), "| with GC and length annotation:",
    length(genes), "(the rest are left out of every analysis)\n")

counts <- count_matrix[genes, ]
gene_info <- annotation %>%
  filter(gene_id %in% genes) %>%
  tibble::column_to_rownames("gene_id")
gene_info <- gene_info[genes, ]
gc <- setNames(gene_info$gc_content, genes)
gene_length <- setNames(gene_info$length, genes)

# ---- EDASeq: GC, length and library size ------------------------------------

# Length first, then GC (see the order note at the top)
normalized <- counts %>%
  withinLaneNormalization(gene_length, which = "full", round = FALSE) %>%
  withinLaneNormalization(gc, which = "full", round = FALSE) %>%
  betweenLaneNormalization(which = "upper", round = FALSE)

# Normalization factors for DESeq2: count / factor = corrected expression.
# Centered per gene (geometric mean 1), as DESeq2 recommends for gene-wise factors.
factors <- (counts + 0.1) / (normalized + 0.1)
factors <- factors / exp(rowMeans(log(factors)))
saveRDS(factors, snakemake@output[["factors"]])

# ---- NOISeq diagnostics -----------------------------------------------------

noiseq_factors <- data.frame(group = samples$group, row.names = samples$specimenID)

# R^2 and p of every sample's spline regression, from a NOISeq bias object
bias_table <- function(bias, bias_name, stage) {
  models <- bias@dat$RegressionModels
  tibble(specimenID = names(models),
         fit = lapply(models, \(m) if (inherits(m, "summary.lm")) m else summary(m))) %>%
    mutate(region = region_name,
           bias = bias_name,
           stage = stage,
           r2 = vapply(fit, \(s) s$r.squared, numeric(1)),
           p = vapply(fit, \(s) {
             f <- s$fstatistic
             pf(f[1], f[2], f[3], lower.tail = FALSE)
           }, numeric(1))) %>%
    dplyr::select(region, specimenID, bias, stage, r2, p)
}

diagnose <- function(matrix, normalized_input, stage) {
  data <- readData(data = matrix, factors = noiseq_factors,
                   gc = gc, length = gene_length)
  gc_bias <- dat(data, type = "GCbias", k = 0, factor = NULL, norm = normalized_input)
  length_bias <- dat(data, type = "lengthbias", k = 0, factor = NULL, norm = normalized_input)
  bind_rows(bias_table(gc_bias, "GC", stage), bias_table(length_bias, "length", stage))
}

bias <- bind_rows(diagnose(counts, FALSE, "before"),
                  diagnose(normalized, TRUE, "after"))

cat("\nMedian R^2 of the biases across samples:\n")
print(as.data.frame(bias %>%
                      group_by(bias, stage) %>%
                      summarise(median_r2 = median(r2), samples_p_below_0.05 = sum(p < 0.05),
                                .groups = "drop")))

write_tsv(bias, snakemake@output[["bias"]])

cat("\nDone\n")
