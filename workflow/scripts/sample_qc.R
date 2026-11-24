# sample_qc.R
# Outlier samples of one brain region. The samples kept here are used by both
# branches: differential expression (deseq2) and networks (network_qc).
#   1. Expression: genes of gc_normalize (Sage's filter, at least one group),
#      VST with the GC/length normalization factors and removal of the
#      technical covariates of the DESeq2 model, protecting the groups.
#   2. Standardized connectivity (Oldham et al. 2012, the usual criterion for
#      co-expression networks): adjacency between samples = ((1 + r) / 2)^2,
#      r = Pearson correlation; connectivity = sum of a sample's adjacencies;
#      z = connectivity standardized. z < max_connectivity_z (e.g. -3): the
#      sample is unlike the rest -> outlier.
#   3. Sex check (optional): predicted sex from XIST vs UTY (raw CPM); samples
#      whose annotated sex does not match are removed.
# Writes every sample with its QC columns (outlier, reason) and the kept ones.
# Run by Snakemake (rule sample_qc).

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
max_z <- snakemake@params[["max_connectivity_z"]]
check_sex <- isTRUE(snakemake@params[["check_sex"]])

samples <- read_tsv(snakemake@input[["samples"]],
                    col_types = cols(individualID = col_character(),
                                     batch = col_character(),
                                     model_batch = col_character(),
                                     apoe4_allele = col_character()))
count_matrix <- readRDS(snakemake@input[["counts"]]) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()
count_matrix <- count_matrix[, samples$specimenID]

cat("Region:", region_name, "| samples:", nrow(samples), "\n")

# ---- Standardized connectivity ----------------------------------------------

# Genes of the GC/length correction (Sage's filter, with annotation)
factors <- readRDS(snakemake@input[["factors"]])
genes <- rownames(factors)
expression <- corrected_expression(count_matrix[genes, ], samples, covariates,
                                   normalization_factors = factors)

adjacency <- ((1 + cor(expression)) / 2)^2
connectivity <- colSums(adjacency) - 1 # without the sample itself
connectivity_z <- setNames(as.numeric(scale(connectivity)), colnames(expression))

# ---- Sex check --------------------------------------------------------------

sex_genes <- c(XIST = "ENSG00000229807", UTY = "ENSG00000183878")
cpm <- t(t(count_matrix[sex_genes, , drop = FALSE]) / colSums(count_matrix)) * 1e6

# Library size and detected genes (CPM > 1) of every sample, for the QC panel
library_size <- colSums(count_matrix)
detected_genes <- colSums(t(t(count_matrix) / library_size) * 1e6 > 1)

qc <- samples %>%
  mutate(library_size = library_size[specimenID],
         detected_genes = detected_genes[specimenID],
         connectivity_z = connectivity_z[specimenID],
         XIST_cpm = cpm[sex_genes[["XIST"]], specimenID],
         UTY_cpm = cpm[sex_genes[["UTY"]], specimenID],
         predicted_sex = if_else(log2(UTY_cpm + 1) > log2(XIST_cpm + 1), "male", "female"),
         low_connectivity = connectivity_z < max_z,
         sex_mismatch = check_sex & predicted_sex != sex,
         outlier = low_connectivity | sex_mismatch,
         reason = case_when(
           low_connectivity & sex_mismatch ~ "low connectivity; sex mismatch",
           low_connectivity ~ "low connectivity",
           sex_mismatch ~ "sex mismatch",
           .default = NA_character_
         ))

cat("\nOutliers (connectivity z <", max_z, if (check_sex) "or sex mismatch", "):\n")
print(as.data.frame(qc %>%
                      filter(outlier) %>%
                      dplyr::select(specimenID, group, connectivity_z, sex, predicted_sex, reason)))
cat("\nSamples kept per group:\n")
print(as.data.frame(qc %>% filter(!outlier) %>% dplyr::count(group)))

write_tsv(qc, snakemake@output[["qc"]])
write_tsv(qc %>% filter(!outlier), snakemake@output[["kept"]])

cat("\nDone\n")
