# split_groups.R
# Splits the corrected expression of one region (network_qc) into one matrix per
# group (diagnosis x sex), the input of ARACNe-AP: tab-separated, first column
# 'gene' (Ensembl ID), one column per sample, log-scale expression.
# All the groups of a region have the same genes.
# Also writes the metadata of the samples of each group.
# Output files end in _<group>.tsv (e.g. matrizQC_MAYO_CRB_AD_female.tsv).
# Run by Snakemake (rule split_groups).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

samples <- read_tsv(snakemake@input[["samples"]], show_col_types = FALSE)
expression <- readRDS(snakemake@input[["expression"]])

# The file of one group among the outputs: the one that ends in _<group>.tsv
output_of <- function(files, group_name) {
  file <- grep(paste0("_", group_name, "\\.tsv$"), files, value = TRUE)
  stopifnot(length(file) == 1)
  file
}

for (group_name in unlist(snakemake@params[["groups"]])) {
  group_samples <- samples %>% filter(group == group_name)
  ids <- group_samples$specimenID
  stopifnot(length(ids) > 0, all(ids %in% colnames(expression)))

  expression[, ids, drop = FALSE] %>%
    as_tibble(rownames = "gene") %>%
    write_tsv(output_of(snakemake@output[["matrices"]], group_name))
  write_tsv(group_samples, output_of(snakemake@output[["metadata"]], group_name))

  cat(group_name, ":", nrow(expression), "genes x", length(ids), "samples\n")
}

cat("\nDone\n")
