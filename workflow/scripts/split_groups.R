# split_groups.R
# Splits the corrected expression of one region (network_qc) into one matrix per
# group (diagnosis x sex), the input of ARACNe-AP: tab-separated, first column
# 'gene' (Ensembl ID), one column per sample, log-scale expression.
# All the groups of a region have the same genes.
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

output_of <- function(group_name) {
  grep(paste0("/", group_name, "\\.tsv$"), snakemake@output[["matrices"]], value = TRUE)
}

for (group_name in unlist(snakemake@params[["groups"]])) {
  ids <- samples %>% filter(group == group_name) %>% pull(specimenID)
  stopifnot(length(ids) > 0, all(ids %in% colnames(expression)))

  expression[, ids, drop = FALSE] %>%
    as_tibble(rownames = "gene") %>%
    write_tsv(output_of(group_name))

  cat(group_name, ":", nrow(expression), "genes x", length(ids), "samples\n")
}

cat("\nDone\n")
