# prepare_counts.R
# Raw count matrix and sample table of one brain region, shared by the
# differential expression and the network branches:
#   - samples of the region with diagnosis AD or CT (Sage's harmonized diagnosis)
#   - one sample per individual (the one with the highest RIN)
#   - raw counts of those samples, from the study's count matrices
# Genes are NOT filtered here: each branch applies its own filter.
# Run by Snakemake (rule prepare_counts).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

region_name <- snakemake@wildcards[["region"]]
diagnoses <- unlist(snakemake@params[["diagnoses"]])

region_info <- read_tsv(snakemake@input[["regions"]], show_col_types = FALSE) %>%
  filter(region == region_name)

# ---- Samples ----------------------------------------------------------------

samples <- read_tsv(snakemake@input[["covariates"]],
                    col_types = cols(.default = col_character())) %>%
  filter(tissue == region_info$covariates_tissue,
         diagnosis %in% diagnoses) %>%
  transmute(
    specimenID,
    individualID,
    region = region_name,
    study = region_info$study,
    diagnosis,
    sex,
    group = paste(diagnosis, sex, sep = "_"),
    batch = .data[[region_info$batch_column]],
    RIN = as.numeric(RIN),
    pmi = as.numeric(pmi),
    age_death = parse_number(age_death), # "90+" (censored) becomes 90
    age_censored = age_death >= 90,
    apoe4_allele,
    pct_aligned = as.numeric(AlignmentSummaryMetrics_PCT_PF_READS_ALIGNED),
    pct_intronic = as.numeric(RnaSeqMetrics_PCT_INTRONIC_BASES),
    pct_intergenic = as.numeric(RnaSeqMetrics_PCT_INTERGENIC_BASES),
    pct_coding = as.numeric(RnaSeqMetrics_PCT_CODING_BASES),
    # Mayo cerebellum is "_CBE" in Sage's covariates but "_CER" in the count matrix
    counts_id = if_else(study == "Mayo", sub("_CBE$", "_CER", specimenID), specimenID)
  )

cat("Samples of", region_name, "with diagnosis", paste(diagnoses, collapse = "/"), ":",
    nrow(samples), "\n")

# One sample per individual: the one with the highest RIN
samples <- samples %>%
  group_by(individualID) %>%
  slice_max(RIN, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("After keeping one sample per individual:", nrow(samples), "\n")

# ---- Sequencing protocol and batch for the model ----------------------------

# Protocol = library prep + kit, from the original metadata. In ROSMAP it
# separates the main clusters of the PCA (e.g. polyA vs rRNA depletion).
# Read length is left out: Kapa samples with 100 and 125 bp reads cluster together.
protocols <- read_csv(snakemake@input[["metadata"]],
                      col_types = cols(.default = col_character())) %>%
  filter(assay == "rnaSeq") %>%
  distinct(specimenID, .keep_all = TRUE) %>%
  transmute(counts_id = specimenID,
            protocol = paste(libraryPrep, libraryPreparationMethod, sep = "_"))

samples <- left_join(samples, protocols, by = "counts_id")
stopifnot(!anyNA(samples$protocol))

# Batches with fewer than min_batch_size samples carry almost no information in
# the model: they are merged into the largest batch of the same protocol
min_batch_size <- snakemake@params[["min_batch_size"]]

samples <- samples %>%
  add_count(batch, name = "batch_size") %>%
  group_by(protocol) %>%
  mutate(model_batch = if_else(batch_size < min_batch_size,
                               batch[which.max(batch_size)],
                               batch)) %>%
  ungroup() %>%
  select(-batch_size)

cat("\nBatches merged into model_batch (fewer than", min_batch_size, "samples):\n")
print(as.data.frame(samples %>% filter(batch != model_batch) %>% count(protocol, batch, model_batch)))
cat("\nSamples per protocol and model_batch:\n")
print(as.data.frame(count(samples, protocol, model_batch)))

# ---- Counts -----------------------------------------------------------------

# Reads only the columns of this region's samples
read_counts <- function(path) {
  read_tsv(path,
           col_select = any_of(c("feature", samples$counts_id)),
           col_types = cols(feature = col_character(), .default = col_double()),
           progress = FALSE)
}

counts <- snakemake@input[["counts"]] %>%
  lapply(read_counts) %>%
  Reduce(function(x, y) inner_join(x, y, by = "feature"), .) %>%
  # Drop STAR summary rows (N_unmapped...) and the PAR copies on chrY
  filter(startsWith(feature, "ENSG"), !endsWith(feature, "_PAR_Y")) %>%
  mutate(gene_id = sub("\\..*$", "", feature), .before = 1) %>%
  select(-feature)

stopifnot(!anyDuplicated(counts$gene_id))

missing <- setdiff(samples$counts_id, names(counts))
if (length(missing) > 0) {
  cat("Samples without counts, removed:", paste(missing, collapse = ", "), "\n")
  samples <- filter(samples, !counts_id %in% missing)
}

# Count columns named by specimenID, in the order of the sample table
counts <- counts %>%
  rename(all_of(setNames(samples$counts_id, samples$specimenID))) %>%
  select(gene_id, all_of(samples$specimenID))

cat("\nFinal samples per group:\n")
print(as.data.frame(count(samples, group)))
cat("\nGenes:", nrow(counts), "\n")

saveRDS(counts, snakemake@output[["counts"]])
write_tsv(samples, snakemake@output[["samples"]])
