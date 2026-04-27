library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

# =============================================================================
# LOGIC OVERVIEW
# =============================================================================
# Goal: find genes with CPM < 1 in ALL women and CPM >= 1 in >= 90% of men,
# per brain region — independently of disease phenotype.
#
# CPM threshold of 1 is the standard "expressed/not expressed" cutoff,
# consistent with the QC filtering applied upstream in this pipeline.
#
# Raw counts are used (not TMM-normalized) because TMM produces continuous
# values with no true zeros, making the absent/present criterion meaningless.
#
# filteredQC metadata is used for all ROSMAP regions to restrict to samples
# that passed quality control, avoiding false "absent in females" calls from
# low-quality samples with globally depressed expression.
#
# Steps per region:
#   1. Load raw count matrix, remove STAR summary rows, clean ENSEMBL suffixes
#   2. Load filteredQC metadata and restrict specimenIDs to those in the matrix
#   3. Split samples by sex, compute CPM per group
#   4. Apply criteria: CPM < 1 in ALL females, CPM >= 1 in >= 90% of males
#   5. Translate ENSEMBL -> HGNC symbol (local annotation first, org.Hs.eg.db fallback)
# =============================================================================

# ---- Paths ----------------------------------------------------------------

base_rosmap  <- "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue"
base_qc      <- "/datos/rosmap/data_by_counts/new_AD_TF_sex"
base_mayo    <- "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue"

# Local annotation: ENSEMBL -> HGNC symbol (from DLPFC QC pipeline)
annotation_path <- file.path(base_rosmap,
                             "metadata/DLPFC/RNA_seq_filteredQC_annotation_DLPFC.txt")

# Region definitions:
#   counts     = raw count matrix (.rds)
#   meta       = filteredQC metadata with specimenID and sex column
#   sex_col    = name of the sex column
#   female_val = value indicating female
#   male_val   = value indicating male
regions <- list(
  
  DLPFC = list(
    counts    = file.path(base_rosmap, "DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds"),
    meta      = file.path(base_rosmap, "metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt"),
    sex_col   = "sex", female_val = 0, male_val = 1
  ),
  
  HCN = list(
    counts    = file.path(base_rosmap, "HCN/ROSMAP_RNAseq_rawcounts_HCN.rds"),
    meta      = file.path(base_qc, "HCN/RNA_seq_metadata_filteredQC_HCN.txt"),
    sex_col   = "sex", female_val = 0, male_val = 1
  ),
  
  PCC = list(
    counts    = file.path(base_rosmap, "PCC/ROSMAP_RNAseq_rawcounts_PCC.rds"),
    meta      = file.path(base_qc, "PCC/RNA_seq_metadata_filteredQC_PCC.txt"),
    sex_col   = "sex", female_val = 0, male_val = 1
  ),
  
  CRB = list(
    counts    = file.path(base_mayo, "CRB/Mayo_RNAseq_rawcounts_CRB.rds"),
    meta      = file.path(base_mayo, "CRB/Mayo_RNAseq_metadata_CRB.txt"),
    sex_col   = "sex", female_val = "female", male_val = "male"
  ),
  
  TC = list(
    counts    = file.path(base_mayo, "TC/Mayo_RNAseq_rawcounts_TC.rds"),
    meta      = file.path(base_mayo, "TC/Mayo_RNAseq_metadata_TC.txt"),
    sex_col   = "sex", female_val = "female", male_val = "male"
  )
)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# ---- Load local ENSEMBL -> HGNC symbol annotation -------------------------

load_annotation <- function(path) {
  vroom::vroom(path, show_col_types = FALSE) %>%
    select(feature, hgnc_symbol) %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
    distinct(feature, .keep_all = TRUE)
}

# ---- Translate ENSEMBL IDs to HGNC symbols --------------------------------
# Primary:  local annotation file from QC pipeline
# Fallback: org.Hs.eg.db for any IDs not covered by the local file

translate_symbols <- function(ensembl_ids, annotation) {
  
  # Step 1: local lookup
  result <- tibble(ensembl_id = ensembl_ids) %>%
    left_join(annotation, by = c("ensembl_id" = "feature"))
  names(result)[names(result) == "hgnc_symbol"] <- "symbol"
  
  n_local <- sum(!is.na(result$symbol))
  
  # Step 2: org.Hs.eg.db fallback for remaining NAs
  missing_ids <- result %>% filter(is.na(symbol)) %>% pull(ensembl_id)
  
  if (length(missing_ids) > 0) {
    fallback <- mapIds(
      org.Hs.eg.db,
      keys      = missing_ids,
      column    = "SYMBOL",
      keytype   = "ENSEMBL",
      multiVals = "first"
    ) %>%
      enframe(name = "ensembl_id", value = "symbol_fallback")
    
    result <- result %>%
      left_join(fallback, by = "ensembl_id") %>%
      mutate(symbol = coalesce(symbol, symbol_fallback)) %>%
      select(-symbol_fallback)
  }
  
  n_fallback <- sum(!is.na(result$symbol)) - n_local
  n_missing  <- sum(is.na(result$symbol))
  
  message("  Symbol mapping — local: ", n_local,
          " | org.Hs.eg.db fallback: ", n_fallback,
          " | unresolved: ", n_missing)
  
  result$symbol
}

# ---- Load and clean raw count matrix --------------------------------------

load_raw_counts <- function(counts_path) {
  
  raw <- readRDS(counts_path)
  
  # Move "feature" column to rownames if present
  if ("feature" %in% colnames(raw)) {
    raw <- raw %>% column_to_rownames("feature")
  }
  
  raw <- as.matrix(raw)
  
  # Remove STAR alignment summary rows
  aligner_rows <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
  raw <- raw[!rownames(raw) %in% aligner_rows, , drop = FALSE]
  
  # Strip ENSEMBL version suffixes: ENSG00000223972.5 -> ENSG00000223972
  rownames(raw) <- str_remove(rownames(raw), "\\.\\d+$")
  
  raw
}

# ---- Split matrix by sex and compute CPM ----------------------------------

split_by_sex <- function(raw, meta_path, sex_col, female_val, male_val) {
  
  meta <- vroom::vroom(meta_path, show_col_types = FALSE) %>%
    filter(!is.na(.data[[sex_col]]))
  
  female_ids <- meta %>% filter(.data[[sex_col]] == female_val) %>% pull(specimenID)
  male_ids   <- meta %>% filter(.data[[sex_col]] == male_val)   %>% pull(specimenID)
  
  # Restrict to samples present in the count matrix
  female_ids <- intersect(female_ids, colnames(raw))
  male_ids   <- intersect(male_ids,   colnames(raw))
  
  message("  Female samples: ", length(female_ids),
          " | Male samples: ", length(male_ids))
  
  if (length(female_ids) == 0 | length(male_ids) == 0) {
    stop("No samples found for one or both sexes — check sex_col and values")
  }
  
  list(
    female_cpm = cpm(raw[, female_ids, drop = FALSE], log = FALSE),
    male_cpm   = cpm(raw[, male_ids,   drop = FALSE], log = FALSE)
  )
}

# ---- Apply absent-in-females / present-in-males criteria ------------------

find_absent_in_females <- function(female_cpm, male_cpm,
                                   cpm_threshold  = 1,
                                   male_threshold = 0.90) {
  
  shared_genes <- intersect(rownames(female_cpm), rownames(male_cpm))
  female_cpm   <- female_cpm[shared_genes, , drop = FALSE]
  male_cpm     <- male_cpm[shared_genes,   , drop = FALSE]
  
  absent_in_all_females <- rowSums(female_cpm >= cpm_threshold) == 0
  present_in_most_males <- rowSums(male_cpm   >= cpm_threshold) / ncol(male_cpm) >= male_threshold
  
  message("  Genes absent in ALL females (CPM < 1):      ", sum(absent_in_all_females))
  message("  Genes present in >=90% of males (CPM >= 1): ", sum(present_in_most_males))
  message("  Genes satisfying BOTH criteria:             ",
          sum(absent_in_all_females & present_in_most_males))
  
  hits <- shared_genes[absent_in_all_females & present_in_most_males]
  
  if (length(hits) == 0) return(tibble())
  
  tibble(
    ensembl_id          = hits,
    n_female_samples    = ncol(female_cpm),
    n_male_samples      = ncol(male_cpm),
    pct_males_expressed = round(
      rowSums(male_cpm[hits, , drop = FALSE] >= cpm_threshold) / ncol(male_cpm) * 100, 1)
  )
}

# =============================================================================
# LOAD ANNOTATION (once, reused for all regions)
# =============================================================================

message("Loading annotation table...")
annotation <- load_annotation(annotation_path)
message("  Annotation entries: ", nrow(annotation))

# =============================================================================
# RUN PER REGION
# =============================================================================

results <- imap(regions, function(r, region_name) {
  
  message("\n========== ", region_name, " ==========")
  
  # Step 1: load and clean raw counts
  raw <- load_raw_counts(r$counts)
  message("  Raw matrix: ", nrow(raw), " genes x ", ncol(raw), " samples")
  message("  Count zeros: ", round(mean(raw == 0) * 100, 1), "%")
  message("  Preview (5x5):")
  print(raw[1:5, 1:5])
  
  # Step 2: split by sex and compute CPM
  sex_data <- split_by_sex(raw, r$meta, r$sex_col, r$female_val, r$male_val)
  
  message("  CPM female preview (5x3):")
  print(round(sex_data$female_cpm[1:5, 1:3], 2))
  
  # Step 3: apply absent/present criteria
  hits <- find_absent_in_females(sex_data$female_cpm, sex_data$male_cpm)
  
  if (nrow(hits) == 0) {
    message("  No genes found for ", region_name)
    return(tibble())
  }
  
  # Step 4: translate ENSEMBL -> SYMBOL
  hits %>%
    mutate(
      symbol = translate_symbols(ensembl_id, annotation),
      region = region_name
    ) %>%
    select(region, ensembl_id, symbol, everything())
})

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

results_df <- bind_rows(results)

write_csv(results_df, "sex_absent_genes_per_region.csv")

# ---- Summary ---------------------------------------------------------------

message("\n===== SUMMARY PER REGION =====")
results_df %>%
  count(region, name = "n_genes") %>%
  print()

message("\nGenes found in more than one region:")
results_df %>%
  count(ensembl_id, symbol, name = "n_regions") %>%
  filter(n_regions > 1) %>%
  arrange(desc(n_regions)) %>%
  print()