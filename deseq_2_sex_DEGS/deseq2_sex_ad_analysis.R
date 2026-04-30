#!/usr/bin/env Rscript
# ==============================================================================
# Differential Expression Analysis: Sex-stratified and Interaction Effects in AD
# ==============================================================================
# Purpose: Run three DESeq2 contrasts:
#   1. AD vs Control in MALES only
#   2. AD vs Control in FEMALES only  
#   3. Sex:Diagnosis INTERACTION (differential AD effect between sexes)
#
# Usage: Rscript deseq2_sex_ad_analysis.R \
#          --counts counts_matrix.rds \
#          --metadata sample_metadata.csv \
#          --outdir results/ \
#          --brain_area DLPFC \
#          --padj 0.05 \
#          --lfc 0.5
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(optparse)
})

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

option_list <- list(
  make_option(c("--counts"), type = "character", default = NULL,
              help = "Path to counts matrix (RDS or CSV)"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Path to sample metadata (CSV/TSV)"),
  make_option(c("--outdir"), type = "character", default = "results",
              help = "Output directory [default: results]"),
  make_option(c("--brain_area"), type = "character", default = "unknown",
              help = "Brain area label for output files"),
  make_option(c("--padj"), type = "numeric", default = 0.05,
              help = "Adjusted p-value threshold [default: 0.05]"),
  make_option(c("--lfc"), type = "numeric", default = 0.5,
              help = "Log2 fold-change threshold [default: 0.5]"),
  make_option(c("--cores"), type = "integer", default = 1,
              help = "Number of cores for parallel processing [default: 1]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$counts) || is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("Both --counts and --metadata are required.", call. = FALSE)
}

# Create output directory
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Set parallel backend if requested
if (opt$cores > 1) {
  BiocParallel::register(BiocParallel::MulticoreParam(opt$cores))
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Load and validate counts matrix
#' @param path Path to counts file (RDS or CSV/TSV)
#' @return Matrix with genes in rows, samples in columns
load_counts <- function(path) {
  message("Loading counts from: ", path)
  
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    counts <- readRDS(path)
  } else if (grepl("\\.(csv|tsv|txt)$", path, ignore.case = TRUE)) {
    sep <- if (grepl("\\.csv$", path)) "," else "\t"
    counts <- read.table(path, header = TRUE, row.names = 1, sep = sep, check.names = FALSE)
    counts <- as.matrix(counts)
  } else {
    stop("Counts file must be .rds, .csv, or .tsv")
  }
  
  # Validate
  if (!is.numeric(counts)) stop("Counts must be numeric")
  if (any(counts < 0, na.rm = TRUE)) stop("Counts contain negative values")
  if (any(is.na(counts))) {
    warning("Counts contain NA values, replacing with 0")
    counts[is.na(counts)] <- 0
  }
  
  message(sprintf("  Loaded %d genes x %d samples", nrow(counts), ncol(counts)))
  return(counts)
}

#' Load and validate metadata
#' @param path Path to metadata file
#' @param required_cols Required column names
#' @return Data frame with validated metadata
load_metadata <- function(path, required_cols = c("specimenID", "sex", "diagnosis")) {
  message("Loading metadata from: ", path)
  
  sep <- if (grepl("\\.csv$", path)) "," else "\t"
  metadata <- read.table(path, header = TRUE, sep = sep, stringsAsFactors = FALSE)
  
  # Check required columns
  missing <- setdiff(required_cols, colnames(metadata))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Standardize sex values
  metadata$sex <- tolower(metadata$sex)
  metadata$sex[metadata$sex %in% c("m", "male")] <- "male"
  metadata$sex[metadata$sex %in% c("f", "female")] <- "female"
  
  # Standardize diagnosis values  
  metadata$diagnosis <- tolower(metadata$diagnosis)
  metadata$diagnosis[metadata$diagnosis %in% c("ctrl", "control", "normal")] <- "control"
  metadata$diagnosis[metadata$diagnosis %in% c("ad", "alzheimer", "alzheimers")] <- "AD"
  
  message(sprintf("  Loaded metadata for %d samples", nrow(metadata)))
  message(sprintf("  Sex distribution: %s", 
                  paste(table(metadata$sex), collapse = " / ")))
  message(sprintf("  Diagnosis distribution: %s", 
                  paste(table(metadata$diagnosis), collapse = " / ")))
  
  return(metadata)
}

#' Scale numeric covariates if mean or SD > 5
#' @param coldata Sample metadata data frame
#' @param vars Variables to check and scale
#' @return Modified coldata with scaled variables
scale_covariates <- function(coldata, vars = c("RIN", "pmi")) {
  for (var in vars) {
    if (!var %in% colnames(coldata)) {
      warning(sprintf("Variable '%s' not found in metadata, skipping", var))
      next
    }
    
    val_mean <- mean(coldata[[var]], na.rm = TRUE)
    val_sd <- sd(coldata[[var]], na.rm = TRUE)
    
    if (val_mean >= 5 || val_sd >= 5) {
      message(sprintf("  Scaling %s (mean=%.2f, sd=%.2f)", var, val_mean, val_sd))
      coldata[[var]] <- as.numeric(scale(coldata[[var]]))
    }
  }
  return(coldata)
}

#' Convert Ensembl IDs to gene symbols
#' @param ensembl_ids Vector of Ensembl gene IDs
#' @return Data frame with ensembl_gene_id and external_gene_name
convert_ens_to_symbol <- function(ensembl_ids) {
  # This is a placeholder - implement with biomaRt or org.Hs.eg.db
  # For now, return a minimal data frame
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys      = ensembl_ids,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  
  data.frame(
    ensembl_gene_id    = names(symbols),
    external_gene_name = unname(symbols),
    stringsAsFactors   = FALSE
  )
}

# ==============================================================================
# CORE DESEQ2 FUNCTIONS
# ==============================================================================

#' Run DESeq2 analysis with specified design
#' @param counts Count matrix (genes x samples)
#' @param metadata Sample metadata
#' @param design_formula Design formula as string
#' @param subset_condition Optional: subset metadata before analysis (e.g., "sex == 'male'")
#' @return List with dds object, coldata, and design formula
run_deseq <- function(counts, metadata, design_formula, subset_condition = NULL) {
  
  # Prepare coldata
  coldata <- data.frame(specimenID = colnames(counts), stringsAsFactors = FALSE) %>%
    left_join(metadata, by = "specimenID") %>%
    column_to_rownames("specimenID")
  
  # Apply subset if specified
  if (!is.null(subset_condition)) {
    message(sprintf("  Subsetting samples: %s", subset_condition))
    n_before <- nrow(coldata)
    coldata <- coldata %>% filter(eval(parse(text = subset_condition)))
    n_after <- nrow(coldata)
    message(sprintf("    Kept %d / %d samples", n_after, n_before))
    
    # Subset counts to match
    counts <- counts[, rownames(coldata), drop = FALSE]
  }
  
  # Convert to factors
  coldata$sex <- factor(coldata$sex, levels = c("female", "male"))
  coldata$diagnosis <- factor(coldata$diagnosis, levels = c("control", "AD"))
  
  # Additional factor conversions if present
  if ("specimenIdSource" %in% colnames(coldata)) {
    coldata$specimenIdSource <- as.factor(coldata$specimenIdSource)
  }
  
  # Drop unused levels
  coldata <- droplevels(coldata)
  
  # Scale covariates
  coldata <- scale_covariates(coldata, vars = c("RIN", "pmi"))
  
  # Build design formula
  design <- as.formula(design_formula)
  message(sprintf("  Design: %s", design_formula))
  message(sprintf("  Final sample size: %d", nrow(coldata)))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = design
  )
  
  # Run DESeq2
  message("  Running DESeq2...")
  dds <- DESeq(dds, parallel = opt$cores > 1)
  
  return(list(dds = dds, coldata = coldata, design_formula = design))
}

#' Extract DEGs from DESeq2 results
#' @param dds DESeqDataSet object
#' @param coldata Sample metadata (for baseMean calculations)
#' @param coef_name Coefficient name (for results extraction)
#' @param contrast_vec Contrast vector (alternative to coef_name)
#' @param group_info Named list describing the groups (for baseMean labeling)
#' @param padj_thresh Adjusted p-value threshold
#' @param lfc_thresh Log2 fold-change threshold
#' @return List with res, res.df, DEGS, and normcounts
extract_degs <- function(dds, coldata, coef_name = NULL, contrast_vec = NULL,
                         group_info = NULL, padj_thresh = 0.05, lfc_thresh = 0.5) {
  
  # Extract results
  if (!is.null(coef_name)) {
    message(sprintf("  Extracting coefficient: %s", coef_name))
    res <- results(dds, name = coef_name, pAdjustMethod = "BH", alpha = padj_thresh)
  } else if (!is.null(contrast_vec)) {
    message(sprintf("  Extracting contrast: %s vs %s", 
                    contrast_vec[3], contrast_vec[2]))
    res <- results(dds, contrast = contrast_vec, pAdjustMethod = "BH", alpha = padj_thresh)
  } else {
    stop("Must provide either coef_name or contrast_vec")
  }
  
  # Convert to data frame
  res.df <- as.data.frame(res) %>%
    filter(!is.na(padj)) %>%
    rownames_to_column("ensembl_gene_id")
  
  # Classify genes
  res.df$diffexpressed <- "NO"
  res.df$diffexpressed[res.df$log2FoldChange >  lfc_thresh & res.df$padj < padj_thresh] <- "UP"
  res.df$diffexpressed[res.df$log2FoldChange < -lfc_thresh & res.df$padj < padj_thresh] <- "DOWN"
  
  # Add gene symbols
  symbols <- convert_ens_to_symbol(res.df$ensembl_gene_id)
  res.df <- res.df %>% left_join(symbols, by = "ensembl_gene_id")
  
  # Calculate baseMean per group
  normcounts <- counts(dds, normalized = TRUE)
  
  if (!is.null(group_info)) {
    # Calculate group-specific baseMeans
    basemean_list <- lapply(names(group_info), function(grp) {
      condition <- group_info[[grp]]
      idx <- eval(parse(text = condition), envir = as.list(coldata))
      if (sum(idx) == 0) {
        warning(sprintf("No samples found for group: %s", grp))
        return(rep(NA, nrow(normcounts)))
      }
      rowMeans(normcounts[, idx, drop = FALSE])
    })
    
    basemean_df <- data.frame(
      ensembl_gene_id = rownames(normcounts),
      do.call(cbind, basemean_list)
    )
    colnames(basemean_df)[-1] <- paste0("baseMean_", names(group_info))
    
    res.df <- res.df %>% left_join(basemean_df, by = "ensembl_gene_id")
  }
  
  # Filter DEGs
  degs <- res.df %>% filter(diffexpressed != "NO")
  message(sprintf("  DEGs found: %d UP, %d DOWN (total: %d)",
                  sum(degs$diffexpressed == "UP"),
                  sum(degs$diffexpressed == "DOWN"),
                  nrow(degs)))
  
  # Add labels
  res.df$delabel <- ifelse(res.df$diffexpressed != "NO", res.df$external_gene_name, NA)
  degs$delabel <- degs$external_gene_name
  
  return(list(
    res = res,
    res.df = res.df,
    degs = degs,
    normcounts = normcounts
  ))
}

#' Save results to files
#' @param results_list Output from extract_degs
#' @param prefix Output file prefix
#' @param outdir Output directory
save_results <- function(results_list, prefix, outdir) {
  # Save full results
  write.csv(results_list$res.df, 
            file.path(outdir, paste0(prefix, "_all_genes.csv")),
            row.names = FALSE)
  
  # Save DEGs only
  write.csv(results_list$degs,
            file.path(outdir, paste0(prefix, "_DEGs.csv")),
            row.names = FALSE)
  
  # Save as RDS for downstream analysis
  saveRDS(results_list, file.path(outdir, paste0(prefix, "_results.rds")))
  
  message(sprintf("  Saved results to: %s", file.path(outdir, prefix)))
}

#' Create volcano plot
#' @param res.df Results data frame
#' @param title Plot title
#' @param outfile Output file path
#' @param padj_thresh Adjusted p-value threshold line
#' @param lfc_thresh Log2 fold-change threshold line
create_volcano <- function(res.df, title, outfile, padj_thresh = 0.05, lfc_thresh = 0.5) {
  
  # Color palette
  colors <- c("UP" = "#E64B35", "DOWN" = "#4DBBD5", "NO" = "grey70")
  
  # Plot
  p <- ggplot(res.df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey30") +
    labs(title = title, x = "log2 Fold Change", y = "-log10(adj. p-value)") +
    theme_bw() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
  message(sprintf("  Saved volcano plot: %s", outfile))
}

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

main <- function() {
  message("\n========================================")
  message("DESeq2 Sex-Stratified Analysis Pipeline")
  message("========================================\n")
  
  # Load data
  counts <- load_counts(opt$counts)
  metadata <- load_metadata(opt$metadata)
  
  # Align samples
  common_samples <- intersect(colnames(counts), metadata$specimenID)
  message(sprintf("\nAligning data: %d common samples", length(common_samples)))
  counts <- counts[, common_samples]
  metadata <- metadata %>% filter(specimenID %in% common_samples)
  
  # ============================================================================
  # CONTRAST 1: AD vs Control in MALES only
  # ============================================================================
  message("\n\n=== CONTRAST 1: AD vs Control (MALES only) ===\n")
  
  male_dds <- run_deseq(
    counts = counts,
    metadata = metadata,
    design_formula = "~ specimenIdSource + pmi + diagnosis",
    subset_condition = "sex == 'male'"
  )
  
  male_results <- extract_degs(
    dds = male_dds$dds,
    coldata = male_dds$coldata,
    coef_name = "diagnosis_AD_vs_control",
    group_info = list(
      control = "diagnosis == 'control'",
      AD = "diagnosis == 'AD'"
    ),
    padj_thresh = opt$padj,
    lfc_thresh = opt$lfc
  )
  
  save_results(male_results,
               prefix = sprintf("%s_males_AD_vs_control", opt$brain_area),
               outdir = opt$outdir)
  
  create_volcano(male_results$res.df,
                 title = sprintf("%s - Males: AD vs Control", opt$brain_area),
                 outfile = file.path(opt$outdir, sprintf("%s_males_volcano.png", opt$brain_area)),
                 padj_thresh = opt$padj,
                 lfc_thresh = opt$lfc)
  
  # ============================================================================
  # CONTRAST 2: AD vs Control in FEMALES only
  # ============================================================================
  message("\n\n=== CONTRAST 2: AD vs Control (FEMALES only) ===\n")
  
  female_dds <- run_deseq(
    counts = counts,
    metadata = metadata,
    design_formula = "~ specimenIdSource + pmi + diagnosis",
    subset_condition = "sex == 'female'"
  )
  
  female_results <- extract_degs(
    dds = female_dds$dds,
    coldata = female_dds$coldata,
    coef_name = "diagnosis_AD_vs_control",
    group_info = list(
      control = "diagnosis == 'control'",
      AD = "diagnosis == 'AD'"
    ),
    padj_thresh = opt$padj,
    lfc_thresh = opt$lfc
  )
  
  save_results(female_results,
               prefix = sprintf("%s_females_AD_vs_control", opt$brain_area),
               outdir = opt$outdir)
  
  create_volcano(female_results$res.df,
                 title = "Females: AD vs Control",
                 outfile = file.path(opt$outdir, sprintf("%s_females_volcano.png", opt$brain_area)),
                 padj_thresh = opt$padj,
                 lfc_thresh = opt$lfc)
  
  # ============================================================================
  # CONTRAST 3: SEX:DIAGNOSIS INTERACTION
  # ============================================================================
  message("\n\n=== CONTRAST 3: Sex × Diagnosis Interaction ===\n")
  message("This tests: Does AD affect males differently than females?\n")
  
  interaction_dds <- run_deseq(
    counts = counts,
    metadata = metadata,
    design_formula = "~ specimenIdSource + pmi + diagnosis + sex + sex:diagnosis",
    subset_condition = NULL  # Use all samples
  )
  
  # Auto-detect interaction coefficient (name varies by factor level order)
  all_coefs <- resultsNames(interaction_dds$dds)
  message("\nAvailable coefficients:")
  print(all_coefs)
  
  interaction_coef <- all_coefs[grepl("diagnosis.*sex|sex.*diagnosis", all_coefs, ignore.case = TRUE)]
  if (length(interaction_coef) == 0) stop("Could not find interaction coefficient. Check resultsNames(dds).")
  message(sprintf("  Using interaction coefficient: %s", interaction_coef))
  
  interaction_results <- extract_degs(
    dds = interaction_dds$dds,
    coldata = interaction_dds$coldata,
    coef_name = interaction_coef,
    group_info = list(
      female_control = "sex == 'female' & diagnosis == 'control'",
      female_AD = "sex == 'female' & diagnosis == 'AD'",
      male_control = "sex == 'male' & diagnosis == 'control'",
      male_AD = "sex == 'male' & diagnosis == 'AD'"
    ),
    padj_thresh = opt$padj,
    lfc_thresh = opt$lfc
  )
  
  save_results(interaction_results,
               prefix = sprintf("%s_sex_diagnosis_interaction", opt$brain_area),
               outdir = opt$outdir)
  
  create_volcano(interaction_results$res.df,
                 title = "Sex × Diagnosis Interaction\n(Differential AD Effect)",
                 outfile = file.path(opt$outdir, sprintf("%s_interaction_volcano.png", opt$brain_area)),
                 padj_thresh = opt$padj,
                 lfc_thresh = opt$lfc)
  
  # ============================================================================
  # SUMMARY REPORT
  # ============================================================================
  message("\n\n========================================")
  message("ANALYSIS COMPLETE")
  message("========================================\n")
  
  summary_df <- data.frame(
    Contrast = c("Males: AD vs Control", 
                 "Females: AD vs Control", 
                 "Sex × Diagnosis Interaction"),
    Total_DEGs = c(nrow(male_results$degs),
                   nrow(female_results$degs),
                   nrow(interaction_results$degs)),
    UP = c(sum(male_results$degs$diffexpressed == "UP"),
           sum(female_results$degs$diffexpressed == "UP"),
           sum(interaction_results$degs$diffexpressed == "UP")),
    DOWN = c(sum(male_results$degs$diffexpressed == "DOWN"),
             sum(female_results$degs$diffexpressed == "DOWN"),
             sum(interaction_results$degs$diffexpressed == "DOWN"))
  )
  
  print(summary_df)
  
  write.csv(summary_df, 
            file.path(opt$outdir, sprintf("%s_summary.csv", opt$brain_area)),
            row.names = FALSE)
  
  message(sprintf("\nAll results saved to: %s\n", opt$outdir))
  
  # Save session info
  writeLines(capture.output(sessionInfo()), 
             file.path(opt$outdir, "sessionInfo.txt"))
}

# ==============================================================================
# RUN
# ==============================================================================

if (!interactive()) {
  main()
}