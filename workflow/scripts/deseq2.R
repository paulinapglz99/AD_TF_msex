# deseq2.R
# Differential expression of one brain region with DESeq2. One model per region,
#   ~ <covariates> + group,   group = diagnosis_sex (CT_female, CT_male, AD_female, AD_male)
# and every contrast of config.yaml is taken from that same fit, so they share
# samples, genes, covariates and dispersions and can be compared.
# Normalization: the GC/length factors of gc_normalize (EDASeq), not DESeq2's
# size factors.
# Each contrast gives a weight per group (numerator +1, denominator -1), e.g.
#   AD_vs_CT_female            AD_female - CT_female
#   sex_diagnosis_interaction  (AD_female - CT_female) - (AD_male - CT_male)
# Contrasts that compare women with men (the interaction) are tested only on
# genes expressed in both sexes: column 'tested', padj recomputed among them.
# Writes one table per contrast (all genes of the model) and a summary of DEGs.
# Run by Snakemake (rule deseq2).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

# Bioconductor first and the tidyverse last: Bioconductor packages mask dplyr
# functions such as count() and select() if they are loaded after it
suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(readr)
  source(snakemake@input[["gene_filter"]])
  source(snakemake@input[["expression_utils"]])
})

region_name <- snakemake@wildcards[["region"]]
covariates <- unlist(snakemake@params[["covariates"]])
contrasts <- snakemake@params[["contrasts"]]
padj_max <- snakemake@params[["padj"]]
log2fc_min <- snakemake@params[["log2fc"]]
threads <- snakemake@threads

# The first level is the reference of the model
group_levels <- c("CT_female", "CT_male", "AD_female", "AD_male")

# ---- Samples and covariates -------------------------------------------------

samples <- read_tsv(snakemake@input[["samples"]],
                    col_types = cols(specimenID = col_character(),
                                     individualID = col_character(),
                                     batch = col_character(),
                                     model_batch = col_character(),
                                     apoe4_allele = col_character())) %>%
  mutate(group = factor(group, levels = group_levels),
         # Categorical covariates as factors, numeric ones centered and scaled
         across(all_of(covariates) & where(is.character), factor),
         across(all_of(covariates) & where(is.numeric), \(x) as.numeric(scale(x))))

incomplete <- samples %>% filter(if_any(all_of(covariates), is.na))
if (nrow(incomplete) > 0) {
  cat("Samples removed for missing covariates:", paste(incomplete$specimenID, collapse = ", "), "\n")
  samples <- anti_join(samples, incomplete, by = "specimenID")
}

design <- reformulate(c(covariates, "group"))
cat("Region:", region_name, "\nDesign:", deparse(design), "\n\nSamples per group:\n")
print(as.data.frame(dplyr::count(samples, group)))

# ---- Genes ------------------------------------------------------------------

count_matrix <- readRDS(snakemake@input[["counts"]]) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()
count_matrix <- round(count_matrix[, samples$specimenID])

min_cpm <- snakemake@params[["min_cpm"]]
min_fraction <- snakemake@params[["min_fraction"]]

# GC/length normalization factors (gc_normalize): only genes with GC and
# length annotation have them, so the model uses those genes only
factors <- readRDS(snakemake@input[["factors"]])

passes <- filter_genes(count_matrix, samples, min_cpm, min_fraction)
keep <- passes & rownames(count_matrix) %in% rownames(factors)
cat("\nGenes passing the filter:", sum(passes), "of", length(passes),
    "| with GC/length correction (kept):", sum(keep), "\n")

# Genes expressed in at least one group of each sex. In contrasts that compare
# women with men (the interaction) only these genes are tested: for genes
# expressed in one sex only (chrY, XIST...) such a contrast is not valid.
# The model itself is fitted with all the genes that pass the filter.
both_sexes <- expressed_in_both_sexes(count_matrix, samples, min_cpm, min_fraction)[keep]

# ---- Model ------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(count_matrix[keep, ],
                              colData = tibble::column_to_rownames(samples, "specimenID"),
                              design = design)
# GC/length correction: gene x sample normalization factors replace DESeq2's
# size factors (DESeq() then skips estimateSizeFactors)
normalizationFactors(dds) <- normalization_factors_for(factors, rownames(dds), colnames(dds))
dds <- DESeq(dds, parallel = threads > 1, BPPARAM = MulticoreParam(threads))
saveRDS(dds, snakemake@output[["dds"]])

cat("\nCoefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

# Weights per group -> numeric contrast over the model coefficients.
# The weights add up to 0, so the intercept and the reference group drop out.
contrast_vector <- function(weights) {
  weights <- unlist(weights)
  stopifnot(all(names(weights) %in% group_levels), abs(sum(weights)) < 1e-8)
  coefficients <- setNames(numeric(length(resultsNames(dds))), resultsNames(dds))
  for (g in setdiff(names(weights), group_levels[1])) {
    coefficients[[paste0("group_", g, "_vs_", group_levels[1])]] <- weights[[g]]
  }
  coefficients
}

# TRUE if the contrast has groups of both sexes (e.g. the interaction)
compares_sexes <- function(weights) {
  sexes <- unique(sub("^[^_]+_", "", names(unlist(weights))))
  length(sexes) > 1
}

# Benjamini-Hochberg adjustment among the eligible genes only; NA for the rest
bh_within <- function(pvalue, eligible) {
  padj <- rep(NA_real_, length(pvalue))
  padj[eligible] <- p.adjust(pvalue[eligible], method = "BH")
  padj
}

# ---- Annotation and mean expression per group -------------------------------

annotation <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(dds), keytype = "ENSEMBL",
                                    columns = c("SYMBOL", "GENETYPE", "CHR")) %>%
  as_tibble() %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  rename(gene_id = ENSEMBL, symbol = SYMBOL, gene_type = GENETYPE, chromosome = CHR)

one_sex_genes <- annotation %>% filter(gene_id %in% names(both_sexes)[!both_sexes])
cat("\nGenes not expressed in both sexes (not tested in contrasts between sexes):",
    sum(!both_sexes), "\n  on chrY:", sum(one_sex_genes$chromosome == "Y", na.rm = TRUE),
    "\n  symbols:", paste(sort(na.omit(one_sex_genes$symbol)), collapse = ", "), "\n")

group_means <- counts(dds, normalized = TRUE) %>%
  as_tibble(rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "specimenID", values_to = "expression") %>%
  left_join(dplyr::select(samples, specimenID, group), by = "specimenID") %>%
  group_by(gene_id, group) %>%
  summarise(mean = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean, names_prefix = "mean_")

# ---- One table per contrast -------------------------------------------------

output_of <- function(name) {
  grep(paste0("/", name, "\\.tsv$"), snakemake@output[["tables"]], value = TRUE)
}

summaries <- lapply(names(contrasts), function(name) {
  result <- results(dds, contrast = contrast_vector(contrasts[[name]]), alpha = padj_max)
  between_sexes <- compares_sexes(contrasts[[name]])

  table <- as.data.frame(result) %>%
    as_tibble(rownames = "gene_id") %>%
    mutate(expressed_both_sexes = unname(both_sexes[gene_id]),
           # Genes this contrast is tested on: all of them, or only those
           # expressed in both sexes if the contrast compares women with men.
           # padj is then recomputed among the tested genes; genes that DESeq2's
           # independent filtering left with padj = NA stay NA.
           tested = if (between_sexes) expressed_both_sexes else TRUE,
           padj = if (between_sexes) bh_within(pvalue, tested & !is.na(padj)) else padj,
           region = region_name,
           contrast = name,
           direction = case_when(
             !tested ~ "not_tested",
             padj < padj_max & log2FoldChange >= log2fc_min ~ "up",
             padj < padj_max & log2FoldChange <= -log2fc_min ~ "down",
             .default = "ns"
           )) %>%
    left_join(annotation, by = "gene_id") %>%
    left_join(group_means, by = "gene_id") %>%
    relocate(region, contrast, gene_id, symbol, chromosome, gene_type) %>%
    arrange(pvalue)

  write_tsv(table, output_of(name))
  dplyr::count(table, region, contrast, direction)
})

summary <- bind_rows(summaries) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)

cat("\nDEGs (padj <", padj_max, ", |log2FC| >=", log2fc_min, "):\n")
print(as.data.frame(summary))
write_tsv(summary, snakemake@output[["summary"]])

cat("\nDone\n")
