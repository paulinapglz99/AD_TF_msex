# gsea.R
# Gene set enrichment analysis (GSEA) of one DESeq2 contrast of one region.
#   - Genes: all genes tested in the contrast (column 'tested'), ranked by the
#     Wald statistic 'stat' (log2FC / its standard error), not by log2FC alone,
#     so noisy genes with large but uncertain fold changes do not dominate.
#   - Collections: GO biological process (clusterProfiler::gseGO, Ensembl IDs)
#     and Reactome (ReactomePA::gsePathway, Entrez IDs).
#   - Every gene set is kept, significant or not, so NES can be compared
#     between contrasts (e.g. AD effect in women vs in men).
# Same parameters for every region and contrast. Everything runs offline.
# Run by Snakemake (rule gsea).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

# Bioconductor first and the tidyverse last: Bioconductor packages mask dplyr
# functions such as count() and select() if they are loaded after it
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(readr)
})

region_name <- snakemake@wildcards[["region"]]
contrast_name <- snakemake@wildcards[["contrast"]]
min_size <- snakemake@params[["min_size"]]
max_size <- snakemake@params[["max_size"]]
seed <- snakemake@params[["seed"]]

# ---- Ranked gene lists ------------------------------------------------------

ranks <- read_tsv(snakemake@input[["table"]],
                  col_types = cols(chromosome = col_character())) %>%
  filter(tested, !is.na(stat)) %>%
  dplyr::select(gene_id, stat)

cat("Region:", region_name, "| contrast:", contrast_name,
    "\nRanked genes:", nrow(ranks), "\n")

ranked_ensembl <- sort(setNames(ranks$stat, ranks$gene_id), decreasing = TRUE)

# Reactome uses Entrez IDs. If several Ensembl genes map to the same Entrez ID,
# keep the one with the largest |stat|.
entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = ranks$gene_id,
                                keytype = "ENSEMBL", columns = "ENTREZID") %>%
  as_tibble() %>%
  filter(!is.na(ENTREZID))

ranks_entrez <- ranks %>%
  inner_join(entrez, by = c("gene_id" = "ENSEMBL"), relationship = "many-to-many") %>%
  group_by(ENTREZID) %>%
  slice_max(abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()

ranked_entrez <- sort(setNames(ranks_entrez$stat, ranks_entrez$ENTREZID), decreasing = TRUE)
cat("Ranked genes with an Entrez ID (Reactome):", length(ranked_entrez), "\n")

# ---- GSEA -------------------------------------------------------------------

set.seed(seed)
go <- gseGO(ranked_ensembl,
            ont = "BP",
            OrgDb = org.Hs.eg.db,
            keyType = "ENSEMBL",
            minGSSize = min_size,
            maxGSSize = max_size,
            pvalueCutoff = 1, # keep every gene set
            pAdjustMethod = "BH",
            eps = 0, # exact p-values, also for very small ones
            seed = TRUE,
            verbose = FALSE)

set.seed(seed)
reactome <- gsePathway(ranked_entrez,
                       organism = "human",
                       minGSSize = min_size,
                       maxGSSize = max_size,
                       pvalueCutoff = 1,
                       pAdjustMethod = "BH",
                       eps = 0,
                       seed = TRUE,
                       verbose = FALSE)

# ---- One table for both collections ------------------------------------------

# Leading edge genes as symbols instead of IDs
tidy_gsea <- function(result, collection, key_type) {
  if (is.null(result) || nrow(result@result) == 0) return(tibble())
  setReadable(result, OrgDb = org.Hs.eg.db, keyType = key_type)@result %>%
    as_tibble() %>%
    transmute(region = region_name,
              contrast = contrast_name,
              collection = collection,
              id = ID,
              description = Description,
              size = setSize,
              NES,
              pvalue,
              padj = p.adjust,
              leading_edge = core_enrichment)
}

table <- bind_rows(tidy_gsea(go, "GO_BP", "ENSEMBL"),
                   tidy_gsea(reactome, "Reactome", "ENTREZID")) %>%
  arrange(collection, pvalue)

cat("\nGene sets tested and significant (padj < 0.05):\n")
print(as.data.frame(table %>%
                      group_by(collection) %>%
                      summarise(tested = n(), significant = sum(padj < 0.05))))

write_tsv(table, snakemake@output[[1]])

cat("\nDone\n")
