# rewiring_enrichment.R
# Function of the genes rewired by AD in one sex only, in one region:
# over-representation of GO biological process and Reactome in the genes
# rewired in women only, in men only and in both (column rewired_in of
# gene_rewiring.R).
# Universe: genes that pass the degree filter in BOTH sexes (rewired_in not NA),
# i.e. the genes that could have been called rewired in either sex. BH within
# each gene set and collection. Only significant terms are kept (padj < padj).
# Run by Snakemake (rule rewiring_enrichment).

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
min_size <- snakemake@params[["min_size"]]
max_size <- snakemake@params[["max_size"]]
padj_max <- snakemake@params[["padj"]]

genes <- read_tsv(snakemake@input[["genes"]], show_col_types = FALSE,
                  col_types = cols(gene = col_character(), rewired_in = col_character())) %>%
  filter(!is.na(rewired_in))
universe <- genes$gene

# Reactome uses Entrez IDs
entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "ENSEMBL",
                                columns = "ENTREZID") %>%
  as_tibble() %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
universe_entrez <- unique(entrez$ENTREZID)

gene_sets <- c("female", "male", "both")
cat(region_name, ": universe", length(universe), "genes\n")
print(as.data.frame(dplyr::count(genes, rewired_in)))

tidy_result <- function(result, collection, set_name) {
  if (is.null(result) || nrow(result@result) == 0) return(tibble())
  as_tibble(result@result) %>%
    filter(p.adjust < padj_max) %>%
    transmute(region = region_name, rewired_in = set_name,
              collection = collection, id = ID, description = Description,
              gene_ratio = GeneRatio, background_ratio = BgRatio,
              pvalue, padj = p.adjust, genes = geneID)
}

# A gene set without enough mappable genes gives NULL instead of stopping the job
safely <- function(expr, label) tryCatch(expr, error = function(e) {
  cat(label, ":", conditionMessage(e), "\n")
  NULL
})

enrichment <- lapply(gene_sets, function(set_name) {
  set_genes <- genes %>% filter(rewired_in == set_name) %>% pull(gene)
  set_entrez <- entrez %>% filter(ENSEMBL %in% set_genes) %>% pull(ENTREZID)

  go <- safely(enrichGO(set_genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                        universe = universe, minGSSize = min_size, maxGSSize = max_size,
                        pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                        readable = TRUE), paste(set_name, "GO_BP"))
  reactome <- safely(enrichPathway(set_entrez, organism = "human", universe = universe_entrez,
                                   minGSSize = min_size, maxGSSize = max_size,
                                   pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                                   readable = TRUE), paste(set_name, "Reactome"))

  bind_rows(tidy_result(go, "GO_BP", set_name), tidy_result(reactome, "Reactome", set_name))
}) %>%
  bind_rows()

cat("\nSignificant terms per gene set and collection:\n")
print(as.data.frame(dplyr::count(enrichment, rewired_in, collection)))

# Empty table with its columns if nothing is significant
if (nrow(enrichment) == 0) {
  enrichment <- tibble(region = character(), rewired_in = character(),
                       collection = character(), id = character(), description = character(),
                       gene_ratio = character(), background_ratio = character(),
                       pvalue = numeric(), padj = numeric(), genes = character())
}
write_tsv(enrichment, snakemake@output[[1]])

cat("\nDone\n")
