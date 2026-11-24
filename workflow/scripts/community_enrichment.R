# community_enrichment.R
# Over-representation of GO biological process and Reactome in every Leiden
# community of one network (communities with >= min_community_size genes).
# Universe: the genes of that network, so a term is enriched compared with
# genes that could have been in the community. BH within each community and
# collection. Only significant terms are kept (padj < padj).
# Run by Snakemake (rule community_enrichment).

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
group_name <- snakemake@wildcards[["group"]]
min_community <- snakemake@params[["min_community_size"]]
min_size <- snakemake@params[["min_size"]]
max_size <- snakemake@params[["max_size"]]
padj_max <- snakemake@params[["padj"]]

nodes <- read_tsv(snakemake@input[["nodes"]], col_types = cols(gene = col_character()))
universe <- nodes$gene

# Reactome uses Entrez IDs
entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = universe, keytype = "ENSEMBL",
                                columns = "ENTREZID") %>%
  as_tibble() %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)
universe_entrez <- unique(entrez$ENTREZID)

communities <- nodes %>%
  dplyr::count(community, name = "size") %>%
  filter(size >= min_community) %>%
  pull(community)

cat(region_name, group_name, ":", length(unique(nodes$community)), "communities,",
    length(communities), "with >=", min_community, "genes\n")

tidy_result <- function(result, collection, community_id) {
  if (is.null(result) || nrow(result@result) == 0) return(tibble())
  as_tibble(result@result) %>%
    filter(p.adjust < padj_max) %>%
    transmute(region = region_name, group = group_name, community = community_id,
              collection = collection, id = ID, description = Description,
              gene_ratio = GeneRatio, background_ratio = BgRatio,
              pvalue, padj = p.adjust, genes = geneID)
}

enrichment <- lapply(communities, function(community_id) {
  genes <- nodes %>% filter(community == community_id) %>% pull(gene)
  genes_entrez <- entrez %>% filter(ENSEMBL %in% genes) %>% pull(ENTREZID)

  # A community without mappable genes gives NULL instead of stopping the job
  safely <- function(expr) tryCatch(expr, error = function(e) {
    cat("Community", community_id, ":", conditionMessage(e), "\n")
    NULL
  })

  go <- safely(enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                        universe = universe, minGSSize = min_size, maxGSSize = max_size,
                        pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                        readable = TRUE))
  reactome <- safely(enrichPathway(genes_entrez, organism = "human",
                                   universe = universe_entrez,
                                   minGSSize = min_size, maxGSSize = max_size,
                                   pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
                                   readable = TRUE))

  bind_rows(tidy_result(go, "GO_BP", community_id),
            tidy_result(reactome, "Reactome", community_id))
}) %>%
  bind_rows()

cat("\nSignificant terms per community and collection:\n")
print(as.data.frame(dplyr::count(enrichment, community, collection)))

# Empty table with its columns if nothing is significant
if (nrow(enrichment) == 0) {
  enrichment <- tibble(region = character(), group = character(), community = integer(),
                       collection = character(), id = character(), description = character(),
                       gene_ratio = character(), background_ratio = character(),
                       pvalue = numeric(), padj = numeric(), genes = character())
}
write_tsv(enrichment, snakemake@output[[1]])

cat("\nDone\n")
