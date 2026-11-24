# gene_annotation.R
# GC content and length of every human gene, from a FIXED Ensembl archive with
# biomaRt (as in the old 2.pre-promRNA script), so the result does not change
# with new Ensembl releases. Used by gc_normalize (EDASeq correction and NOISeq
# diagnostics of the GC and length biases).
#   gc_content  percentage_gene_gc_content
#   length      end - start + 1 (genomic span of the gene, as in the old script)
# Needs internet only when this rule runs (once).
# Run by Snakemake (rule gene_annotation).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

# biomaRt first and the tidyverse last: Bioconductor packages mask dplyr
# functions such as select() if they are loaded after it
suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
})

host <- snakemake@params[["host"]]
cat("Ensembl archive:", host, "\n")

mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = host)
cat("Mart version:", listMarts(host = host)$version[1], "\n")

annotation <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype",
                                   "chromosome_name", "start_position", "end_position",
                                   "percentage_gene_gc_content"),
                    mart = mart) %>%
  as_tibble() %>%
  transmute(gene_id = ensembl_gene_id,
            symbol = na_if(external_gene_name, ""),
            gene_biotype,
            chromosome = chromosome_name,
            gc_content = percentage_gene_gc_content,
            length = end_position - start_position + 1) %>%
  # Genes on patches/haplotypes have extra rows; keep the primary assembly one
  arrange(gene_id, nchar(chromosome)) %>%
  distinct(gene_id, .keep_all = TRUE)

cat("Genes annotated:", nrow(annotation),
    "| with GC and length:", sum(!is.na(annotation$gc_content) & !is.na(annotation$length)), "\n")

write_tsv(annotation, snakemake@output[[1]])

cat("\nDone\n")
