# gene_annotation.R
# GC content and length of every human gene, from a FIXED Ensembl archive, so
# the result does not change with new Ensembl releases. Used by gc_normalize
# (EDASeq correction and NOISeq diagnostics of the GC and length biases).
#   gc_content  percentage_gene_gc_content
#   length      end - start + 1 (genomic span of the gene, as in the old
#               2.pre-promRNA script)
# The query goes straight to the archive's BioMart service (what biomaRt does
# internally). biomaRt itself is not used: before querying an archive it checks
# the list of archives on www.ensembl.org, which can redirect to a mirror with
# an expired certificate and fail.
# Needs internet only when this rule runs (once).
# Run by Snakemake (rule gene_annotation).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

host <- snakemake@params[["host"]]
cat("Ensembl archive:", host, "\n")

attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "chromosome_name",
                "start_position", "end_position", "percentage_gene_gc_content")

# BioMart XML query; completionStamp adds a final "[success]" line, so a
# truncated download can be detected
query <- paste0(
  '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
  '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" completionStamp="1">',
  '<Dataset name="hsapiens_gene_ensembl" interface="default">',
  paste0('<Attribute name="', attributes, '"/>', collapse = ""),
  "</Dataset></Query>"
)
url <- paste0(host, "/biomart/martservice?query=", URLencode(query, reserved = TRUE))

download <- tempfile(fileext = ".tsv")
download.file(url, download, method = "libcurl", quiet = TRUE)

lines <- readLines(download)
if (length(lines) == 0 || tail(lines, 1) != "[success]") {
  stop("Incomplete BioMart download; last line: ", tail(lines, 1))
}

annotation <- read_tsv(I(head(lines, -1)), col_names = attributes,
                       col_types = cols(.default = col_character(),
                                        start_position = col_double(),
                                        end_position = col_double(),
                                        percentage_gene_gc_content = col_double())) %>%
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
