# Gene filter shared by the PCA and the differential expression, so both use
# exactly the same genes. Sage's criterion: keep a gene if its CPM is above
# min_cpm in at least min_fraction of the samples of at least one group.
# Scripts load it with:
#   source(snakemake@input[["gene_filter"]])

library(dplyr)

# count_matrix: genes x samples (raw counts), column names = specimenID
# samples: table with columns specimenID, group and sex
# Returns a named list, one logical vector per group: TRUE if the gene passes
# the filter in that group
passes_by_group <- function(count_matrix, samples, min_cpm, min_fraction) {
  cpm <- t(t(count_matrix) / colSums(count_matrix)) * 1e6
  grouped <- group_by(samples, group)

  grouped %>%
    group_map(\(g, key) {
      rowMeans(cpm[, g$specimenID, drop = FALSE] > min_cpm) >= min_fraction
    }) %>%
    setNames(as.character(group_keys(grouped)$group))
}

# TRUE for the genes that pass the filter in at least one group
filter_genes <- function(count_matrix, samples, min_cpm, min_fraction) {
  passes_by_group(count_matrix, samples, min_cpm, min_fraction) %>%
    Reduce(`|`, .)
}

# TRUE for the genes that pass the filter in at least one group of EACH sex.
# Genes expressed in one sex only (chrY genes, XIST...) are FALSE: a contrast
# between sexes, such as the sex x diagnosis interaction, is not valid for them.
expressed_in_both_sexes <- function(count_matrix, samples, min_cpm, min_fraction) {
  passes <- passes_by_group(count_matrix, samples, min_cpm, min_fraction)
  groups_of_sex <- \(s) samples %>% filter(sex == s) %>% distinct(group) %>% pull(group) %>% as.character()

  Reduce(`|`, passes[groups_of_sex("female")]) & Reduce(`|`, passes[groups_of_sex("male")])
}
