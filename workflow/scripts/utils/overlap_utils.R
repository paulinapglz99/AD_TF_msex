# Overlap test between a gene set and a gene label (e.g. rewired genes x DEGs),
# with a stratified permutation null. Copied from de_utils.R of the
# Region-specific_co-expression-_network_rewiring_AD repository (without the
# functions that read Sage's DE tables).
# Scripts load it with:
#   source(snakemake@input[["overlap_utils"]])

library(dplyr)
library(purrr)

# Stratum of each gene: bins (deciles by default) of one or more variables
add_strata <- function(df, vars, bins) {
  df |>
    mutate(across(all_of(vars), \(v) ntile(v, bins), .names = "bin_{.col}")) |>
    tidyr::unite("stratum", starts_with("bin_"), sep = "_", remove = TRUE)
}

# Function that returns a row order with the rows shuffled only within each
# stratum. The strata are split once, because the shuffle runs once per
# permutation (thousands of times per test).
shuffler_within <- function(stratum) {
  groups <- unname(split(seq_along(stratum), stratum))
  positions <- unlist(groups, use.names = FALSE)
  function() {
    order <- integer(length(stratum))
    order[positions] <- unlist(map(groups, \(g) g[sample.int(length(g))]), use.names = FALSE)
    order
  }
}

# Two-sided permutation p-value: how extreme the observed value is in the null
perm_p <- function(observed, null) {
  null <- null[!is.na(null)]
  upper <- (1 + sum(null >= observed)) / (1 + length(null))
  lower <- (1 + sum(null <= observed)) / (1 + length(null))
  min(1, 2 * min(upper, lower))
}

# Odds ratio of the 2 x 2 table of two logical vectors (0 or Inf when a cell is empty)
odds_ratio <- function(x, y) {
  (sum(x & y) * sum(!x & !y)) / (sum(x & !y) * sum(!x & y))
}

# Overlap test shared by analyses 4, 5 and 7: does gene set `in_set` contain
# more genes with `label` than expected among genes of the same stratum?
#   null:      `label` shuffled among genes of the same stratum, n_perm times
#   statistic: genes in both. The shuffle keeps the size of the set and the
#              number of labels, so the odds ratio only depends on this count
#   expected:  mean of the genes in both in the null
#   effect:    odds ratio of the 2 x 2 table in_set x label
# All three vectors have one element per gene of the universe.
overlap_test <- function(in_set, label, stratum, n_perm) {
  observed <- sum(in_set & label)
  shuffle <- shuffler_within(stratum)
  null <- map_int(seq_len(n_perm), \(i) sum(in_set & label[shuffle()]))
  tibble(
    n_universe = length(in_set),
    n_set = sum(in_set),
    n_label = sum(label),
    n_both = observed,
    expected_both = mean(null),
    odds_ratio = odds_ratio(in_set, label),
    p_perm = if (sum(in_set) == 0 || sum(label) == 0) NA_real_ else perm_p(observed, null),
    n_perm = n_perm
  )
}
