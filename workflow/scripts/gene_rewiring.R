# gene_rewiring.R
# Rewiring of every gene by AD, in women and in men, in one region; and whether
# the rewired genes are the DEGs. Adapted from gene_rewiring.R and
# rewiring_vs_de.R of the Region-specific_co-expression-_network_rewiring_AD
# repository, with one AD and one CT network per sex.
#   1. neighbor_change = 1 - Jaccard of the gene's neighbors in AD and in CT
#      (universe: every gene of the region's networks; degree 0 if absent)
#   2. Degree filter: genes with >= min_degree neighbors in AD or in CT
#   3. Adjusted rewiring = residual of neighbor_change ~ log(mean degree) +
#      log(mean expression): genes with few neighbors change more by chance,
#      and expression drives both degree and DE power
#   4. Rewired genes = top fraction by adjusted rewiring (main: top_rewired)
#   5. Rewired in women only / men only / both / neither
#   6. Overlap test rewired x DEG of the same sex (overlap_utils.R): DEG labels
#      shuffled within strata of expression x mean degree; odds ratio and
#      permutation p, for the main fraction and the sensitivity ones.
# WARNING: the AD and CT networks of each sex have different numbers of samples.
# Run by Snakemake (rule gene_rewiring).

source(snakemake@input[["network_utils"]])
source(snakemake@input[["overlap_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region_name <- snakemake@wildcards[["region"]]
groups <- unlist(snakemake@params[["groups"]])
min_degree <- snakemake@params[["min_degree"]]
top_rewired <- snakemake@params[["top_rewired"]]
fractions <- unique(c(top_rewired, unlist(snakemake@params[["sensitivity"]])))
n_perm <- snakemake@params[["n_perm"]]
bins <- snakemake@params[["strata"]]
set.seed(snakemake@params[["seed"]])

# Inputs are in the order of `groups` (expand in the Snakefile)
networks <- set_names(map(unlist(snakemake@input[["networks"]]), read_network), groups)
universe <- readLines(snakemake@input[["regulators"]]) # every gene of the networks

read_de <- function(file) {
  vroom(file, col_types = cols(gene_id = "c", symbol = "c", chromosome = "c")) |>
    transmute(gene = gene_id, symbol, log2FoldChange, padj, baseMean,
              is_de = direction %in% c("up", "down"))
}
de <- list(female = read_de(snakemake@input[["de_female"]]),
           male = read_de(snakemake@input[["de_male"]]))
interaction <- read_de(snakemake@input[["de_interaction"]]) |>
  transmute(gene, interaction_log2fc = log2FoldChange, interaction_padj = padj,
            interaction_de = is_de)

# Each edge seen from both ends: one row per gene and neighbor
neighbor_table <- function(edges) {
  bind_rows(transmute(edges, gene = a, neighbor = b), transmute(edges, gene = b, neighbor = a))
}

# ---- 1-4. Rewiring of every gene in one sex ---------------------------------

rewiring_in <- function(sex) {
  ad <- neighbor_table(networks[[paste0("AD_", sex)]])
  ct <- neighbor_table(networks[[paste0("CT_", sex)]])
  shared <- inner_join(ad, ct, by = c("gene", "neighbor")) |> dplyr::count(gene, name = "shared")

  tibble(gene = universe) |>
    left_join(dplyr::count(ad, gene, name = "degree_ad"), by = "gene") |>
    left_join(dplyr::count(ct, gene, name = "degree_ct"), by = "gene") |>
    left_join(shared, by = "gene") |>
    replace_na(list(degree_ad = 0L, degree_ct = 0L, shared = 0L)) |>
    mutate(union = degree_ad + degree_ct - shared,
           neighbor_change = if_else(union > 0, 1 - shared / union, NA_real_),
           mean_degree = (degree_ad + degree_ct) / 2) |>
    left_join(de[[sex]] |> dplyr::select(gene, baseMean, is_de), by = "gene") |>
    filter(degree_ad >= min_degree | degree_ct >= min_degree, !is.na(baseMean)) |>
    mutate(log_expression = log(baseMean + 1)) |>
    (\(df) mutate(df, rewiring_adj = residuals(
      lm(neighbor_change ~ log(mean_degree) + log_expression, data = df))))() |>
    mutate(rewired = rewiring_adj >= quantile(rewiring_adj, 1 - top_rewired),
           sex = sex)
}

by_sex <- map(c("female", "male"), rewiring_in) |> set_names(c("female", "male"))

# ---- 5. One row per gene: women, men and the comparison ---------------------

genes <- bind_rows(by_sex) |>
  dplyr::select(gene, sex, degree_ad, degree_ct, neighbor_change, rewiring_adj, rewired, is_de) |>
  pivot_wider(names_from = sex,
              values_from = c(degree_ad, degree_ct, neighbor_change, rewiring_adj, rewired, is_de)) |>
  mutate(rewired_in = case_when(
           rewired_female & rewired_male ~ "both",
           rewired_female ~ "female",
           rewired_male ~ "male",
           !is.na(rewired_female) & !is.na(rewired_male) ~ "none",
           .default = NA_character_ # passes the degree filter in one sex only
         ),
         # > 0: the gene is rewired more by AD in women than in men
         rewiring_difference = rewiring_adj_female - rewiring_adj_male) |>
  left_join(de$female |> dplyr::select(gene, symbol), by = "gene") |>
  left_join(interaction, by = "gene") |>
  mutate(region = region_name, .before = 1) |>
  relocate(symbol, .after = gene) |>
  arrange(desc(abs(rewiring_difference)))

# ---- 6. Are the rewired genes the DEGs of the same sex? ---------------------

tests <- expand_grid(sex = c("female", "male"), top_fraction = fractions) |>
  pmap(\(sex, top_fraction) {
    df <- add_strata(by_sex[[sex]], c("log_expression", "mean_degree"), bins)
    in_set <- df$rewiring_adj >= quantile(df$rewiring_adj, 1 - top_fraction)
    overlap_test(in_set, df$is_de, df$stratum, n_perm) |>
      mutate(region = region_name, sex = sex, top_fraction = top_fraction,
             main = top_fraction == top_rewired, .before = 1)
  }) |>
  list_rbind()

vroom_write(genes, snakemake@output[["genes"]])
vroom_write(tests, snakemake@output[["tests"]])

message(region_name, ": genes passing the degree filter - women ", nrow(by_sex$female),
        ", men ", nrow(by_sex$male))
print(as.data.frame(dplyr::count(genes, rewired_in)))
print(as.data.frame(tests))
