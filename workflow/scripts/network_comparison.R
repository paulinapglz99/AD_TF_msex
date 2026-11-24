# network_comparison.R
# Comparison of the 4 networks of one region (diagnosis x sex), every pair:
#   edge_jaccard     edges in both / edges in at least one
#   node_jaccard     genes with edges in both / in at least one
#   degree_spearman  Spearman rho of the degree, genes with edges in both
#   community_nmi    NMI of the Leiden communities, genes with edges in both
# The pairs that answer the question:
#   AD_female vs CT_female   rewiring by AD in women
#   AD_male   vs CT_male     rewiring by AD in men
#   CT_female vs CT_male     normal difference between the sexes
#   AD_female vs AD_male     difference between the sexes with AD
# Convergence/divergence index (per measure) = similarity of the sexes with AD
# minus their similarity in controls: > 0 convergence (with AD the networks of
# women and men are more alike), < 0 divergence.
# WARNING: each network has its own number of samples (columns n_x, n_y). With
# more samples ARACNe finds more edges, so these measures mix biology with n.
# Run by Snakemake (rule network_comparison).

source(snakemake@input[["network_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region_name <- snakemake@wildcards[["region"]]
groups <- unlist(snakemake@params[["groups"]])

# Inputs are in the order of `groups` (expand in the Snakefile)
networks <- set_names(map(unlist(snakemake@input[["networks"]]), read_network), groups)
communities <- set_names(
  map(unlist(snakemake@input[["nodes"]]),
      \(f) vroom(f, col_types = cols(gene = "c")) |> dplyr::select(gene, community)),
  groups
)
n_samples <- vroom(snakemake@input[["samples"]], col_types = cols(.default = "c")) |>
  dplyr::count(group)
n_of <- \(g) n_samples$n[n_samples$group == g]

node_jaccard <- function(x, y) {
  gx <- unique(c(x$a, x$b))
  gy <- unique(c(y$a, y$b))
  length(intersect(gx, gy)) / length(union(gx, gy))
}

pairs <- combn(groups, 2, simplify = FALSE) |>
  map(\(p) {
    x <- p[1]
    y <- p[2]
    tibble(region = region_name, x = x, y = y,
           n_x = n_of(x), n_y = n_of(y),
           edges_x = nrow(networks[[x]]), edges_y = nrow(networks[[y]]),
           edge_jaccard = jaccard(networks[[x]], networks[[y]]),
           node_jaccard = node_jaccard(networks[[x]], networks[[y]]),
           degree_spearman = degree_spearman(networks[[x]], networks[[y]]),
           community_nmi = communities_nmi(communities[[x]], communities[[y]]))
  }) |>
  list_rbind()

# Similarity of the two sexes within one diagnosis
between_sexes <- function(dx) {
  pairs |>
    filter((x == paste0(dx, "_female") & y == paste0(dx, "_male")) |
             (x == paste0(dx, "_male") & y == paste0(dx, "_female")))
}

measures <- c("edge_jaccard", "node_jaccard", "degree_spearman", "community_nmi")
convergence <- tibble(measure = measures) |>
  mutate(region = region_name,
         between_sexes_AD = map_dbl(measure, \(m) between_sexes("AD")[[m]]),
         between_sexes_CT = map_dbl(measure, \(m) between_sexes("CT")[[m]]),
         index = between_sexes_AD - between_sexes_CT,
         direction = if_else(index > 0, "convergence", "divergence"),
         .before = 1)

vroom_write(pairs, snakemake@output[["pairs"]])
vroom_write(convergence, snakemake@output[["convergence"]])

print(as.data.frame(pairs))
print(as.data.frame(convergence))
