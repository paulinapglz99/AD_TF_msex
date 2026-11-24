# network_metrics.R
# Topology of one ARACNe-AP network (one region x group).
# Undirected edges: with all genes as regulators ARACNe can report A->B and
# B->A, and both count as one edge (network_utils.R).
#   metrics.tsv  one row: nodes, edges, density, components, giant component,
#                clustering (global and mean local), degree assortativity,
#                degree (mean, median, sd, max), max k-core, Leiden modularity,
#                number of communities and size of the largest one
#   nodes.tsv    one row per gene with edges: degree, PageRank, coreness and
#                Leiden community (1 = largest)
# Path lengths, diameter and percolation are left out: with ~25,000 nodes and
# hundreds of thousands of edges they take very long and add little.
# Run by Snakemake (rule network_metrics).

source(snakemake@input[["network_utils"]])
library(tidyr)

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

region_name <- snakemake@wildcards[["region"]]
group_name <- snakemake@wildcards[["group"]]

edges <- read_network(snakemake@input[["network"]])
graph <- graph_from_data_frame(edges, directed = FALSE)

# Leiden communities, renumbered by size (1 = largest)
communities <- leiden_communities(edges,
                                  resolution = snakemake@params[["resolution"]],
                                  seed = snakemake@params[["seed"]])
community_sizes <- dplyr::count(communities, community, sort = TRUE) |>
  mutate(rank = row_number())
communities <- communities |>
  left_join(community_sizes |> dplyr::select(community, rank), by = "community") |>
  transmute(gene, community = rank)

membership <- communities$community[match(V(graph)$name, communities$gene)]
component_sizes <- components(graph)$csize
degree_values <- degree(graph)

metrics <- tibble(
  region = region_name,
  group = group_name,
  n_nodes = vcount(graph),
  n_edges = ecount(graph),
  density = edge_density(graph),
  n_components = length(component_sizes),
  giant_component_fraction = max(component_sizes) / vcount(graph),
  clustering_global = transitivity(graph, type = "global"),
  clustering_local_mean = transitivity(graph, type = "localaverage"),
  assortativity = assortativity_degree(graph, directed = FALSE),
  degree_mean = mean(degree_values),
  degree_median = median(degree_values),
  degree_sd = sd(degree_values),
  degree_max = max(degree_values),
  kcore_max = max(coreness(graph)),
  modularity = modularity(graph, membership),
  n_communities = max(communities$community),
  largest_community = max(community_sizes$n)
)

nodes <- tibble(gene = V(graph)$name,
                degree = degree_values,
                pagerank = page_rank(graph)$vector,
                coreness = coreness(graph)) |>
  left_join(communities, by = "gene") |>
  mutate(region = region_name, group = group_name, .before = 1)

vroom_write(metrics, snakemake@output[["metrics"]])
vroom_write(nodes, snakemake@output[["nodes"]])

message(region_name, " ", group_name, ": ", metrics$n_nodes, " nodes, ", metrics$n_edges,
        " edges, ", metrics$n_communities, " communities, modularity ",
        round(metrics$modularity, 3))
