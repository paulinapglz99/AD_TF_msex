#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pacman::p_load(
    data.table,
    optparse,
    tools,
    igraph,
    tidyverse, 
    vroom, 
    cowplot
  )
})

#Parser  ============================================================

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character",
              help = "Directory with network files"),
  make_option(c("-o", "--out_dir"), type = "character",
              default = "jaccard_output",
              help = "Output directory"),
  make_option(c("-p", "--pattern"), type = "character",
              default = "\\.txt$",
              help = "File name regex pattern"),
  make_option(c("-t", "--type"), type = "character",
              default = "auto",
              help = "Network type: auto | edgelist | adjacency")
)

opt <- parse_args(OptionParser(option_list = option_list))

# opt$out_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/results_conn_comparisons"
# opt$input_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex"
# opt$pattern <- "\\.txt$"
# opt$type <- "edgelist"

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

setwd(opt$out_dir)

#Functions ============================================================

#Read networks
read_network <- function(path, type = c("auto", "edgelist", "adjacency")) {
  message("Reading: ", path)
  type <- match.arg(type)
  ext  <- tolower(file_ext(path))
  
  if (ext == "graphml") {
    
    g <- read_graph(path, format = "graphml")
    
  } else if (ext %in% c("csv", "tsv", "txt")) {
    
    df <- vroom(path, col_names = TRUE)
    
    if (nrow(df) == 0)
      stop("Empty file: ", path)
    
    if (type == "edgelist") {
      
      g <- graph_from_data_frame(df, directed = FALSE)
      
    } else if (type == "adjacency") {
      
      mat <- as.matrix(df)
      storage.mode(mat) <- "numeric"
      g <- graph_from_adjacency_matrix(mat,
                                       mode = "undirected",
                                       weighted = TRUE)
      
    } else {  # auto
      
      first_two_char <- all(sapply(df[, 1:2, drop = FALSE], is.character))
      
      if (first_two_char) {
        g <- graph_from_data_frame(df, directed = FALSE)
      } else {
        if (!all(sapply(df, is.numeric)))
          stop("Auto detection failed for: ", path)
        
        mat <- as.matrix(df)
        storage.mode(mat) <- "numeric"
        g <- graph_from_adjacency_matrix(mat,
                                         mode = "undirected",
                                         weighted = TRUE)
      }
    }
    
  } else {
    stop("Unsupported format: ", ext)
  }
  
  g <- as.undirected(g, mode = "collapse")
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  return(g)
}

#Jaccard

jaccard_index <- function(set1, set2) {
  length(intersect(set1, set2)) /
    length(union(set1, set2))
}

get_edge_set <- function(g) {
  edgelist <- as_edgelist(g)
  edgelist_sorted <- t(apply(edgelist, 1, sort))
  apply(edgelist_sorted, 1, paste, collapse = "__")
}

compare_networks <- function(g1, g2) {
  
  tibble(
    jaccard_nodos =
      jaccard_index(V(g1)$name, V(g2)$name),
    
    jaccard_enlaces =
      jaccard_index(get_edge_set(g1),
                    get_edge_set(g2))
  )
}

#Read data ============================================================

files <- list.files(opt$input_dir,
                    pattern = opt$pattern,
                    full.names = TRUE)

networks <- lapply(files, read_network, type = opt$type)
names(networks) <- basename(files)

#Build metadata ============================================================

meta <- tibble(
  file = names(networks),
  sex = str_extract(file, "female|male"),
  condition = str_extract(file, "AD|control"),
  region = str_extract(file, "CRB|DLPFC|HCN|PCC|TC"))

#Comparisons per region ============================================================

regions <- unique(meta$region)

results <- list()

for (r in regions) {
  
  sub <- meta %>% filter(region == r)
  
  g_f_AD  <- networks[[ sub$file[sub$sex=="female" & sub$condition=="AD"] ]]
  g_f_C   <- networks[[ sub$file[sub$sex=="female" & sub$condition=="control"] ]]
  g_m_AD  <- networks[[ sub$file[sub$sex=="male"   & sub$condition=="AD"] ]]
  g_m_C   <- networks[[ sub$file[sub$sex=="male"   & sub$condition=="control"] ]]
  
  comp1 <- compare_networks(g_f_AD, g_f_C) %>%
    mutate(region = r,
           comparison = "female_AD_vs_female_control")
  
  comp2 <- compare_networks(g_f_AD, g_m_AD) %>%
    mutate(region = r,
           comparison = "female_AD_vs_male_AD")
  
  comp3 <- compare_networks(g_m_AD, g_m_C) %>%
    mutate(region = r,
           comparison = "male_AD_vs_male_control")
  
  results[[r]] <- bind_rows(comp1, comp2, comp3)
}

final_table <- bind_rows(results) %>%
  select(region, comparison,
         jaccard_enlaces, jaccard_nodos)

#Save results ============================================================

write_csv(final_table, file.path(opt$out_dir, "jaccard_summary.csv"))

print(final_table)

#Plot results ============================================================

#Heatmap by region × comparison (Jaccard links)

plot_data <- final_table %>%
  mutate(
    comparison = str_replace_all(comparison, "_", " "),
    region = str_replace_all(region, "_", " ")) 

plot_data <- plot_data %>%
  mutate(
    comparison = case_when(
      comparison == "female AD vs female control" ~ "F-AD vs F-CTL",
      comparison == "female AD vs male AD"        ~ "F-AD vs M-AD",
      comparison == "male AD vs male control"     ~ "M-AD vs M-CTL"))

plot_data$region <- factor(
  plot_data$region,
  levels = c("CRB", "DLPFC", "HCN", "PCC", "TC")
)

plot_data$comparison <- factor(
  plot_data$comparison,
  levels = c(
    "F-AD vs F-CTL",
    "F-AD vs M-AD",
    "M-AD vs M-CTL"))

#Plot edges --- ---

ht_edges <- ggplot(plot_data,
       aes(x = comparison,
           y = region,
           fill = jaccard_enlaces)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = round(jaccard_enlaces, 3)),
            size = 4,
            fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#B2182B", "#F7F7F7", "#2166AC"),
    name = "Jaccard\n(edges)"
  )+
  labs(
    x = "Comparison",
    y = "Region"
  )+
  theme(
  legend.position = "bottom",
  legend.title = element_text(face = "bold", size = 14),
  legend.text = element_text(size = 11),
  legend.margin = margin(t = 5),
  legend.box.margin = margin(t = -5)
) +
  theme_minimal(base_size = 14) 


#Vis
ht_edges

#Save plot
ggsave(filename = "heatmap_edges_AD_TF_msex.jpeg",
        plot = ht_edges, 
        device = "jpg",
        width = 20,
        height = 20,
        units = "cm",
        dpi = 300)

#Plot nodes --- --- 

ht_nodes <- ggplot(plot_data,
                   aes(x = comparison,
                       y = region,
                       fill = jaccard_nodos)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = round(jaccard_enlaces, 3)),
            size = 4,
            fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#B2182B", "#F7F7F7", "#2166AC"),
    name = "Jaccard\n(genes)"
  )+
  labs(
    x = "Comparison",
    y = "Region"
  )+
theme(
  legend.position = "bottom",
  legend.title = element_text(face = "bold", size = 14),
  legend.text = element_text(size = 11),
  legend.margin = margin(t = 5),
  legend.box.margin = margin(t = -5)
) +
theme_minimal(base_size = 14)

#Vis
ht_nodes

#Save plot
# ggsave(filename = "heatmap_nodes_AD_TF_msex.jpeg",
#        plot = ht_nodes, 
#        device = "jpg",
#        width = 20,
#        height = 20,
#        units = "cm",
#        dpi = 300)

#Plot panel --- ---

panel <- cowplot::plot_grid(
  ht_edges,
  ht_nodes,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 2,
  align = "hv"
)

#Vis
panel

#Save panel
ggsave(filename = "panel_nodes_edges_heatmap_AD_TF_msex.jpeg",
       plot = panel, 
       device = "jpg",
       width = 40,
       height = 20,
       units = "cm",
       dpi = 300)

#THE END