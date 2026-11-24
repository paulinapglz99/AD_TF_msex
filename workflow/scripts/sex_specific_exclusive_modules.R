#!/usr/bin/env Rscript

#Libraries
pacman::p_load("optparse", 
               "data.table", 
               "dplyr", 
               "stringr", 
               "purrr", 
               "ggplot2", 
               "cowplot", 
               "ggsci")

#Parser

option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character",
              help = "Directory with *_nodes_summary.csv files"),
  make_option(c("-o", "--out_dir"), type = "character",
              default = "sex_exclusivity_output",
              help = "Output directory"),
  make_option(c("-p", "--pattern"), type = "character",
              default = "_nodes_summary\\.csv$",
              help = "Regex pattern for files"),
  make_option(c("-u", "--upper"), type = "double",
              default = 0.7,
              help = "Upper Jaccard threshold"),
  make_option(c("-l", "--lower"), type = "double",
              default = 0.4,
              help = "Lower Jaccard threshold"),
  make_option(c("-m", "--min_genes"), type = "integer",
              default = 3,
              help = "Minimum genes per module")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$input_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/topos_results"
opt$out_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/exclusive_modules_results"
opt$upper <- 0.7
opt$lower <- 0.3
opt$min_genes <- 3
setwd(opt$out_dir)

if (is.null(opt$input_dir)) {
  stop("ERROR: --input_dir must be provided")
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Input directory:", opt$input_dir, "\n")
cat("Output directory:", opt$out_dir, "\n")
cat("Upper threshold:", opt$upper, "\n")
cat("Lower threshold:", opt$lower, "\n")
cat("Min genes:", opt$min_genes, "\n")

#Functions --- ---

#Jaccard 
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#Filename parser to extract metadata

parse_filename <- function(filepath) {
  
  fname <- basename(filepath)
  
  sex <- case_when(
    str_detect(fname, "female") ~ "F",
    str_detect(fname, "male")   ~ "M",
    TRUE ~ NA_character_
  )
  
  phenotype <- case_when(
    str_detect(fname, "_AD_") ~ "AD",
    str_detect(fname, "_control_") ~ "CTL",
    TRUE ~ NA_character_
  )
  
  region <- fname %>%
    str_remove("_nodes_summary\\.csv") %>%
    str_extract("counts_[A-Za-z0-9]+") %>%
    str_remove("counts_")
  
  tibble(
    file = filepath,
    region = region,
    sex = sex,
    phenotype = phenotype
  )
}

#Get files

files <- list.files(
  opt$input_dir,
  pattern = opt$pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files found with pattern: ", opt$pattern)
}

meta <- map_dfr(files, parse_filename)

#Get module data

load_modules <- function(file, region, sex, phenotype) {
  
  df <- fread(file)
  df <- df[, .(node, membership)]
  
  keep_modules <- df[, .N, by = membership][N >= opt$min_genes]$membership
  df <- df[membership %in% keep_modules]
  
  df[, module_id := paste(region, sex, phenotype, membership, sep = "_")]
  
  split(df$node, df$module_id)
}

modules_list <- list()

for (i in seq_len(nrow(meta))) {
  
  mods <- load_modules(
    meta$file[i],
    meta$region[i],
    meta$sex[i],
    meta$phenotype[i]
  )
  
  modules_list <- c(modules_list, mods)
}

#Pairwise Jaccard

regions <- unique(meta$region)
jaccard_results <- list()

for (reg in regions) {
  
  cat("Processing region:", reg, "\n")
  
  region_modules <- modules_list[str_detect(names(modules_list), paste0("^", reg, "_"))]
  mod_names <- names(region_modules)
  
  if (length(mod_names) < 2) next
  
  for (i in 1:(length(mod_names)-1)) {
    for (j in (i+1):length(mod_names)) {
      
      m1 <- mod_names[i]
      m2 <- mod_names[j]
      
      g1 <- region_modules[[m1]]
      g2 <- region_modules[[m2]]
      
      jval <- jaccard_index(g1, g2)
      
      classif <- case_when(
        jval >= opt$upper ~ "Similar",
        jval <= opt$lower ~ "Dissimilar",
        TRUE ~ "Intermediate"
      )
      
      jaccard_results[[length(jaccard_results)+1]] <- tibble(
        Region = reg,
        Module_1 = m1,
        Module_2 = m2,
        Jaccard = jval,
        Classification = classif
      )
    }
  }
}

jaccard_df <- bind_rows(jaccard_results)

#Save data

vroom::vroom_write(jaccard_df, file.path(opt$out_dir, "jaccard_all_pairs.csv"))

#Structural exclusivity

exclusive_list <- list()

for (reg in regions) {
  
  region_modules <- modules_list[str_detect(names(modules_list), paste0("^", reg, "_"))]
  region_jacc <- jaccard_df %>% filter(Region == reg)
  
  for (mod in names(region_modules)) {
    
    comps <- region_jacc %>%
      filter(Module_1 == mod | Module_2 == mod)
    
    has_similar <- any(comps$Jaccard >= opt$upper)
    
    if (!has_similar) {
      
      exclusive_list[[length(exclusive_list)+1]] <- tibble(
        Region = reg,
        Module = mod,
        N_genes = length(region_modules[[mod]])
      )
    }
  }
}

exclusive_df <- bind_rows(exclusive_list)

vroom::vroom_write(exclusive_df, file.path(opt$out_dir, "exclusive_modules.csv"))

#Summary 

summary_df <- exclusive_df %>%
  mutate(
    parts = str_split(Module, "_"),
    Region = map_chr(parts, 1),
    Sex = map_chr(parts, 2),
    Phenotype = map_chr(parts, 3)
  ) %>%
  select(-parts) %>%
  group_by(Region, Sex, Phenotype) %>%
  summarise(N_exclusive = n(), .groups = "drop")

#Save data

#vroom::vroom_write(summary_df, file.path(opt$out_dir, "exclusive_summary.csv"))

#BUILD COMPLETE MODULE CLASSIFICATION ============================================================

module_class_list <- list()

for (reg in regions) {
  
  region_modules <- modules_list[str_detect(names(modules_list),
                                            paste0("^", reg, "_"))]
  
  region_jacc <- jaccard_df %>% filter(Region == reg)
  
  for (mod in names(region_modules)) {
    
    comps <- region_jacc %>%
      filter(Module_1 == mod | Module_2 == mod)
    
    if (nrow(comps) == 0) {
      classif <- "Exclusive"
    } else {
      maxJ <- max(comps$Jaccard)
      
      if (maxJ >= opt$upper) {
        classif <- "Similar"
      } else if (maxJ <= opt$lower) {
        classif <- "Exclusive"
      } else {
        classif <- "Intermediate"
      }
    }
    
    module_class_list[[length(module_class_list)+1]] <- data.frame(
      Region = reg,
      Module = mod,
      Classification = classif,
      stringsAsFactors = FALSE
    )
  }
}

module_class_df <- bind_rows(module_class_list)

#Clean metadata from Module
module_class_df <- module_class_df %>%
  mutate(
    parts = str_split(Module, "_"),
    Region = map_chr(parts, 1),
    Sex = map_chr(parts, 2),
    Phenotype = map_chr(parts, 3)
  ) %>%
  select(-parts)

#Module size dataframe

module_size_df <- tibble(
  Module = names(modules_list),
  N_genes = map_int(modules_list, length)
) %>%
  mutate(
    parts = str_split(Module, "_"),
    Region = map_chr(parts, 1),
    Sex = map_chr(parts, 2),
    Phenotype = map_chr(parts, 3)
  ) %>%
  select(-parts)

#Check normalidad 
#Test Shapiro-Wilk

#The Shapiro–Wilk test evaluates:
#H0 (null hypothesis): the data follow a normal distribution
#H1: the data do NOT follow a normal distributionl

#Within each sex, does the distribution of module size follow a normal distribution?
by(module_size_df$N_genes,
   module_size_df$Sex,
   shapiro.test)

# module_size_df$Sex: F
# 
# Shapiro-Wilk normality test
# 
# data:  dd[x, ]
# W = 0.91656, p-value = 2.168e-06
# 
# --------------------------------------------------------------------------------------------------------------- 
#   module_size_df$Sex: M
# 
# Shapiro-Wilk normality test
# 
# data:  dd[x, ]
# W = 0.95153, p-value = 0.001048

#Similar shapes, different centers

wilcox.test(N_genes ~ Sex,
            data = module_size_df)

# Wilcoxon rank sum test with continuity correction
# 
# data:  N_genes by Sex
# W = 5132.5, p-value = 0.1453
# alternative hypothesis: true location shift is not equal to 0

#No statistical evidence that the size of the modules differs between men and women.

#Do the men’s and women’s sections vary in size within each region?
module_size_region_test <- module_size_df %>%
  group_by(Region) %>%
  summarise(
    p_value = wilcox.test(N_genes ~ Sex)$p.value
  )

module_size_region_test$p_adj <- p.adjust(
  module_size_region_test$p_value,
  method = "BH"
)
# 
# Region p_value p_adj
# <chr>    <dbl> <dbl>
#   1 CRB      0.522 0.653
# 2 DLPFC    0.427 0.653
# 3 HCN      0.232 0.653
# 4 PCC      0.933 0.933
# 5 TC       0.392 0.653

module_size_region_test

#Plot data ============================================================

#Plot LOCAL proportions (Region - Sex - Phenotype)

proportions.df <- module_class_df %>%
  group_by(Region, Sex, Phenotype, Classification) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Region, Sex, Phenotype) %>%
  mutate(Percentage = N / sum(N))


proportions.p <- ggplot(proportions.df,
                        aes(x = Region,
                            y = Percentage,
                            fill = Classification)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Sex ~ Phenotype) +
  scale_fill_aaas(name = "Module classification") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_blank(),
    legend.position = "none"
  )

#Save plot

ggsave("proportions_barplot.pdf",
       plot = proportions.p,
       width = 10,
       height = 10,
       dpi = 300)

ggsave("proportions_barplot.jpg",
       plot = proportions.p,
       width = 10,
       height = 10,
       dpi = 300)

#Plot GLOBAL proportions (Region - Sex)

global.df <- module_class_df %>%
  group_by(Region, Sex, Classification) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Region, Sex) %>%
  mutate(Percentage = N / sum(N))

#Set factors

global.df$Classification <- factor(global.df$Classification,
                                   levels = c("Exclusive",
                                              "Intermediate",
                                              "Similar"))

global.df$Sex <- factor(global.df$Sex,
                        levels = c("F", "M"))

global.p <- ggplot(global.df,
                   aes(x = Region,
                       y = Percentage,
                       fill = Classification)) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8) +
  facet_wrap(~ Sex, nrow = 1) +
  scale_fill_aaas(name = "Classification  ") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(y = element_blank()) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    #axis.title = "Region",
    legend.title = element_text(size = 12),
    legend.position = "bottom")

#Save plots

ggsave("proportions_global_barplot.pdf",
       plot = global.p,
       width = 10,
       height = 10,
       dpi = 300)

ggsave("proportions_global_barplot.jpg",
       plot = global.p,
       width = 10,
       height = 10,
       dpi = 300)

#Build panel
panel_combined <- plot_grid(
  proportions.p,
  global.p,
  labels = c("A", "B"),
  ncol = 1,
  rel_heights = c(1.2, 0.8))

#Save panel

ggsave("Figure_Modular_Reorganization_Sex_AD_PANEL.jpg",
       panel_combined,
       width = 10,
       height = 14,
       dpi = 300)

ggsave("Figure_Modular_Reorganization_Sex_AD_PANEL.pdf",
       panel_combined,
       width = 10,
       height = 14,
       dpi = 300)

#Plot module size 
#Si usas conteos, el histograma estará dominado por el grupo que tenga más módulos.

module_size_hist <- ggplot(module_size_df,
                           aes(x = N_genes)) +
  geom_histogram(binwidth = 120,
                 color = "black",
                 fill = "#3B4992FF") +
  theme_classic(base_size = 14) +
  labs(x = element_blank(),
       y = "Count")

#Save plot

ggsave("module_size_distribution_hist.jpg",
       module_size_hist,
       width = 8,
       height = 8,
       dpi = 300)

ggsave("module_size_distribution_hist.pdf",
       module_size_hist,
       width = 8,
       height = 8,
       dpi = 300)

#Module size per sex

module_size_sex <- ggplot(module_size_df,
                          aes(x = N_genes,
                              fill = Sex)) +
  geom_histogram(binwidth = 100,
                 color = "black",
                 alpha = 0.9,
                 position = "identity") +
  scale_fill_aaas() +
  facet_wrap(~ Sex) +
  theme_classic(base_size = 14) +
  labs(x = "Module size",
       y = "Count") +
  theme(legend.position = "bottom")

#Save plot

ggsave("module_size_distribution_hist_sex.jpg",
       module_size_sex,
       width = 13,
       height = 10,
       dpi = 300)

ggsave("module_size_distribution_hist_sex.pdf",
       module_size_sex,
       width = 13,
       height = 10,
       dpi = 300)

#Panel de module size

module_size_panel <- plot_grid(
  module_size_hist,
  module_size_sex,
  labels = c("A", "B"),
  ncol = 1)

#Save panel

ggsave("module_size_distribution_panel.jpg",
       module_size_panel,
       width = 8,
       height = 10,
       dpi = 300)

ggsave("module_size_distribution_panel.pdf",
       module_size_panel,
       width = 8,
       height = 10,
       dpi = 300)

#Si usas conteos, el histograma estará dominado por el grupo que tenga más módulos, entonces
#Para comparar distribuciones entre grupos con tamaños distintos, es mejor usar probability density o proportion
#¿Cuál es la probabilidad de observar un módulo de tamaño X?
#Plot module size (probability)

prob_module_size_hist <- ggplot(module_size_df,
                           aes(x = N_genes,
                               y = after_stat(density))) +
  geom_histogram(binwidth = 120,
                 color = "black",
                 fill = "#3B4992FF") +
  theme_classic(base_size = 14) +
  labs(x = "Module size",
       y = "Probability density")

#Save plot

ggsave("prob_module_size_distribution_hist.jpg",
       prob_module_size_hist,
       width = 8,
       height = 8,
       dpi = 300)

ggsave("prob_module_size_distribution_hist.pdf",
       prob_module_size_hist,
       width = 8,
       height = 8,
       dpi = 300)

#Module size per sex (probability)

prob_module_size_sex <- ggplot(module_size_df,
                          aes(x = N_genes,
                              y = after_stat(density),
                              fill = Sex)) +
  geom_histogram(binwidth = 100,
                 color = "black",
                 alpha = 0.9,
                 position = "identity") +
  scale_fill_aaas() +
  facet_wrap(~ Sex) +
  theme_classic(base_size = 14) +
  labs(x = "Module size",
       y = "Probability density") +
  theme(legend.position = "bottom")

#Save plot

ggsave("prob_module_size_distribution_hist_sex.jpg",
       prob_module_size_sex,
       width = 13,
       height = 10,
       dpi = 300)

ggsave("prob_module_size_distribution_hist_sex.pdf",
       prob_module_size_sex,
       width = 13,
       height = 10,
       dpi = 300)

#Panel de module size

prob_module_size_panel <- plot_grid(
  module_size_hist,
  module_size_sex,
  labels = c("A", "B"),
  ncol = 1)

#Save panel

ggsave("prob_module_size_distribution_panel.jpg",
       prob_module_size_panel,
       width = 9,
       height = 11,
       dpi = 300)

ggsave("prob_module_size_distribution_panel.pdf",
       prob_module_size_panel,
       width = 9,
       height = 11,
       dpi = 300)

#And by region?

module_size_sex_region <- ggplot(module_size_df,
                                 aes(x = N_genes,
                                     y = after_stat(density),
                                     fill = Sex)) +
  geom_histogram(binwidth = 100,
                 color = "black",
                 alpha = 0.8,
                 position = "identity") +
  scale_fill_aaas() +
  facet_grid(Sex ~ Region) +
  theme_classic(base_size = 14) +
  labs(x = "Module size (genes)",
       y = "Probability density") +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave("module_size_distribution_sex_region.jpg",
       module_size_sex_region,
       width = 16,
       height = 10,
       dpi = 300)

ggsave("module_size_distribution_sex_region.pdf",
       module_size_sex_region,
       width = 16,
       height = 10,
       dpi = 300)

#Diferent
module_size_region_test <- module_size_region_test %>%
  mutate(label = paste0("p adj= ", signif(p_adj, 2))) %>%
  left_join(label_positions, by = "Region")

label_positions <- module_size_df %>%
  group_by(Region) %>%
  summarise(
    x = 2000,
    y = 0.00035   # altura fija en la densidad
  )

module_size_region_test <- left_join(module_size_region_test,
                                     label_positions,
                                     by = "Region")

module_size_sex_region2 <- ggplot(module_size_df,
                                  aes(x = N_genes,
                                      y = after_stat(density),
                                      fill = Sex)) +
  geom_density(color = "black",
               alpha = 0.7) +
  scale_fill_aaas() +
  facet_wrap(~ Region, nrow = 5) +
  geom_text(data = module_size_region_test,
            aes(x = x,
                y = y,
                label = label),
            inherit.aes = FALSE,
            hjust = 0,
            vjust = 1.5,
            size = 4) +
  theme_classic(base_size = 14) +
  labs(x = "Module size (genes)",
       y = "Probability density") +
  theme(legend.position = "bottom")

#Save plots

ggsave("2-module_size_distribution_sex_region.jpg",
       module_size_sex_region2,
       width = 5,
       height = 16,
       dpi = 300)

ggsave("2-module_size_distribution_sex_region.pdf",
       module_size_sex_region2,
       width = 5,
       height = 16,
       dpi = 300)

#Panel

panel2 <- plot_grid(
  module_size_sex_region2,
  global.p,
    labels = c("A", "B"),
    ncol = 2)

#Save panel

ggsave("panel2.jpg",
       panel2,
       width = 9,
       height = 10,
       dpi = 300)

ggsave("panel2.pdf",
       panel2,
       width = 9,
       height = 10,
       dpi = 300)

#Condenado heatmap --- ---

jaccard_df <- jaccard_df %>%
  mutate(
    sex1 = str_extract(Module_1, "_[MF]_") %>% str_remove_all("_"),
    sex2 = str_extract(Module_2, "_[MF]_") %>% str_remove_all("_"),
    
    pheno1 = str_extract(Module_1, "_(AD|CTL)_") %>% str_remove_all("_"),
    pheno2 = str_extract(Module_2, "_(AD|CTL)_") %>% str_remove_all("_"),
    
    group1 = paste0(sex1, pheno1),
    group2 = paste0(sex2, pheno2)
  )

plot_heatmap <- function(df, title) {
  
  ggplot(df, aes(Module_1, Module_2, fill = Jaccard)) +
    
    geom_tile(color = "grey90") +
    
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0.5,
      limits = c(0,1),
      name = "Jaccard"
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      plot.title = element_text(size = 9)
    ) +
    
    labs(
      title = title,
      x = NULL,
      y = NULL
    )
}

make_region_plots <- function(df, label) {
  
  regions <- unique(df$Region)
  
  plots <- lapply(regions, function(r) {
    
    pdat <- df %>% filter(Region == r)
    
    plot_heatmap(pdat, r)
    
  })
  
  names(plots) <- regions
  plots
}

plots_mc_mad  <- make_region_plots(mc_vs_mad, "MC vs MAD")
plots_fc_fad  <- make_region_plots(fc_vs_fad, "FC vs FAD")
plots_mad_fad <- make_region_plots(mad_vs_fad, "MAD vs FAD")


library(cowplot)

row1 <- plot_grid(plotlist = plots_mc_mad, nrow = 1)
row2 <- plot_grid(plotlist = plots_fc_fad, nrow = 1)
row3 <- plot_grid(plotlist = plots_mad_fad, nrow = 1)

panel <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  labels = c("MC vs MAD", "FC vs FAD", "MAD vs FAD"),
  label_size = 14
)

ggsave(
  "module_overlap_panel.png",
  panel,
  width = 14,
  height = 10
)

###################################


library(dplyr)
library(ggplot2)

# ---- 1. Añadir etiqueta de comparación ----

mc_vs_mad$Comparison  <- "MC_vs_MAD"
fc_vs_fad$Comparison  <- "FC_vs_FAD"
mad_vs_fad$Comparison <- "MAD_vs_FAD"

# ---- 2. Unir todos los resultados ----

all_df <- bind_rows(
  mc_vs_mad,
  fc_vs_fad,
  mad_vs_fad
)

# ---- 3. Crear etiquetas únicas para cada eje ----

all_df <- all_df %>%
  mutate(
    x_id = paste(Comparison, Region, Module_1, sep = "_"),
    y_id = paste(Comparison, Region, Module_2, sep = "_")
  )

# ---- 4. Ordenar para que los bloques queden agrupados ----

all_df <- all_df %>%
  arrange(Comparison, Region)

# ---- 5. Convertir a factor para fijar orden en el heatmap ----

all_df$x_id <- factor(all_df$x_id, levels = unique(all_df$x_id))
all_df$y_id <- factor(all_df$y_id, levels = unique(all_df$y_id))

# ---- 6. Calcular posiciones de separación entre bloques ----

block_sizes <- all_df %>%
  group_by(Comparison, Region) %>%
  summarise(n = n_distinct(x_id), .groups="drop")

block_lines <- cumsum(block_sizes$n)

# ---- 7. Dibujar heatmap global ----

p <- ggplot(all_df, aes(x_id, y_id, fill = Jaccard)) +
  
  geom_tile(color = "grey85", size = 0.2) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.5,
    limits = c(0,1),
    name = "Jaccard"
  ) +
  
  geom_vline(xintercept = block_lines + 0.5, size = 0.4) +
  geom_hline(yintercept = block_lines + 0.5, size = 0.4) +
  
  coord_fixed() +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 90, size = 4),
    axis.text.y = element_text(size = 4),
    panel.grid = element_blank()
  ) +
  
  labs(
    x = "",
    y = "",
    title = "Global module overlap across regions, sex and phenotype"
  )

# ---- 8. Guardar figura ----

ggsave(
  "global_module_overlap_heatmap.png",
  p,
  width = 14,
  height = 14
)

################


# Filtrar edges relevantes
edges <- jaccard_df %>%
  select(Module_1, Module_2, Jaccard, Region)

nodes <- tibble(name = V(g)$name) %>%
  mutate(
    sex = ifelse(grepl("_M_", name), "Male", "Female"),
    phenotype = ifelse(grepl("_AD_", name), "AD", "Control"),
    region = sub("^(.*?)_.*", "\\1", name)
  )

V(g)$sex <- nodes$sex
V(g)$phenotype <- nodes$phenotype
V(g)$region <- nodes$region

# Crear grafo
g <- igraph::graph_from_data_frame(edges, directed = FALSE)

png()
ggraph(g, layout = "fr") +
  
  geom_edge_link(aes(width = Jaccard),
                 alpha = 0.6,
                 colour = "grey50") +
  
  geom_node_point(aes(color = phenotype,
                      shape = sex),
                  size = 4) +
  
  theme_void()
dev.off()
