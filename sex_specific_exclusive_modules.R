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

vroom::vroom_write(summary_df, file.path(opt$out_dir, "exclusive_summary.csv"))

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
    #axis.text.x = element_text(angle = 45, hjust = 1), 
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

#END