#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ============================================================
# ARGUMENT PARSER
# ============================================================

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

if (is.null(opt$input_dir)) {
  stop("ERROR: --input_dir must be provided")
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Input directory:", opt$input_dir, "\n")
cat("Output directory:", opt$out_dir, "\n")
cat("Upper threshold:", opt$upper, "\n")
cat("Lower threshold:", opt$lower, "\n")
cat("Min genes:", opt$min_genes, "\n")

# ============================================================
# FILENAME PARSER
# ============================================================

parse_filename <- function(filepath) {
  
  fname <- basename(filepath)
  
  sex <- case_when(
    str_detect(fname, "male") ~ "M",
    str_detect(fname, "female") ~ "F",
    TRUE ~ NA_character_
  )
  
  phenotype <- case_when(
    str_detect(fname, "_AD_") ~ "AD",
    str_detect(fname, "_control_") ~ "CTL",
    TRUE ~ NA_character_
  )
  
  region <- fname %>%
    str_remove("_nodes_summary\\.csv") %>%
    str_extract("counts_.*") %>%
    str_remove("counts_")
  
  tibble(
    file = filepath,
    region = region,
    sex = sex,
    phenotype = phenotype
  )
}

# ============================================================
# LOAD FILES
# ============================================================

files <- list.files(
  opt$input_dir,
  pattern = opt$pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files found with pattern: ", opt$pattern)
}

meta <- map_dfr(files, parse_filename)

# ============================================================
# LOAD MODULES
# ============================================================

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

# ============================================================
# JACCARD FUNCTION
# ============================================================

jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

# ============================================================
# PAIRWISE JACCARD (WITHIN REGION)
# ============================================================

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

fwrite(jaccard_df,
       file.path(opt$out_dir, "jaccard_all_pairs.csv"))

# ============================================================
# STRUCTURAL EXCLUSIVITY
# ============================================================

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

fwrite(exclusive_df,
       file.path(opt$out_dir, "exclusive_modules.csv"))

# ============================================================
# SUMMARY
# ============================================================

summary_df <- exclusive_df %>%
  mutate(
    Sex = str_extract(Module, "_F_|_M_"),
    Phenotype = str_extract(Module, "_AD_|_CTL_")
  ) %>%
  group_by(Region, Sex, Phenotype) %>%
  summarise(N_exclusive = n(), .groups = "drop")

fwrite(summary_df,
       file.path(opt$out_dir, "exclusive_summary.csv"))

cat("DONE.\n")
