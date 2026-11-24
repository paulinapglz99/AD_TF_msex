#!/usr/bin/env Rscript

# =============================================================================
# Sex-exclusive gene module analysis in Alzheimer's disease
# =============================================================================

# Libraries -------------------------------------------------------------------

pacman::p_load(
  "optparse",
  "data.table",
  "dplyr",
  "stringr",
  "purrr",
  "ggplot2",
  "cowplot",
  "ggsci",
  "vroom",
  "forcats"
)

# Parameters ------------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input_dir"),  type = "character",
              help = "Directory with *_nodes_summary.csv files"),
  make_option(c("-o", "--out_dir"),    type = "character",
              default = "sex_exclusivity_output",
              help = "Output directory"),
  make_option(c("-p", "--pattern"),    type = "character",
              default = "_nodes_summary\\.csv$",
              help = "Regex pattern for input files"),
  make_option(c("-u", "--upper"),      type = "double",  default = 0.7,
              help = "Upper Jaccard threshold (similar)"),
  make_option(c("-l", "--lower"),      type = "double",  default = 0.3,
              help = "Lower Jaccard threshold (exclusive)"),
  make_option(c("-m", "--min_genes"),  type = "integer", default = 3,
              help = "Minimum genes per module")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Override defaults for local run
opt$input_dir  <- "/STORAGE/genut/alzheimer/networks_ad_sex/topos_results"
opt$out_dir    <- "/STORAGE/genut/alzheimer/networks_ad_sex/results_conn_comparisons_renewed"
opt$upper      <- 0.7
opt$lower      <- 0.3
opt$min_genes  <- 3

if (is.null(opt$input_dir)) stop("ERROR: --input_dir must be provided")

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$out_dir)

cat("Input directory :", opt$input_dir, "\n")
cat("Output directory:", opt$out_dir,   "\n")
cat("Upper threshold :", opt$upper,     "\n")
cat("Lower threshold :", opt$lower,     "\n")
cat("Min genes       :", opt$min_genes, "\n")

# =============================================================================
# Functions
# =============================================================================

#' Compute Jaccard index between two gene sets
jaccard_index <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

#' Extract region, sex and phenotype metadata from a file path
parse_filename <- function(filepath) {
  fname <- basename(filepath)
  tibble(
    file      = filepath,
    region    = fname %>%
      str_remove("_nodes_summary\\.csv") %>%
      str_extract("counts_[A-Za-z0-9]+") %>%
      str_remove("counts_"),
    sex       = case_when(
      str_detect(fname, "female") ~ "F",
      str_detect(fname, "male")   ~ "M",
      TRUE ~ NA_character_),
    phenotype = case_when(
      str_detect(fname, "_AD_")      ~ "AD",
      str_detect(fname, "_control_") ~ "CTL",
      TRUE ~ NA_character_)
  )
}

#' Load modules from a nodes_summary CSV, applying a minimum gene filter
load_modules <- function(file, region, sex, phenotype, min_genes) {
  df <- fread(file)[, .(node, membership)]
  keep <- df[, .N, by = membership][N >= min_genes]$membership
  df   <- df[membership %in% keep]
  df[, module_id := paste(region, sex, phenotype, membership, sep = "_")]
  split(df$node, df$module_id)
}

#' Compute all pairwise Jaccard indices for modules within each region
compute_jaccard <- function(modules_list, regions, upper, lower) {
  results <- list()
  for (reg in regions) {
    cat("Processing region:", reg, "\n")
    reg_mods  <- modules_list[str_detect(names(modules_list), paste0("^", reg, "_"))]
    mod_names <- names(reg_mods)
    if (length(mod_names) < 2) next
    for (i in seq_len(length(mod_names) - 1)) {
      for (j in (i + 1):length(mod_names)) {
        jval <- jaccard_index(reg_mods[[mod_names[i]]], reg_mods[[mod_names[j]]])
        results[[length(results) + 1]] <- tibble(
          Region         = reg,
          Module_1       = mod_names[i],
          Module_2       = mod_names[j],
          Jaccard        = jval,
          Classification = case_when(
            jval >= upper ~ "Similar",
            jval <= lower ~ "Dissimilar",
            TRUE          ~ "Intermediate"
          )
        )
      }
    }
  }
  bind_rows(results)
}

#' Classify every module as Similar / Intermediate / Exclusive
classify_modules <- function(modules_list, regions, jaccard_df, upper, lower) {
  results <- list()
  for (reg in regions) {
    reg_mods  <- modules_list[str_detect(names(modules_list), paste0("^", reg, "_"))]
    reg_jacc  <- filter(jaccard_df, Region == reg)
    for (mod in names(reg_mods)) {
      comps <- filter(reg_jacc, Module_1 == mod | Module_2 == mod)
      classif <- if (nrow(comps) == 0) {
        "Exclusive"
      } else {
        maxJ <- max(comps$Jaccard)
        if      (maxJ >= upper) "Similar"
        else if (maxJ <= lower) "Exclusive"
        else                    "Intermediate"
      }
      results[[length(results) + 1]] <- tibble(Region = reg, Module = mod,
                                               Classification = classif)
    }
  }
  bind_rows(results) %>%
    mutate(parts     = str_split(Module, "_"),
           Region    = map_chr(parts, 1),
           Sex       = map_chr(parts, 2),
           Phenotype = map_chr(parts, 3)) %>%
    select(-parts)
}

#' Add Region / Sex / Phenotype columns parsed from the Module ID
parse_module_metadata <- function(df) {
  df %>%
    mutate(parts     = str_split(Module, "_"),
           Region    = map_chr(parts, 1),
           Sex       = map_chr(parts, 2),
           Phenotype = map_chr(parts, 3)) %>%
    select(-parts)
}

#' Save a ggplot in both PDF and JPG
save_plot <- function(plot, filename, width, height, dpi = 300) {
  for (ext in c("pdf", "jpg")) {
    ggsave(paste0(filename, ".", ext), plot = plot,
           width = width, height = height, dpi = dpi)
  }
}

#' Annotate jaccard_df rows with group labels and canonical pair name
annotate_jaccard_pairs <- function(jaccard_df, comparison_levels) {
  jaccard_df %>%
    mutate(
      sex1   = str_extract(Module_1, "(?<=_)[MF](?=_)"),
      pheno1 = str_extract(Module_1, "(?<=_)(AD|CTL)(?=_)"),
      sex2   = str_extract(Module_2, "(?<=_)[MF](?=_)"),
      pheno2 = str_extract(Module_2, "(?<=_)(AD|CTL)(?=_)"),
      group1 = paste(sex1, pheno1, sep = "-"),
      group2 = paste(sex2, pheno2, sep = "-"),
      comparison = map2_chr(group1, group2, ~ {
        sorted <- sort(c(.x, .y))
        paste(sorted[1], "vs", sorted[2])
      })
    ) %>%
    filter(comparison %in% comparison_levels) %>%
    mutate(comparison = factor(comparison, levels = comparison_levels))
}

#' Sort module IDs by region → phenotype → sex → module number
sort_module_id <- function(ids) {
  parts     <- str_split_fixed(ids, "_", 4)
  order_idx <- order(parts[, 1], parts[, 3], parts[, 2], as.integer(parts[, 4]))
  ids[order_idx]
}

#' Build a blocked Jaccard heatmap for a given subset of data
#' block_lines and block_mids are pre-computed from the full dataset
#' so that axis positions stay consistent across region plots
plot_blocked_heatmap <- function(df, block_lines, block_mids,
                                 comparison_levels, fill_scale,
                                 title = NULL) {
  
  mod_levels  <- levels(df$Module_1)
  n_modules_y <- nlevels(df$Module_2)
  
  # geom_text() with explicit dataframe avoids annotate() length mismatch
  block_labels_df <- tibble(
    Module_1 = factor(mod_levels[round(block_mids)], levels = mod_levels),
    Module_2 = factor(mod_levels[n_modules_y],       levels = mod_levels),
    label    = comparison_levels
  )
  
  p <- ggplot(df, aes(x = Module_1, y = Module_2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_vline(xintercept = block_lines + 0.5, color = "grey20", linewidth = 0.5) +
    geom_hline(yintercept = block_lines + 0.5, color = "grey20", linewidth = 0.5) +
    geom_text(data        = block_labels_df,
              aes(x = Module_1, y = Module_2, label = label),
              inherit.aes = FALSE,
              size        = 2.5,
              hjust       = 0.5,
              vjust       = -1,
              color       = "grey20") +
    fill_scale +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          panel.grid       = element_blank(),
          plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
          legend.position  = "right",
          legend.title     = element_text(size = 9),
          plot.margin      = margin(t = 25, r = 10, b = 10, l = 10))
  
  return(p)
}

# =============================================================================
# Data loading
# =============================================================================

files <- list.files(opt$input_dir, pattern = opt$pattern, full.names = TRUE)
if (length(files) == 0) stop("No files found with pattern: ", opt$pattern)

meta         <- map_dfr(files, parse_filename)
modules_list <- pmap(list(meta$file, meta$region, meta$sex, meta$phenotype),
                     load_modules, min_genes = opt$min_genes) %>%
  flatten()

regions <- unique(meta$region)

# =============================================================================
# Pairwise Jaccard
# =============================================================================

jaccard_df <- compute_jaccard(modules_list, regions, opt$upper, opt$lower)
vroom_write(jaccard_df, file.path(opt$out_dir, "jaccard_all_pairs.csv"))

# =============================================================================
# Module classification
# =============================================================================

module_class_df <- classify_modules(modules_list, regions, jaccard_df,
                                    opt$upper, opt$lower)

# Exclusive modules table
exclusive_df <- module_class_df %>%
  filter(Classification == "Exclusive") %>%
  mutate(N_genes = map_int(Module, ~ length(modules_list[[.x]])))

vroom_write(exclusive_df, file.path(opt$out_dir, "exclusive_modules.csv"))

# Module size table
module_size_df <- tibble(
  Module  = names(modules_list),
  N_genes = map_int(modules_list, length)
) %>% parse_module_metadata()

# =============================================================================
# Statistical tests: module size by sex
# =============================================================================

# Shapiro-Wilk normality test within each sex
by(module_size_df$N_genes, module_size_df$Sex, shapiro.test)
# Both groups deviate from normality → use non-parametric tests

# Global Wilcoxon test
wilcox.test(N_genes ~ Sex, data = module_size_df)
# W = 5132.5, p = 0.145 → no significant difference in module size between sexes

# Wilcoxon per region with BH correction
module_size_region_test <- module_size_df %>%
  group_by(Region) %>%
  summarise(p_value = wilcox.test(N_genes ~ Sex)$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))
# All adjusted p-values > 0.6 → no regional differences either

# =============================================================================
# Plots
# =============================================================================

# ---- 1. Local proportions (Region × Sex × Phenotype) -----------------------

proportions_df <- module_class_df %>%
  count(Region, Sex, Phenotype, Classification) %>%
  group_by(Region, Sex, Phenotype) %>%
  mutate(Percentage = n / sum(n))

proportions_p <- ggplot(proportions_df,
                        aes(x = Region, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(Sex ~ Phenotype) +
  scale_fill_aaas(name = "Module classification") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme_classic() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1),
        axis.title      = element_blank(),
        legend.position = "none")

save_plot(proportions_p, "proportions_barplot", width = 10, height = 10)

# ---- 2. Global proportions (Region × Sex) -----------------------------------

global_df <- module_class_df %>%
  count(Region, Sex, Classification) %>%
  group_by(Region, Sex) %>%
  mutate(Percentage     = n / sum(n),
         Classification = factor(Classification,
                                 levels = c("Exclusive", "Intermediate", "Similar")),
         Sex            = factor(Sex, levels = c("F", "M")))

global_p <- ggplot(global_df,
                   aes(x = Region, y = Percentage, fill = Classification)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8) +
  facet_wrap(~ Sex, nrow = 1) +
  scale_fill_aaas(name = "Classification  ") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  labs(y = element_blank()) +
  theme_classic(base_size = 14) +
  theme(strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.8),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        legend.title     = element_text(size = 12),
        legend.position  = "bottom")

save_plot(global_p, "proportions_global_barplot", width = 10, height = 10)

# ---- 3. Combined panel A/B --------------------------------------------------

panel_AB <- plot_grid(proportions_p, global_p,
                      labels = c("A", "B"), ncol = 1, rel_heights = c(1.2, 0.8))

save_plot(panel_AB, "Figure_Modular_Reorganization_Sex_AD_PANEL", width = 10, height = 14)

# ---- 4. Module size distribution --------------------------------------------

# 4a. Overall histogram (counts)
size_hist <- ggplot(module_size_df, aes(x = N_genes)) +
  geom_histogram(binwidth = 120, color = "black", fill = "#3B4992FF") +
  theme_classic(base_size = 14) +
  labs(x = "Module size (genes)", y = "Count")

save_plot(size_hist, "module_size_distribution_hist", width = 8, height = 8)

# 4b. By sex — probability density
size_sex_dens <- ggplot(module_size_df,
                        aes(x = N_genes, y = after_stat(density), fill = Sex)) +
  geom_histogram(binwidth = 100, color = "black", alpha = 0.9, position = "identity") +
  scale_fill_aaas() +
  facet_wrap(~ Sex) +
  theme_classic(base_size = 14) +
  labs(x = "Module size (genes)", y = "Probability density") +
  theme(legend.position = "bottom")

save_plot(size_sex_dens, "prob_module_size_distribution_hist_sex", width = 13, height = 10)

# 4c. By sex × region — density curves with p-adjusted labels
label_positions <- module_size_df %>%
  group_by(Region) %>%
  summarise(x = 2000, y = 0.00035, .groups = "drop")

region_labels <- module_size_region_test %>%
  mutate(label = paste0("p adj = ", signif(p_adj, 2))) %>%
  left_join(label_positions, by = "Region")

size_sex_region <- ggplot(module_size_df,
                          aes(x = N_genes, y = after_stat(density), fill = Sex)) +
  geom_density(color = "black", alpha = 0.7) +
  scale_fill_aaas() +
  facet_wrap(~ Region, nrow = 5) +
  geom_text(data = region_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.5, size = 4) +
  theme_classic(base_size = 14) +
  labs(x = "Module size (genes)", y = "Probability density") +
  theme(legend.position = "bottom")

save_plot(size_sex_region, "module_size_distribution_sex_region", width = 5, height = 16)

# =============================================================================
# Panel C — Jaccard heatmaps: M-AD vs F-AD, M-AD vs M-CTL, F-AD vs F-CTL
# =============================================================================

# AAAS discrete palette repurposed as a continuous gradient
# (anchored at the 3 most informative stops: low, mid, high)
aaas_low  <- "#3B4992FF"   # blue
aaas_mid  <- "#f7f7f7"     # white
aaas_high <- "#EE0000FF"   # red

jaccard_fill_c <- scale_fill_gradient2(
  low      = aaas_low,
  mid      = aaas_mid,
  high     = aaas_high,
  midpoint = 0.5,
  limits   = c(0, 1),
  name     = "Jaccard"
)

# -----------------------------------------------------------------------------
# Helper: one heatmap per comparison pair
# -----------------------------------------------------------------------------

make_comparison_heatmap <- function(comp_pair, group_left, group_right,
                                    jaccard_df_full, fill_scale) {
  
  df <- jaccard_df_full %>%
    mutate(
      sex1   = str_extract(Module_1, "(?<=_)[MF](?=_)"),
      pheno1 = str_extract(Module_1, "(?<=_)(AD|CTL)(?=_)"),
      sex2   = str_extract(Module_2, "(?<=_)[MF](?=_)"),
      pheno2 = str_extract(Module_2, "(?<=_)(AD|CTL)(?=_)"),
      group1 = paste(sex1, pheno1, sep = "-"),
      group2 = paste(sex2, pheno2, sep = "-")
    ) %>%
    filter(
      (group1 == group_left  & group2 == group_right) |
        (group1 == group_right & group2 == group_left)
    ) %>%
    mutate(
      m1 = if_else(group1 == group_left, Module_1, Module_2),
      m2 = if_else(group1 == group_left, Module_2, Module_1)
    )
  
  order_within_group <- function(ids) {
    parts <- str_split_fixed(ids, "_", 4)
    ids[order(parts[, 1], as.integer(parts[, 4]))]
  }
  
  x_levels <- order_within_group(unique(df$m1))
  y_levels <- order_within_group(unique(df$m2))
  
  df <- df %>%
    mutate(m1 = factor(m1, levels = x_levels),
           m2 = factor(m2, levels = y_levels))
  
  # Axis breaks: one label per region block, showing first and last module number
  make_axis_breaks <- function(ids, group_label) {
    parts <- str_split_fixed(as.character(ids), "_", 4)
    tbl <- tibble(
      id     = as.character(ids),
      region = parts[, 1],
      num    = as.integer(parts[, 4]),
      pos    = as.integer(factor(ids, levels = levels(ids)))
    ) %>%
      group_by(region) %>%
      summarise(
        mid   = mean(pos),
        first = min(num),
        last  = max(num),
        .groups = "drop"
      ) %>%
      mutate(label = paste0(region, "\n", group_label, "\n[", first, "-", last, "]"))
    
    # Find the actual module ID closest to the midpoint position
    mid_ids <- map_chr(tbl$mid, function(m) {
      levels(ids)[round(m)]
    })
    
    list(breaks = mid_ids, labels = tbl$label)
  }
  
  x_axis <- make_axis_breaks(df$m1, group_left)
  y_axis <- make_axis_breaks(df$m2, group_right)
  
  # Region separator lines
  x_sep <- df %>%
    mutate(region = str_extract(as.character(m1), "^[^_]+")) %>%
    group_by(region) %>%
    summarise(last_pos = max(as.integer(m1)), .groups = "drop") %>%
    filter(last_pos < nlevels(df$m1)) %>%
    pull(last_pos)
  
  y_sep <- df %>%
    mutate(region = str_extract(as.character(m2), "^[^_]+")) %>%
    group_by(region) %>%
    summarise(last_pos = max(as.integer(m2)), .groups = "drop") %>%
    filter(last_pos < nlevels(df$m2)) %>%
    pull(last_pos)
  
  ggplot(df, aes(x = m1, y = m2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_vline(xintercept = x_sep + 0.5, color = "grey20", linewidth = 0.5) +
    geom_hline(yintercept = y_sep + 0.5, color = "grey20", linewidth = 0.5) +
    fill_scale +
    scale_x_discrete(breaks = x_axis$breaks, labels = x_axis$labels) +
    scale_y_discrete(breaks = y_axis$breaks, labels = y_axis$labels) +
    labs(title = comp_pair, x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_text(size = 6.5, color = "black",
                                     face = "bold", lineheight = 1.2,
                                     hjust = 0.5),
      axis.text.y     = element_text(size = 6.5, color = "black",
                                     face = "bold", lineheight = 1.2),
      axis.ticks      = element_line(linewidth = 0.3),
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 15),
      legend.position = "none",
      plot.margin     = margin(t = 8, r = 4, b = 8, l = 4)
    )
}

# -----------------------------------------------------------------------------
# Build the three heatmaps
# -----------------------------------------------------------------------------

heatmap_mad_fad <- make_comparison_heatmap(
  "M-AD vs F-AD", "M-AD", "F-AD", jaccard_df, jaccard_fill_c
)

heatmap_mad_mc <- make_comparison_heatmap(
  "M-AD vs M-CTL", "M-AD", "M-CTL", jaccard_df, jaccard_fill_c
)

heatmap_fad_fc <- make_comparison_heatmap(
  "F-AD vs F-CTL", "F-AD", "F-CTL", jaccard_df, jaccard_fill_c
)

# -----------------------------------------------------------------------------
# Shared legend — extracted from a temporary full-legend version
# -----------------------------------------------------------------------------

legend_source <- make_comparison_heatmap(
  "M-AD vs F-AD", "M-AD", "F-AD", jaccard_df,
  scale_fill_gradient2(
    low      = aaas_low,
    mid      = aaas_mid,
    high     = aaas_high,
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Jaccard"
  )
) +
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust    = 0.5,
    barwidth       = unit(0.4, "cm"),
    barheight      = unit(5,   "cm"),
    ticks          = TRUE,
    ticks.linewidth = 0.5
  )) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 8)
  )

shared_legend <- get_legend(legend_source)

# -----------------------------------------------------------------------------
# Assemble panel C: heatmaps left, legend centered on the right
# -----------------------------------------------------------------------------

heatmap_stack <- plot_grid(
  heatmap_mad_fad,
  heatmap_mad_mc,
  heatmap_fad_fc,
  ncol  = 1,
  align = "v",
  axis  = "lr"
)

panel_C <- plot_grid(
  heatmap_stack,
  shared_legend,
  ncol       = 2,
  rel_widths = c(0.88, 0.12)
)

save_plot(panel_C, "panel_C_jaccard_heatmaps", width = 7, height = 14)

# -----------------------------------------------------------------------------
# Final panel_2: A + B + C
# -----------------------------------------------------------------------------

panel_2 <- plot_grid(
  plot_grid(size_sex_region, global_p,
            labels     = c("A", "B"),
            ncol       = 2,
            rel_widths = c(0.45, 0.55)),
  panel_C,
  labels     = c("", "C"),
  ncol       = 2,
  rel_widths = c(0.6, 0.4)
)

save_plot(panel_2, "panel2", width = 16, height = 14)

#END

# =============================================================================
# Lollipop plot — exclusive modules by region, sex and phenotype
# Requires: exclusive_df (output of sex_exclusive_modules.R)
# =============================================================================

# ---- Data prep --------------------------------------------------------------

lollipop_df <- exclusive_df %>%
  count(Region, Sex, Phenotype) %>%
  mutate(
    group      = paste(Sex, Phenotype, sep = "-"),
    group      = factor(group, levels = c("F-AD", "F-CTL", "M-AD", "M-CTL")),
    Region     = factor(Region, levels = c("CRB", "TC", "PCC", "HCN", "DLPFC")),
    Sex        = factor(Sex, levels = c("F", "M"))
  )

# Total per region (for ordering annotation)
region_totals <- lollipop_df %>%
  group_by(Region) %>%
  summarise(total = sum(n), .groups = "drop")

# ---- Colors (AAAS palette) --------------------------------------------------

group_colors <- c(
  "F-AD"  = "#3B4992FF",
  "F-CTL" = "#008B45FF",
  "M-AD"  = "#EE0000FF",
  "M-CTL" = "#631879FF"
)

# ---- Plot -------------------------------------------------------------------

lollipop_p <- ggplot(lollipop_df,
                     aes(x = n, y = Region, color = group, group = group)) +
  
  # Segment from 0 to the point
  geom_segment(aes(x = 0, xend = n, y = Region, yend = Region),
               linewidth  = 0.7,
               alpha      = 0.6,
               position   = position_dodge(width = 0.7)) +
  
  # Point at the tip
  geom_point(size     = 3.5,
             alpha    = 0.95,
             position = position_dodge(width = 0.7)) +
  
  # n label next to each point
  geom_text(aes(label = n),
            size     = 3,
            hjust    = -0.5,
            position = position_dodge(width = 0.7),
            show.legend = FALSE) +
  
  # Region total annotation on the right margin
  geom_text(data = region_totals,
            aes(x = total + 2, y = Region,
                label = paste0("n = ", total)),
            inherit.aes = FALSE,
            hjust  = 0,
            size   = 3,
            color  = "grey40",
            fontface = "italic") +
  
  scale_color_manual(values = group_colors, name = "Group") +
  
  # Give space for the n= labels
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  
  facet_wrap(~ Sex,
             ncol    = 1,
             labeller = labeller(Sex = c(F = "Female", M = "Male"))) +
  
  labs(
    x = "Number of exclusive modules",
    y = "Brain region"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black",
                                    linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 12),
    axis.text.y      = element_text(size = 11),
    axis.text.x      = element_text(size = 10),
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    plot.margin      = margin(t = 10, r = 40, b = 10, l = 10)
  )

save_plot(lollipop_p, "lollipop_exclusive_modules", width = 7, height = 7)

# =============================================================================
# UpSet plot — exclusive modules by group (Sex × Phenotype) and region
# Requires: exclusive_df, modules_list (output of sex_exclusive_modules.R)
# Packages: ComplexHeatmap (Bioconductor), ggplot2, dplyr
# =============================================================================

# Install if needed:
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(purrr)

# =============================================================================
# Option A — UpSet by group (F-AD, F-CTL, M-AD, M-CTL)
# Shows which exclusive modules are shared across Sex × Phenotype combinations
# =============================================================================

# Build a named list: one character vector of module IDs per group
exclusive_by_group <- exclusive_df %>%
  mutate(group = paste(Sex, Phenotype, sep = "-")) %>%
  split(.$group) %>%
  map(~ .x$Module)

# Consistent group order and colors (AAAS palette)
group_order  <- c("F-AD", "F-CTL", "M-AD", "M-CTL")
group_colors <- c(
  "F-AD"  = "#3B4992FF",
  "F-CTL" = "#008B45FF",
  "M-AD"  = "#EE0000FF",
  "M-CTL" = "#631879FF"
)

exclusive_by_group <- exclusive_by_group[group_order]

# Create the combination matrix
m_group <- make_comb_mat(exclusive_by_group, mode = "distinct")

# UpSet plot — by group
pdf("upset_exclusive_by_group.pdf", width = 8, height = 5)
UpSet(
  m_group,
  set_order        = group_order,
  comb_order       = order(comb_size(m_group), decreasing = TRUE),
  top_annotation   = HeatmapAnnotation(
    "Intersection\nsize" = anno_barplot(
      comb_size(m_group),
      ylim     = c(0, max(comb_size(m_group)) * 1.15),
      border   = FALSE,
      gp       = gpar(fill = "grey30", col = NA),
      axis_param = list(gp = gpar(fontsize = 9))
    ),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 9)
  ),
  left_annotation  = rowAnnotation(
    "Set size" = anno_barplot(
      set_size(m_group),
      border     = FALSE,
      gp         = gpar(fill = unname(group_colors[group_order]), col = NA),
      axis_param = list(gp = gpar(fontsize = 9)),
      width      = unit(2.5, "cm")
    ),
    annotation_name_gp = gpar(fontsize = 9)
  ),
  pt_size          = unit(4, "mm"),
  lwd              = 2,
  row_names_gp     = gpar(fontsize = 11, fontface = "bold",
                          col = unname(group_colors[group_order])),
  column_title     = "Exclusive modules — intersections by group",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# Also save as PNG for quick inspection
png("upset_exclusive_by_group.png", width = 8, height = 5,
    units = "in", res = 300)
UpSet(
  m_group,
  set_order        = group_order,
  comb_order       = order(comb_size(m_group), decreasing = TRUE),
  top_annotation   = HeatmapAnnotation(
    "Intersection\nsize" = anno_barplot(
      comb_size(m_group),
      ylim     = c(0, max(comb_size(m_group)) * 1.15),
      border   = FALSE,
      gp       = gpar(fill = "grey30", col = NA),
      axis_param = list(gp = gpar(fontsize = 9))
    ),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 9)
  ),
  left_annotation  = rowAnnotation(
    "Set size" = anno_barplot(
      set_size(m_group),
      border     = FALSE,
      gp         = gpar(fill = unname(group_colors[group_order]), col = NA),
      axis_param = list(gp = gpar(fontsize = 9)),
      width      = unit(2.5, "cm")
    ),
    annotation_name_gp = gpar(fontsize = 9)
  ),
  pt_size          = unit(4, "mm"),
  lwd              = 2,
  row_names_gp     = gpar(fontsize = 11, fontface = "bold",
                          col = unname(group_colors[group_order])),
  column_title     = "Exclusive modules — intersections by group",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# =============================================================================
# UpSet plot — genes in exclusive modules, by group (Sex × Phenotype)
# Built with ggplot2 + ggupset (no ComplexHeatmap)
# Requires: exclusive_df, modules_list (output of sex_exclusive_modules.R)
# install.packages("ggupset")
# =============================================================================

pacman::p_load("dplyr", "purrr", "tidyr", "ggplot2", "ggupset", "cowplot")

# AAAS colors
group_colors <- c(
  "F-AD"  = "#3B4992FF",
  "F-CTL" = "#008B45FF",
  "M-AD"  = "#EE0000FF",
  "M-CTL" = "#631879FF"
)

group_order <- c("F-AD", "F-CTL", "M-AD", "M-CTL")

# =============================================================================
# 1. Build gene-level membership table
# =============================================================================

# For each exclusive module, get its genes from modules_list
gene_membership <- exclusive_df %>%
  mutate(group = paste(Sex, Phenotype, sep = "-")) %>%
  select(Module, group) %>%
  mutate(genes = map(Module, ~ modules_list[[.x]])) %>%
  unnest(genes) %>%
  select(gene = genes, group)

# For each gene, collect the set of groups it appears in
gene_groups <- gene_membership %>%
  group_by(gene) %>%
  summarise(groups = list(sort(unique(group))), .groups = "drop")

cat("Total unique genes:", nrow(gene_groups), "\n")
cat("Distribution of membership size:\n")
gene_groups %>%
  mutate(n_groups = lengths(groups)) %>%
  count(n_groups) %>%
  print()

# =============================================================================
# 2. UpSet plot — by group
# =============================================================================

# Intersection counts (for the bar chart on top)
intersection_counts <- gene_groups %>%
  count(groups) %>%
  arrange(desc(n))

# Set sizes (for the bar chart on the left)
set_sizes <- gene_membership %>%
  distinct(gene, group) %>%
  count(group) %>%
  mutate(group = factor(group, levels = rev(group_order)))

# --- Top bar: intersection sizes ---
p_bars <- ggplot(gene_groups, aes(x = groups)) +
  geom_bar(fill = "grey25", width = 0.6) +
  scale_x_upset(order_by = "freq", sets = group_order) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Genes") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.line.x     = element_blank(),
    plot.margin     = margin(t = 10, r = 10, b = 0, l = 10)
  )

# --- Bottom: intersection matrix (ggupset handles this via theme_combmatrix) ---
p_matrix <- ggplot(gene_groups, aes(x = groups)) +
  geom_bar(fill = "grey25", width = 0.6) +
  scale_x_upset(order_by = "freq", sets = group_order) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Genes") +
  theme_combmatrix(
    combmatrix.label.text         = element_text(
      size = 11, face = "bold",
      colour = group_colors[group_order]
    ),
    combmatrix.label.extra_spacing = 8,
    combmatrix.panel.point.size    = 4,
    combmatrix.panel.line.size     = 0.8,
    combmatrix.panel.point.color.fill   = "grey20",
    combmatrix.panel.point.color.empty  = "grey88"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
  )

# --- Left bar: set sizes ---
p_sets <- ggplot(set_sizes,
                 aes(x = n, y = group, fill = group)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  scale_y_discrete(limits = rev(group_order)) +
  labs(x = "Set size", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_text(face = "bold",
                                colour = group_colors[rev(group_order)]),
    plot.margin  = margin(t = 0, r = 0, b = 10, l = 10)
  )

# --- Combine ---
upset_group <- p_matrix

save_plot(upset_group, "upset_genes_exclusive_groups",
          width = 9, height = 5)

# =============================================================================
# 3. UpSet plot — by region, one per group
# =============================================================================

region_colors <- c(
  CRB   = "#3B4992FF",
  TC    = "#EE0000FF",
  PCC   = "#008B45FF",
  HCN   = "#631879FF",
  DLPFC = "#A20056FF"
)

make_region_upset <- function(sex_filter, pheno_filter) {
  
  grp_label <- paste0(sex_filter, "-", pheno_filter)
  
  # Genes per region for this group's exclusive modules
  gene_region <- exclusive_df %>%
    filter(Sex == sex_filter, Phenotype == pheno_filter) %>%
    mutate(genes = map(Module, ~ modules_list[[.x]])) %>%
    unnest(genes) %>%
    select(gene = genes, Region)
  
  regions_present <- sort(unique(gene_region$Region))
  
  if (length(regions_present) < 2) {
    cat("Skipping", grp_label, "— fewer than 2 regions\n")
    return(invisible(NULL))
  }
  
  # One list entry per gene: which regions contain it
  gene_regions_list <- gene_region %>%
    group_by(gene) %>%
    summarise(regions = list(sort(unique(Region))), .groups = "drop")
  
  p <- ggplot(gene_regions_list, aes(x = regions)) +
    geom_bar(fill = "grey25", width = 0.6) +
    scale_x_upset(order_by = "freq", sets = regions_present) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = paste0("Gene overlap by region — ", grp_label,
                     " exclusive modules"),
      y     = "Genes"
    ) +
    theme_combmatrix(
      combmatrix.label.text = element_text(
        size   = 10,
        face   = "bold",
        colour = region_colors[regions_present]
      ),
      combmatrix.label.extra_spacing  = 8,
      combmatrix.panel.point.size     = 4,
      combmatrix.panel.line.size      = 0.8,
      combmatrix.panel.point.color.fill  = "grey20",
      combmatrix.panel.point.color.empty = "grey88"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 11),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )
  
  fname <- paste0("upset_genes_", sex_filter, "_", pheno_filter, "_byregion")
  save_plot(p, fname, width = 8, height = 5)
  cat("Saved:", fname, "\n")
}

walk2(
  c("F", "F", "M", "M"),
  c("AD", "CTL", "AD", "CTL"),
  make_region_upset
)

cat("\nDone.\n")
