#Script to identify DEGs intersection between sex and conditions

#Libraries --- ---
library(tidyverse)

# ==========================================================
# 1. Definir paths
# ==========================================================

base_path <- "~/Desktop/local_work/AD_TF_sex/DEGS_new/full_DEGs"
output_path <- file.path(base_path, "output_overlap")

# Crear carpeta si no existe
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# ==========================================================
# 2. Listar y leer archivos
# ==========================================================

files <- list.files(path = base_path,
                    pattern = "\\.txt$",
                    full.names = TRUE)

deg_list <- list()

for (f in files) {
  
  name <- basename(f)
  
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (!"delabel" %in% colnames(df)) {
    stop(paste("El archivo", name, "no tiene columna 'delabel'"))
  }
  
  deg_list[[name]] <- df
}

# ==========================================================
# 3. Función para extraer sets
# ==========================================================

get_region_sets <- function(region_name) {
  
  region_files <- names(deg_list)[grepl(region_name, names(deg_list))]
  
  file_A <- region_files[grepl("control_female_vs_AD_female", region_files)]
  file_B <- region_files[grepl("AD_female_vs_AD_male", region_files)]
  file_C <- region_files[grepl("control_male_vs_AD_male", region_files)]
  
  if (length(file_A) == 0 |
      length(file_B) == 0 |
      length(file_C) == 0) {
    
    stop(paste("Archivos faltantes para región", region_name))
  }
  
  A <- unique(deg_list[[file_A]]$delabel)
  B <- unique(deg_list[[file_B]]$delabel)
  C <- unique(deg_list[[file_C]]$delabel)
  
  list(A = A,
       B = B,
       C = C)
}

# ==========================================================
# 4. Calcular intersecciones y guardar resultados
# ==========================================================

regions <- c("DLPFC",
             "HCN",
             "Mayo_CRB",
             "Mayo_TC",
             "PCC")

summary_results <- data.frame()

for (r in regions) {
  
  sets <- get_region_sets(r)
  
  I_AB  <- intersect(sets$A, sets$B)
  I_ABC <- Reduce(intersect, sets)
  
  # Guardar A ∩ B
  write.table(I_AB,
              file = file.path(output_path,
                               paste0(r, "_overlap_A_B.txt")),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  # Guardar A ∩ B ∩ C
  write.table(I_ABC,
              file = file.path(output_path,
                               paste0(r, "_overlap_A_B_C.txt")),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  tmp <- data.frame(
    region = r,
    n_A = length(sets$A),
    n_B = length(sets$B),
    n_C = length(sets$C),
    n_AB = length(I_AB),
    n_ABC = length(I_ABC)
  )
  
  summary_results <- rbind(summary_results, tmp)
}

# Guardar tabla resumen
write.csv(summary_results,
          file = file.path(output_path, "summary_overlap_counts.csv"),
          row.names = FALSE)

summary_results

install.packages("ComplexUpset")
library(ComplexUpset)
library(ggplot2)

for (r in regions) {
  
  sets <- get_region_sets(r)
  
  # Crear dataframe de presencia/ausencia
  all_genes <- unique(c(sets$A, sets$B, sets$C))
  
  presence_df <- data.frame(
    gene = all_genes,
    female_AD_vs_control = all_genes %in% sets$A,
    female_AD_vs_male_AD = all_genes %in% sets$B,
    male_AD_vs_control   = all_genes %in% sets$C
  )
  
  # Generar UpSet plot
  p <- upset(
    presence_df,
    c("female_AD_vs_control",
      "female_AD_vs_male_AD",
      "male_AD_vs_control"),
    name = "DEG sets"
  ) +
    ggtitle(paste("UpSet Plot -", r))
  
  # Guardar en PDF
  ggsave(
    filename = file.path(output_path,
                         paste0("UpSet_", r, ".pdf")),
    plot = p,
    width = 8,
    height = 6
  )
}


plot_list <- list()

for (r in regions) {
  
  sets <- get_region_sets(r)
  all_genes <- unique(c(sets$A, sets$B, sets$C))
  
  presence_df <- data.frame(
    gene = all_genes,
    female_AD_vs_control = all_genes %in% sets$A,
    female_AD_vs_male_AD = all_genes %in% sets$B,
    male_AD_vs_control   = all_genes %in% sets$C
  )
  
  p <- upset(
    presence_df,
    c("female_AD_vs_control",
      "female_AD_vs_male_AD",
      "male_AD_vs_control"),
    name = "DEGs",
    width_ratio = 0.2
  ) +
    ggtitle(r) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      text = element_text(size = 10)
    )
  
  plot_list[[r]] <- p
}

# ----------------------------------------------------------
# Construir panel (3 arriba, 2 abajo centrados)
# ----------------------------------------------------------

panel_plot <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) /
  (plot_list[[4]] | plot_list[[5]])

panel_plot <- panel_plot +
  plot_annotation(
    title = "Overlap of Differentially Expressed Genes Across Sex-Specific AD Comparisons",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

panel_plot

# ----------------------------------------------------------
# Guardar figura final
# ----------------------------------------------------------

ggsave(
  filename = file.path(output_path, "UpSet_panel_all_regions.pdf"),
  plot = panel_plot,
  width = 18,
  height = 12
)