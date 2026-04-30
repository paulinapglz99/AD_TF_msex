#!/usr/bin/env Rscript
# ==============================================================================
# Análisis Comparativo: Efecto del AD entre Sexos
# ==============================================================================
# Pregunta: ¿Qué genes se expresan diferencialmente en AD de manera
#           diferente entre hombres y mujeres?
#
# Estrategia:
#   1. Actualizar símbolos en todos los CSVs existentes
#   2. Análisis de interacción (estadísticamente riguroso)
#   3. Comparación de listas (intersección/exclusión)
#   4. Visualizaciones: Venn, UpSet, scatter LFC, heatmap, volcano conjunto
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ==============================================================================
# CONFIGURACIÓN
# ==============================================================================

TISSUES   <- c("CRB", "TC")
RESULTS_DIR <- "results"
OUT_DIR     <- "results/comparative"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5

# ==============================================================================
# 1. FUNCIÓN: CONVERSIÓN ENSEMBL → SYMBOL
# ==============================================================================

ensembl_to_symbol <- function(ensembl_ids) {
  symbols <- mapIds(
    org.Hs.eg.db,
    keys      = ensembl_ids,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  data.frame(
    ensembl_gene_id    = names(symbols),
    external_gene_name = unname(symbols),
    stringsAsFactors   = FALSE
  )
}

# ==============================================================================
# 2. ACTUALIZAR SÍMBOLOS EN TODOS LOS CSVs EXISTENTES
# ==============================================================================

message("═══════════════════════════════════════════")
message("Actualizando símbolos en resultados existentes")
message("═══════════════════════════════════════════\n")

update_symbols <- function(csv_path) {
  df <- read.csv(csv_path)
  
  # Convertir IDs
  sym <- ensembl_to_symbol(df$ensembl_gene_id)
  
  # Reemplazar external_gene_name y delabel
  df <- df %>%
    dplyr::select(-any_of(c("external_gene_name", "delabel"))) %>%
    left_join(sym, by = "ensembl_gene_id") %>%
    mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA))
  
  write.csv(df, csv_path, row.names = FALSE)
  message(sprintf("  ✅ Actualizado: %s", basename(csv_path)))
  return(df)
}

# Actualizar todos los CSVs
all_csvs <- list.files(RESULTS_DIR, pattern = "\\.csv$", 
                       recursive = TRUE, full.names = TRUE)
all_csvs <- all_csvs[!grepl("summary|comparative|sessionInfo", all_csvs)]

for (csv in all_csvs) {
  tryCatch(update_symbols(csv), error = function(e) {
    message(sprintf("  ⚠️  Error en %s: %s", basename(csv), e$message))
  })
}

# ==============================================================================
# 3. CARGAR DATOS ACTUALIZADOS
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("Cargando resultados")
message("═══════════════════════════════════════════\n")

load_results <- function(tissue) {
  list(
    males       = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_males_AD_vs_control_DEGs.csv", tissue))),
    females     = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_females_AD_vs_control_DEGs.csv", tissue))),
    interaction = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_sex_diagnosis_interaction_DEGs.csv", tissue))),
    males_all   = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_males_AD_vs_control_all_genes.csv", tissue))),
    females_all = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_females_AD_vs_control_all_genes.csv", tissue))),
    interaction_all = read.csv(file.path(RESULTS_DIR, tissue, sprintf("%s_sex_diagnosis_interaction_all_genes.csv", tissue)))
  )
}

data <- lapply(setNames(TISSUES, TISSUES), load_results)

# Resumen
for (t in TISSUES) {
  cat(sprintf("\n%s:\n", t))
  cat(sprintf("  Males DEGs:       %d (UP: %d | DOWN: %d)\n",
              nrow(data[[t]]$males),
              sum(data[[t]]$males$diffexpressed == "UP"),
              sum(data[[t]]$males$diffexpressed == "DOWN")))
  cat(sprintf("  Females DEGs:     %d (UP: %d | DOWN: %d)\n",
              nrow(data[[t]]$females),
              sum(data[[t]]$females$diffexpressed == "UP"),
              sum(data[[t]]$females$diffexpressed == "DOWN")))
  cat(sprintf("  Interaction DEGs: %d (UP: %d | DOWN: %d)\n",
              nrow(data[[t]]$interaction),
              sum(data[[t]]$interaction$diffexpressed == "UP"),
              sum(data[[t]]$interaction$diffexpressed == "DOWN")))
}

# ==============================================================================
# 4. CATEGORIZACIÓN DE GENES
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("Categorizando genes por patrón sex-specific")
message("═══════════════════════════════════════════\n")

categorize_genes <- function(tissue) {
  m_ids <- data[[tissue]]$males$ensembl_gene_id
  f_ids <- data[[tissue]]$females$ensembl_gene_id
  
  # Categorías
  male_only   <- setdiff(m_ids, f_ids)
  female_only <- setdiff(f_ids, m_ids)
  shared      <- intersect(m_ids, f_ids)
  
  # Para genes compartidos: ¿mismo sentido?
  shared_m <- data[[tissue]]$males %>% 
    filter(ensembl_gene_id %in% shared) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, lfc_male = log2FoldChange, dir_male = diffexpressed)
  
  shared_f <- data[[tissue]]$females %>%
    filter(ensembl_gene_id %in% shared) %>%
    dplyr::select(ensembl_gene_id, lfc_female = log2FoldChange, dir_female = diffexpressed)
  
  shared_df <- shared_m %>% left_join(shared_f, by = "ensembl_gene_id") %>%
    mutate(concordant = dir_male == dir_female)
  
  concordant   <- shared_df %>% filter(concordant) %>% pull(ensembl_gene_id)
  discordant   <- shared_df %>% filter(!concordant) %>% pull(ensembl_gene_id)
  
  cat(sprintf("\n%s:\n", tissue))
  cat(sprintf("  Male-only DEGs:          %d\n", length(male_only)))
  cat(sprintf("  Female-only DEGs:        %d\n", length(female_only)))
  cat(sprintf("  Shared (same direction): %d\n", length(concordant)))
  cat(sprintf("  Shared (opposite dir):   %d  ← más interesantes biológicamente\n", length(discordant)))
  cat(sprintf("  Interaction DEGs:        %d  ← estadísticamente rigurosos\n", nrow(data[[tissue]]$interaction)))
  
  list(
    male_only   = male_only,
    female_only = female_only,
    concordant  = concordant,
    discordant  = discordant,
    shared_df   = shared_df,
    interaction = data[[tissue]]$interaction$ensembl_gene_id
  )
}

cats <- lapply(setNames(TISSUES, TISSUES), categorize_genes)

# ==============================================================================
# 5. PLOTS
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("Generando visualizaciones")
message("═══════════════════════════════════════════\n")

# ── PLOT 1: Scatter LFC Males vs Females ────────────────────────────────────
# El más informativo: muestra la relación entre los efectos de AD en cada sexo

plot_lfc_scatter <- function(tissue) {
  
  # Usar all_genes para tener todos los puntos
  lfc_m <- data[[tissue]]$males_all %>%
    dplyr::select(ensembl_gene_id, external_gene_name, lfc_male = log2FoldChange, 
           padj_male = padj)
  
  lfc_f <- data[[tissue]]$females_all %>%
    dplyr::select(ensembl_gene_id, lfc_female = log2FoldChange, padj_female = padj)
  
  df <- inner_join(lfc_m, lfc_f, by = "ensembl_gene_id") %>%
    filter(!is.na(lfc_male), !is.na(lfc_female)) %>%
    mutate(
      category = case_when(
        padj_male < PADJ_THRESH & abs(lfc_male) > LFC_THRESH &
          padj_female < PADJ_THRESH & abs(lfc_female) > LFC_THRESH ~ "Both sexes",
        padj_male < PADJ_THRESH & abs(lfc_male) > LFC_THRESH       ~ "Male only",
        padj_female < PADJ_THRESH & abs(lfc_female) > LFC_THRESH   ~ "Female only",
        TRUE                                                          ~ "Not DEG"
      ),
      category = factor(category, levels = c("Not DEG", "Both sexes", 
                                             "Male only", "Female only")),
      # Label solo genes de interacción o discordantes interesantes
      label = ifelse(
        ensembl_gene_id %in% data[[tissue]]$interaction$ensembl_gene_id,
        external_gene_name, NA
      )
    )
  
  colors <- c("Not DEG" = "grey85", "Both sexes" = "#7B2D8B",
              "Male only" = "#2196F3", "Female only" = "#E91E8C")
  
  # Correlación
  r <- cor(df$lfc_male, df$lfc_female, use = "complete.obs")
  
  p <- ggplot(df, aes(x = lfc_male, y = lfc_female, color = category)) +
    geom_point(data = df %>% filter(category == "Not DEG"), 
               alpha = 0.15, size = 0.4) +
    geom_point(data = df %>% filter(category != "Not DEG"), 
               alpha = 0.7, size = 1.2) +
    # Líneas de referencia
    geom_hline(yintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                color = "grey40", linewidth = 0.5) +
    # Labels de genes de interacción
    ggrepel::geom_text_repel(
      data = df %>% filter(!is.na(label)),
      aes(label = label), size = 3, max.overlaps = 20,
      fontface = "italic", color = "black"
    ) +
    scale_color_manual(values = colors) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.3f", r), size = 3.5, color = "grey30") +
    labs(
      title    = sprintf("%s: AD Effect in Males vs Females", tissue),
      subtitle = "Diagonal = same effect in both sexes | Off-diagonal = sex-specific",
      x        = "log2FC (Males: AD vs Control)",
      y        = "log2FC (Females: AD vs Control)",
      color    = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))
  
  out <- file.path(OUT_DIR, sprintf("%s_lfc_scatter_sex.png", tissue))
  ggsave(out, p, width = 7, height = 6, dpi = 300)
  message(sprintf("  ✅ %s", basename(out)))
  
  return(invisible(df))
}

# ── PLOT 2: Volcano conjunto (Males y Females lado a lado) ───────────────────

plot_volcano_combined <- function(tissue) {
  
  m <- data[[tissue]]$males_all %>% mutate(sex = "Males")
  f <- data[[tissue]]$females_all %>% mutate(sex = "Females")
  
  df <- bind_rows(m, f) %>%
    filter(!is.na(padj)) %>%
    mutate(
      diffexpressed = case_when(
        log2FoldChange >  LFC_THRESH & padj < PADJ_THRESH ~ "UP",
        log2FoldChange < -LFC_THRESH & padj < PADJ_THRESH ~ "DOWN",
        TRUE ~ "NO"
      ),
      sex = factor(sex, levels = c("Males", "Females"))
    )
  
  # Top 10 genes por sexo para labelear
  top_genes <- df %>%
    filter(diffexpressed != "NO") %>%
    group_by(sex) %>%
    slice_min(padj, n = 10) %>%
    pull(ensembl_gene_id)
  
  df <- df %>%
    mutate(label = ifelse(ensembl_gene_id %in% top_genes, external_gene_name, NA))
  
  colors <- c("UP" = "#E64B35", "DOWN" = "#4DBBD5", "NO" = "grey80")
  
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
    geom_point(data = df %>% filter(diffexpressed == "NO"), 
               alpha = 0.2, size = 0.4) +
    geom_point(data = df %>% filter(diffexpressed != "NO"), 
               alpha = 0.7, size = 0.8) +
    ggrepel::geom_text_repel(
      aes(label = label), size = 2.5, max.overlaps = 15,
      fontface = "italic", na.rm = TRUE
    ) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH), 
               linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_hline(yintercept = -log10(PADJ_THRESH), 
               linetype = "dashed", color = "grey30", linewidth = 0.4) +
    scale_color_manual(values = colors) +
    facet_wrap(~sex, ncol = 2) +
    labs(
      title = sprintf("%s: AD vs Control by Sex", tissue),
      x = "log2 Fold Change", y = "-log10(adj. p-value)", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  
  out <- file.path(OUT_DIR, sprintf("%s_volcano_combined.png", tissue))
  ggsave(out, p, width = 12, height = 6, dpi = 300)
  message(sprintf("  ✅ %s", basename(out)))
}

# ── PLOT 3: UpSet plot (requiere UpSetR) ─────────────────────────────────────

plot_upset <- function(tissue) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("  ⚠️  UpSetR no instalado: install.packages('UpSetR')")
    return(invisible(NULL))
  }
  library(UpSetR)
  
  all_genes <- unique(c(
    data[[tissue]]$males$ensembl_gene_id,
    data[[tissue]]$females$ensembl_gene_id,
    data[[tissue]]$interaction$ensembl_gene_id
  ))
  
  upset_df <- data.frame(
    gene       = all_genes,
    Males      = as.integer(all_genes %in% data[[tissue]]$males$ensembl_gene_id),
    Females    = as.integer(all_genes %in% data[[tissue]]$females$ensembl_gene_id),
    Interaction = as.integer(all_genes %in% data[[tissue]]$interaction$ensembl_gene_id)
  )
  
  out <- file.path(OUT_DIR, sprintf("%s_upset.png", tissue))
  png(out, width = 10, height = 6, units = "in", res = 300)
  print(upset(
    upset_df,
    sets           = c("Males", "Females", "Interaction"),
    order.by       = "freq",
    sets.bar.color = c("#2196F3", "#E91E8C", "#FF9800"),
    text.scale     = 1.4,
    mainbar.y.label = "DEG Intersections",
    sets.x.label    = "DEGs per Group"
  ))
  dev.off()
  message(sprintf("  ✅ %s", basename(out)))
}

# ── PLOT 4: Heatmap de genes discordantes ────────────────────────────────────
# Los más interesantes: misma dirección opuesta en AD según sexo

plot_discordant_heatmap <- function(tissue) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("  ⚠️  pheatmap no instalado: install.packages('pheatmap')")
    return(invisible(NULL))
  }
  library(pheatmap)
  
  disc_ids <- cats[[tissue]]$discordant
  
  if (length(disc_ids) == 0) {
    message(sprintf("  ⚠️  %s: sin genes discordantes para heatmap", tissue))
    return(invisible(NULL))
  }
  
  # LFC de males y females para genes discordantes
  lfc_m <- data[[tissue]]$males_all %>%
    filter(ensembl_gene_id %in% disc_ids) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, lfc_male = log2FoldChange)
  
  lfc_f <- data[[tissue]]$females_all %>%
    filter(ensembl_gene_id %in% disc_ids) %>%
    dplyr::select(ensembl_gene_id, lfc_female = log2FoldChange)
  
  mat <- lfc_m %>%
    left_join(lfc_f, by = "ensembl_gene_id") %>%
    mutate(label = ifelse(!is.na(external_gene_name), external_gene_name, ensembl_gene_id)) %>%
    column_to_rownames("label") %>%
    dplyr::select(lfc_male, lfc_female) %>%
    as.matrix()
  
  colnames(mat) <- c("Males (AD vs Ctrl)", "Females (AD vs Ctrl)")
  
  # Limitar a top 50 por magnitud
  if (nrow(mat) > 50) {
    top_idx <- order(rowSums(abs(mat)), decreasing = TRUE)[1:50]
    mat <- mat[top_idx, ]
  }
  
  out <- file.path(OUT_DIR, sprintf("%s_discordant_heatmap.png", tissue))
  
  lim <- max(abs(mat), na.rm = TRUE)
  
  pheatmap(
    mat,
    color            = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    breaks           = seq(-lim, lim, length.out = 101),
    cluster_cols     = FALSE,
    cluster_rows     = TRUE,
    fontsize_row     = 7,
    fontsize_col     = 11,
    border_color     = NA,
    main             = sprintf("%s: Discordant DEGs\n(Opposite AD effect by sex)", tissue),
    filename         = out,
    width            = 6,
    height           = max(4, nrow(mat) * 0.2 + 2)
  )
  message(sprintf("  ✅ %s", basename(out)))
}

# ── PLOT 5: Bar plot resumen ──────────────────────────────────────────────────

plot_summary_bars <- function() {
  
  summary_df <- bind_rows(lapply(TISSUES, function(t) {
    data.frame(
      Tissue   = t,
      Category = c("Male-only", "Female-only", "Shared concordant", 
                   "Shared discordant", "Interaction"),
      Count    = c(
        length(cats[[t]]$male_only),
        length(cats[[t]]$female_only),
        length(cats[[t]]$concordant),
        length(cats[[t]]$discordant),
        length(cats[[t]]$interaction)
      )
    )
  }))
  
  summary_df$Category <- factor(summary_df$Category,
                                levels = c("Male-only", "Female-only", "Shared concordant",
                                           "Shared discordant", "Interaction"))
  
  colors <- c("Male-only"          = "#2196F3",
              "Female-only"        = "#E91E8C",
              "Shared concordant"  = "#7B2D8B",
              "Shared discordant"  = "#FF5722",
              "Interaction"        = "#FF9800")
  
  p <- ggplot(summary_df, aes(x = Tissue, y = Count, fill = Category)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9),
              vjust = -0.3, size = 3.2) +
    scale_fill_manual(values = colors) +
    labs(
      title = "DEG Categories by Sex and Brain Region",
      x = NULL, y = "Number of Genes", fill = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  out <- file.path(OUT_DIR, "summary_bars.png")
  ggsave(out, p, width = 9, height = 5, dpi = 300)
  message(sprintf("  ✅ %s", basename(out)))
}

# ── GENERAR TODOS LOS PLOTS ───────────────────────────────────────────────────

# Instalar ggrepel si falta
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggrepel)

for (t in TISSUES) {
  message(sprintf("\n--- %s ---", t))
  plot_lfc_scatter(t)
  plot_volcano_combined(t)
  plot_upset(t)
  plot_discordant_heatmap(t)
}
plot_summary_bars()

# ==============================================================================
# 6. TABLAS DE SALIDA
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("Guardando tablas de genes de interés")
message("═══════════════════════════════════════════\n")

for (t in TISSUES) {
  
  # Tabla: genes discordantes (mismo gen, dirección opuesta en AD según sexo)
  if (length(cats[[t]]$discordant) > 0) {
    disc_table <- cats[[t]]$shared_df %>%
      filter(ensembl_gene_id %in% cats[[t]]$discordant) %>%
      arrange(desc(abs(lfc_male - lfc_female))) %>%
      mutate(tissue = t,
             interpretation = ifelse(lfc_male > 0,
                                     "UP in male AD, DOWN in female AD",
                                     "DOWN in male AD, UP in female AD"))
    
    write.csv(disc_table,
              file.path(OUT_DIR, sprintf("%s_discordant_genes.csv", t)),
              row.names = FALSE)
    message(sprintf("  ✅ %s: %d discordant genes", t, nrow(disc_table)))
  }
  
  # Tabla: genes de interacción con contexto de ambos sexos
  if (nrow(data[[t]]$interaction) > 0) {
    
    int_context <- data[[t]]$interaction %>%
      dplyr::select(ensembl_gene_id, external_gene_name,
             lfc_interaction = log2FoldChange,
             padj_interaction = padj, diffexpressed) %>%
      left_join(
        data[[t]]$males_all %>% dplyr::select(ensembl_gene_id,
                                       lfc_male = log2FoldChange, padj_male = padj),
        by = "ensembl_gene_id"
      ) %>%
      left_join(
        data[[t]]$females_all %>% dplyr::select(ensembl_gene_id,
                                         lfc_female = log2FoldChange, padj_female = padj),
        by = "ensembl_gene_id"
      ) %>%
      mutate(
        tissue = t,
        interpretation = case_when(
          lfc_interaction > 0 & lfc_male > 0  ~ "Stronger UP in males",
          lfc_interaction > 0 & lfc_male < 0  ~ "UP in males, DOWN in females",
          lfc_interaction < 0 & lfc_female > 0 ~ "Stronger UP in females",
          lfc_interaction < 0 & lfc_female < 0 ~ "UP in females, DOWN in males",
          TRUE ~ "Complex"
        )
      ) %>%
      arrange(padj_interaction)
    
    write.csv(int_context,
              file.path(OUT_DIR, sprintf("%s_interaction_genes_context.csv", t)),
              row.names = FALSE)
    message(sprintf("  ✅ %s: %d interaction genes with context", t, nrow(int_context)))
    
    message(sprintf("\n  Top interaction genes in %s:", t))
    print(int_context %>% 
            dplyr::select(external_gene_name, lfc_male, lfc_female, 
                   lfc_interaction, padj_interaction, interpretation) %>%
            head(10),
          row.names = FALSE)
  }
}

# ==============================================================================
# RESUMEN FINAL
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("PLOTS GENERADOS en results/comparative/")
message("═══════════════════════════════════════════\n")

message("Por tejido (CRB y TC):")
message("  *_lfc_scatter_sex.png      ← LFC males vs females (MÁS INFORMATIVO)")
message("  *_volcano_combined.png     ← Volcanos lado a lado")
message("  *_upset.png                ← Intersecciones entre grupos")
message("  *_discordant_heatmap.png   ← Genes con efecto opuesto por sexo")
message("")
message("Global:")
message("  summary_bars.png           ← Resumen comparativo CRB vs TC")
message("")
message("Tablas:")
message("  *_discordant_genes.csv     ← Genes con efecto AD opuesto por sexo")
message("  *_interaction_genes_context.csv ← Genes de interacción con LFC por sexo")