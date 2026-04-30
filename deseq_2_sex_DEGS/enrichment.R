#!/usr/bin/env Rscript
# ==============================================================================
# GSEA y ORA: Análisis funcional por categoría sex-specific en AD
# ==============================================================================
# Estrategia por categoría:
#   Male-only, Female-only, Concordant → GSEA (ranked por LFC × -log10 padj)
#   Discordant, Interaction            → ORA  (lista fija, n < 100)
#
# Bases de datos: GO:MF + KEGG
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ==============================================================================
# CONFIGURACIÓN
# ==============================================================================

TISSUES     <- c("CRB", "TC")
RESULTS_DIR <- "results"
OUT_DIR     <- "results/enrichment"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PADJ_THRESH <- 0.05
LFC_THRESH  <- 0.5

# ==============================================================================
# HELPERS
# ==============================================================================

#' Convertir ENSEMBL → ENTREZ (requerido por KEGG)
ensembl_to_entrez <- function(ensembl_ids) {
  mapIds(org.Hs.eg.db,
         keys      = ensembl_ids,
         column    = "ENTREZID",
         keytype   = "ENSEMBL",
         multiVals = "first")
}

#' Construir ranked list para GSEA
#' Métrica: sign(LFC) × -log10(padj)
#' Captura dirección + significancia
make_ranked_list <- function(all_genes_df) {
  all_genes_df %>%
    filter(!is.na(padj), padj > 0) %>%
    mutate(
      entrez = ensembl_to_entrez(ensembl_gene_id),
      rank   = sign(log2FoldChange) * -log10(padj)
    ) %>%
    filter(!is.na(entrez)) %>%
    # Si hay duplicados de entrez, quedarse con el de mayor |rank|
    group_by(entrez) %>%
    slice_max(abs(rank), n = 1) %>%
    ungroup() %>%
    arrange(desc(rank)) %>%
    { setNames(.$rank, .$entrez) }
}

#' Convertir lista de ENSEMBL a ENTREZ limpia
ensembl_list_to_entrez <- function(ensembl_ids) {
  entrez <- ensembl_to_entrez(ensembl_ids)
  entrez[!is.na(entrez)]
}

#' Guardar plot con manejo de error
save_plot <- function(plot, path, width = 10, height = 7) {
  tryCatch({
    ggsave(path, plot, width = width, height = height, dpi = 300)
    message(sprintf("    ✅ %s", basename(path)))
  }, error = function(e) {
    message(sprintf("    ⚠️  No se pudo guardar %s: %s", basename(path), e$message))
  })
}

# ==============================================================================
# FUNCIÓN PRINCIPAL: GSEA (para categorías con muchos genes)
# ==============================================================================

run_gsea <- function(ranked_list, label, tissue, db = c("GO_MF", "KEGG")) {
  
  results <- list()
  
  # ── GO Molecular Function ──────────────────────────────────────────────────
  if ("GO_MF" %in% db) {
    message(sprintf("    Running GSEA GO:MF for %s...", label))
    
    gsea_go <- tryCatch(
      gseGO(
        geneList     = ranked_list,
        OrgDb        = org.Hs.eg.db,
        ont          = "MF",
        keyType      = "ENTREZID",
        minGSSize    = 15,
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose      = FALSE,
        eps          = 0
      ),
      error = function(e) { message("      ⚠️  GO:MF falló: ", e$message); NULL }
    )
    
    if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
      results$go_mf <- gsea_go
      
      # Dotplot
      p <- dotplot(gsea_go, showCategory = 20, split = ".sign") +
        facet_grid(~.sign) +
        labs(title = sprintf("%s %s — GSEA GO:MF", tissue, label)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot(p,
                file.path(OUT_DIR, sprintf("%s_%s_gsea_gomf_dotplot.png", tissue, label)),
                width = 14, height = 8)
      
      # Ridge plot (distribución de rankings)
      p2 <- ridgeplot(gsea_go, showCategory = 20, fill = "pvalue") +
        labs(title = sprintf("%s %s — GSEA GO:MF ranking distribution", tissue, label)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot(p2,
                file.path(OUT_DIR, sprintf("%s_%s_gsea_gomf_ridgeplot.png", tissue, label)),
                width = 12, height = 9)
      
      # Guardar tabla
      write.csv(as.data.frame(gsea_go),
                file.path(OUT_DIR, sprintf("%s_%s_gsea_gomf.csv", tissue, label)),
                row.names = FALSE)
      
      message(sprintf("      %d significant GO:MF terms", nrow(gsea_go@result)))
    } else {
      message("      No significant GO:MF terms found")
    }
  }
  
  # ── KEGG ──────────────────────────────────────────────────────────────────
  if ("KEGG" %in% db) {
    message(sprintf("    Running GSEA KEGG for %s...", label))
    
    gsea_kegg <- tryCatch(
      gseKEGG(
        geneList      = ranked_list,
        organism      = "hsa",
        minGSSize     = 15,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        pAdjustMethod = "BH",
        verbose       = FALSE,
        eps           = 0
      ),
      error = function(e) { message("      ⚠️  KEGG falló: ", e$message); NULL }
    )
    
    if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
      results$kegg <- gsea_kegg
      
      # Dotplot
      p <- dotplot(gsea_kegg, showCategory = 20, split = ".sign") +
        facet_grid(~.sign) +
        labs(title = sprintf("%s %s — GSEA KEGG", tissue, label)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot(p,
                file.path(OUT_DIR, sprintf("%s_%s_gsea_kegg_dotplot.png", tissue, label)),
                width = 14, height = 8)
      
      # GSEA enrichment plot para top 3 pathways
      top_paths <- head(gsea_kegg@result$ID[order(gsea_kegg@result$p.adjust)], 3)
      for (pid in top_paths) {
        p3 <- gseaplot2(gsea_kegg, geneSetID = pid,
                        title = gsea_kegg@result$Description[gsea_kegg@result$ID == pid])
        save_plot(p3,
                  file.path(OUT_DIR, sprintf("%s_%s_kegg_%s_gseaplot.png", tissue, label, pid)),
                  width = 10, height = 6)
      }
      
      write.csv(as.data.frame(gsea_kegg),
                file.path(OUT_DIR, sprintf("%s_%s_gsea_kegg.csv", tissue, label)),
                row.names = FALSE)
      
      message(sprintf("      %d significant KEGG pathways", nrow(gsea_kegg@result)))
    } else {
      message("      No significant KEGG pathways found")
    }
  }
  
  return(results)
}

# ==============================================================================
# FUNCIÓN PRINCIPAL: ORA (para categorías con pocos genes)
# ==============================================================================

run_ora <- function(gene_entrez, background_entrez, label, tissue,
                    db = c("GO_MF", "KEGG")) {
  
  results <- list()
  
  if (length(gene_entrez) < 5) {
    message(sprintf("    ⚠️  %s: muy pocos genes (%d) para ORA, saltando",
                    label, length(gene_entrez)))
    return(invisible(NULL))
  }
  
  # ── GO Molecular Function ──────────────────────────────────────────────────
  if ("GO_MF" %in% db) {
    message(sprintf("    Running ORA GO:MF for %s (%d genes)...",
                    label, length(gene_entrez)))
    
    ora_go <- tryCatch(
      enrichGO(
        gene          = gene_entrez,
        universe      = background_entrez,
        OrgDb         = org.Hs.eg.db,
        ont           = "MF",
        keyType       = "ENTREZID",
        minGSSize     = 10,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        pAdjustMethod = "BH",
        readable      = TRUE
      ),
      error = function(e) { message("      ⚠️  GO:MF ORA falló: ", e$message); NULL }
    )
    
    if (!is.null(ora_go) && nrow(ora_go@result) > 0) {
      results$go_mf <- ora_go
      
      p <- dotplot(ora_go, showCategory = 20) +
        labs(title = sprintf("%s %s — ORA GO:MF\n(%d genes)", tissue, label, length(gene_entrez))) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot(p,
                file.path(OUT_DIR, sprintf("%s_%s_ora_gomf_dotplot.png", tissue, label)),
                width = 10, height = 8)
      
      write.csv(as.data.frame(ora_go),
                file.path(OUT_DIR, sprintf("%s_%s_ora_gomf.csv", tissue, label)),
                row.names = FALSE)
      
      message(sprintf("      %d significant GO:MF terms", nrow(ora_go@result)))
    } else {
      message("      No significant GO:MF terms found")
    }
  }
  
  # ── KEGG ──────────────────────────────────────────────────────────────────
  if ("KEGG" %in% db) {
    message(sprintf("    Running ORA KEGG for %s (%d genes)...",
                    label, length(gene_entrez)))
    
    ora_kegg <- tryCatch(
      enrichKEGG(
        gene          = gene_entrez,
        universe      = background_entrez,
        organism      = "hsa",
        minGSSize     = 10,
        maxGSSize     = 500,
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        pAdjustMethod = "BH"
      ),
      error = function(e) { message("      ⚠️  KEGG ORA falló: ", e$message); NULL }
    )
    
    if (!is.null(ora_kegg) && nrow(ora_kegg@result) > 0) {
      results$kegg <- ora_kegg
      
      p <- dotplot(ora_kegg, showCategory = 20) +
        labs(title = sprintf("%s %s — ORA KEGG\n(%d genes)", tissue, label, length(gene_entrez))) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      save_plot(p,
                file.path(OUT_DIR, sprintf("%s_%s_ora_kegg_dotplot.png", tissue, label)),
                width = 10, height = 8)
      
      write.csv(as.data.frame(ora_kegg),
                file.path(OUT_DIR, sprintf("%s_%s_ora_kegg.csv", tissue, label)),
                row.names = FALSE)
      
      message(sprintf("      %d significant KEGG pathways", nrow(ora_kegg@result)))
    } else {
      message("      No significant KEGG pathways found")
    }
  }
  
  return(results)
}

# ==============================================================================
# PLOT COMPARATIVO: dotplot de múltiples grupos en un solo panel
# ==============================================================================

plot_comparative_dotplot <- function(enrichment_list, tissue, db_label) {
  # enrichment_list: lista nombrada de objetos enrichResult/gseaResult
  
  # Extraer top 10 terms por grupo
  combined <- bind_rows(lapply(names(enrichment_list), function(grp) {
    obj <- enrichment_list[[grp]]
    if (is.null(obj)) return(NULL)
    
    df <- as.data.frame(obj) %>%
      arrange(p.adjust) %>%
      head(10) %>%
      mutate(
        group     = grp,
        log_padj  = -log10(p.adjust),
        # GeneRatio como numérico
        gene_ratio = if ("GeneRatio" %in% colnames(.)) {
          sapply(GeneRatio, function(x) eval(parse(text = x)))
        } else {
          NES  # para GSEA
        }
      )
    df
  }))
  
  if (is.null(combined) || nrow(combined) == 0) {
    message(sprintf("  ⚠️  No hay términos para el dotplot comparativo %s", db_label))
    return(invisible(NULL))
  }
  
  # Top 20 términos más frecuentes entre grupos
  top_terms <- combined %>%
    count(Description) %>%
    slice_max(n, n = 25) %>%
    pull(Description)
  
  plot_df <- combined %>%
    filter(Description %in% top_terms) %>%
    mutate(
      group = factor(group, levels = names(enrichment_list)),
      Description = stringr::str_wrap(Description, width = 40)
    )
  
  p <- ggplot(plot_df, aes(x = group, y = Description, size = gene_ratio,
                           color = log_padj)) +
    geom_point() +
    scale_color_gradient(low = "#74C0FC", high = "#E03131",
                         name = "-log10(padj)") +
    scale_size_continuous(name = "Gene ratio", range = c(2, 8)) +
    labs(
      title = sprintf("%s — %s enrichment by category", tissue, db_label),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold"),
      axis.text.y  = element_text(size = 9),
      plot.title   = element_text(face = "bold", hjust = 0.5),
      panel.grid.major = element_line(color = "grey92")
    )
  
  save_plot(p,
            file.path(OUT_DIR, sprintf("%s_comparative_%s_dotplot.png", tissue, db_label)),
            width = 14, height = max(8, length(top_terms) * 0.35 + 3))
}

# ==============================================================================
# MAIN: CORRER TODO POR TEJIDO
# ==============================================================================

message("═══════════════════════════════════════════")
message("Análisis Funcional Sex-Specific en AD")
message("═══════════════════════════════════════════\n")

for (tissue in TISSUES) {
  
  message(sprintf("\n╔══════════════════════════════╗"))
  message(sprintf("║  %s                           ", tissue))
  message(sprintf("╚══════════════════════════════╝\n"))
  
  # Cargar datos
  males_all   <- read.csv(file.path(RESULTS_DIR, tissue,
                                    sprintf("%s_males_AD_vs_control_all_genes.csv", tissue)))
  females_all <- read.csv(file.path(RESULTS_DIR, tissue,
                                    sprintf("%s_females_AD_vs_control_all_genes.csv", tissue)))
  interaction_degs <- read.csv(file.path(RESULTS_DIR, tissue,
                                         sprintf("%s_sex_diagnosis_interaction_DEGs.csv", tissue)))
  
  # Background: todos los genes testeados (unión de males y females)
  background_ensembl <- union(males_all$ensembl_gene_id, females_all$ensembl_gene_id)
  background_entrez  <- ensembl_list_to_entrez(background_ensembl)
  message(sprintf("  Background: %d genes → %d con ENTREZ",
                  length(background_ensembl), length(background_entrez)))
  
  # Categorías desde cats (generado en sex_comparative_analysis.R)
  male_only_ids   <- cats[[tissue]]$male_only
  female_only_ids <- cats[[tissue]]$female_only
  concordant_ids  <- cats[[tissue]]$concordant
  discordant_ids  <- cats[[tissue]]$discordant
  interaction_ids <- interaction_degs$ensembl_gene_id
  
  # ── GSEA: Male-only ─────────────────────────────────────────────────────────
  message("\n  [GSEA] Male-only genes")
  # Ranked list solo con genes que son DEG en males
  male_ranked <- make_ranked_list(
    males_all %>% filter(ensembl_gene_id %in% male_only_ids)
  )
  message(sprintf("    Ranked list: %d genes", length(male_ranked)))
  gsea_male <- run_gsea(male_ranked, "male_only", tissue)
  
  # ── GSEA: Female-only ───────────────────────────────────────────────────────
  message("\n  [GSEA] Female-only genes")
  female_ranked <- make_ranked_list(
    females_all %>% filter(ensembl_gene_id %in% female_only_ids)
  )
  message(sprintf("    Ranked list: %d genes", length(female_ranked)))
  gsea_female <- run_gsea(female_ranked, "female_only", tissue)
  
  # ── GSEA: Concordant (shared, same direction) ────────────────────────────────
  message("\n  [GSEA] Concordant genes (shared AD effect)")
  # Usar el promedio del LFC de males y females como métrica de ranking
  conc_males   <- males_all %>%
    filter(ensembl_gene_id %in% concordant_ids) %>%
    select(ensembl_gene_id, lfc_m = log2FoldChange, padj_m = padj)
  conc_females <- females_all %>%
    filter(ensembl_gene_id %in% concordant_ids) %>%
    select(ensembl_gene_id, lfc_f = log2FoldChange, padj_f = padj)
  
  concordant_df <- inner_join(conc_males, conc_females, by = "ensembl_gene_id") %>%
    mutate(
      # Métrica: promedio de los dos ranks
      mean_lfc  = (lfc_m + lfc_f) / 2,
      mean_padj = (padj_m + padj_f) / 2,
      entrez    = ensembl_to_entrez(ensembl_gene_id),
      rank      = sign(mean_lfc) * -log10(mean_padj)
    ) %>%
    filter(!is.na(entrez)) %>%
    group_by(entrez) %>%
    slice_max(abs(rank), n = 1) %>%
    ungroup() %>%
    arrange(desc(rank))
  
  concordant_ranked <- setNames(concordant_df$rank, concordant_df$entrez)
  message(sprintf("    Ranked list: %d genes", length(concordant_ranked)))
  gsea_concordant <- run_gsea(concordant_ranked, "concordant", tissue)
  
  # ── ORA: Discordant ──────────────────────────────────────────────────────────
  message("\n  [ORA] Discordant genes (opposite AD effect by sex)")
  disc_entrez <- ensembl_list_to_entrez(discordant_ids)
  ora_discordant <- run_ora(disc_entrez, background_entrez, "discordant", tissue)
  
  # ── ORA: Interaction (solo CRB tiene suficientes) ───────────────────────────
  if (length(interaction_ids) >= 5) {
    message("\n  [ORA] Interaction DEGs")
    int_entrez <- ensembl_list_to_entrez(interaction_ids)
    ora_interaction <- run_ora(int_entrez, background_entrez, "interaction", tissue)
  } else {
    message(sprintf("\n  [ORA] Interaction: solo %d genes, saltando",
                    length(interaction_ids)))
    ora_interaction <- NULL
  }
  
  # ── PLOT COMPARATIVO: todos los grupos en un panel ───────────────────────────
  message("\n  Generando plots comparativos...")
  
  # GO:MF comparativo
  go_list <- list(
    male_only   = if (!is.null(gsea_male$go_mf))       gsea_male$go_mf       else NULL,
    female_only = if (!is.null(gsea_female$go_mf))     gsea_female$go_mf     else NULL,
    concordant  = if (!is.null(gsea_concordant$go_mf)) gsea_concordant$go_mf else NULL,
    discordant  = if (!is.null(ora_discordant$go_mf))  ora_discordant$go_mf  else NULL
  )
  go_list <- Filter(Negate(is.null), go_list)
  if (length(go_list) >= 2) plot_comparative_dotplot(go_list, tissue, "GO_MF")
  
  # KEGG comparativo
  kegg_list <- list(
    male_only   = if (!is.null(gsea_male$kegg))       gsea_male$kegg       else NULL,
    female_only = if (!is.null(gsea_female$kegg))     gsea_female$kegg     else NULL,
    concordant  = if (!is.null(gsea_concordant$kegg)) gsea_concordant$kegg else NULL,
    discordant  = if (!is.null(ora_discordant$kegg))  ora_discordant$kegg  else NULL
  )
  kegg_list <- Filter(Negate(is.null), kegg_list)
  if (length(kegg_list) >= 2) plot_comparative_dotplot(kegg_list, tissue, "KEGG")
}

# ==============================================================================
# RESUMEN DE OUTPUTS
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("OUTPUTS en results/enrichment/")
message("═══════════════════════════════════════════\n")

all_outputs <- list.files("results/enrichment", full.names = FALSE)
csvs  <- all_outputs[grepl("\\.csv$", all_outputs)]
plots <- all_outputs[grepl("\\.png$", all_outputs)]

message(sprintf("Tablas CSV:  %d archivos", length(csvs)))
message(sprintf("Plots PNG:   %d archivos\n", length(plots)))

message("Plots por tipo:")
message("  *_gsea_gomf_dotplot.png    ← GSEA GO:MF activado/reprimido")
message("  *_gsea_gomf_ridgeplot.png  ← distribución de rankings")
message("  *_gsea_kegg_dotplot.png    ← GSEA KEGG activado/reprimido")
message("  *_kegg_*_gseaplot.png      ← curva de enriquecimiento top pathways")
message("  *_ora_*_dotplot.png        ← ORA para discordant/interaction")
message("  *_comparative_*_dotplot.png ← TODOS los grupos en un panel ← MÁS ÚTIL")