#!/usr/bin/env Rscript
# ==============================================================================
# Preprocesar datos Mayo RNA-seq para el pipeline DESeq2
# ==============================================================================
# Problemas que resuelve:
#   1. Counts tiene columna 'feature' en lugar de rownames
#   2. Primeras 4 filas son stats de STAR (N_unmapped, etc.) - no genes
#   3. IDs de genes tienen versión (ENSG00000223972.5 → ENSG00000223972)
#   4. Metadata necesita columnas estándar
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(vroom)
})

# ==============================================================================
# CONFIGURACIÓN - EDITA AQUÍ
# ==============================================================================

tissues <- list(
  CRB = list(
    counts   = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_rawcounts_CRB.rds",
    metadata = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_metadata_CRB.txt"
  ),
  TC = list(
    counts   = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_rawcounts_TC.rds",
    metadata = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_metadata_TC.txt"
  )
)

# Filas de STAR stats que hay que eliminar
STAR_STATS_ROWS <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

# ==============================================================================
# FUNCIÓN PRINCIPAL
# ==============================================================================

preprocess_tissue <- function(tissue_name, counts_path, metadata_path) {
  
  message(sprintf("\n─── Procesando: %s ───", tissue_name))
  
  # ==================
  # COUNTS
  # ==================
  message("  Cargando counts...")
  raw <- readRDS(counts_path)
  message(sprintf("    Dimensiones originales: %d rows × %d cols", nrow(raw), ncol(raw)))
  
  # Identificar columna de genes
  # Puede llamarse 'feature', 'gene_id', 'GeneID', etc.
  gene_col <- NULL
  for (col in c("feature", "gene_id", "GeneID", "Gene", "gene")) {
    if (col %in% colnames(raw)) {
      gene_col <- col
      break
    }
  }
  
  if (is.null(gene_col)) {
    # Si no hay columna obvia, asumir que la primera es genes
    if (!is.numeric(raw[, 1])) {
      gene_col <- colnames(raw)[1]
      message(sprintf("    Usando primera columna como genes: '%s'", gene_col))
    } else {
      stop("No se encontró columna de IDs de genes")
    }
  }
  
  message(sprintf("    Columna de genes: '%s'", gene_col))
  
  # Mover columna de genes a rownames
  counts <- raw %>%
    column_to_rownames(gene_col) %>%
    as.matrix()
  
  # Eliminar filas de STAR stats
  star_rows <- rownames(counts)[rownames(counts) %in% STAR_STATS_ROWS]
  if (length(star_rows) > 0) {
    message(sprintf("    Eliminando %d filas de STAR stats: %s", 
                    length(star_rows), paste(star_rows, collapse = ", ")))
    counts <- counts[!rownames(counts) %in% STAR_STATS_ROWS, ]
  }
  
  # Eliminar versión de ENSEMBL IDs (ENSG00000223972.5 → ENSG00000223972)
  if (any(grepl("\\.", rownames(counts)))) {
    message("    Eliminando versiones de Ensembl IDs (ENSG.X → ENSG)...")
    rownames(counts) <- sub("\\..*", "", rownames(counts))
  }
  
  # Verificar que sea numérico
  if (!is.numeric(counts)) {
    counts <- apply(counts, 2, as.numeric)
    rownames(counts) <- sub("\\..*", "", rownames(raw %>% column_to_rownames(gene_col) %>% rownames()))
  }
  
  # Remover genes con todos ceros
  keep <- rowSums(counts) > 0
  message(sprintf("    Genes con counts > 0: %d / %d", sum(keep), nrow(counts)))
  counts <- counts[keep, ]
  
  # Remover rownames duplicados (suma sus counts)
  if (any(duplicated(rownames(counts)))) {
    n_dup <- sum(duplicated(rownames(counts)))
    message(sprintf("    Colapsando %d rownames duplicados (sumando counts)...", n_dup))
    counts <- counts %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      group_by(gene) %>%
      summarise(across(everything(), sum)) %>%
      column_to_rownames("gene") %>%
      as.matrix()
  }
  
  message(sprintf("    Dimensiones finales: %d genes × %d samples", nrow(counts), ncol(counts)))
  
  # ==================
  # METADATA
  # ==================
  message("  Cargando metadata...")
  metadata <- vroom(metadata_path, show_col_types = FALSE) %>% as.data.frame()
  message(sprintf("    Cargada: %d samples × %d columnas", nrow(metadata), ncol(metadata)))
  
  # Mostrar columnas disponibles para diagnóstico
  message("    Columnas disponibles:")
  message(sprintf("    %s", paste(colnames(metadata), collapse = ", ")))
  
  # Identificar columna de diagnosis
  diag_col <- NULL
  for (col in c("diagnosis", "Diagnosis", "disease", "Disease", "condition", 
                "Condition", "diagnosis_label", "DX")) {
    if (col %in% colnames(metadata)) {
      diag_col <- col
      break
    }
  }
  
  # Identificar columna de sexo
  sex_col <- NULL
  for (col in c("sex", "Sex", "GENDER", "gender", "Gender", "msex")) {
    if (col %in% colnames(metadata)) {
      sex_col <- col
      break
    }
  }
  
  if (is.null(diag_col)) {
    message("    ⚠️  No se encontró columna de diagnóstico automáticamente")
    message("    Columnas disponibles: ", paste(colnames(metadata), collapse = ", "))
    message("    ACCIÓN REQUERIDA: Edita DIAGNOSIS_COL en este script")
  } else {
    message(sprintf("    Columna de diagnóstico: '%s'", diag_col))
    message(sprintf("    Valores únicos: %s", 
                    paste(unique(metadata[[diag_col]]), collapse = ", ")))
  }
  
  if (is.null(sex_col)) {
    message("    ⚠️  No se encontró columna de sexo automáticamente")
  } else {
    message(sprintf("    Columna de sexo: '%s'", sex_col))
    message(sprintf("    Valores únicos: %s", 
                    paste(unique(metadata[[sex_col]]), collapse = ", ")))
  }
  
  # Estandarizar columnas
  std_meta <- data.frame(specimenID = metadata$specimenID)
  
  if (!is.null(sex_col)) {
    std_meta$sex <- tolower(as.character(metadata[[sex_col]]))
    std_meta$sex[std_meta$sex %in% c("m", "male", "1")] <- "male"
    std_meta$sex[std_meta$sex %in% c("f", "female", "0", "2")] <- "female"
  }
  
  if (!is.null(diag_col)) {
    std_meta$diagnosis <- tolower(as.character(metadata[[diag_col]]))
    std_meta$diagnosis[std_meta$diagnosis %in% c("control", "ctrl", "ct", "normal", "0")] <- "control"
    std_meta$diagnosis[std_meta$diagnosis %in% c("ad", "alzheimer", "alzheimers", 
                                                 "alzheimer disease", "1", "case")] <- "AD"
    # Excluir diagnósticos no relevantes para la comparación AD vs control
    excluded <- unique(std_meta$diagnosis[!std_meta$diagnosis %in% c("AD", "control")])
    if (length(excluded) > 0) {
      n_excluded <- sum(std_meta$diagnosis %in% excluded)
      message(sprintf("    Excluyendo %d muestras con diagnóstico: %s",
                      n_excluded, paste(excluded, collapse = ", ")))
      std_meta <- std_meta %>% filter(diagnosis %in% c("AD", "control"))
    }
  }
  
  # Agregar covariables via join (respeta el filtro previo de diagnosis)
  extra_cols <- intersect(
    c("pmi", "RIN", "Braak", "thal", "apoeGenotype", "platform",
      "libraryPrep", "individualID", "specimenIdSource"),
    colnames(metadata)
  )
  if (length(extra_cols) > 0) {
    extra <- metadata %>% select(specimenID, all_of(extra_cols))
    std_meta <- std_meta %>% left_join(extra, by = "specimenID")
  }
  
  # Imputar NAs en pmi con la mediana
  if ("pmi" %in% colnames(std_meta)) {
    n_na_pmi <- sum(is.na(std_meta$pmi))
    if (n_na_pmi > 0) {
      pmi_median <- median(std_meta$pmi, na.rm = TRUE)
      message(sprintf("    Imputando %d NAs en pmi con mediana (%.2f)", n_na_pmi, pmi_median))
      std_meta$pmi[is.na(std_meta$pmi)] <- pmi_median
    }
  }
  
  # specimenIdSource como variable de batch
  if (!"specimenIdSource" %in% colnames(std_meta)) {
    if ("platform" %in% colnames(std_meta)) {
      std_meta$specimenIdSource <- std_meta$platform
      message("    Usando 'platform' como specimenIdSource (batch)")
    } else if ("libraryPrep" %in% colnames(std_meta)) {
      std_meta$specimenIdSource <- std_meta$libraryPrep
      message("    Usando 'libraryPrep' como specimenIdSource (batch)")
    } else {
      std_meta$specimenIdSource <- "batch1"
      message("    Sin batch info, usando 'batch1' para todos")
    }
  }
  
  # ==================
  # ALINEAR MUESTRAS
  # ==================
  message("  Alineando muestras entre counts y metadata...")
  
  common <- intersect(colnames(counts), std_meta$specimenID)
  message(sprintf("    Muestras en counts: %d", ncol(counts)))
  message(sprintf("    Muestras en metadata: %d", nrow(std_meta)))
  message(sprintf("    Muestras comunes: %d", length(common)))
  
  if (length(common) == 0) {
    stop("¡No hay muestras en común! Verifica los IDs de muestras")
  }
  
  counts <- counts[, common]
  std_meta <- std_meta %>% filter(specimenID %in% common)
  
  # Eliminar muestras con NAs en cualquier variable del diseño
  design_vars <- c("sex", "diagnosis", "specimenIdSource", "pmi", "RIN")
  design_vars <- intersect(design_vars, colnames(std_meta))
  
  for (var in design_vars) {
    n_na <- sum(is.na(std_meta[[var]]))
    if (n_na > 0) {
      message(sprintf("    ⚠️  %d NAs en '%s' - eliminando esas muestras", n_na, var))
      std_meta <- std_meta %>% filter(!is.na(.data[[var]]))
      counts   <- counts[, std_meta$specimenID, drop = FALSE]
    }
  }
  
  # Mostrar distribución final
  if (!is.null(sex_col) && !is.null(diag_col)) {
    message("    Distribución final:")
    print(table(std_meta$sex, std_meta$diagnosis))
  }
  
  # ==================
  # GUARDAR
  # ==================
  dir.create("data/counts", recursive = TRUE, showWarnings = FALSE)
  dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
  
  counts_out <- sprintf("data/counts/%s_counts.rds", tissue_name)
  meta_out <- sprintf("data/metadata/%s_metadata.csv", tissue_name)
  
  saveRDS(counts, counts_out)
  write.csv(std_meta, meta_out, row.names = FALSE)
  
  message(sprintf("  ✅ Guardado: %s", counts_out))
  message(sprintf("  ✅ Guardado: %s", meta_out))
  
  return(list(counts = counts, metadata = std_meta))
}

# ==============================================================================
# CORRER PARA CADA TEJIDO
# ==============================================================================

message("═══════════════════════════════════════════")
message("Preprocesando datos Mayo RNA-seq")
message("═══════════════════════════════════════════")

results <- list()

for (tissue_name in names(tissues)) {
  tryCatch({
    results[[tissue_name]] <- preprocess_tissue(
      tissue_name   = tissue_name,
      counts_path   = tissues[[tissue_name]]$counts,
      metadata_path = tissues[[tissue_name]]$metadata
    )
  }, error = function(e) {
    message(sprintf("  ❌ Error en %s: %s", tissue_name, e$message))
  })
}

# ==============================================================================
# RESUMEN FINAL
# ==============================================================================

message("\n═══════════════════════════════════════════")
message("RESUMEN")
message("═══════════════════════════════════════════")

for (tissue_name in names(results)) {
  res <- results[[tissue_name]]
  message(sprintf("\n%s:", tissue_name))
  message(sprintf("  Genes:   %d", nrow(res$counts)))
  message(sprintf("  Samples: %d", ncol(res$counts)))
  if (!is.null(res$metadata$sex) && !is.null(res$metadata$diagnosis)) {
    message("  Grupos:")
    print(table(res$metadata$sex, res$metadata$diagnosis))
  }
}

message("\n✅ Preprocesamiento completo!")
message("Siguiente paso:")
message("  Rscript deseq2_sex_ad_analysis.R --counts data/counts/cerebellum_counts.rds --metadata data/metadata/cerebellum_metadata.csv --outdir results/cerebellum --brain_area cerebellum --cores 4")