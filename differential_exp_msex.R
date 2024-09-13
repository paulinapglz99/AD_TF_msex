#differential_exp_msex.R
#This script performs Differential expression analyis for 3 types of comparisons
#AD women vs control women
#AD men vs control men
#AD men vs AD women

# Load necessary libraries
pacman::p_load("tidyverse", 
               "DESeq2", 
               "pheatmap", 
               "ggplot2", 
               "ggrepel", 
               "stringr",
               "biomaRt", 
               "RColorBrewer", 
               "gridExtra",
               "clusterProfiler",
               "org.Hs.eg.db",
               "fgsea", 
               "enrichplot", 
               "igraph")

# Set seed
set.seed(10)

# Define function to convert from ENSMBL to SYMBOL
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

convert_ens_to_symbol <- function(ensembl_ids) {
  trad <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_gene_id",
                values = ensembl_ids,
                mart = ensembl)
  trad$external_gene_name <- ifelse(trad$external_gene_name == "", trad$ensembl_gene_id, trad$external_gene_name)
  return(trad)
}

# Function to delete duplicates
del_dupl <- function(counts) {
    repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)
  
  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  repeated_rows <- repeated_rows[order(repeated_rows$feature),]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), median, na.rm = TRUE))
  
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#Load data

load_data <- function() {
  counts <- readRDS(file = "~/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/data/ROSMAP_RNAseq_rawcounts_DLPFC.rds") %>% as.data.frame()
  metadata <- vroom::vroom(file = "~/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/metadata/RNA_seq_metadata_DLPFC.txt")
  
  return(list(counts = counts, metadata = metadata))
}

#Pre-process data

preprocess_data <- function(counts, metadata, condition_column) {
  counts <- counts[-c(1:4),]  # Remove the first four rows
  counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))
  counts <- del_dupl(counts)
  rownames(counts) <- NULL
  counts <- counts %>% column_to_rownames(var = "feature")
  metadata <- metadata %>% filter(!is.na(!!sym(condition_column)))
  counts <- counts %>% dplyr::select(all_of(metadata$specimenID))
  
  return(list(counts = counts, metadata = metadata))
}

#Run deseq

run_deseq <- function(counts, metadata, condition_column) {
  coldata <- as.data.frame(colnames(counts))
  colnames(coldata) <- "specimenID"
  
  coldata <- coldata %>%
    left_join(metadata, by = "specimenID") %>%
    dplyr::select("specimenID", "sequencingBatch","msex","cogdx", "ceradsc", condition_column)
  
  coldata[[condition_column]] <- as.factor(coldata[[condition_column]])  # Convert to factor
  
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = coldata,
                                design = as.formula(paste("~ sequencingBatch +", condition_column)))  
  dds <- DESeq(dds)
  
  return(list(dds = dds, 
              coldata = coldata))
}

#Extract DEGS

extract_DEGS <- function(dds, condition_column, coldata){
  
  res <- results(dds, contrast = c(condition_column, "1", "4"), #condition, problem, control 
                 pAdjustMethod = 'BH', alpha = 0.05)
  res.df <- res %>% as.data.frame() %>% filter(!is.na(padj))
  
  # add a column of NAs
  res.df$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  res.df$diffexpressed[res.df$log2FoldChange < -0.5 & res.df$padj < 0.05] <- "DOWN"
  
  #Add gene names in SYMBOL --- ---
  
  #Create dictionary
  symbol <- convert_ens_to_symbol(rownames(res))
  
  # Make rownames columns
  res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
  res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge
  
  # BaseMean by condition --- ---
  
  #Extract counts normalized 
  
  normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()
  
  # Calcula los valores de baseMean por condición
  baseMean_4 <- rowMeans(normcounts[, coldata$ceradsc == "4"]) %>% as.data.frame()
  colnames(baseMean_4) <- "baseMean_ceradsc_4"
  
  baseMean_1 <- rowMeans(normcounts[, coldata$ceradsc == "1"])%>% as.data.frame()
  colnames(baseMean_1) <- "baseMean_ceradsc_1"
  
  baseMean_genes <- cbind(baseMean_4, baseMean_1)
  baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")
  
  #DEGS are...
  DEGS <- res.df %>% filter(diffexpressed != "NO")
  dim(DEGS)
  
  baseMean_DEGS <- baseMean_genes %>% filter(feature %in% DEGS$ensembl_gene_id)
  DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))
  
  
  return(list(res.df = res.df, DEGS = DEGS, normcounts = normcounts))
  
}


extract_DEGS <- function(dds, condition_column, coldata){
  
  res <- results(dds, contrast = c(condition_column, "1", "4"), #condition, problem, control 
                 pAdjustMethod = 'BH', alpha = 0.05)
  res.df <- res %>% as.data.frame() %>% filter(!is.na(padj))
  
  # add a column of NAs
  res.df$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  res.df$diffexpressed[res.df$log2FoldChange < -0.5 & res.df$padj < 0.05] <- "DOWN"
  
  # Add gene names in SYMBOL --- ---
  # Create dictionary
  symbol <- convert_ens_to_symbol(rownames(res))
  
  # Make rownames columns
  res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
  res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge
  
  # BaseMean by condition --- ---
  # Extract counts normalized 
  normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()
  
  # Calcula los valores de baseMean por condición
  baseMean_4 <- rowMeans(normcounts[, coldata$ceradsc == "4"]) %>% as.data.frame()
  colnames(baseMean_4) <- "baseMean_ceradsc_4"
  
  baseMean_1 <- rowMeans(normcounts[, coldata$ceradsc == "1"]) %>% as.data.frame()
  colnames(baseMean_1) <- "baseMean_ceradsc_1"
  
  baseMean_genes <- cbind(baseMean_4, baseMean_1)
  baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")
  
  # DEGs are...
  DEGS <- res.df %>% filter(diffexpressed != "NO")
  dim(DEGS)
  
  baseMean_DEGS <- baseMean_genes %>% filter(feature %in% DEGS$ensembl_gene_id)
  DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))
  
  # Add delabel column with gene names for DEGs
  DEGS <- DEGS %>% mutate(delabel = external_gene_name)
  
  # For non-DEGs, delabel will be NA or empty
  res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA))
  
  return(list(res.df = res.df, DEGS = DEGS, normcounts = normcounts))
  
}


#Visualize results

create_volcano_plot <- function(res.df, condition_compared ,output_file) {

  vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj),  col = diffexpressed, label = delabel)) +
    geom_point() +
    geom_vline(xintercept=c(0.5, -0.5), col="red") + 
    geom_hline(yintercept=-log10(0.05), col="red") + 
    scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 25) +
    labs(title = "Differential expression", 
         subtitle = paste0(condition_compared, " using CERAD score as criteria"))
  
  ggsave(output_file, plot = vplot, width = 11, height = 10, units = "in", dpi = 300)
}

#run GSEA
run_gsea <- function(DEGS.g, output_file) {
  gse_BP <- gseGO(geneList = DEGS.g, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "BH")
  
  dotplot_BP <- dotplot(gse_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  ggsave(output_file, plot = dotplot_BP, height = 15, width = 20, dpi = 300)
}

#Main -----

####### hombres sanos vs. enfermos ##########

# Load data
data <- load_data()

# Filtrar hombres únicamente
data$metadata <- data$metadata %>% filter(msex == 1)

# Preprocesamiento de datos (filtrar los hombres)
preprocessed_1 <- preprocess_data(data$counts, data$metadata, condition_column = "ceradsc")

# Análisis diferencial de expresión entre hombres sanos (0) y enfermos (1)
de_results_1 <- run_deseq(preprocessed_1$counts, preprocessed_1$metadata, condition_column = "ceradsc")

DEGS_1 <- extract_DEGS(de_results_1$dds, condition_column = "ceradsc", coldata = de_results_1$coldata)

# Guardar los resultados del volcano plot
output_file_1 <- "DEG_control_male_vs_AD_male_volcano_plot.jpg"
create_volcano_plot(DEGS_1$res.df, condition_compared = "Control males vs AD males", output_file_1)

# Ejecutar GSEA para los genes diferencialmente expresados
DEGS.g_1 <- setNames(de_results_1$res$log2FoldChange, rownames(de_results_1$res))
DEGS.g_1 <- sort(DEGS.g, decreasing = TRUE)

run_gsea(DEGS.g, "GSEA_control_male_vs_AD_male.jpg")

########## Mujeres sanas vs mujeres enfermas  ##########

# Filtrar hombres únicamente
data$metadata <- data$metadata %>% filter(msex == 0)

# Preprocesamiento de datos (filtrar los hombres)
preprocessed_0 <- preprocess_data(data$counts, data$metadata, condition_column = "ceradsc")

# Análisis diferencial de expresión entre hombres sanos (0) y enfermos (1)
de_results_0 <- run_deseq(preprocessed$counts, preprocessed$metadata, condition_column = "ceradsc")

# Guardar los resultados del volcano plot
output_file_0 <- "DEG_women_control_vs_AD_volcano_plot.png"
create_volcano_plot(de_results_0$res, output_file)

# Ejecutar GSEA para los genes diferencialmente expresados
DEGS.g <- setNames(de_results_0$res$log2FoldChange, rownames(de_results_0$res))
DEGS.g <- sort(DEGS.g, decreasing = TRUE)

run_gsea(DEGS.g, "GSEA_women_control_vs_AD.png")

########## Mujeres enfermas contra hombres enfemos #######

# Filtrar mujeres y hombres enfermos únicamente
data$metadata <- data$metadata %>% filter(ceradsc == 1)

# Preprocesamiento de datos (filtrar por condición enferma)
preprocessed_AD <- preprocess_data(data$counts, data$metadata, condition_column = "msex")

# Análisis diferencial de expresión entre mujeres enfermas (0) y hombres enfermos (1), corrigiendo por sequencing batch
de_results_AD <- run_deseq(preprocessed$counts, preprocessed$metadata, condition_column = "msex")

# Guardar los resultados del volcano plot
output_file_AD <- "DEG_mujeres_enfermas_vs_hombres_enfermos_volcano_plot.png"
create_volcano_plot(de_results$res, output_file_AD)

# Ejecutar GSEA para los genes diferencialmente expresados
DEGS.g <- setNames(de_results$res$log2FoldChange, rownames(de_results$res))
DEGS.g <- sort(DEGS.g, decreasing = TRUE)

run_gsea(DEGS.g, "GSEA_AD_women_vs_AD_males.png")