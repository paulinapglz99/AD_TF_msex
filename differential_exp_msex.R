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

preprocess_data <- function(counts, metadata, condition_column, sex, ceradsc_condition) {
  counts <- counts[-c(1:4),]  # Remove the first four rows
  counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))
  counts <- del_dupl(counts)
  rownames(counts) <- NULL
  counts <- counts %>% column_to_rownames(var = "feature")
  metadata <- metadata %>% filter(!is.na(!!sym(condition_column)))
  metadata <- metadata %>% filter(msex %in% sex)
  metadata <- metadata %>% filter(ceradsc %in% ceradsc_condition)
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

extract_DEGS <- function(dds, condition_column, coldata, contrast_problem, contrast_control){
  
  res <- results(dds, contrast = c(condition_column, contrast_problem, contrast_control), #condition, problem, control 
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
  baseMean_4 <- rowMeans(normcounts[, coldata$ceradsc == contrast_control]) %>% as.data.frame()
  colnames(baseMean_4) <- paste0("baseMean_", contrast_control)
  
  baseMean_1 <- rowMeans(normcounts[, coldata$ceradsc == contrast_problem])%>% as.data.frame()
  colnames(baseMean_1) <- paste0("baseMean_", contrast_problem)
  
  baseMean_genes <- cbind(baseMean_4, baseMean_1)
  baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")
  
  #DEGS are...
  DEGS <- res.df %>% filter(diffexpressed != "NO")
  dim(DEGS)
  
  baseMean_DEGS <- baseMean_genes %>% filter(feature %in% DEGS$ensembl_gene_id)
  DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))
  
  
  return(list(res.df = res.df, DEGS = DEGS, normcounts = normcounts))
  
}


extract_DEGS <- function(dds, condition_column, coldata, padj.v, log2FoldChange.v, contrast_problem, contrast_control){
  
  res <- results(dds, contrast = c(condition_column, contrast_problem, contrast_control), #condition, problem, control 
                 pAdjustMethod = 'BH', alpha = 0.05)
  res.df <- res %>% as.data.frame() %>% filter(!is.na(padj))
  
  # add a column of NAs
  res.df$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  res.df$diffexpressed[res.df$log2FoldChange > log2FoldChange.v & res.df$padj < padj.v] <- "UP"
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  res.df$diffexpressed[res.df$log2FoldChange < -(log2FoldChange) & res.df$padj < padj.v] <- "DOWN"
  
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
  baseMean_4 <- rowMeans(normcounts[, coldata$ceradsc == contrast_control]) %>% as.data.frame()
  colnames(baseMean_4) <- paste0("baseMean_", contrast_control)
  
  baseMean_1 <- rowMeans(normcounts[, coldata$ceradsc == contrast_problem]) %>% as.data.frame()
  colnames(baseMean_1) <- paste0("baseMean_", contrast_problem)
  
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
  
  return(list(dds = dds, res = res, res.df = res.df, DEGS = DEGS, normcounts = normcounts))
  
}

#Visualize results

create_volcano_plot <- function(res.df, condition_compared, brain_region, output_file) {

  vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj),  col = diffexpressed, label = delabel)) +
    geom_point() +
    geom_vline(xintercept=c(0.5, -0.5), col="red") + 
    geom_hline(yintercept=-log10(0.05), col="red") + 
    scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 25) +
    labs(title = "Differential expression", 
         subtitle = paste0(condition_compared, " using CERAD score as criteria \n in the ", brain_region))
  
  ggsave(output_file, plot = vplot, width = 11, height = 10, units = "in", dpi = 300)
}

#Create heatmap

create_heatmap_vst <- function(dds, geneNames, output_file) {
  # Aplicar la transformación de variancia estabilizada (vst)
  vst_data <- vst(dds, blind = FALSE)
  
  # Normalizar las cuentas regularizadas por filas
  Z <- t(scale(t(assay(vst_data))))
  
  # Filtrar la matriz Z solo para los genes que están en geneNames
  Z_filtered <- Z[rownames(Z) %in% geneNames, ]
  
  # Asegurarse de que los nombres de los genes coinciden
  geneNames_filtered <- geneNames[geneNames %in% rownames(Z_filtered)]
  
  # Convertir la matriz a formato apropiado para heatmap
  mat_gn <- Z_filtered[geneNames_filtered, , drop = FALSE]  # Filtrar solo los genes de interés
  rownames(mat_gn) <- geneNames_filtered  # Usar los nombres de los genes para las filas
  
  # Crear y guardar el heatmap directamente con pheatmap
  pheatmap(mat_gn, 
           name = "Row Z-Score", 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
           clustering_distance_rows = "euclidean", 
           clustering_method = "ward.D2", 
           fontsize = 10, 
           angle_col = 45,  # Ángulo de las etiquetas de las columnas
           show_rownames = TRUE,  # Cambia a TRUE si quieres mostrar los nombres de los genes
           show_colnames = FALSE, 
           filename = output_file)  # Guardar directamente el heatmap en el archivo
}


#Plot MA 

plot_MA <- function(res, output_MAplot){
  # Abrir dispositivo gráfico según la extensión del archivo
  jpeg(output_MAplot, width = 11, height = 10, units = "in", res = 300)  # Guardar en formato PNG

  # Generar el MA plot
  plotMA(res, ylim=c(-2,2))
  
  # Cerrar el dispositivo gráfico
  dev.off()
}


#run GSEA
run_gsea <- function(DEGS.g, name_BP_dotplot, name_MF_dotplot, name_simmatrix_BP, name_simmatrix_MF) {
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
  
  gse_MF <-  gseGO(geneList = DEGS.g, 
                   ont ="MF", 
                   keyType = "ENSEMBL", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = "org.Hs.eg.db", 
                   pAdjustMethod = "BH")
  
  
  dotplot_MF <- dotplot(gse_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  
  #emmap plot BP
  simatrix_BP <- pairwise_termsim(gse_BP)
  simatrix_BP.m <- as.matrix(simatrix_BP@termsim)
  simatrix_BP.g <- graph_from_adjacency_matrix(simatrix_BP.m)
  degree_BP <- degree(simatrix_BP.g) %>% as.data.frame()
  degree_BP$. <- sort(degree_BP$., decreasing = TRUE)
  
  emmap_BP <- emapplot(simatrix_BP, showCategory = 25, group_legend = T,
                       layout = "kk", repel = T) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  #emmap plot BP
  simatrix_MF <- pairwise_termsim(gse_MF)
  simatrix_MF.m <- as.matrix(simatrix_MF@termsim)
  simatrix_MF.g <- graph_from_adjacency_matrix(simatrix_MF.m)
  degree_MF <- degree(simatrix_MF.g) %>% as.data.frame()
  degree_MF$. <- sort(degree_MF$., decreasing = TRUE)
  
  emmap_MF <- emapplot(simatrix_BP, showCategory = 25, group_legend = T,
                       layout = "kk", repel = T) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  
  ggsave(name_BP_dotplot, plot = dotplot_BP, height = 15, width = 20, dpi = 300)
  ggsave(name_MF_dotplot, plot = dotplot_MF, height = 15, width = 20, dpi = 300)
  ggsave(name_simmatrix_BP, plot = emmap_BP, height = 15, width = 20, dpi = 300)
  ggsave(name_simmatrix_MF, plot = emmap_MF, height = 15, width = 20, dpi = 300)
  
  return(list(simatrix_MF = simatrix_MF,
              simatrix_BP = simatrix_BP)
  )
}

#Main -----

# Load data
data <- load_data()

####### Control male vs AD males ##########

# Pre-proces data, filter by condition and sex
preprocessed_1 <- preprocess_data(data$counts, data$metadata, 
                                  condition_column = "ceradsc", 
                                  sex = "1", ceradsc_condition = c("1", "4"))

# DEA control males vs AD males
de_results_1 <- run_deseq(preprocessed_1$counts, preprocessed_1$metadata, condition_column = "ceradsc")

DEGS_1 <- extract_DEGS(de_results_1$dds, condition_column = "ceradsc", 
                       coldata = de_results_1$coldata, 
                       contrast_problem = "1",
                       contrast_control = "4")

#Save log2fold change full list 

vroom::vroom_write(DEGS_1$res.df, file = "ROSMAP_DLFPC_DE_full_gene_list_control_male_vs_AD_male.txt")

#Save DEGs list 

vroom::vroom_write(DEGS_1$DEGS, file = "ROSMAP_DLFPC_DEGS_control_male_vs_AD_male.txt")

#Save volcano plot

create_volcano_plot(DEGS_1$res.df, condition_compared = "Control males vs AD males",
                    brain_region = "DLPFC - ROSMAP", 
                    output_file = "ROSMAP_DLFPC_DEG_control_male_vs_AD_male_volcano_plot.jpg")

#Save heatmap

create_heatmap_vst(dds = de_results_1$dds, geneNames = DEGS_1$DEGS$ensembl_gene_id, 
                   output_file = "ROSMAP_DLFPC_DEG_control_male_vs_AD_male_heatmap.jpg")

#Save MA plot

plot_MA(DEGS_1$res, output_MAplot = "ROSMAP_DLFPC_DEG_control_male_vs_AD_male_MA_plot.jpg")

# Run GSEA for DEGS
DEGS.g_1 <- setNames(DEGS_1$res$log2FoldChange, rownames(DEGS_1$res))
DEGS.g_1 <- sort(DEGS.g_1, decreasing = TRUE)

GSEA_1 <- run_gsea(DEGS.g_1,
         name_BP_dotplot = "ROSMAP_DLFPC_GSEA_control_male_vs_AD_male_BP_dotplot.jpg",
         name_MF_dotplot  = "ROSMAP_DLFPC_GSEA_control_male_vs_AD_male_MF_dotplot.jpg", 
         name_simmatrix_BP = "ROSMAP_DLFPC_GSEA_control_male_vs_AD_male_BP_emmapplot.jpg",
         name_simmatrix_MF = "ROSMAP_DLFPC_GSEA_control_male_vs_AD_male_MF_emmapplot.jpg")

#Save enrichments

vroom::vroom_write(GSEA_1$simatrix_BP@result, file ="ROSMAP_DLFPC_GO_BP_results_control_male_vs_AD_male.txt")
vroom::vroom_write(GSEA_1$simatrix_MF@result, file ="ROSMAP_DLFPC_GO_MF_results_control_male_vs_AD_male.txt")

########## Control female vs AD females  ##########

#Pre-processing data

preprocessed_0 <- preprocess_data(data$counts, data$metadata, 
                                  condition_column = "ceradsc", 
                                  sex = "0", ceradsc_condition = c("1", "4"))

# DEA Control female vs AD females
de_results_0 <- run_deseq(preprocessed_0$counts, preprocessed_0$metadata, condition_column = "ceradsc")

DEGS_0 <- extract_DEGS(de_results_0$dds, condition_column = "ceradsc", 
                       coldata = de_results_0$coldata, 
                       contrast_problem = "1",
                       contrast_control = "4")

#Save log2fold change full list 

vroom::vroom_write(DEGS_0$res.df, file = "ROSMAP_DLFPC_DE_full_gene_list_control_female_vs_AD_female.txt")

#Save DEGs list 

vroom::vroom_write(DEGS_0$DEGS, file = "ROSMAP_DLFPC_DEGS_control_female_vs_AD_female.txt")

# Save volcano plot
 
create_volcano_plot(DEGS_0$res.df, condition_compared = "Control females vs AD females", 
                    brain_region = "DLPFC - ROSMAP",
                    output_file ="ROSMAP_DLFPC_DEG_control_female_vs_AD_female_volcano_plot.jpg")

#Save heatmap

create_heatmap_vst(dds = de_results_0$dds, geneNames = DEGS_0$DEGS$ensembl_gene_id, 
                   output_file = "ROSMAP_DLFPC_DEG_control_female_vs_AD_female_heatmap.jpg")

#Save MA plot

plot_MA(DEGS_0$res, output_MAplot = "ROSMAP_DLFPC_DEG_control_female_vs_AD_female_MA_plot.jpg")

# Ejecutar GSEA para los genes diferencialmente expresados
DEGS.g_0 <- setNames(DEGS_0$res$log2FoldChange, rownames(DEGS_0$res))
DEGS.g_0 <- sort(DEGS.g_0, decreasing = TRUE)

GSEA_0 <- run_gsea(DEGS.g_0,
         name_BP_dotplot = "ROSMAP_DLFPC_GSEA_control_female_vs_AD_female_BP_dotplot.jpg",
         name_MF_dotplot  = "ROSMAP_DLFPC_GSEA_control_female_vs_AD_female_MF_dotplot.jpg", 
         name_simmatrix_BP = "ROSMAP_DLFPC_GSEA_control_female_vs_AD_female_BP_emmapplot.jpg",
         name_simmatrix_MF = "ROSMAP_DLFPC_GSEA_control_female_vs_AD_female_MF_emmapplot.jpg")

#Save enrichments

vroom::vroom_write(GSEA_0$simatrix_BP@result, file ="ROSMAP_DLFPC_GO_BP_results_control_female_vs_AD_female.txt")
vroom::vroom_write(GSEA_0$simatrix_MF@result, file ="ROSMAP_DLFPC_GO_MF_results_control_female_vs_AD_female.txt")

########## AD females vs AD males  #######

# Data pre-process
preprocessed_AD <- preprocess_data(data$counts, data$metadata,
                                   condition_column = "msex", 
                                   sex = c("0", "1"), ceradsc_condition = "4")

# DEA Control female vs AD females
de_results_AD <- run_deseq(preprocessed_AD$counts, preprocessed_AD$metadata, condition_column = "msex")

DEGS_AD <- extract_DEGS(de_results_AD$dds, condition_column = "msex",
                        coldata = de_results_AD$coldata,
                        contrast_problem = "0",
                        contrast_control = "1")

#Save log2fold change full list 

vroom::vroom_write(DEGS_AD$res.df, file = "ROSMAP_DLFPC_DE_full_gene_list_AD_female_vs_AD_male.txt")

#Save DEGs list 

vroom::vroom_write(DEGS_AD$DEGS, file = "ROSMAP_DLFPC_DEGS_AD_female_vs_AD_male.txt")

# Save volcano plot

create_volcano_plot(DEGS_AD$res.df, condition_compared = "Control AD males vs AD females",
                    brain_region = "DLPFC - ROSMAP",
                    output_file = "ROSMAP_DLFPC_DEG_AD_female_vs_AD_male_volcano_plot.jpg")

#Save heatmap

create_heatmap_vst(dds = de_results_AD$dds, geneNames = DEGS_AD$DEGS$ensembl_gene_id, 
                   output_file = "ROSMAP_DLFPC_DEG_AD_female_vs_AD_male_heatmap.jpg")

#Save MA plot

plot_MA(DEGS_0$res, output_MAplot = "ROSMAP_DLFPC_DEG_AD_female_vs_AD_male__MA_plot.jpg")


# Ejecutar GSEA para los genes diferencialmente expresados
DEGS.g_AD <- setNames(DEGS_AD$res$log2FoldChange, rownames(DEGS_AD$res))
DEGS.g_AD <- sort(DEGS.g_AD, decreasing = TRUE)

GSEA_AD <- run_gsea(DEGS.g_AD,
         name_BP_dotplot = "ROSMAP_DLFPC_GSEA_control_AD_female_vs_AD_male_BP_dotplot.jpg",
         name_MF_dotplot  = "ROSMAP_DLFPC_GSEA_control_AD_female_vs_AD_male_MF_dotplot.jpg", 
         name_simmatrix_BP = "ROSMAP_DLFPC_GSEA_control_AD_female_vs_AD_male_BP_emmapplot.jpg",
         name_simmatrix_MF = "ROSMAP_DLFPC_GSEA_control_AD_female_vs_AD_male_MF_emmapplot.jpg")

#END ----