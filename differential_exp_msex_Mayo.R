#differential_exp_sex_MSBB.R
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
               "enrichplot")
time <- Sys.time()
#Args

args <- commandArgs(trailingOnly = TRUE)

counts <- args[1]
metadata <- args[2]
brain_area <- args[3]
brain_area <- as.character(brain_area)

cat("Interpreted args:\n")
cat("Counts file: ", args[1], "\n")
cat("Metadata file: ", args[2], "\n")
cat("Brain area: ", args[3], "\n")
cat("Batch correction: ", args[4], "\n")

counts <- readRDS("/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_rawcounts_CRB.rds") %>% as.data.frame()
metadata <- vroom::vroom("/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_metadata_CRB.txt")
brain_area <- "TC_Mayo"
tissue_batch_correction <- FALSE

# Set seed
set.seed(10)

#Load data

load_data <- function(counts, metadata) {
  counts <- readRDS(file = args[1]) %>% as.data.frame()
  metadata <- vroom::vroom(file = args[2])
  
  return(list(counts = counts, metadata = metadata))
}

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
    summarize(across(everything(), \(x) median(x, na.rm = TRUE)))
  
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#Pre-process data

preprocess_data <- function(counts, metadata) {
  print(dim(metadata))
  print(dim(counts))
  counts <- counts[-c(1:4), ]
  print(dim(counts))
  
  counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))
  counts <- del_dupl(counts)
  rownames(counts) <- NULL
  counts <- counts %>% column_to_rownames(var = "feature")
  
  #Ensure diagnosis is treated as a factor
  metadata <- metadata %>%
    filter(diagnosis %in% c("Alzheimer Disease", "control")) %>% 
    mutate(diagnosis = factor(diagnosis))
  
  #Filter metadata based on conditions
  if ("exclude" %in% colnames(metadata)) {
    metadata <- metadata[metadata$exclude != TRUE | is.na(metadata$exclude), ]
  }
  
  #Filtering metadata before applying batch correction
  if ("sequencingBatch" %in% colnames(metadata)) {
    metadata <- metadata[!is.na(metadata$sequencingBatch), ]
  }
 
  print(dim(metadata))
  
  specimenIDs_in_counts <- colnames(counts) 
  metadata <- metadata %>%
    dplyr::filter(specimenID %in% specimenIDs_in_counts)
  print(dim(metadata))
  
  specimenIDs <- metadata$specimenID
  counts <- counts %>% 
    dplyr::select(c(1, all_of(specimenIDs)))
  counts <- counts[rowSums(counts) != 0, ]
  print(dim(counts))
  
  specimenIDs_in_counts <- colnames(counts)
  metadata <- metadata %>% dplyr::filter(specimenID %in% specimenIDs_in_counts)
  print(dim(metadata))
  return(list(counts = counts, metadata = metadata))
}

#Process data 

process_data <- function(counts, metadata, condition_column, sex_condition, diagnosis_condition) {
  metadata <- metadata %>% filter(!is.na(!!sym(condition_column))) %>% 
    filter(sex %in% sex_condition) %>% 
    filter(diagnosis %in% diagnosis_condition) %>% 
    filter(!is.na(specimenIdSource)) %>% 
    droplevels()  #Eliminate empty factor levels
  
  #Select columns in counts that match specimenID in metadata
  counts <- counts %>% dplyr::select(all_of(metadata$specimenID))
  
  return(list(counts = counts, metadata = metadata))
}

#Run deseq
run_deseq <- function(counts, metadata, condition_column, design_vars) {
  #Prepare coldata
  coldata <- data.frame(specimenID = colnames(counts)) %>%
    left_join(metadata, by = "specimenID")
  #Ensure that msex is a factor if present
  if (!is.factor(coldata$sex)) {
    coldata$sex <- as.factor(coldata$sex)
  }
  #Convert column_condition to factor
  coldata[[condition_column]] <- as.factor(coldata[[condition_column]])
  coldata <- droplevels(coldata)
  #Reorder levels if ‘control’ is present
  if ("control" %in% levels(coldata[[condition_column]])) {
    coldata[[condition_column]] <- relevel(coldata[[condition_column]], ref = "control")
  }
  #
  if ("specimenIdSource" %in% colnames(coldata) && nlevels(coldata$specimenIdSource) > 2 ) {
    coldata[["specimenIdSource"]] <- as.factor(coldata[["specimenIdSource"]])
  }
  #Verify sufficient levels in condition
  if (nlevels(coldata[[condition_column]]) < 2) {
    stop(paste("Variable", condition_column, "has less than two levels."))
  }
 
  #Scaling RIN and pmi if mean or sd > 5
  for (var in c("RIN", "pmi")) {
    val_mean <- mean(coldata[[var]], na.rm = TRUE)
    val_sd <- sd(coldata[[var]], na.rm = TRUE)
    if (val_mean >= 5 || val_sd >= 5) {
      message(paste("Scaling", var, "for having mean or sd > 5"))
      coldata[[var]] <- scale(coldata[[var]])
    }
  }
  #Build design formula
  
  if (nlevels(coldata$sex) > 1 && condition_column != "sex") {
    design_vars <- c(design_vars, "sex")
  }
  design_vars <- c(design_vars, condition_column)
  design_formula <- as.formula(paste("~", paste(design_vars, collapse = " + ")))
  
  #Execute DESeq2
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = coldata,
                                design = design_formula)
  dds <- DESeq(dds)
  
  return(list(dds = dds, 
              coldata = coldata))
}

#Extract DEGS

extract_DEGS <- function(dds, condition_column, coldata, padj.v, log2FoldChange.v, 
                         contrast_problem, contrast_control){
  
  res <- results(dds, contrast = c(condition_column, contrast_problem, contrast_control), #condition, problem, control 
                 pAdjustMethod = 'BH', alpha = 0.05)
  res.df <- res %>% as.data.frame() %>% filter(!is.na(padj))
  
  # add a column of NAs
  res.df$diffexpressed <- "NO"
  # if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
  res.df$diffexpressed[res.df$log2FoldChange > log2FoldChange.v & res.df$padj < padj.v] <- "UP"
  # if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
  res.df$diffexpressed[res.df$log2FoldChange < -(log2FoldChange.v) & res.df$padj < padj.v] <- "DOWN"
  
  # Add gene names in SYMBOL --- ---
  # Create dictionary
  symbol <- convert_ens_to_symbol(rownames(res))
  
  # Make rownames columns
  res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
  res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge
  
  # Extract counts normalized 
  normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()
  
  #BaseMean per condition
  baseMean_4 <- rowMeans(normcounts[, coldata$diagnosis == contrast_control]) %>% as.data.frame()
  colnames(baseMean_4) <- paste0("baseMean_", contrast_control)
  
  baseMean_1 <- rowMeans(normcounts[, coldata$diagnosis == contrast_problem]) %>% as.data.frame()
  colnames(baseMean_1) <- paste0("baseMean_", contrast_problem)
  
  baseMean_genes <- cbind(baseMean_4, baseMean_1)
  baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")
  
  #DEGs are...
  DEGS <- res.df %>% filter(diffexpressed != "NO")
  dim(DEGS)
  
  baseMean_DEGS <- baseMean_genes %>% filter(feature %in% DEGS$ensembl_gene_id)
  DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))
  
  #Add delabel column with gene names for DEGs
  DEGS <- DEGS %>% mutate(delabel = external_gene_name)
  
  #For non-DEGs, delabel will be NA or empty
  res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA))
  
  return(list(res = res, res.df = res.df, DEGS = DEGS, normcounts = normcounts))
  
}

#Visualize results

create_volcano_plot <- function(res.df, condition_compared, brain_region, output_file, log2fc) {

  vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj),  col = diffexpressed, label = delabel)) +
    geom_point() +
    geom_vline(xintercept=c(log2fc, -log2fc), col="red") + 
    geom_hline(yintercept=-log10(0.05), col="red") + 
    scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 25) +
    labs(title = "Differential expression", 
         subtitle = paste0(condition_compared, " using CERAD score as criteria \n in ", brain_region))
  
  ggsave(output_file, plot = vplot, width = 11, height = 10, units = "in", dpi = 300)
}

#Create heatmap
create_heatmap_vst <- function(dds, geneNames, output_file, contrast) {
  #Check if there are enough SDRs to generate the heatmap
  if (length(geneNames) < 2) {
    message("There are not enough DEGs to generate a heatmap. The process will not run.")
    return(NULL)
  } else{
  #Applying the variance stabilised transformation (vst)
  vst_data <- vst(dds, blind = FALSE)
  #Normalise the accounts regularised by rows
  Z <- t(scale(t(assay(vst_data))))
  #Filter the Z-array only for genes that are in geneNames
  Z_filtered <- Z[rownames(Z) %in% geneNames, ]
  #Make sure that the gene names match
  geneNames_filtered <- geneNames[geneNames %in% rownames(Z_filtered)]
  #Convert matrix to heatmap-appropriate format
  mat_gn <- Z_filtered[geneNames_filtered, , drop = FALSE]   #Filter out only genes of interest
  rownames(mat_gn) <- geneNames_filtered #Use gene names for the rows
  #Extract annotations from columns
  #my_sample_col <- as.data.frame(colData(dds))
  #Create and save the heatmap directly with pheatmap
  pheatmap(mat_gn, 
           main = paste0(brain_area, " ", contrast," DEGs"), 
           name = "Row Z-Score", 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
           clustering_distance_rows = "euclidean", 
           clustering_method = "ward.D2", 
           #annotation_col = my_sample_col,
           fontsize = 10, 
           angle_col = 45,  
           show_rownames = TRUE,
           show_colnames = FALSE, 
           filename = output_file)
} }

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
  #Run GSE BP
  gse_BP <- gseGO(geneList = DEGS.g, 
                  ont = "BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "BH")
  
  if (nrow(gse_BP) > 0) {
    # BP dotplot
    dotplot_BP <- dotplot(gse_BP, showCategory = 10, split = ".sign") + 
      facet_grid(.~.sign) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    # emmap plot BP
    simatrix_BP <- pairwise_termsim(gse_BP)
    
    emmap_BP <- emapplot(simatrix_BP, showCategory = 25, group_legend = TRUE,
                         layout = "kk", repel = TRUE) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    ggsave(name_BP_dotplot, plot = dotplot_BP, height = 15, width = 20, dpi = 300)
    ggsave(name_simmatrix_BP, plot = emmap_BP, height = 15, width = 20, dpi = 300)
  } else {
    warning("No enriched terms found for BP")
  }
  
  #Run GSE MF
  gse_MF <- gseGO(geneList = DEGS.g, 
                  ont = "MF", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = "org.Hs.eg.db", 
                  pAdjustMethod = "BH")
  
  if (nrow(gse_MF) > 0) {
    # MF dotplot
    dotplot_MF <- dotplot(gse_MF, showCategory = 10, split = ".sign") + 
      facet_grid(.~.sign) + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    # emmap plot MF
    simatrix_MF <- pairwise_termsim(gse_MF)
    
    emmap_MF <- emapplot(simatrix_MF, showCategory = 25, group_legend = TRUE,
                         layout = "kk", repel = TRUE) +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    ggsave(name_MF_dotplot, plot = dotplot_MF, height = 15, width = 20, dpi = 300)
    ggsave(name_simmatrix_MF, plot = emmap_MF, height = 15, width = 20, dpi = 300)
  } else {
    warning("No enriched terms found for MF")
  }
  
  return(list(simatrix_MF = if (exists("simatrix_MF")) simatrix_MF else NULL,
              simatrix_BP = if (exists("simatrix_BP")) simatrix_BP else NULL))
}

###### MAIN #####

# Load data
data <- load_data()

data <- preprocess_data(data$counts, metadata = data$metadata)

####### Control male vs AD males ##########

cat("####### Control male vs AD males ##########")

#Pre-proces data, filter by condition and sex
preprocessed_1 <- process_data(counts = data$counts, 
                               metadata = data$metadata, 
                                  condition_column = "diagnosis", 
                                  sex = "male", 
                                  diagnosis_condition = c("Alzheimer Disease", "control"))
table(preprocessed_1$metadata$sex)
table(preprocessed_1$metadata$diagnosis)

#DEA control males vs AD males
de_results_1 <- run_deseq(preprocessed_1$counts,
                          preprocessed_1$metadata, 
                          condition_column = "diagnosis", 
                          design_vars = c("RIN", "specimenIdSource"))

DEGS_1 <- extract_DEGS(de_results_1$dds, condition_column = "diagnosis", 
                       coldata = de_results_1$coldata, 
                       padj.v = 0.05,
                       log2FoldChange.v = 0.5,
                       contrast_problem = "Alzheimer Disease",
                       contrast_control = "control")

#Save log2fold change full list 

vroom::vroom_write(DEGS_1$res.df, file = paste0(brain_area,"_DE_full_gene_list_control_male_vs_AD_male.txt"))

#Save DEGs list 

vroom::vroom_write(DEGS_1$DEGS, file = paste0(brain_area, "_DEGS_control_male_vs_AD_male.txt"))

#Save volcano plot

create_volcano_plot(DEGS_1$res.df, condition_compared = "Control males vs AD males",
                    log2fc = 0.5,
                    brain_region = brain_area, 
                    output_file = paste0(brain_area, "_DEG_control_male_vs_AD_male_volcano_plot.jpg"))

#Save heatmap

create_heatmap_vst(dds = de_results_1$dds, geneNames = DEGS_1$DEGS$external_gene_name, 
                   contrast = "Control_male_vs_AD_male", 
                   output_file = paste0(brain_area, "_DEG_control_male_vs_AD_male_heatmap.jpg"))

#Save MA plot

plot_MA(DEGS_1$res, output_MAplot = paste0(brain_area, "_DEG_control_male_vs_AD_male_MA_plot.jpg"))

#Run GSEA for DEGS
DEGS.g_1 <- setNames(DEGS_1$res$log2FoldChange, rownames(DEGS_1$res))
DEGS.g_1 <- sort(DEGS.g_1, decreasing = TRUE)

GSEA_1 <- run_gsea(DEGS.g_1,
         name_BP_dotplot = paste0(brain_area, "_GSEA_control_male_vs_AD_male_BP_dotplot.jpg"),
         name_MF_dotplot  = paste0(brain_area, "_GSEA_control_male_vs_AD_male_MF_dotplot.jpg"), 
         name_simmatrix_BP = paste0(brain_area, "_GSEA_control_male_vs_AD_male_BP_emmapplot.jpg"),
         name_simmatrix_MF = paste0(brain_area, "_GSEA_control_male_vs_AD_male_MF_emmapplot.jpg"))

#Save enrichments

if (!is.null(GSEA_1$simatrix_BP)) {
  vroom::vroom_write(GSEA_1$simatrix_BP@result, 
                     file = paste0(brain_area, "_GO_BP_results_control_male_vs_AD_male.txt"))
} else {
  message("No enriched terms found for BP.")
}

# Verifica si hay términos enriquecidos en MF antes de escribir los resultados
if (!is.null(GSEA_1$simatrix_MF)) {
  vroom::vroom_write(GSEA_1$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_male_vs_AD_male.txt"))
} else {
  message("No enriched terms found for MF.")
}

########## Control female vs AD females  ##########
cat("########## Control female vs AD females  ##########")

#Pre-processing data

preprocessed_0 <- process_data(data$counts, data$metadata, 
                                  condition_column = "diagnosis", 
                                  sex = "female", diagnosis_condition = c("Alzheimer Disease", "control"))
table(preprocessed_0$metadata$sex)
table(preprocessed_0$metadata$diagnosis)

# DEA Control female vs AD females
de_results_0 <- run_deseq(preprocessed_0$counts, 
                          preprocessed_0$metadata,
                          condition_column = "diagnosis", 
                          design_vars = c("RIN", "specimenIdSource"))

DEGS_0 <- extract_DEGS(de_results_0$dds, 
                       condition_column = "diagnosis", 
                       coldata = de_results_0$coldata, 
                       padj.v = 0.05,
                       log2FoldChange.v = 0.5, 
                       contrast_problem = "Alzheimer Disease",
                       contrast_control = "control")

#Save log2fold change full list 

vroom::vroom_write(DEGS_0$res.df, file = paste0(brain_area, "_DE_full_gene_list_control_female_vs_AD_female.txt"))

#Save DEGs list 

vroom::vroom_write(DEGS_0$DEGS, file = paste0(brain_area, "_DEGS_control_female_vs_AD_female.txt"))

#Save volcano plot
 
create_volcano_plot(DEGS_0$res.df, condition_compared = "Control females vs AD females", 
                    log2fc = 0.5,
                    brain_region = brain_area, 
                    output_file = paste0(brain_area, "_DEG_control_female_vs_AD_female_volcano_plot.jpg"))

#Save heatmap

create_heatmap_vst(dds = de_results_0$dds, geneNames = DEGS_0$DEGS$external_gene_name, 
                   contrast = "Control_female_vs_AD_female",
                   output_file = paste0(brain_area, "_DEG_control_female_vs_AD_female_heatmap.jpg"))

#Save MA plot

plot_MA(DEGS_0$res, output_MAplot = paste0(brain_area, "_DEG_control_female_vs_AD_female_MA_plot.jpg"))

#Running GSEA for differentially expressed genes
DEGS.g_0 <- setNames(DEGS_0$res$log2FoldChange, rownames(DEGS_0$res))
DEGS.g_0 <- sort(DEGS.g_0, decreasing = TRUE)

GSEA_0 <- run_gsea(DEGS.g_0,
         name_BP_dotplot = paste0(brain_area, "_GSEA_control_female_vs_AD_female_BP_dotplot.jpg"),
         name_MF_dotplot  = paste0(brain_area, "_GSEA_control_female_vs_AD_female_MF_dotplot.jpg"), 
         name_simmatrix_BP = paste0(brain_area, "_GSEA_control_female_vs_AD_female_BP_emmapplot.jpg"),
         name_simmatrix_MF = paste0(brain_area, "_GSEA_control_female_vs_AD_female_MF_emmapplot.jpg"))

#Save enrichments

if (!is.null(GSEA_0$simatrix_BP)) {
  vroom::vroom_write(GSEA_0$simatrix_BP@result, 
                     file = paste0(brain_area, "_GO_BP_results_control_female_vs_AD_female.txt"))
} else {
  message("No enriched terms found for BP.")
}

# Verifica si hay términos enriquecidos en MF antes de escribir los resultados
if (!is.null(GSEA_0$simatrix_MF)) {
  vroom::vroom_write(GSEA_0$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_female_vs_AD_female.txt"))
} else {
  message("No enriched terms found for MF.")
}

########## AD females vs AD males  #######

cat("########## AD females vs AD males  #######\n")

#Data pre-process
preprocessed_AD <- process_data(data$counts, data$metadata,
                                   condition_column = "sex", 
                                   sex = c("female", "male"), 
                                   diagnosis_condition = "Alzheimer Disease")
table(preprocessed_AD$metadata$sex)
table(preprocessed_AD$metadata$diagnosis)

#DEA Control female vs AD females
de_results_AD <- run_deseq(preprocessed_AD$counts,
                           preprocessed_AD$metadata,
                           condition_column = "sex", 
                           design_vars = c("RIN"))

DEGS_AD <- extract_DEGS(de_results_AD$dds,
                        condition_column = "sex",
                        coldata = de_results_AD$coldata,
                        padj.v = 0.05,
                        log2FoldChange.v = 0.5,
                        contrast_problem = "female",
                        contrast_control = "male")

#Save log2fold change full list 

vroom::vroom_write(DEGS_AD$res.df, file = paste0(brain_area, "_DE_full_gene_list_AD_female_vs_AD_male.txt"))

#Save DEGs list 

vroom::vroom_write(DEGS_AD$DEGS, file = paste0(brain_area, "_DEGS_AD_female_vs_AD_male.txt"))

#Save volcano plot

create_volcano_plot(DEGS_AD$res.df, condition_compared = "Control AD males vs AD females",
                    brain_region = brain_area,
                    log2fc = 0.5, 
                    output_file = paste0(brain_area, "_DEG_AD_female_vs_AD_male_volcano_plot.jpg"))

#Save heatmap

create_heatmap_vst(dds = de_results_AD$dds, geneNames = DEGS_AD$DEGS$external_gene_name, 
                   contrast = "AD_female_vs_AD_male", 
                   output_file = paste0(brain_area, "_DEG_AD_female_vs_AD_male_heatmap.jpg"))

#Save MA plot

plot_MA(DEGS_0$res, output_MAplot = paste0(brain_area, "_DEG_AD_female_vs_AD_male_MA_plot.jpg"))

#Running GSEA for AD-AD differentially expressed genes
DEGS.g_AD <- setNames(DEGS_AD$res$log2FoldChange, rownames(DEGS_AD$res))
DEGS.g_AD <- sort(DEGS.g_AD, decreasing = TRUE)

GSEA_AD <- run_gsea(DEGS.g_AD,
         name_BP_dotplot = paste0(brain_area, "_GSEA_control_AD_female_vs_AD_male_BP_dotplot.jpg"),
         name_MF_dotplot  = paste0(brain_area, "_GSEA_control_AD_female_vs_AD_male_MF_dotplot.jpg"), 
         name_simmatrix_BP = paste0(brain_area, "_GSEA_control_AD_female_vs_AD_male_BP_emmapplot.jpg"),
         name_simmatrix_MF = paste0(brain_area, "_GSEA_control_AD_female_vs_AD_male_MF_emmapplot.jpg"))

#Save enrichments

if (!is.null(GSEA_AD$simatrix_BP)) {
  vroom::vroom_write(GSEA_AD$simatrix_BP@result, 
                     file = paste0(brain_area, "_GO_BP_results_control_AD_female_vs_AD_male.txt"))
} else {
  message("No enriched terms found for BP.")
}

#Check for rich terms in MF before writing results
if (!is.null(GSEA_AD$simatrix_MF)) {
  vroom::vroom_write(GSEA_AD$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_AD_female_vs_AD_male.txt"))
} else {
  message("No enriched terms found for MF.")
}

#Table for IPA

#Function to rename the columns except ‘ensembl_gene_id’.
rename_columns <- function(df, suffix) {
  colnames(df)[-1] <- paste0(colnames(df)[-1], "_", suffix)
  return(df)
}

#Rename the columns for each table
IPA_0 <- rename_columns(DEGS_0$res.df, paste0(brain_area, "_0"))
IPA_1 <- rename_columns(DEGS_1$res.df,  paste0(brain_area, "_1"))
IPA_AD <- rename_columns(DEGS_AD$res.df,  paste0(brain_area, "_AD"))

IPA <- IPA_0 %>% 
  left_join(IPA_1) %>% 
  left_join(IPA_AD)

dim(IPA)

IPA <- IPA %>%
  filter(complete.cases(dplyr::select(., starts_with("log2FoldChange"), starts_with("padj"))))
dim(IPA)

#Save table for IPA

vroom::vroom_write(IPA,
                   file = paste0(brain_area, "_IPA_table.txt") )

cat("Elapsed time:", Sys.time() - time, "minutes")

cat("END ----")