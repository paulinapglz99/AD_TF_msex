#!/usr/bin/env Rscript

#differential_exp_msex.R
#This script performs Differential expression analyis for 3 types of comparisons
#AD women vs control women
#AD men vs control men
#AD men vs AD women

time <- Sys.time()

cat("# Args:\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Error: Expected 3 arguments:\n
       1. Path to counts file (.rds)\n 
       2. Path to metadata file (.txt or .csv)\n
       3. Brain area name\n",
       call. = FALSE)
}

cat("Interpreted args:\n")
cat("Counts file: ", args[1], "\n")
cat("Metadata file: ", args[2], "\n")
cat("Brain area: ", args[3], "\n")

counts <- args[1]
metadata <- args[2]
brain_area <- args[3]

cat("# Load libraries\n")
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

# Set seed
set.seed(10)

counts<- readRDS(file = "/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.rds") %>% as.data.frame()
metadata <- vroom::vroom("/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.csv")
brain_area <- "PCC_ROSMAP"

cat("# Define functions \n")

#Load data

load_data <- function(args) {
  counts <- readRDS(file = args[1]) %>% as.data.frame()
  metadata <- vroom::vroom(file = args[2])
  brain_area <- as.character(args[3])  #Text in ""
  
  return(list(counts = counts, metadata = metadata))
}


#Define function to convert from ENSMBL to SYMBOL

ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")

convert_ens_to_symbol <- function(ensembl_ids) {
  trad <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                filters = "ensembl_gene_id",
                values = ensembl_ids,
                mart = ensembl)
  trad$external_gene_name <- ifelse(trad$external_gene_name == "", trad$ensembl_gene_id, trad$external_gene_name)
  return(trad)
}

#Function to delete duplicates
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
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")
  
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#Pre-process data

preprocess_data <- function(counts, metadata, condition_column) {
  print(dim(metadata))
  print(dim(counts))
  counts <- counts[-c(1:4), ]
  print(dim(counts))
  
  counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))
  counts <- del_dupl(counts)
  rownames(counts) <- NULL
  counts <- counts %>% column_to_rownames(var = "feature")
  
  #Ensure condition_column is treated as a factor
  metadata <- metadata %>% mutate(condition_column = factor(condition_column))
  
  cat("# Filter metadata based on conditions\n")
  if ("exclude" %in% colnames(metadata)) {
    metadata <- metadata[metadata$exclude != TRUE | is.na(metadata$exclude), ]
  }
  metadata <- metadata %>% filter(!is.na(condition_column))
  print(dim(metadata))
  
  cat("#Filter NAs in sequencingBatch\n")
  if ("sequencingBatch" %in% colnames(metadata)) {
    metadata <- metadata[!is.na(metadata$sequencingBatch), ]
  }
  print(dim(metadata))
  
  #Remove value ‘0, 6, 7’ if it exists
  if ("0, 6, 7" %in% metadata$sequencingBatch) {
    metadata <- metadata %>% filter(sequencingBatch != "0, 6, 7")
  }
  
  specimenIDs_in_counts <- colnames(counts)  #Assuming that the first column is not a specimenID
  metadata <- metadata %>%
    dplyr::filter(specimenID %in% specimenIDs_in_counts)
  print(dim(metadata))
  
  specimenIDs <- metadata$specimenID
  counts <- counts %>% 
    dplyr::select(c(1, any_of(specimenIDs)))
  print(dim(counts))
  
  specimenIDs_in_counts <- colnames(counts)  #Assuming that the first column is not a specimenID
  metadata <- metadata %>% dplyr::filter(specimenID %in% specimenIDs_in_counts)
  
  #Drop NA in metadata
  metadata <- metadata %>% drop_na(pmi)
  
  print(dim(metadata))
  return(list(counts = counts, metadata = metadata))
  
}

counts = preprocessed_data$counts
metadata = preprocessed_data$metadata
condition_column = "msex"
design_vars = "~ sequencingBatch + pmi + RIN + msex + ceradsc + msex:ceradsc"


#Run deseq

run_deseq <- function(counts, metadata, condition_column, design_vars) {
  
  # Verificar alineamiento counts y metadata
  if (!all(colnames(counts) %in% metadata$specimenID)) {
    stop("Not all columns in counts have corresponding specimenIDs in metadata!")
  }
  
  # Filtrar y alinear metadata
  metadata <- metadata %>%
    filter(specimenID %in% colnames(counts)) %>%
    arrange(match(specimenID, colnames(counts)))
  
  # Extraer variables de la fórmula
  design_formula <- as.formula(design_vars)
  vars <- all.vars(design_formula)
  
  # Construir coldata alineado
  coldata <- metadata %>%
    select(specimenID, all_of(vars))
  
  # Convertir variables categóricas a factor
  for (v in vars) {
    if (is.character(coldata[[v]])) {
      coldata[[v]] <- as.factor(coldata[[v]])
    }
  }
  
  # Convertir condition_column a factor (por seguridad)
  coldata[[condition_column]] <- as.factor(coldata[[condition_column]])
  coldata <- droplevels(coldata)
  
  # Verificar niveles en condition_column
  if (nlevels(coldata[[condition_column]]) < 2) {
    stop(paste("Variable", condition_column, "has less than two levels."))
  }
  
  # Relevel si existe "control"
  if ("control" %in% levels(coldata[[condition_column]])) {
    coldata[[condition_column]] <- relevel(coldata[[condition_column]], ref = "control")
  }
  
  # Escalamiento RIN y PMI
  for (var in c("RIN", "pmi")) {
    if (var %in% vars) {
      val_mean <- mean(coldata[[var]], na.rm = TRUE)
      val_sd <- sd(coldata[[var]], na.rm = TRUE)
      if (val_mean >= 5 || val_sd >= 5) {
        message(paste("Scaling", var, "for having mean or sd > 5"))
        coldata[[var]] <- scale(coldata[[var]])
      }
    }
  }
  
  # Finalmente correr DESeq2
  dds <- DESeqDataSetFromMatrix(countData = counts, 
                                colData = coldata,
                                design = design_formula)
  dds <- DESeq(dds)
  
  return(list(dds = dds, 
              coldata = coldata,
              design_formula = design_formula))
}



run_deseq <- function(counts, metadata, condition_column, design_vars) {
  #Build design formula
  design_formula <- as.formula(design_vars)
  vars <- all.vars(design_formula)
  
  #Prepare coldata
  coldata <- data.frame(specimenID = colnames(counts)) %>%
    left_join(metadata, by = "specimenID") %>% 
    dplyr::select(specimenID, all_of(vars)) %>% 
    drop_na()
  #Ensure that msex is a factor if present
  if (!is.factor(coldata$msex)) {
    coldata$msex <- as.factor(coldata$msex)
  }
  #Convert column_condition to factor
  coldata[[condition_column]] <- as.factor(coldata[[condition_column]])
  coldata <- droplevels(coldata)
  #Reorder levels if ‘control’ is present
  if ("control" %in% levels(coldata[[condition_column]])) {
    coldata[[condition_column]] <- relevel(coldata[[condition_column]], ref = "control")
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

  #Run DESeq2
   dds <- DESeqDataSetFromMatrix(countData = counts, 
                                 colData = coldata,
                                 design = design_formula)
   dds <- DESeq(dds)
  
  return(list(dds = dds, 
              coldata = coldata,
              design_formula = design_formula))
}

de_results_AD <- run_deseq(counts = preprocessed_data$counts, 
                           metadata = preprocessed_data$metadata,
                           condition_column = "msex", 
                           design_vars = "~ sequencingBatch + pmi + RIN + msex + ceradsc + msex:ceradsc") #I'm expecting interaction between sex and condition


#Extract DEGS

extract_DEGS <- function(dds, condition_column, coldata, padj.v, log2FoldChange.v, 
                         contrast_problem, contrast_control, name = NULL, label, brain_area){
  
  print(resultsNames(dds))
  
  #Extract results according to whether a 'name' or a 'contrast' is given
  if (!is.null(name)) {
    message(paste("Using specific coefficient:", name))
    res <- results(dds, name = name, pAdjustMethod = 'BH', alpha = padj.v)
  } else {
    if (is.null(condition_column) || is.null(contrast_problem) || is.null(contrast_control)) {
      stop("If 'name' is not defined, you must provide 'condition_column',
           'contrast_problem' and 'contrast_control'. 
           and 'contrast_control'.")
    }
    message(paste("Using contrast: ", contrast_problem, " vs ", contrast_control))
    res <- results(dds, 
                   contrast = c(condition_column, contrast_problem, contrast_control), 
                   pAdjustMethod = 'BH', alpha = padj.v)
  }
  
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
  
  # BaseMean by condition --- ---
  # Extract counts normalized 
  normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()
  
  #baseMean per condition
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
  
  # Return the result
  return(list(res = res, res.df = res.df, DEGS = DEGS, normcounts = normcounts))
}

#Visualize results

create_volcano_plot <- function(res.df, condition_compared, brain_region, output_file, log2fc) {
  
  vplot <- ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj), 
                                   col = diffexpressed, label = delabel)) +
    geom_point() +
    geom_vline(xintercept=c(log2fc, -log2fc), col="red") + 
    geom_hline(yintercept=-log10(1), col="red") + 
    scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
    theme_minimal() +
    geom_text_repel(max.overlaps = 25) +
    labs(title = "Differential expression", 
         subtitle = paste0(condition_compared, " using CERAD score as criteria \n in ", brain_region))
  
  ggsave(output_file, plot = vplot, width = 11, height = 15, units = "in", dpi = 300)
}

#Create heatmap
# create_heatmap_vst <- function(dds, geneNames, output_file, contrast) {
#   #Verify if there are enough DEGS to plot the heatmap
#   if (length(geneNames) < 2) {
#     message("There are not enough DEGs to generate a heatmap. The process will not run.\n")
#     return(NULL)
#   } else{
#     
#     #Apply the vst transformation
#     vst_data <- vst(dds, blind = FALSE)
#     #Normalize per row
#     Z <- t(scale(t(assay(vst_data))))
#     #Filter matrix Z
#     Z_filtered <- Z[rownames(Z) %in% geneNames, ]
#     geneNames_filtered <- geneNames[geneNames %in% rownames(Z_filtered)]
#     #Convert matrix to appropiate format for the heatmap
#     mat_gn <- Z_filtered[geneNames_filtered, , drop = FALSE]  #Only genes of interest
#     rownames(mat_gn) <- geneNames_filtered  #Use the names of genes for the rows
#     #Extract annotatoins from columns
#     my_sample_col <- as.data.frame(colData(dds))
#     my_sample_col <- my_sample_col[, c("ceradsc", "msex")]
#     #Create and save the heatmap directly with pheatmap
#     pheatmap(mat_gn, 
#              main = paste0(brain_area, " ", contrast," DEGs"), 
#              name = "Row Z-Score", 
#              color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
#              clustering_distance_rows = "euclidean", 
#              clustering_method = "ward.D2", 
#              annotation_col = my_sample_col,
#              fontsize = 10, 
#              angle_col = 45, 
#              show_rownames = TRUE,
#              show_colnames = FALSE, 
#              filename = output_file)  
#   } }

heatmap.f <- function(dds, DEGS_names, metadata, filename){
 
   #Apply VST
  vst_counts <- vst(dds, blind = TRUE)
  vst_mat <- assay(vst_counts)  #Transformed matrix
  
  normcounts <- vst_mat[rownames(vst_mat) %in% DEGS_names, ]
  normcounts_degs_scaled <- t(scale(t(normcounts)))

  annotation_col <- metadata %>%
    dplyr::select(specimenID, ceradsc, msex)%>%
    filter(specimenID %in% colnames(normcounts_degs_scaled)) %>%
    column_to_rownames("specimenID")
  
  gene_names <- convert_ens_to_symbol(rownames(normcounts_degs_scaled))
  gene_names <- setNames(gene_names$external_gene_name, gene_names$ensembl_gene_id)
  
  rownames(normcounts_degs_scaled) <- gene_names[rownames(normcounts_degs_scaled)]
  
  pheatmap(normcounts_degs_scaled,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           #scale = "row",
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100), 
           filename = filename)
  
}

#Plot MA 

plot_MA <- function(res, output_MAplot){
  jpeg(output_MAplot, width = 11, height = 15, units = "in", res = 300)
  
  #MA plot
  plotMA(res, ylim=c(-2,2))
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

cat("# Load data \n")
data <- load_data(args)
# data <- list()
# data$counts <- counts
# data$metadata <- metadata
preprocessed_data <- preprocess_data(data$counts, 
                                     metadata = data$metadata, 
                                     condition_column = "ceradsc")

####### Control male vs AD males ##########

cat("####### Control male vs AD males ##########\n")

male_metadata <- preprocessed_data$metadata %>% filter(msex == "1")
male_counts <- preprocessed_data$counts[, male_metadata$specimenID]

#DEA control males vs AD males
de_results_1 <- run_deseq(male_counts,
                          male_metadata,
                          condition_column = "ceradsc",
                          design_vars = "~ sequencingBatch + pmi + RIN + ceradsc")

#Results males
DEGS_1 <- extract_DEGS(de_results_1$dds, 
                      condition_column = "ceradsc", 
                      coldata = de_results_1$coldata, 
                      padj.v = 0.05,
                      log2FoldChange.v = 0.5, 
                      contrast_problem = "1",
                      contrast_control = "4",
                      #   label = "c_male_vs_AD_male",
                      brain_area = brain_area)

cat("#Save log2fold change full list\n")

vroom::vroom_write(DEGS_1$res.df,
                   file = paste0(brain_area,"_DE_full_gene_list_control_male_vs_AD_male.txt"))

cat("#Save DEGs list\n")

vroom::vroom_write(DEGS_1$DEGS, 
                   file = paste0(brain_area, "_DEGS_control_male_vs_AD_male.txt"))

cat("#Save volcano plot\n")

create_volcano_plot(DEGS_1$res.df,
                    condition_compared = "Control males vs AD males",
                    log2fc = 0.5,
                    brain_region = brain_area, 
                    output_file = paste0(brain_area, "_DEG_control_male_vs_AD_male_volcano_plot.jpg"))

cat("#Save heatmap")

heatmap.f(dds = de_results_1$dds,
          DEGS_names = DEGS_1$DEGS$ensembl_gene_id,
          metadata = de_results_1$coldata, 
          filename = paste0(brain_area, "_DEG_control_male_vs_AD_male_heatmap.jpg"))

# create_heatmap_vst(dds = de_results_1$dds, geneNames = DEGS_1$DEGS$external_gene_name, 
#                    contrast = "Control male vs AD male" ,
#                    output_file = paste0(brain_area, "_DEG_control_male_vs_AD_male_heatmap.jpg"))

cat("#Save MA plot \n")

plot_MA(DEGS_1$res, output_MAplot = paste0(brain_area, "_DEG_control_male_vs_AD_male_MA_plot.jpg"))

cat("#Run GSEA for male DEGS \n")
DEGS.g_1 <- setNames(DEGS_1$res$log2FoldChange, rownames(DEGS_1$res))
DEGS.g_1 <- sort(DEGS.g_1, decreasing = TRUE)

GSEA_1 <- run_gsea(DEGS.g_1,
                   name_BP_dotplot = paste0(brain_area, "_GSEA_control_male_vs_AD_male_BP_dotplot.jpg"),
                   name_MF_dotplot  = paste0(brain_area, "_GSEA_control_male_vs_AD_male_MF_dotplot.jpg"), 
                   name_simmatrix_BP = paste0(brain_area, "_GSEA_control_male_vs_AD_male_BP_emmapplot.jpg"),
                   name_simmatrix_MF = paste0(brain_area, "_GSEA_control_male_vs_AD_male_MF_emmapplot.jpg"))

cat("# Save enrichments \n")

if (!is.null(GSEA_1$simatrix_BP)) {
  vroom::vroom_write(GSEA_1$simatrix_BP@result, 
                     file = paste0(brain_area, "_GO_BP_results_control_male_vs_AD_male.txt"))
} else {
  message("No enriched terms found for BP.")
}

#Check for enriched terms in MF before writing the results
if (!is.null(GSEA_1$simatrix_MF)) {
  vroom::vroom_write(GSEA_1$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_male_vs_AD_male.txt"))
} else {
  message("No enriched terms found for MF.")
}

########## Control female vs AD females  ##########
cat("########## Control female vs AD females  ##########\n")

#Processing data

female_metadata <- preprocessed_data$metadata %>% filter(msex == "0")
female_counts <- preprocessed_data$counts[, female_metadata$specimenID]

# DEA Control female vs AD females
de_results_0 <- run_deseq(female_counts,
                          female_metadata, 
                          condition_column = "ceradsc", 
                          design_vars = "~ sequencingBatch  + RIN + pmi + ceradsc")

DEGS_0 <- extract_DEGS(de_results_0$dds, 
                       condition_column = "ceradsc", 
                       coldata = de_results_0$coldata, 
                       padj.v = 0.05,
                       log2FoldChange.v = 0.5, 
                       contrast_problem = "1",
                       contrast_control = "4",
                       #label = "c_female_vs_AD_female",
                       brain_area = brain_area)

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

# create_heatmap_vst(dds = de_results_0$dds, geneNames = DEGS_0$DEGS$external_gene_name, 
#                    contrast = "Control female vs AD female",
#                    output_file = paste0(brain_area, "_DEG_control_female_vs_AD_female_heatmap.jpg"))

heatmap.f(dds = de_results_0$dds,
          DEGS_names = DEGS_0$DEGS$ensembl_gene_id,
          metadata = de_results_0$coldata, 
          filename = paste0(brain_area, "_DEG_control_female_vs_AD_female_heatmap.jpg"))

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

#Check for rich terms in MF before writing results
if (!is.null(GSEA_0$simatrix_MF)) {
  vroom::vroom_write(GSEA_0$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_female_vs_AD_female.txt"))
} else {
  message("No enriched terms found for MF.")
}

########## AD females vs AD males  #######

cat("########## AD females vs AD males  #######\n")

#DEA Control female vs AD females
de_results_AD <- run_deseq(counts = preprocessed_data$counts, 
                           metadata = preprocessed_data$metadata,
                           condition_column = "msex", 
                           design_vars = "~ sequencingBatch + pmi + RIN + msex + ceradsc + msex:ceradsc") #I'm expecting interaction between sex and condition

DEGS_AD <- extract_DEGS(de_results_AD$dds, 
                        condition_column = "msex",
                        name = "msex1.ceradsc4", 
                        coldata = de_results_AD$coldata,
                        padj.v = 0.05,
                        log2FoldChange.v = 0.5,
                        contrast_problem = "0", #female
                        contrast_control = "1", #male
                        #  label = "AD_male_vs_AD_female",
                        brain_area = brain_area)

#res_interaction <- results(de_results_AD$dds, name = "msex1.ceradsc4")

#Save log2fold change full list 

vroom::vroom_write(DEGS_AD$res.df, file = paste0(brain_area, "_DE_full_gene_list_AD_female_vs_AD_male.txt"))

#Save DEGs list 

vroom::vroom_write(DEGS_AD$DEGS, file = paste0(brain_area, "_DEGS_AD_female_vs_AD_male.txt"))

# Save volcano plot

create_volcano_plot(DEGS_AD$res.df, 
                    condition_compared = "AD males vs AD females",
                    log2fc = 0.5,
                    brain_region = brain_area,
                    output_file = paste0(brain_area, "_DEG_AD_female_vs_AD_male_volcano_plot.jpg"))

#Save heatmap

# create_heatmap_vst(dds = de_results_AD$dds,
#                    geneNames = DEGS_AD$DEGS$ensembl_gene_id, 
#                    contrast = "AD female vs AD male", 
#                    output_file = paste0(brain_area, "_DEG_AD_female_vs_AD_male_heatmap.jpg"))

heatmap.f(
  dds = de_results_AD$dds,
  DEGS_names = DEGS_AD$DEGS$ensembl_gene_id,
  metadata = de_results_AD$coldata, 
  filename = paste0(brain_area, "_DEG_AD_female_vs_AD_male_heatmap.jpg"))

#Save MA plot

plot_MA(DEGS_AD$res, output_MAplot = paste0(brain_area, "_DEG_AD_female_vs_AD_male_MA_plot.jpg"))

#Run GSEA for DEGS
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

#Verify if there are enriched terms in MF before write results
if (!is.null(GSEA_AD$simatrix_MF)) {
  vroom::vroom_write(GSEA_AD$simatrix_MF@result, 
                     file = paste0(brain_area, "_GO_MF_results_control_AD_female_vs_AD_male.txt"))
} else {
  message("No enriched terms found for MF.")
}

cat("#Table for IPA")

#Rename cols in data
DEGS_0_renamed <- DEGS_0$res.df %>% 
  rename_with(~ paste0("DEGS_0_", .x), -ensembl_gene_id)

DEGS_1_renamed <- DEGS_1$res.df %>% 
  rename_with(~ paste0("DEGS_1_", .x), -ensembl_gene_id)

DEGS_AD_renamed <- DEGS_AD$res.df %>% 
  rename_with(~ paste0("DEGS_AD_", .x), -ensembl_gene_id)

#Join
IPA <- DEGS_0_renamed %>%
  left_join(DEGS_1_renamed, by = "ensembl_gene_id") %>%
  left_join(DEGS_AD_renamed, by = "ensembl_gene_id")

#Save table for IPA

vroom::vroom_write(IPA,
                   file = paste0(brain_area, "_IPA_table.txt") )

cat("END ----")
cat("Elapsed time:", Sys.time() - time, "mins")