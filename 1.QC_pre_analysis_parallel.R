#
#1.QC_pre_analysis_parallel_ROSMAP.R
#Script to inspect RNA-seq data counts 
#paulinapglz.99@gmail.com

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

###Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2", 
               "gridExtra", 
               "furrr")

#FUNCTIONS --- ---

#Function to filter columns of a count matrix according to specimenIDs in the metadata

filter_by_metadata <- function(counts, metadata) {
  specimenIDs <- metadata$specimenID
  filtered_counts <- counts %>%
  dplyr::select(c(1, any_of(specimenIDs)))
  
  return(filtered_counts)
}

# Function to filter metadata according to the specimenIDs in the count matrix
filter_metadata_by_counts <- function(metadata, counts) {
  specimenIDs_in_counts <- colnames(counts)[-1]  # Asumiendo que la primera columna no es un specimenID
  filtered_metadata <- metadata %>%
             dplyr::filter(specimenID %in% specimenIDs_in_counts)
  
  return(filtered_metadata)
}

# Delete duplicates

del_dupl <- function(counts) {
  # Verificar genes repetidos
  repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)
  
  # Ver las filas duplicadas
  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  
  # Ordenar y calcular la mediana de los valores duplicados
  repeated_rows <- repeated_rows[order(repeated_rows$feature),]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), median, na.rm = TRUE))
  
  # Eliminar las filas duplicadas y agregar las filas con la mediana calculada
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  
  # Convertir las columnas seleccionadas a enteros
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#Get metadata --- ---

#Obtain factors from metadata
metadata_files <- list(
  #FC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FC/RNA_seq_metadata_FC.txt",
  #TC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/TC/RNA_seq_metadata_TC.txt",
  #DLPFC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt",
  HCN = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN/RNA_seq_metadata_HCN.txt",
  PCC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.txt",
  FRONTAL_ROSMAP = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FRONTAL/RNA_seq_metadata_FRONTAL_ROSMAP.txt",
  #FP_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FP/MSBB_RNAseq_metadata_FP.txt",
  #IFG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/IFG/MSBB_RNAseq_metadata_IFG.txt",
  #PFC_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PFC/MSBB_RNAseq_metadata_PFC.txt",
  PHG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PHG/MSBB_RNAseq_metadata_PHG.txt",
  STG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/STG/metadata/MSBB_RNAseq_metadata_STG.txt",
  FRONTAL_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FRONTAL/MSBB_RNAseq_metadata_FRONTAL_MSBB.txt",
  CRB_Mayo = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_metadata_CRB.txt",
  TC_Mayo = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_metadata_TC.txt"
)

# Leer todos los archivos de metadatos en una lista
metadata <- lapply(metadata_files,  function(file) {
  vroom::vroom(file, delim = "\t")
})
length(metadata)
lapply(metadata, dim)

#Get data --- ---

files <- list(
  #FC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FC/ROSMAP_RNAseq_rawcounts_FC.rds",
  #TC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/TC/ROSMAP_RNAseq_rawcounts_TC.rds",
  #DLPFC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds",
  HCN = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/HCN/ROSMAP_RNAseq_rawcounts_HCN.rds",
  PCC = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.rds",
  FRONTAL_ROSMAP = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FRONTAL/ROSMAP_RNAseq_rawcounts_FRONTAL_ROSMAP.rds",
  #FP_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FP/MSBB_RNAseq_rawcounts_FP.rds",
  #IFG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/IFG/MSBB_RNAseq_rawcounts_IFG.rds",
  #PFC_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PFC/MSBB_RNAseq_rawcounts_counts_PFC_MSBB.rds",
  PHG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PHG/MSBB_RNAseq_rawcounts_PHG.rds",
  STG_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/STG/MSBB_RNAseq_rawcounts_STG.rds",
  FRONTAL_MSBB = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FRONTAL/MSBB_RNAseq_rawcounts_counts_FRONTAL_MSBB.rds",
  CRB_Mayo = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_rawcounts_CRB.rds",
  TC_Mayo = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_rawcounts_TC.rds"
)

# Leer todos los archivos en una lista
counts <- lapply(files, readRDS)
counts <- lapply(counts, as.data.frame)
names(counts) <- names(files)
lapply(counts, dim)

#Prepare counts --- ---

#Change name of the first column to 

counts <- lapply(counts, function(counts) {
  colnames(counts)[1] <- "feature"
  return(counts)
})

#Trim STAR alignment stats

counts <- lapply(counts, function(x) x[-c(1:4), ])
lapply(counts, dim)

# Delete the "_PAR_Y" elements

counts <- lapply(counts, function(x) x[-grep("_PAR_Y", x$feature), ])
lapply(counts, dim)

#Prepare feature column
counts <- lapply(counts, function(x) {
  x <- x %>% mutate(feature = str_remove(feature, "\\..*$"))
  return(x)
})

#Create sub metadata

metadata <- lapply(metadata, function(x) {
  if ("sequencingBatch" %in% colnames(x)) {
    x <- x[!is.na(x$sequencingBatch), ]
  }
  if ("exclude" %in% colnames(x)) {
    x <- x[x$exclude != TRUE | is.na(x$exclude), ]
  }
  return(x)
})

lapply(metadata, dim)

#Filter counts to have only specimenes with metadata, just to be sure

# Aplicar la función a cada par de matrices de conteos y metadatos
counts <- lapply(seq_along(counts), function(i) {
  filter_by_metadata(counts[[i]], metadata[[i]])
})

lapply(counts, dim)

# Aplicar la función a tu lista de metadata y counts
metadata <- lapply(seq_along(metadata), function(i) {
  filter_metadata_by_counts(metadata = metadata[[i]], counts = counts[[i]])
})

lapply(metadata, dim)

counts <- lapply(seq_along(counts), function(i) {
  filter_by_metadata(counts[[i]], metadata[[i]])
})


# Aplicar la función a cada dataframe en la lista
counts <- lapply(counts, del_dupl)
lapply(counts, dim)

#Set rownames

counts <- lapply(counts, function(x) {
  x <- column_to_rownames(x, var = colnames(x)[1])  # Usa la primera columna como rownames
  return(x)
})

# PCA --- ---

#Transpose matrices 

matrix.l <- map(counts, t) #rows must be = samples and columns = variables
lapply(matrix.l, dim)

#control

class(matrix.l[[1]])
#print
lapply(matrix.l, function(x) x[1:5, 1:5])

# PCA --- ---

matrix.l <- lapply(matrix.l, function(x) {
  x <- x[, apply(x, 2, var) != 0]  # Eliminar columnas con varianza 0
  return(x)
})

pca.l <- lapply(matrix.l, function(x) prcomp(x, retx = TRUE, center = TRUE, scale. = TRUE)) #slow

#PCA to table
pca_df.l <- lapply(pca.l, function(x) x$x %>% as.data.frame() %>%  rownames_to_column(var = 'specimenID'))
names(pca_df.l) <- names(files)

lapply(pca_df.l, function(x) x[1:5, 1:5])

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

var_table.df <- lapply(pca.l, function(pca) {
  # Calculate the variance of each component
  variance_percentage <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  # Calculate the cumulative variance
  cumulative_percentage <- cumsum(variance_percentage)
  # Create a dataframe with the information
  data.frame(
    PC = 1:length(pca$sdev),
    Variance_Percentage = variance_percentage,
    Cumulative_Percentage = cumulative_percentage
  )
})

#Plot results --- ---

# Replace msex with sex
metadata <- lapply(metadata, function(df) {
  names(df) <- ifelse(names(df) == "msex", "sex", names(df))
  return(df)
})

# List names
names_list <- names(pca_df.l)

#Create elbow plots

elbow_plots.l <- lapply(seq_along(var_table.df), function(i) {
  ggplot(var_table.df[[i]], aes(x = PC, y = Variance_Percentage)) +
    geom_bar(stat = 'identity', fill = var_table.df[[i]]$PC, color = 'black', position = 'dodge') +
    labs(title = paste('Scree plot -', names(var_table.df)[i]),
         x = 'Principal Components',
         y = 'Variance percentage') +
    scale_x_discrete() +
    theme_minimal()
})

##### PCA for SEX ####

#PCA scatterplot coloured by msex function 
scatterplot_msex.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$sex)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "Sex") +
    theme_minimal()
}

scatterplot_msex.f(pca_df.l[[4]], var_table = var_table.df[[4]], metadata = metadata[[4]], name = "PHG - MSBB" )

#Plot

pca_plots_msex.l <- lapply(seq_along(pca_df.l), function(i) {
  scatterplot_msex.f(pca_df.l[[i]], var_table.df[[i]], metadata[[i]], names_list[[i]])
})

pca_plots_msex_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_msex.l[1:3], ncol = 3))
pca_plots_msex_MSBB.p <- do.call(grid.arrange, c(pca_plots_msex.l[4:6], ncol = 3))
pca_plots_msex_Mayo.p <- do.call(grid.arrange, c(pca_plots_msex.l[7:9], ncol = 2))

#### PCA for exclude ####

scatterplot_excl.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$exclude)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "Sex") +
    theme_minimal()
}

#Plot

pca_plots_excl.l <-  mapply(scatterplot_excl.f, 
                            pca_df.l, 
                            var_table.df, 
                            metadata, 
                            names_list,
                            SIMPLIFY = FALSE)

pca_plots_excl_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_excl.l[1:3], ncol = 3))
pca_plots_excl_MSBB.p <- do.call(grid.arrange, c(pca_plots_excl.l[4:6], ncol = 3))
pca_plots_excl_Mayo.p <- do.call(grid.arrange, c(pca_plots_excl.l[7:9], ncol = 2))

##### PCA for TISSUE (especially for FRONTAL) ####

scatterplot_tissue.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$tissue)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "Sex") +
    theme_minimal()
}

#Plot

pca_plots_tissue.l <-  mapply(scatterplot_tissue.f, 
                            pca_df.l, 
                            var_table.df, 
                            metadata, 
                            names_list,
                            SIMPLIFY = FALSE)

pca_plots_tissue_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_tissue.l[1:3], ncol = 3))
pca_plots_tissue_MSBB.p <- do.call(grid.arrange, c(pca_plots_tissue.l[4:6], ncol = 3))
pca_plots_tissue_Mayo.p <- do.call(grid.arrange, c(pca_plots_tissue.l[7:9], ncol = 2))

#Save
# Define la carpeta donde se guardarán los archivos
output_dir <- "/datos/rosmap/data_by_counts/AD_TF_sex/PCAs_ori"

# Lista de gráficos y nombres de archivos con la ruta completa
plot_list <- list(pca_plots_msex_ROSMAP.p, pca_plots_msex_MSBB.p, pca_plots_msex_Mayo.p, 
                  pca_plots_excl_ROSMAP.p, pca_plots_excl_MSBB.p,  pca_plots_excl_Mayo.p, 
                  pca_plots_tissue_ROSMAP.p, pca_plots_tissue_MSBB.p, pca_plots_tissue_Mayo.p)
file_names <- c("pca_plots_msex_ROSMAP_FRONTAL.png", "pca_plots_msex_MSBB.png", "pca_plots_msex_Mayo.png", 
                "pca_plots_excl_ROSMAP.png", "pca_plots_excl_MSBB.png",  "pca_plots_excl_Mayo.png", 
                "pca_plots_tissue_ROSMAP.png", "pca_plots_tissue_MSBB.png", "pca_plots_tissue_Mayo.png")
file_paths <- file.path(output_dir, file_names)

# Guardar cada gráfico en un archivo
mapply(ggsave, filename = file_paths, plot = plot_list, width = 25, height = 10, dpi = 300)

#Here something happens, and that is that the Mayo metadata is different and does not have batch sequencing,
#so I will do the scatterplots separately for them.

pca_df.mayo <- pca_df.l[(length(pca_df.l)-1):length(pca_df.l)]
var_table.df.mayo <- var_table.df[(length(var_table.df)-1):length(var_table.df)]
metadata.mayo <- metadata[(length(metadata)-1):length(metadata)]
names_list.mayo <- names_list[(length(names_list)-1):length(names_list)]

#ROSMAP and MSBB 

pca_df.l <- pca_df.l[1:(length(pca_df.l)-2)]
var_table.df <- var_table.df[1:(length(var_table.df)-2)]
metadata <- metadata[1:(length(metadata)-2)]
names_list <-  names_list[1:(length(names_list)-2)]

#PCA scatterplot coloured by library batch function for ROSMAP and MSBB

scatterplot_seqbatch.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$sequencingBatch)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "Sequencing batch") +
    theme_minimal()
}

#PCA scatterplot coloured by NIA Reagan function for ROSMAP and MSBB

scatterplot_NIA_Reagan.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$dicho_NIA_reagan)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "NIA Reagan") +
    theme_minimal()
}

#PCA scatterplot coloured by CERAD function for ROSMAP and MSBB

scatterplot_CERAD.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$ceradsc)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),  # Incluye el nombre en el título
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "CERAD score") +
    theme_minimal()
}

#Plot and save

##### PCA for Seq batch ####

pca_plots_seqbatch.l <- mapply(scatterplot_seqbatch.f, 
                               pca_df.l, 
                               var_table.df, 
                               metadata, 
                               names_list,
                               SIMPLIFY = FALSE)

pca_plots_seqbatch_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_seqbatch.l[1:3], ncol = 3))
pca_plots_seqbatch_MSBB.p <- do.call(grid.arrange, c(pca_plots_seqbatch.l[4:6], ncol = 3))

#NIA Reagan
pca_plots_NIAR.l <-  mapply(scatterplot_NIA_Reagan.f, 
                            pca_df.l, 
                            var_table.df, 
                            metadata,
                            names_list,
                            SIMPLIFY = FALSE)

pca_plots_NIAR_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_NIAR.l[1:3], ncol = 3))
pca_plots_NIAR_MSBB.p <- do.call(grid.arrange, c(pca_plots_NIAR.l[4:6], ncol = 3))

##### PCA CERAD ###

pca_plots_cerad.l <-  mapply(scatterplot_CERAD.f, 
                            pca_df.l, 
                            var_table.df, 
                            metadata, 
                            names_list,
                            SIMPLIFY = FALSE)

pca_plots_cerad_ROSMAP.p <- do.call(grid.arrange, c(pca_plots_cerad.l[1:3], ncol = 3))
pca_plots_cerad_MSBB.p <- do.call(grid.arrange, c(pca_plots_cerad.l[4:6], ncol = 3))

#### PCA DIAGNOSIS MAYO ###

scatterplot_diagnosis_mayo.f <- function(pca_df, var_table, metadata, name) {
  ggplot(pca_df) +
    aes(x = PC1, y = PC2, colour = as.factor(metadata$diagnosis)) +
    geom_point() +
    geom_text(mapping = aes(label = specimenID), hjust = 1, vjust = 1) +
    labs(title = paste("PCA Scatterplot -", name),
         subtitle = "PC1 vs PC2",
         x = paste("PC1 (", sprintf("%.2f", var_table$Variance_Percentage[1]), "%)"),
         y = paste("PC2 (", sprintf("%.2f", var_table$Variance_Percentage[2]), "%)"),
         colour = "Diagnosis") +
    theme_minimal()
}

#plot

pca_plots_dx_mayo.l <-  mapply(scatterplot_diagnosis_mayo.f, 
                             pca_df.mayo, 
                             var_table.df.mayo, 
                             metadata.mayo, 
                             names_list.mayo,
                             SIMPLIFY = FALSE)

pca_plots_dx_mayo.p <- do.call(grid.arrange, c(pca_plots_dx_mayo.l, ncol = 2))

#Save

# Define la carpeta donde se guardarán los archivos
output_dir <-  "/datos/rosmap/data_by_counts/AD_TF_sex/PCAs_ori"

# Lista de gráficos y nombres de archivos con la ruta completa
plot_list <- list(pca_plots_seqbatch_ROSMAP.p, pca_plots_seqbatch_MSBB.p, pca_plots_NIAR_ROSMAP.p,
                  pca_plots_NIAR_MSBB.p, pca_plots_cerad_ROSMAP.p, pca_plots_cerad_MSBB.p, pca_plots_dx_mayo.p)
file_names <- c("pca_plots_seqbatch_ROSMAP.png", "pca_plots_seqbatch_MSBB.png", "pca_plots_NIAR_ROSMAP.png",
                "pca_plots_NIAR_MSBB.png", "pca_plots_cerad_ROSMAP.png","pca_plots_cerad_MSBB.png", "pca_plots_dx_mayo.png")
file_paths <- file.path(output_dir, file_names)

# Guardar cada gráfico en un archivo
mapply(ggsave, filename = file_paths, plot = plot_list, width = 25, height = 10, dpi = 300)

#END