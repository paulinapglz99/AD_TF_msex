#
#Script to plot some metadata exploration

#Libraries --- ---

pacman::p_load("tidyverse",
               "ggplot2", 
               "viridis", 
               "gridExtra", 
               "vroom", 
               "cowplot")

#Functions --- ---

#Filter and summarize data by brain region

summarize_by_tissue <- function(metadata, tissue_name) {
  tissue_data <- metadata %>% filter(tissue == tissue_name)
  cat(paste0("Metadata dim ", tissue_name, ": "), dim(tissue_data), "\n")
  return(tissue_data)
}

##################################### ROSMAP ##################################### 

#Read metadata --- ---

metadata_ROSMAP <- vroom::vroom(file = '/STORAGE/genut/alzheimer/data_by_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv')
dim(metadata_ROSMAP)
#[1] 3400   38

#Filter only regions of interest

colnames(metadata_ROSMAP)

table(metadata_ROSMAP$tissue)
table(metadata_ROSMAP$exclude)

metadata_ROSMAP <- metadata_ROSMAP %>% 
  filter(exclude != TRUE | is.na(exclude)) %>% 
  filter(tissue %in% c("dorsolateral prefrontal cortex", "Head of caudate nucleus", "posterior cingulate cortex")) %>% 
  filter(ceradsc %in% c("1", "4"))
dim(metadata_ROSMAP)
metadata_ROSMAP$age_death <- as.numeric(gsub("\\+", "", metadata_ROSMAP$age_death))  #no +

#Violin plots for

counts <- metadata_ROSMAP %>%
  group_by(ceradsc) %>%
  summarise(n = n())

pmi_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = pmi, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("PMI") +
  theme_minimal() +
  theme(legend.position = "none")

age_d_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = age_death, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("Age of death") +
  theme_minimal() +
  theme(legend.position = "none")

rin_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = RIN, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("RIN") +
  theme_minimal() +
  theme(legend.position = "none")

grid_rosmap <- cowplot::plot_grid(pmi_rosmap.p, age_d_rosmap.p, rin_rosmap.p, ncol = 3, labels = c("a", "b", "b"))

RINvsPMI_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = RIN, y = pmi)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  xlab("RIN") +
  ylab("PMI") +
  theme_minimal()

##################################### Mayo Clinic ##################################### 

metadata_Mayo <- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/metadata/RNAseq_Harmonization_Mayo_combined_metadata.csv")
