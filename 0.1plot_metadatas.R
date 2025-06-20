#
#Script to plot some metadata exploration

#Libraries --- ---

pacman::p_load("tidyverse",
               "ggplot2", 
               "viridis", 
               "gridExtra", 
               "vroom",  'biomaRt',
               "cowplot", 
               "corrplot")

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
  filter(assay == "rnaSeq") %>% 
  filter(tissue %in% c("dorsolateral prefrontal cortex", "Head of caudate nucleus", "posterior cingulate cortex")) %>% 
  filter(ceradsc %in% c("1", "4"))
dim(metadata_ROSMAP)
metadata_ROSMAP$age_death <- as.numeric(gsub("\\+", "", metadata_ROSMAP$age_death))  #no +

#Demographics 

summary <- metadata_ROSMAP %>%
  group_by(ceradsc, msex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(age_death, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Education = paste0(round(mean(educ, na.rm = TRUE), 1), " ± ", round(sd(educ, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Males = as.character(sum(msex == 1, na.rm = TRUE)),
    Females = as.character(sum(msex == 0, na.rm = TRUE))
  ) %>%
  mutate(Sex = ifelse(msex == 1, "Male", "Female"))

#Violin plots for

counts <- metadata_ROSMAP %>%
  group_by(ceradsc) %>%
  summarise(n = n())

pmi_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = pmi, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  #ggtitle(label = "PMI - ROSMAP") +
  xlab("CERAD Score") + 
  ylab("PMI") +
  theme_minimal() +
  theme(legend.position = "none")

age_d_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = age_death, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  #ggtitle(label = "Age of Death - ROSMAP") +
  xlab("CERAD Score") + 
  ylab("Age of death") +
  theme_minimal() +
  theme(legend.position = "none")

rin_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = factor(ceradsc), y = RIN, fill = factor(ceradsc))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("1" = "firebrick", "4" = "navyblue")) +
  #ggtitle(label = "CERAD Score - ROSMAP") +
  xlab("CERAD Score") + 
  ylab("RIN") +
  theme_minimal() +
  theme(legend.position = "none")

#Grid violins
title_rosmap <- ggdraw() + draw_label("ROSMAP- metadata", fontface='bold')
grid_rosmap <- cowplot::plot_grid(pmi_rosmap.p, age_d_rosmap.p, rin_rosmap.p, ncol = 3, labels = c("a", "b", "c"))
grid_rosmap <- plot_grid(title_rosmap, grid_rosmap, 
                       ncol = 1, 
                       rel_heights = c(0.1, 1))

ggsave("grid_metadatas_rosmap.jpeg",
       grid_rosmap,
       device = "jpeg")

#RIN vs PMI

RINvsPMI_rosmap.p <- ggplot(metadata_ROSMAP, aes(x = RIN, y = pmi)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  ggtitle(label = "RIN vs PMI - ROSMAP") +
  xlab("RIN") +
  ylab("PMI") +
  theme_minimal()

ggsave("RINvsPMI_rosmap.jpeg", 
       RINvsPMI_rosmap.p,
       device = "jpeg")

#Relationships between variables

#Only numeric values
metadata_cor_rosmap <- metadata_ROSMAP %>% 
  dplyr::select(where(is.numeric)) %>% 
  dplyr::select(-BrodmannArea, -readLength)

cor_matrix <- cor(metadata_cor_rosmap, use = "pairwise.complete.obs", method = "pearson")

png("corplot_metadata_rosmap.png", width = 1200, height = 1000, res = 150)
corrplot(cor_matrix, 
         method = "ellipse", 
         type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.cex = 0.8)
title(main = "Metadata correlation ROSMAP", cex.main = 1.5)
dev.off()

#Check for balance in covariates

#Continous covariates

numeric_rosmap <- c("age_death", "RIN", "pmi")
categorical_rosmap <- c("msex", "braaksc", "ceradsc", "cogdx","race", "educ")

covariates <- data.frame(Variable = character(),
                         Test = character(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

for (var in numeric_rosmap) {
  x <- metadata_ROSMAP[[var]]
  
  #Test normality (Shapiro-Wilk) (Shapiro needs at least 3 values)
  if (sum(!is.na(x)) >= 3) {
    shapiro_result <- shapiro.test(x)
    
    if (shapiro_result$p.value > 0.05) {
      #Normal: ANOVA
      formula <- as.formula(paste(var, "~ ceradsc"))
      aov_result <- aov(formula, data = metadata_ROSMAP_DLPFC)
      p_val <- summary(aov_result)[[1]][["Pr(>F)"]][1]
      test_used <- "ANOVA"
    } else {#No normal: Kruskal-Wallis
      formula <- as.formula(paste(var, "~ ceradsc"))
      kruskal_result <- kruskal.test(formula, data = metadata_ROSMAP_DLPFC)
      p_val <- kruskal_result$p.value
      test_used <- "Kruskal-Wallis"
    }
    
    covariates <- rbind(covariates, data.frame(Variable = var, Test = test_used, p_value = p_val))
    
  } else {
    cat(paste("Variable", var, "has too little data to perform the normality test.\n"))
  }
}

#Evaluate categorial variables
for (var in categorical_rosmap) {
  tabla <- table(metadata_ROSMAP[[var]], metadata_ROSMAP$ceradsc)
  chisq_result <- chisq.test(tabla)
  covariates <- rbind(covariates, data.frame(Variable = var, Test = "Chi-square", p_value = chisq_result$p.value))
}

#See expression 

counts_ROSMAP1 <- readRDS(file = '/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/data/ROSMAP_RNAseq_rawcounts_DLPFC.rds')%>% 
            as.data.frame() %>%
            mutate(feature = sub("\\..*", "", feature))
counts_ROSMAP2 <- readRDS(file = '/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/data/ROSMAP_RNAseq_rawcounts_HCN.rds')%>% 
  as.data.frame() %>%
  mutate(feature = sub("\\..*", "", feature))
counts_ROSMAP3 <- readRDS(file = '/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/data/ROSMAP_RNAseq_rawcounts_PCC.rds')%>% 
  as.data.frame() %>%
  mutate(feature = sub("\\..*", "", feature))

counts_ROSMAP <- dplyr::left_join(counts_ROSMAP1, counts_ROSMAP2, by = 'feature') %>%
  dplyr::left_join(counts_ROSMAP3, by = 'feature')
counts_ROSMAP <- counts_ROSMAP[-c(1:4), ]

#XIST and UTY genes

XIST_UTY <- counts_ROSMAP %>% 
  filter(str_detect(feature, "^ENSG00000229807|^ENSG00000183878")) %>% 
  column_to_rownames(var = "feature") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "specimenID")%>% 
  left_join(metadata_ROSMAP,  by = "specimenID") %>% 
  filter(!is.na(msex))

colnames(XIST_UTY)[2] <- "UTY"
colnames(XIST_UTY)[3] <- "XIST"

XIST_UTY.p  <-ggplot(XIST_UTY, aes(x = log(UTY), y = log(XIST), color = as.factor(msex))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("0" = "#4B0082", "1" = "#FFA500")) + 
  labs(title = "Expression of XIST and UTY - ROSMAP",
       x = "XIST",
       y = "UTY")

ggsave("XIST_UTY_rosmap.jpeg", 
       XIST_UTY.p,
       device = "jpeg")

#Transcript composition 

library(org.Hs.eg.db)
annot_rosmap <- AnnotationDbi::select(org.Hs.eg.db,
                               keys = unique(counts_ROSMAP$feature),
                               keytype = "ENSEMBL",
                               columns = c("ENSEMBL", "SYMBOL", "GENENAME", "GENETYPE"))

annot_rosmap_n <- as.data.frame(table(annot_rosmap$GENETYPE, useNA = "ifany"))

annot_rosmap_n <- annot_rosmap_n %>% rename(group = Var1, value = Freq)

#Demographics per brain region 

#DLPFC

metadata_ROSMAP_DLPFC <- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
metadata_ROSMAP_DLPFC$age_death <- as.numeric(gsub("\\+", "", metadata_ROSMAP_DLPFC$age_death))  #no +

summary_DLPFC <- metadata_ROSMAP_DLPFC %>%
  group_by(ceradsc, msex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(age_death, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Education = paste0(round(mean(educ, na.rm = TRUE), 1), " ± ", round(sd(educ, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(RIN, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(pmi, na.rm = TRUE), 2)),
    braak = paste(paste0("Braak", names(table(braaksc)), ": ", as.numeric(table(braaksc))), collapse = ", "),
    Cogdx  = paste(paste0("cogdx", names(table(cogdx)), ": ", as.numeric(table(cogdx))), collapse = ", "),  
    APOE = paste(paste0("APOE", names(table(apoe_genotype)), ": ", as.numeric(table(apoe_genotype))), collapse = ", ")) %>% 
  ungroup()

#HCN

metadata_ROSMAP_HCN <-vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN/RNA_seq_metadata_HCN.csv") %>% as.data.frame()
metadata_ROSMAP_HCN$age_death <- as.numeric(gsub("\\+", "", metadata_ROSMAP_HCN$age_death))  #no +

summary_HCN <- metadata_ROSMAP_HCN %>%
  group_by(ceradsc, msex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(age_death, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Education = paste0(round(mean(educ, na.rm = TRUE), 1), " ± ", round(sd(educ, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(RIN, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(pmi, na.rm = TRUE), 2)),
    braak = paste(paste0("Braak", names(table(braaksc)), ": ", as.numeric(table(braaksc))), collapse = ", "),
    Cogdx  = paste(paste0("cogdx", names(table(cogdx)), ": ", as.numeric(table(cogdx))), collapse = ", "),  
    APOE = paste(paste0("APOE", names(table(apoe_genotype)), ": ", as.numeric(table(apoe_genotype))), collapse = ", ")) %>% 
  ungroup()

#PCC

metadata_ROSMAP_PCC<- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.csv") %>% as.data.frame()
metadata_ROSMAP_PCC$age_death <- as.numeric(gsub("\\+", "", metadata_ROSMAP_PCC$age_death))  #no +

summary_PCC <- metadata_ROSMAP_PCC %>%
  group_by(ceradsc, msex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(age_death, na.rm = TRUE), 1), " ± ", round(sd(age_death, na.rm = TRUE), 2)),
    Education = paste0(round(mean(educ, na.rm = TRUE), 1), " ± ", round(sd(educ, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(RIN, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(pmi, na.rm = TRUE), 2)),
    braak = paste(paste0("Braak", names(table(braaksc)), ": ", as.numeric(table(braaksc))), collapse = ", "),
    Cogdx  = paste(paste0("cogdx", names(table(cogdx)), ": ", as.numeric(table(cogdx))), collapse = ", "),  
    APOE = paste(paste0("APOE", names(table(apoe_genotype)), ": ", as.numeric(table(apoe_genotype))), collapse = ", ")) %>% 
  ungroup()

##################################### Mayo Clinic ##################################### 

metadata_Mayo <- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/metadata/RNAseq_Harmonization_Mayo_combined_metadata.csv")

colnames(metadata_Mayo)

table(metadata_Mayo$diagnosis)

metadata_Mayo <- metadata_Mayo %>% 
  filter(exclude != TRUE | is.na(exclude)) %>% 
  #filter(tissue %in% c("", "Head of caudate nucleus", "posterior cingulate cortex")) %>% 
  filter(diagnosis %in% c("Alzheimer Disease", "control"))
dim(metadata_Mayo)
metadata_Mayo$ageDeath <- as.numeric(gsub("\\+", "", metadata_Mayo$ageDeath))  #no +

#Violin plots for

#PMI
pmi_mayo.p <- ggplot(metadata_Mayo, aes(x = factor(diagnosis), y = pmi, fill = factor(diagnosis))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("Alzheimer Disease" = "firebrick", "control" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("PMI") +
  theme_minimal() +
  theme(legend.position = "none")
  
#Age of death
age_d_mayo.p <-  ggplot(metadata_Mayo, aes(x = factor(diagnosis), y = ageDeath, fill = factor(diagnosis))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("Alzheimer Disease" = "firebrick", "control" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("Age of death") +
  theme_minimal() +
  theme(legend.position = "none")

#RIN
rin_mayo.p <-  ggplot(metadata_Mayo, aes(x = factor(diagnosis), y = RIN, fill = factor(diagnosis))) +
  geom_violin(trim = FALSE, color = "black", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("Alzheimer Disease" = "firebrick", "control" = "navyblue")) +
  xlab("CERAD Score") + 
  ylab("RIN") +
  theme_minimal() +
  theme(legend.position = "none")

#Grid violins
title_Mayo <- ggdraw() + draw_label("Mayo - metadata", fontface='bold')
grid_mayo <- plot_grid(pmi_mayo.p, age_d_mayo.p, rin_mayo.p, 
                       ncol = 3, labels = c("a", "b", "c"))
grid_mayo <- plot_grid(title_Mayo, grid_mayo, 
          ncol = 1, 
          rel_heights = c(0.1, 1))  

ggsave("grid_metadata_mayo.jpeg", 
       grid_mayo,
       device = "jpeg")

#RIN vs PMI
RINvsPMI_mayo.p <- ggplot(metadata_Mayo, aes(x = RIN, y = pmi)) +
  geom_point(alpha = 0.7, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "red", se = F) +
  ggtitle(label = "RIN vs PMI - Mayo") +
  xlab("RIN") +
  ylab("PMI") +
  theme_minimal()

ggsave("RINvsPMI_mayo.jpeg", 
       RINvsPMI_mayo.p,
       device = "jpeg")

#Only numeric values
metadata_cor_Mayo <- metadata_Mayo %>% 
  dplyr::select(where(is.numeric))

cor_matrix_mayo <- cor(metadata_cor_Mayo, use = "pairwise.complete.obs", method = "pearson")

png("corplot_metadata_Mayo.png", width = 1200, height = 1000, res = 150)
corrplot(cor_matrix_mayo, 
         method = "ellipse", 
         type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.cex = 0.8)
title(main = "Metadata correlation Mayo", cex.main = 1.5, line = 2)
dev.off()

# Check for balance in covariates

#Continous covariates

numeric_mayo <- c("ageDeath", "RIN", "pmi")

covariates_mayo <- data.frame(Variable = character(),
                         Test = character(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

for (var in numeric_mayo) {
  x <- metadata_Mayo[[var]]
  
  #Test normality (Shapiro-Wilk) (Shapiro needs at least 3 values)
  if (sum(!is.na(x)) >= 3) {
    shapiro_result <- shapiro.test(x)
    
    if (shapiro_result$p.value > 0.05) {
      #Normal: ANOVA
      formula <- as.formula(paste(var, "~ diagnosis"))
      aov_result <- aov(formula, data = metadata_Mayo)
      p_val <- summary(aov_result)[[1]][["Pr(>F)"]][1]
      test_used <- "ANOVA"
    } else {#No normal: Kruskal-Wallis
      formula <- as.formula(paste(var, "~ diagnosis"))
      kruskal_result <- kruskal.test(formula, data = metadata_Mayo)
      p_val <- kruskal_result$p.value
      test_used <- "Kruskal-Wallis"
    }
    
    covariates_mayo <- rbind(covariates_mayo, data.frame(Variable = var, Test = test_used, p_value = p_val))
    
  } else {
    cat(paste("Variable", var, "has too little data to perform the normality test.\n"))
  }
}

#Categorical variables
categorical_mayo <- c("sex", "Braak", "thal", "race")

#Evaluate categorial variables
for (var in categorical_mayo) {
  tabla <- table(metadata_Mayo[[var]], metadata_Mayo$diagnosis)
  chisq_result <- chisq.test(tabla)
  covariates_mayo <- rbind(covariates_mayo, data.frame(Variable = var, Test = "Chi-square", p_value = chisq_result$p.value))
}

#See expression 

counts_Mayo <-  vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/Mayo_gene_all_counts_matrix_clean.txt")

#XIST and UTY genes

XIST_UTY_mayo <- counts_Mayo %>% 
  filter(str_detect(feature, "^ENSG00000229807|^ENSG00000183878")) %>% 
  column_to_rownames(var = "feature") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "specimenID")%>% 
  left_join(metadata_Mayo,  by = "specimenID") %>% 
  filter(!is.na(sex))

colnames(XIST_UTY_mayo)[2] <- "UTY"
colnames(XIST_UTY_mayo)[3] <- "XIST"

XIST_UTY_mayo.p  <-ggplot(XIST_UTY_mayo, aes(x = log(UTY), y = log(XIST), color = as.factor(sex))) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("female" = "#4B0082", "male" = "#FFA500")) + 
  labs(title = "Expression of XIST and UTY - Mayo",
       x = "XIST",
       y = "UTY")

ggsave("XIST_UTY_mayo.jpeg", 
       XIST_UTY_mayo.p,
       device = "jpeg")

#Demographics per brain region 

metadata_Mayo_CRB <- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_metadata_CRB.txt")
metadata_Mayo_CRB$ageDeath <- as.numeric(gsub("\\+", "", metadata_Mayo_CRB$ageDeath))  #no +

summary_mayo_CRB <- metadata_Mayo_CRB %>%
  group_by(diagnosis, sex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(ageDeath, na.rm = TRUE), 1), " ± ", round(sd(ageDeath, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(RIN, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(pmi, na.rm = TRUE), 2)),
    braak = paste(paste0("Braak", names(table(Braak)), ": ", as.numeric(table(Braak))), collapse = ", "),
    Thal = paste(paste0("Thal", names(table(thal)), ": ", as.numeric(table(thal))), collapse = ", "),
    APOE = paste(paste0("APOE", names(table(apoeGenotype)), ": ", as.numeric(table(apoeGenotype))), collapse = ", ")
  ) 

#TC

metadata_Mayo_TC <- vroom::vroom(file = "/STORAGE/genut/alzheimer/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_metadata_TC.txt")
metadata_Mayo_TC$ageDeath <- as.numeric(gsub("\\+", "", metadata_Mayo_TC$ageDeath))  #no +
summary_mayo_TC <- metadata_Mayo_TC %>%
  group_by(diagnosis, sex) %>%
  summarise(
    n = as.character(n()),
    Age_Death = paste0(round(mean(ageDeath, na.rm = TRUE), 1), " ± ", round(sd(ageDeath, na.rm = TRUE), 2)),
    RIN = paste0(round(mean(RIN, na.rm = TRUE), 1), " ± ", round(sd(RIN, na.rm = TRUE), 2)),
    PMI = paste0(round(mean(pmi, na.rm = TRUE), 1), " ± ", round(sd(pmi, na.rm = TRUE), 2)),
    braak = paste(paste0("Braak", names(table(Braak)), ": ", as.numeric(table(Braak))), collapse = ", "),
    Thal = paste(paste0("Thal", names(table(thal)), ": ", as.numeric(table(thal))), collapse = ", "),
    APOE = paste(paste0("APOE", names(table(apoeGenotype)), ": ", as.numeric(table(apoeGenotype))), collapse = ", ")
  ) 



#END