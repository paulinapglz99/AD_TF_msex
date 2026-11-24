#Libraries

library(vroom)
library(dplyr)

#Set functions

# Crear una función para renombrar las columnas excepto "ensembl_gene_id"
rename_columns <- function(df, suffix) {
  colnames(df)[-1] <- paste0(colnames(df)[-1], "_", suffix)
  return(df)
}

#Read data

#ROSMAP

FRONTAL_ROSMAP <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/FRONTAL_ROSMAP_IPA_table.txt") %>% as.data.frame()

#HCN

HCN_AD <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/HCN/ROSMAP_HCN_DE_full_gene_list_AD_female_vs_AD_male.txt")
HCN_0 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/HCN/ROSMAP_HCN_DE_full_gene_list_control_female_vs_AD_female.txt")
HCN_1 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/HCN/ROSMAP_HCN_DE_full_gene_list_control_male_vs_AD_male.txt")

#PCC

PCC_AD <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PCC/ROSMAP_PCC_DE_full_gene_list_AD_female_vs_AD_male.txt")
PCC_0 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PCC/ROSMAP_PCC_DE_full_gene_list_AD_female_vs_AD_male.txt")
PCC_1 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PCC/ROSMAP_PCC_DE_full_gene_list_control_male_vs_AD_male.txt")

#MSBB

#PHG

PHG_AD <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PHG_MSBB/PHG_MSBB_DE_full_gene_list_AD_female_vs_AD_male.txt")
PHG_0 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PHG_MSBB/PHG_MSBB_DE_full_gene_list_control_female_vs_AD_female.txt")
PHG_1 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/PHG_MSBB/PHG_MSBB_DE_full_gene_list_control_male_vs_AD_male.txt")

#STG

STG_AD <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/STG_MSBB/STG_MSBB_DE_full_gene_list_AD_female_vs_AD_male.txt")
STG_0 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/STG_MSBB/STG_MSBB_DE_full_gene_list_control_female_vs_AD_female.txt")
STG_1 <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/STG_MSBB/STG_MSBB_DE_full_gene_list_control_male_vs_AD_male.txt")

#Mayo

#CRB
CRB_AD <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/CRB_Mayo/Mayo_CRB_DE_full_gene_list_AD_female_vs_AD_male.txt")
CRB_0 <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/CRB_Mayo/Mayo_CRB_DE_full_gene_list_control_female_vs_AD_female.txt")
CRB_1 <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/CRB_Mayo/Mayo_CRB_DE_full_gene_list_control_male_vs_AD_male.txt")

#TC

TC_AD <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/TC_Mayo/Mayo_TC_DE_full_gene_list_AD_female_vs_AD_male.txt")
TC_0 <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/TC_Mayo/Mayo_TC_DE_full_gene_list_control_female_vs_AD_female.txt")
TC_1 <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/TC_Mayo/Mayo_TC_DE_full_gene_list_control_male_vs_AD_male.txt")

#IPA tables

#ROSMAP

#HCN

# Renombrar las columnas para cada tabla
HCN_AD <- rename_columns(HCN_AD, "HCN_AD")
HCN_0 <- rename_columns(HCN_0, "HCN_0")
HCN_1 <- rename_columns(HCN_1, "HCN_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
HCN_IPA <- HCN_AD %>%
  left_join(HCN_0, by = "ensembl_gene_id") %>%
  left_join(HCN_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
HCN_IPA <- HCN_IPA %>%
  as.data.frame() %>%
  filter(complete.cases(.[, c(grep("^log2FoldChange", names(.)), grep("^padj", names(.)))]))

#Save HCN_ROSMAP

vroom::vroom_write(x = HCN_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/HCN_ROSMAP_IPA_table.txt")

#PCC

# Renombrar las columnas para cada tabla
PCC_AD <- rename_columns(PCC_AD, "PCC_AD")
PCC_0 <- rename_columns(PCC_0, "PCC_0")
PCC_1 <- rename_columns(PCC_1, "PCC_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
PCC_IPA <- PCC_AD %>%
  left_join(PCC_0, by = "ensembl_gene_id") %>%
  left_join(PCC_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
PCC_IPA <- PCC_IPA %>%
  filter(complete.cases(.[, c(grep("^log2FoldChange", names(.)), grep("^padj", names(.)))]))

dim(PCC_IPA)
table(is.na(PCC_IPA$log2FoldChange_PCC_AD))
table(is.na(PCC_IPA$log2FoldChange_PCC_1))
table(is.na(PCC_IPA$log2FoldChange_PCC_0))

# Guardar PCC_ROSMAP
vroom::vroom_write(x = PCC_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/PCC_ROSMAP_IPA_table.txt")


#PHG

# Renombrar las columnas para cada tabla
PHG_AD <- rename_columns(PHG_AD, "PHG_AD")
PHG_0 <- rename_columns(PHG_0, "PHG_0")
PHG_1 <- rename_columns(PHG_1, "PHG_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
PHG_IPA <- PHG_AD %>%
  left_join(PHG_0, by = "ensembl_gene_id") %>%
  left_join(PHG_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
PHG_IPA <- PHG_IPA %>%
  filter(complete.cases(select(., starts_with("log2FoldChange"), starts_with("padj"))))

# Guardar PHG_MSBB
vroom::vroom_write(x = PHG_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/PHG_MSBB_IPA_table.txt")

#STG

# Renombrar las columnas para cada tabla
STG_AD <- rename_columns(STG_AD, "STG_AD")
STG_0 <- rename_columns(STG_0, "STG_0")
STG_1 <- rename_columns(STG_1, "STG_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
STG_IPA <- STG_AD %>%
  left_join(STG_0, by = "ensembl_gene_id") %>%
  left_join(STG_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
STG_IPA <- STG_IPA %>%
  filter(complete.cases(select(., starts_with("log2FoldChange"), starts_with("padj"))))

# Guardar STG_MSBB
vroom::vroom_write(x = STG_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/STG_MSBB_IPA_table.txt")

#CRB

# Renombrar las columnas para cada tabla
CRB_AD <- rename_columns(CRB_AD, "CRB_AD")
CRB_0 <- rename_columns(CRB_0, "CRB_0")
CRB_1 <- rename_columns(CRB_1, "CRB_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
CRB_IPA <- CRB_AD %>%
  left_join(CRB_0, by = "ensembl_gene_id") %>%
  left_join(CRB_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
CRB_IPA <- CRB_IPA %>%
  filter(complete.cases(select(., starts_with("log2FoldChange"), starts_with("padj"))))

# Guardar CRB_Mayo
vroom::vroom_write(x = CRB_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/CRB_Mayo_IPA_table.txt")

#TC

# Renombrar las columnas para cada tabla
TC_AD <- rename_columns(TC_AD, "TC_AD")
TC_0 <- rename_columns(TC_0, "TC_0")
TC_1 <- rename_columns(TC_1, "TC_1")

# Realizar los left joins usando "ensembl_gene_id" como clave
TC_IPA <- TC_AD %>%
  left_join(TC_0, by = "ensembl_gene_id") %>%
  left_join(TC_1, by = "ensembl_gene_id")

# Filtrar las filas que no tengan NAs en las columnas que comienzan con log2FoldChange o padj
TC_IPA <- TC_IPA %>%
  filter(complete.cases(select(., starts_with("log2FoldChange"), starts_with("padj"))))

FRONTAL_ROSMAP <- FRONTAL_ROSMAP %>% 
  filter(complete.cases(.[, c(grep("^log2FoldChange", names(.)), grep("^padj", names(.)))]))

vroom::vroom_write(x = FRONTAL_ROSMAP, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/FRONTAL_ROSMAP_IPA_table.txt")

# Guardar TC_Mayo
vroom::vroom_write(x = TC_IPA, file = "/datos/rosmap/data_by_counts/AD_TF_sex/strat_by_sex/DEGS/IPA_tables/TC_Mayo_IPA_table.txt")

#END