
DLPFC <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_DLPFC_DEGS/DLPFC_ROSMAP_IPA_table.txt")
HCN =vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_HCN_DEGS/HCN_ROSMAP_IPA_table.txt")
PCC = vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_PCC_DEGS/PCC_ROSMAP_IPA_table.txt")
TC = vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/Mayo_TC/Mayo_TC_IPA_table.txt")
CRB = vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/Mayo_CRB/Mayo_CRB_IPA_table.txt")

#
tablas <- list(
  DLPFC = DLPFC,
  HCN   = HCN,
  PCC   = PCC,
  TC    = TC,
  CRB   = CRB
)


# Renombrar columnas excepto la primera (ensembl_gene_id)
tablas_renombradas <- lapply(names(tablas), function(nombre) {
  tablas[[nombre]] %>%
    rename_with(~ paste0(.x, "_", nombre), -ensembl_gene_id)
}) %>%
  setNames(names(tablas))

# Unir todas las tablas por 'ensembl_gene_id'
tabla_final <- reduce(tablas_renombradas, full_join, by = "ensembl_gene_id")

vroom::vroom_write(tabla_final, "tabla_IPA_AD_TF_sex.xlsx")


