#Comparisons of DEGS 
#This script compares DEGS within an anatomic region 

#Libraries --- ---
#install.packages("ggvenn")
library(ggvenn)
library(clusterProfiler)
library(org.Hs.eg.db)

#Análisis de Superposición: ¿Qué genes son únicos y cuáles son compartidos?

#Core AD genes:
#These are the genes that appear in both the control_female_vs_AD_female and control_male_vs_AD_male lists. 
#These represent the molecular signature of AD that is independent of sex.

#Genes with specific effect in females: 
#These are genes that appear only in the control_female_vs_AD_female list.
#They could be involved in disease mechanisms exclusive to females.

#Genes with specific effect in males:
##genes that are only in the control_male_vs_AD_male list.

#Sex-dimorphic genes in AD
#These are the most interesting ones. They are the genes on your AD_female_vs_AD_male list.
#Do these genes, which differ between sexes in the disease, also appear on the disease vs. control lists?
#The intersection of this list with the other two will give you the genes whose alteration in AD is also sexually dimorphic.

#Read DEGS

degs_sex <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_DLPFC_DEGS/DLPFC_ROSMAP_DEGS_AD_female_vs_AD_male.txt")
degs_female <-vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_DLPFC_DEGS/DLPFC_ROSMAP_DEGS_control_female_vs_AD_female.txt")
degs_male <- vroom::vroom(file = "/datos/rosmap/data_by_counts/AD_TF_sex/DEGS_new/ROSMAP_DLPFC_DEGS/DLPFC_ROSMAP_DEGS_control_male_vs_AD_male.txt")

DEGS <- list(`AD female vs AD male` = degs_sex$external_gene_name,
             `AD vs Control (Female)` = degs_female$external_gene_name,
             `AD vs Control (Male)` = degs_male$external_gene_name)

#Venn diagram
ggvenn(
  DEGS, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, 
  set_name_size = 4
)

#Functions of female-specific genes.
specific_fem <-setdiff(degs_female$ensembl_gene_id, union(degs_male$ensembl_gene_id, degs_sex))
  
specific_fem <- bitr(
  geneID     = specific_fem,
  fromType   = "ENSEMBL",    
  toType     = "ENTREZID",   
  OrgDb      = org.Hs.eg.db  
)

fem_enr <- enrichKEGG(
  gene = specific_fem$ENTREZID,
  organism = 'hsa', #Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

#Dotplot
cnetplot(fem_enr, categorySize="geneNum", showCategory = 10)

# Prioridad 2: Analiza los 40 genes específicos de machos y compara sus funciones con las del grupo de hembras para identificar las diferencias biológicas clave. 🧠♂️

specific_masc <-


#AD Core genes
genes_core_ad <- intersect(degs_female$ensembl_gene_id, degs_male$ensembl_gene_id)
genes_core_ad <- bitr(
  geneID     = genes_core_ad,
  fromType   = "ENSEMBL",    
  toType     = "ENTREZID",   
  OrgDb      = org.Hs.eg.db  
)

#Enrichment with KEGG
core_enr <- enrichKEGG(
  gene = genes_core_ad$ENTREZID,
  organism = 'hsa', #Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

#Dotplot
dotplot(core_enr, showCategory = 15)

#All columns identifier 

# Prioridad 3: Analiza los 25 genes "Core" para entender la patología compartida y usa el scatter plot de log2FC para buscar diferencias cuantitativas.

#Genes "Core" de AD, crea un gráfico de dispersión (scatter plot)

degs_female_core <- degs_female %>% 
  filter(ensembl_gene_id %in% genes_core_ad$ENSEMBL)

degs_male_core <- degs_male %>% 
  filter(ensembl_gene_id %in% genes_core_ad$ENSEMBL)

core_join <- inner_join(degs_female_core, degs_male_core, by = "ensembl_gene_id", suffix = c("_fem", "_male"))

#Scatter plot
ggplot(core_join, aes(x = log2FoldChange_male, y = log2FoldChange_fem)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + # Línea de identidad (y=x)
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(
    title = "Comparación del Fold Change en genes 'Core' de AD",
    x = "log2FC (Control Male vs AD Male)",
    y = "log2FC (Control Female vs AD Female)"
  ) +
  theme_bw()

# Puntos en la diagonal (y=x): El gen cambia de manera muy similar en ambos sexos.
# Puntos por encima/debajo de la diagonal: El cambio es más pronunciado en un sexo que en el otro.
# Puntos en cuadrantes opuestos (ej. arriba-izquierda): 
#¡Estos son los más interesantes! Indican una regulación antagónica: el gen está sobre-expresado en un sexo y sub-expresado en el otro, a pesar de ser significativo en ambos.

# Prioridad 4: Investiga individualmente los 10 genes de dimorfismo sexual para entender la naturaleza exacta de su diferencia.





