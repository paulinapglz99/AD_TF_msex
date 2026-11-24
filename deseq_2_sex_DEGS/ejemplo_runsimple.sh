# 1. Generar datos de prueba
#Rscript generate_test_data.R

# 2. Ejecutar pipeline completo
snakemake --cores 8 --use-conda

# 3. O ejecutar standalone
Rscript deseq2_sex_ad_analysis.R \
  --counts data/counts/DLPFC_counts.rds \
  --metadata data/metadata/DLPFC_metadata.csv \
  --outdir results/DLPFC \
  --brain_area DLPFC \
  --padj 0.05 \
  --lfc 0.5