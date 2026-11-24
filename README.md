# AD_TF_msex: diferencias por sexo en el transcriptoma del Alzheimer

## Objetivo

Comparar el panorama transcriptómico de la enfermedad de Alzheimer (AD) entre hombres y mujeres
en cinco regiones del cerebro, usando RNA-seq *bulk* post mortem de los estudios ROSMAP y Mayo
(AMP-AD RNAseq Harmonization Study, Synapse `syn26720675`).

| Región | Estudio | Tejido (columna `tissue` de la metadata) |
|---|---|---|
| DLPFC | ROSMAP | dorsolateral prefrontal cortex |
| HCN | ROSMAP | Head of caudate nucleus |
| PCC | ROSMAP | posterior cingulate cortex |
| CRB | Mayo | cerebellum |
| TC | Mayo | temporal cortex |

## Pipeline

1. **Descarga** de las matrices de conteos y la metadata desde Synapse.
2. **Gráficas de la metadata** descargada.
3. **Separación de las matrices** entre hombres y mujeres, para cada una de las 5 regiones.
4. **Expresión diferencial** entre hombres con AD y mujeres con AD, para identificar los DEGs
   que distinguen la enfermedad entre sexos, con sus gráficas.
5. **Control de calidad** de las matrices de RNA-seq.
6. **Redes de co-expresión con ARACNe-AP** para cada matriz que pasó el control de calidad.
7. **Análisis de topología** de las redes.

### Estado

| Paso | Regla de Snakemake | Scripts en `workflow/scripts/` |
|---|---|---|
| 1. Descarga | `download_synapse` ✅ | — |
| 2. Gráficas de metadata | `plot_metadata` ✅ | `0.1plot_metadatas.R` |
| 3. Matriz de cada región (AD/CT, una muestra por individuo) | `prepare_counts` ✅ | `prepare_counts.R` |
| 4. Expresión diferencial | `pca_pre`, `plot_pca_pre`, `deseq2`, `gsea` ✅ | `pca_pre.R`, `plot_pca.R`, `deseq2.R`, `gsea.R` |
| 5. Control de calidad (muestras atípicas, matrices para redes, separación por grupo) | `gene_annotation`, `gc_normalize`, `sample_qc`, `network_qc`, `plot_pca_post`, `plot_sample_qc`, `split_groups` ✅ | `sample_qc.R`, `network_qc.R`, `split_groups.R` |
| 6. Redes ARACNe-AP (todos los genes como reguladores) | `install_aracne`, `aracne_regulators`, `aracne_threshold`, `aracne_bootstrap`, `aracne_consolidate` ✅ | — |
| 7. Análisis de redes | `network_metrics`, `community_enrichment`, `network_comparison`, `gene_rewiring` ✅ | `network_metrics.R`, `community_enrichment.R`, `network_comparison.R`, `gene_rewiring.R` |

Los scripts de `workflow/scripts/` que no aparecen en esta tabla son análisis anteriores (módulos y genes
exclusivos por sexo, intersección de DEGs, gráficas globales, preprocesamiento viejo); todavía leen y
escriben rutas fijas de servidor (`/STORAGE`, `/datos`) o usan argumentos de línea de comandos.
`deseq_2_sex_DEGS/` es un análisis de expresión diferencial de Mayo hecho por separado, con su propio
Snakefile. Los scripts viejos de expresión diferencial se borraron; siguen en el historial de git.

## Datos

Los archivos a descargar están en [`config/synapse_files.tsv`](config/synapse_files.tsv), con su
ID y versión de Synapse fijados:

- `resources/data/`: matrices de conteos crudos (`*_gene_all_counts_matrix_clean.txt`),
  4 lotes de ROSMAP y 1 de Mayo.
- `resources/metadata/`: metadata combinada de cada estudio (`RNAseq_Harmonization_*_combined_metadata.csv`)
  y covariables de Sage (`*_Covariates_ages_censored.tsv`).

### Diagnóstico

En todo el pipeline se usa la columna `diagnosis` de las covariables de Sage, armonizada entre estudios
([syn26720676](https://www.synapse.org/Synapse:syn26720676)). Se conservan AD y CT (`diagnoses` en
`config/config.yaml`):

- **ROSMAP:** AD = Braak ≥ 4, CERAD ≤ 2 y cogdx = 4; CT = Braak ≤ 3, CERAD ≥ 3 y cogdx = 1. El resto es `OTHER`.
- **Mayo:** solo neuropatología. AD = Braak ≥ 4; CT = Braak ≤ 3. También hay `PSP` y `PATH_AGE`.

En las covariables de ROSMAP, Sage llama **`ACC`** a la cabeza del núcleo caudado (HCN): las 717 muestras
con esa etiqueta son *Head of caudate nucleus* en la metadata original. No es corteza cingulada anterior,
que ROSMAP no tiene. La correspondencia entre nombres está en la columna `covariates_tissue` de
`config/samples.tsv`.

En Mayo, el cerebelo aparece como `_CBE` en los IDs de las covariables y como `_CER` en la matriz de conteos.

## Resultados

- `results/plots/supplementary/metadata_panel.{pdf,jpeg}`: muestras por región, diagnóstico y sexo;
  proporción de mujeres y hombres; edad de muerte, PMI y RIN; y verificación del sexo con XIST y UTY.
- `resources/annotation/ensembl106_genes.tsv`: % de GC y longitud de cada gen, de un *archive* fijo de
  Ensembl (106, abril 2022) con biomaRt.
- `results/gc_normalize/{region}/`: corrección de longitud y luego de GC con EDASeq (dentro de cada muestra,
  *full quantile*) y de tamaño de librería (cuartil superior), guardada como factores de normalización
  que usan `sample_qc`, `deseq2` y `network_qc`. `bias.tsv` tiene el diagnóstico de NOISeq (R² del sesgo
  de GC y de longitud de cada muestra) antes y después.
- `results/sample_qc/{region}/`: muestras atípicas de cada región (conectividad estandarizada z < −3 sobre
  la expresión corregida, y sexo por XIST/UTY distinto al anotado). `sample_qc.tsv` tiene todas las muestras
  con la razón de exclusión; `samples_kept.tsv`, las que usan DESeq2 y las redes.
- `results/network_qc/{region}/expression.rds`: expresión para las redes. Muestras conservadas, genes
  expresados en los 4 grupos, VST y `limma::removeBatchEffect` con las mismas covariables que DESeq2
  (protegiendo diagnóstico × sexo).
- `results/plots/supplementary/pca_post_panel.{pdf,jpeg}`: PCA de cada región después de corregir.
- `results/plots/supplementary/sample_qc_panel.{pdf,jpeg}`: tamaño de librería, genes detectados y
  conectividad de cada muestra (atípicas en rojo), distribución de cada muestra antes y después de corregir,
  y sesgos de GC y de longitud antes y después.
- `results/network_input/{region}/{group}.tsv`: una matriz por grupo (`AD_female`, `AD_male`, `CT_female`,
  `CT_male`), entrada de ARACNe-AP.
- `results/plots/supplementary/pca_pre_panel.{pdf,jpeg}`: PCA de cada región antes de corregir, y
  varianza de cada PC explicada por cada covariable. Con él se eligieron las covariables de DESeq2.
- `results/deseq2/{region}/`: un modelo de DESeq2 por región,
  `~ model_batch + RIN + pct_intronic + pmi + age_death + group` (sin `model_batch` en Mayo), con
  `group` = diagnóstico × sexo. Una tabla por contraste (`AD_vs_CT_female`, `AD_vs_CT_male`,
  `sex_diagnosis_interaction`) y `summary.tsv` con los DEGs (padj < 0.05, |log2FC| ≥ log2(1.2)).
  En la interacción solo se prueban los genes expresados en ambos sexos (columna `tested`); los genes
  del cromosoma Y, XIST, etc. quedan con `padj = NA`.
- `results/gsea/{region}/{contrast}.tsv`: GSEA de cada contraste (GO BP y Reactome), con los genes
  ordenados por el estadístico de Wald. Se guardan todas las vías, significativas o no.
- `results/aracne/{region}/{group}/network.txt`: red de ARACNe-AP de cada región × grupo (20 redes),
  con todos los genes como reguladores. No están en `rule all`: se piden explícitamente (ver abajo).
- `results/network_analysis/{region}/{group}/`: topología de cada red (`metrics.tsv`), grado, PageRank y
  comunidad Leiden de cada gen (`nodes.tsv`), y enriquecimiento de GO BP y Reactome de cada comunidad.
- `results/network_analysis/{region}/comparison_pairs.tsv`: las 4 redes de la región comparadas por pares
  (Jaccard de aristas y de nodos, Spearman del grado, NMI de comunidades), y `convergence_index.tsv`:
  similitud entre sexos con AD menos en controles (> 0 convergencia, < 0 divergencia).
- `results/network_analysis/{region}/gene_rewiring.tsv`: *rewiring* de cada gen por AD en mujeres y en
  hombres (ajustado por grado y expresión), y `rewiring_vs_de.tsv`: si los genes recableados son los DEGs.
- `results/network_analysis/{region}/rewiring_enrichment.tsv`: GO BP y Reactome de los genes recableados
  solo en mujeres, solo en hombres y en ambos (universo: genes que pasan el filtro de grado en los dos sexos).
- `results/plots/main/networks_ad_by_sex.{pdf,jpeg}`: figura principal de redes. (a) *rewiring* en mujeres
  contra hombres, (b) recableados × DEGs, (c) función de los recableados en un solo sexo, (d) índice de
  convergencia/divergencia, (e) Jaccard de aristas AD vs CT en cada sexo.
- `results/plots/supplementary/network_topology_panel.{pdf,jpeg}`: topología de las 20 redes (métricas,
  distribución de grado y aristas contra número de muestras).
- `results/plots/supplementary/networks_panel.{pdf,jpeg}`: los 6 pares de cada región con todas las medidas y
  sensibilidad de recableados × DEGs (top 5%, 10% y 20%).
- `results/plots/supplementary/community_enrichment_panel.{pdf,jpeg}`: función principal de las comunidades
  más grandes de cada red.

**Advertencia:** cada red tiene su propio número de muestras (por ejemplo, 212 mujeres AD contra 77 hombres AD
en DLPFC). Con más muestras ARACNe encuentra más aristas, así que las comparaciones entre redes mezclan biología
con tamaño de muestra.

Todas las gráficas toman colores, fuentes, tamaños y orden de los grupos de
[`workflow/scripts/R/theme_pub.R`](workflow/scripts/R/theme_pub.R). Para cambiar el estilo de todas las
figuras se edita solo ese archivo. Como es un `input` de cada regla de gráficas, Snakemake rehace las
gráficas cuando cambia.

Se necesita una cuenta de [Synapse](https://www.synapse.org) con acceso a los datos de AMP-AD y un
*personal access token* en la variable de ambiente `SYNAPSE_AUTH_TOKEN`. El token nunca se escribe
en el repositorio.

## Uso

```bash
# 1. Entorno que ejecuta Snakemake (una sola vez)
mamba env create -f workflow/envs/ad-tf-msex-smk.yaml
conda activate ad-tf-msex-smk

# 2. Token de Synapse
export SYNAPSE_AUTH_TOKEN=<tu token>

# 3. Revisar qué se va a correr, y correr
snakemake -n
snakemake --cores 4
```

Las redes, su análisis y sus figuras se piden aparte. Todo junto, con 1 *bootstrap* por red (prueba) en una
computadora de 15 GB:

```bash
snakemake --cores 4 all \
  results/plots/main/networks_ad_by_sex.pdf \
  results/plots/supplementary/networks_panel.pdf \
  results/plots/supplementary/community_enrichment_panel.pdf \
  --config 'aracne={"bootstrap": {"n_boot": 1}}' \
  --set-resources aracne_bootstrap:mem_gb=12 aracne_threshold:mem_gb=8 aracne_consolidate:mem_gb=4 \
  --resources mem_gb=13 --keep-going
```

Los objetivos van antes de `--config`, que acepta varios valores y se traga lo que venga después.

El perfil `workflow/profiles/default/profile.yaml` activa los entornos conda de cada regla
(`--sdm conda`), así que no hace falta pasarlo en la línea de comandos.

## Estructura

```
├── config/
│   ├── config.yaml            # parámetros del análisis
│   ├── samples.tsv            # regiones del cerebro: estudio y tejido
│   └── synapse_files.tsv      # archivos de Synapse a descargar
├── workflow/
│   ├── Snakefile              # punto de entrada: configuración, tablas, rule all y todas las reglas
│   ├── rules/                 # vacía por ahora; para separar reglas en .smk más adelante
│   ├── envs/                  # entornos conda con versiones fijadas
│   ├── scripts/               # scripts de R llamados con script:
│   │   └── R/theme_pub.R      # estilo común de todas las figuras
│   └── profiles/default/      # opciones por defecto de Snakemake
├── resources/                 # datos descargados (no se versionan)
└── results/                   # resultados (no se versionan)
```
