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

Los pasos están numerados igual que sus scripts (`workflow/scripts/<paso>/`) y sus resultados
(`results/<paso>/`):

0. **Descarga** de las matrices de conteos y la metadata desde Synapse.
1. **Metadata:** gráficas de la metadata descargada.
2. **Preprocesamiento:** matriz de cada región (AD/CT, una muestra por individuo), corrección de GC y
   longitud, PCA y control de calidad de las muestras.
3. **Expresión diferencial** de AD en mujeres y en hombres, e interacción sexo × diagnóstico, con GSEA
   y sus gráficas.
4. **Inferencia de redes:** matrices corregidas por grupo (diagnóstico × sexo) y redes de co-expresión con
   ARACNe-AP, todos los genes como reguladores.
5. **Topología** de cada red: métricas, comunidades y su función.
6. **Comparación de redes** entre mujeres y hombres: similitud, convergencia y recableado de genes.

### Estado

Nombre de los resultados: `results/<paso>/<región>/<experimento>_<cohorte>_<región>_<diagnóstico>_<sexo>.<ext>`,
por ejemplo `results/04_network_inference/CRB/coexprnet_MAYO_CRB_AD_female.tsv`. Cohorte: `ROSMAP` o
`MAYO`. Grupo (`<diagnóstico>_<sexo>`): `AD_female`, `AD_male`, `CT_female` o `CT_male`; `all` cuando el
archivo tiene todos los grupos de la región. En la expresión diferencial, el contraste va en lugar del
grupo: `ADvsCT_female`, `ADvsCT_male` o `sexXdx` (interacción). En la tabla, `<C>_<R>` = `<cohorte>_<región>`.

| Paso | Reglas de Snakemake | Scripts en `workflow/scripts/` | Resultados (ruta) |
|---|---|---|---|
| 00. Descarga | `download_synapse` ✅ | — | `resources/data/*_gene_all_counts_matrix_clean.txt`<br>`resources/metadata/*` |
| 01. Metadata | `plot_metadata` ✅ | `01_metadata/plot_metadata.R` | `results/plots/supplementary/metadata_panel.{pdf,jpeg}` |
| 02. Preprocesamiento | `prepare_counts`, `gene_annotation`, `gc_normalize`, `pca_pre`, `plot_pca_pre`, `sample_qc` ✅ | `02_preprocessing/`: `prepare_counts.R`, `gene_annotation.R`, `gc_normalize.R`, `pca_pre.R`, `plot_pca.R`, `sample_qc.R` | `results/02_preprocessing/<región>/`:<br>`rawcounts_<C>_<R>.rds`, `rawmetadata_<C>_<R>.tsv`<br>`gcFactors_<C>_<R>.rds`, `gcBias_<C>_<R>.tsv`<br>`pcaPreScores_`, `pcaPreVariance_`, `pcaPreCovariates_<C>_<R>.tsv`<br>`sampleQC_<C>_<R>_all.tsv`, `metadata_<C>_<R>_all.tsv`<br>`resources/annotation/ensembl106_genes.tsv`<br>`results/plots/supplementary/pca_pre_panel.{pdf,jpeg}` |
| 03. Expresión diferencial | `deseq2`, `gsea`, `plot_degs_ad_by_sex`, `plot_de_supplementary`, `plot_pathways_by_sex` ✅ | `03_differential_expression/`: `deseq2.R`, `gsea.R`, `plot_degs_ad_by_sex.R`, `plot_de_supplementary.R`, `plot_pathways_by_sex.R` | `results/03_differential_expression/<región>/`:<br>`dds_<C>_<R>.rds`, `degSummary_<C>_<R>.tsv`<br>`degs_<C>_<R>_<contraste>.tsv`, `gsea_<C>_<R>_<contraste>.tsv`<br>`results/plots/main/degs_ad_by_sex.{pdf,jpeg}`<br>`results/plots/supplementary/de_panel.{pdf,jpeg}`, `pathways_by_sex.{pdf,jpeg}` |
| 04. Inferencia de redes | `network_qc`, `plot_pca_post`, `plot_sample_qc`, `split_groups`, `install_aracne`, `aracne_regulators`, `aracne_threshold`, `aracne_bootstrap`, `aracne_consolidate` ✅ | `04_network_inference/`: `network_qc.R`, `plot_sample_qc.R`, `split_groups.R` (ARACNe-AP son comandos del Snakefile) | `results/04_network_inference/<región>/`:<br>`matrizQC_<C>_<R>_all.rds`, `sampleDist_<C>_<R>.tsv`<br>`pcaPostScores_`, `pcaPostVariance_`, `pcaPostCovariates_<C>_<R>.tsv`<br>`matrizQC_<C>_<R>_<grupo>.tsv`, `metadata_<C>_<R>_<grupo>.tsv`<br>`regulators_<C>_<R>.txt`, `coexprnet_<C>_<R>_<grupo>.tsv`<br>`results/plots/supplementary/pca_post_panel.{pdf,jpeg}`, `sample_qc_panel.{pdf,jpeg}` |
| 05. Topología de redes | `network_metrics`, `community_enrichment`, `plot_network_topology`, `plot_community_enrichment` ✅ | `05_network_topology/`: `network_metrics.R`, `community_enrichment.R`, `plot_network_topology.R`, `plot_community_enrichment.R` | `results/05_network_topology/<región>/`:<br>`netMetrics_`, `netNodes_`, `commEnrich_<C>_<R>_<grupo>.tsv`<br>`results/plots/supplementary/network_topology_panel.{pdf,jpeg}`, `community_enrichment_panel.{pdf,jpeg}` |
| 06. Comparación de redes | `network_comparison`, `gene_rewiring`, `rewiring_enrichment`, `plot_networks_by_sex`, `plot_networks_supplementary` ✅ | `06_network_comparison/`: `network_comparison.R`, `gene_rewiring.R`, `rewiring_enrichment.R`, `plot_networks_by_sex.R`, `plot_networks_supplementary.R` | `results/06_network_comparison/<región>/`:<br>`netPairs_`, `convergence_`, `rewiring_`, `rewiringVsDEG_`, `rewiringEnrich_<C>_<R>.tsv`<br>`results/plots/main/networks_ad_by_sex.{pdf,jpeg}`<br>`results/plots/supplementary/networks_panel.{pdf,jpeg}` |

Las funciones que comparten varios scripts están en `workflow/scripts/utils/`. El diseño de cada figura y
una propuesta de pie de figura están en [`documents/plot_design.md`](documents/plot_design.md).

Los análisis anteriores al pipeline (preprocesamiento viejo, módulos y genes exclusivos por sexo,
intersección de DEGs, gráficas globales, tablas de IPA y la carpeta `deseq_2_sex_DEGS/` con el análisis de
expresión diferencial de Mayo) se borraron del repositorio; siguen en el historial de git, por ejemplo en el
commit `fe9b8a4`.

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

Qué contiene cada archivo de la tabla de Estado (`<C>_<R>` = `<cohorte>_<región>`):

**02. Preprocesamiento** (`results/02_preprocessing/<región>/`)

- `rawcounts_<C>_<R>.rds` y `rawmetadata_<C>_<R>.tsv`: conteos y metadata de las muestras AD/CT de la región,
  una por individuo.
- `resources/annotation/ensembl106_genes.tsv`: % de GC y longitud de cada gen, de un *archive* fijo de
  Ensembl (106, abril 2022).
- `gcFactors_<C>_<R>.rds`: corrección de longitud y luego de GC con EDASeq (dentro de cada muestra,
  *full quantile*) y de tamaño de librería (cuartil superior), guardada como factores de normalización
  que usan `sample_qc`, `deseq2` y `network_qc`. `gcBias_<C>_<R>.tsv`: diagnóstico de NOISeq (R² del sesgo
  de GC y de longitud de cada muestra) antes y después.
- `pcaPre*_<C>_<R>.tsv`: PCA de cada región antes de corregir y varianza de cada PC explicada por cada
  covariable (figura `pca_pre_panel`). Con ella se eligieron las covariables de DESeq2.
- `sampleQC_<C>_<R>_all.tsv`: todas las muestras con su razón de exclusión (conectividad estandarizada
  z < −3 sobre la expresión corregida, o sexo por XIST/UTY distinto al anotado).
- `metadata_<C>_<R>_all.tsv`: las muestras que pasan el control de calidad, las que usan DESeq2 y las redes.

**03. Expresión diferencial** (`results/03_differential_expression/<región>/`)

- `dds_<C>_<R>.rds`: un modelo de DESeq2 por región,
  `~ model_batch + RIN + pct_intronic + pmi + age_death + group` (sin `model_batch` en Mayo), con
  `group` = diagnóstico × sexo.
- `degs_<C>_<R>_<contraste>.tsv`: una tabla por contraste (`ADvsCT_female`, `ADvsCT_male`, `sexXdx`), y
  `degSummary_<C>_<R>.tsv` con el número de DEGs (padj < 0.05, |log2FC| ≥ log2(1.2)). En la interacción solo
  se prueban los genes expresados en ambos sexos (columna `tested`); los genes del cromosoma Y, XIST, etc.
  quedan con `padj = NA`.
- `gsea_<C>_<R>_<contraste>.tsv`: GSEA de cada contraste (GO BP y Reactome), con los genes ordenados por el
  estadístico de Wald. Se guardan todas las vías, significativas o no.

**04. Inferencia de redes** (`results/04_network_inference/<región>/`)

- `matrizQC_<C>_<R>_all.rds`: expresión para las redes. Muestras conservadas, genes expresados en los
  4 grupos, VST y `limma::removeBatchEffect` con las mismas covariables que DESeq2 (protegiendo
  diagnóstico × sexo). `pcaPost*_<C>_<R>.tsv` y `sampleDist_<C>_<R>.tsv`: PCA y distribución de cada muestra
  después de corregir (figuras `pca_post_panel` y `sample_qc_panel`).
- `matrizQC_<C>_<R>_<grupo>.tsv` y `metadata_<C>_<R>_<grupo>.tsv`: matriz y metadata de cada grupo; la
  matriz es la entrada de ARACNe-AP.
- `coexprnet_<C>_<R>_<grupo>.tsv`: red de ARACNe-AP de cada región × grupo (20 redes, 150 *bootstraps*),
  con todos los genes como reguladores (`regulators_<C>_<R>.txt`). Columnas `Regulator`, `Target`, `MI`,
  `pvalue`. Los archivos intermedios de ARACNe-AP quedan en `aracne_<C>_<R>_<grupo>/`.

**05. Topología de redes** (`results/05_network_topology/<región>/`)

- `netMetrics_<C>_<R>_<grupo>.tsv`: topología de cada red. `netNodes_`: grado, PageRank y comunidad Leiden de
  cada gen. `commEnrich_`: GO BP y Reactome de cada comunidad.

**06. Comparación de redes** (`results/06_network_comparison/<región>/`)

- `netPairs_<C>_<R>.tsv`: las 4 redes de la región comparadas por pares (Jaccard de aristas y de nodos,
  Spearman del grado, NMI de comunidades). `convergence_<C>_<R>.tsv`: similitud entre sexos con AD menos en
  controles (> 0 convergencia, < 0 divergencia).
- `rewiring_<C>_<R>.tsv`: *rewiring* de cada gen por AD en mujeres y en hombres (ajustado por grado y
  expresión). `rewiringVsDEG_<C>_<R>.tsv`: si los genes recableados son los DEGs del mismo sexo.
- `rewiringEnrich_<C>_<R>.tsv`: GO BP y Reactome de los genes recableados solo en mujeres, solo en hombres y
  en ambos (universo: genes que pasan el filtro de grado en los dos sexos).

**Figuras** (`results/plots/main/` y `results/plots/supplementary/`): el diseño de cada una está en
[`documents/plot_design.md`](documents/plot_design.md).

**Advertencia:** cada red tiene su propio número de muestras (por ejemplo, 212 mujeres AD contra 77 hombres AD
en DLPFC). Con más muestras ARACNe encuentra más aristas, así que las comparaciones entre redes mezclan biología
con tamaño de muestra.

Todas las gráficas toman colores, fuentes, tamaños y orden de los grupos de
[`workflow/scripts/utils/theme_pub.R`](workflow/scripts/utils/theme_pub.R). Para cambiar el estilo de todas las
figuras se edita solo ese archivo. Como es un `input` de cada regla de gráficas, Snakemake rehace las
gráficas cuando cambia.

Se necesita una cuenta de [Synapse](https://www.synapse.org) con acceso a los datos de AMP-AD y un
*personal access token* en la variable de ambiente `SYNAPSE_AUTH_TOKEN`. El token nunca se escribe
en el repositorio.

## Uso

```bash
# 1. Entorno que ejecuta Snakemake (una sola vez)
conda env create -f workflow/envs/ad-tf-msex-smk.yaml
conda activate ad-tf-msex-smk

# 2. Token de Synapse
export SYNAPSE_AUTH_TOKEN=<tu token>

# 3. Entornos conda de todas las reglas, sin correr nada (opcional)
snakemake --conda-create-envs-only --cores 1 all results/plots/main/networks_ad_by_sex.pdf
```

`snakemake all` corre todo hasta las matrices de entrada de ARACNe-AP (pasos 00–04, sin las redes). Las redes,
su análisis y sus figuras se piden aparte, porque ARACNe-AP tarda días:

```bash
# Hasta antes de las redes
snakemake all --workflow-profile workflow/profiles/server

# Todo, con las redes y sus figuras
snakemake --workflow-profile workflow/profiles/server all \
  results/plots/main/networks_ad_by_sex.pdf \
  results/plots/supplementary/networks_panel.pdf \
  results/plots/supplementary/network_topology_panel.pdf \
  results/plots/supplementary/community_enrichment_panel.pdf
```

Agrega `-n` para ver qué se va a correr sin ejecutar nada. Perfiles (`workflow/profiles/`):

- `default`: se usa solo si no se pide otro; activa los entornos conda y la memoria de cada regla.
- `server`: 16 núcleos y 120 GB.
- `laptop`: 4 núcleos y 13 GB, ARACNe-AP con menos memoria.

Para una prueba rápida con 1 *bootstrap* por red (la red queda vacía; solo prueba que corre), agrega al final
`--config 'aracne={"bootstrap": {"n_boot": 1}}'`. Los objetivos van antes de `--config`, que acepta varios
valores y se traga lo que venga después.

## Estructura

```
├── config/
│   ├── config.yaml              # parámetros del análisis
│   ├── samples.tsv              # regiones del cerebro: estudio y tejido
│   └── synapse_files.tsv        # archivos de Synapse a descargar
├── documents/                   # documentos explicativos (diseño de figuras, etc.)
├── workflow/
│   ├── Snakefile                # configuración, rutas de resultados, rule all y todas las reglas
│   ├── envs/                    # entornos conda con versiones fijadas
│   ├── scripts/                 # scripts de R llamados con script:, uno por carpeta de paso
│   │   ├── 01_metadata/
│   │   ├── 02_preprocessing/
│   │   ├── 03_differential_expression/
│   │   ├── 04_network_inference/
│   │   ├── 05_network_topology/
│   │   ├── 06_network_comparison/
│   │   └── utils/               # funciones comunes y estilo de las figuras (theme_pub.R)
│   └── profiles/                # opciones de Snakemake: default, server, laptop
├── resources/                   # datos descargados (no se versionan)
└── results/                     # resultados (no se versionan)
    ├── 02_preprocessing/<región>/
    ├── 03_differential_expression/<región>/
    ├── 04_network_inference/<región>/
    ├── 05_network_topology/<región>/
    ├── 06_network_comparison/<región>/
    └── plots/{main,supplementary}/
```
