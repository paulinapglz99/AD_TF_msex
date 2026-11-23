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
| 3. Separación por región y sexo | pendiente | `0.handle_data_and_metadata.R`, `3.stratify_by_sex.R` |
| 4. Expresión diferencial | pendiente | `DEGS_local.R` (ROSMAP), `differential_exp_msex_Mayo.R` (Mayo), `make_IPA_tables.R` |
| 5. Control de calidad | pendiente | `1.QC_pre_analysis_parallel.R`, `2.pre-promRNA_parallel.R` |
| 6. Redes ARACNe-AP | pendiente | — |
| 7. Topología de redes | pendiente | — |

`diff_exp_rosmap.R`, `differential_exp_msex.R`, `differential_exp_msex_REDO.R` y
`differential_exp_msex_MSBB.R` son versiones anteriores del análisis de expresión diferencial.
Los scripts todavía leen y escriben rutas fijas de servidor (`/STORAGE`, `/datos`) o usan
argumentos de línea de comandos; se irán adaptando a la directiva `script:` uno por uno.

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
