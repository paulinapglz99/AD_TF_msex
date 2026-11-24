# Diseño de las figuras

Este documento explica qué muestra cada figura del pipeline, cómo está construido cada panel y propone un pie de figura
para el manuscrito. Los pies están en inglés, como el texto de las figuras, y no dependen de los resultados: describen
qué se grafica y cómo leerlo, no qué se encontró.

Todas las figuras se guardan en `results/plots/main/` (figuras principales) o `results/plots/supplementary/`
(suplementarias), en PDF (vectorial) y JPEG a 600 dpi, al ancho de doble columna. Los colores, tipografía, tamaños y el
orden de regiones y grupos vienen de [workflow/scripts/R/theme_pub.R](workflow/scripts/R/theme_pub.R).

## Convenciones comunes

| Elemento | Codificación |
|---|---|
| Sexo | Color: mujeres verde (`#00A087`), hombres naranja (`#D9772B`) |
| Diagnóstico, cuando el color ya es el sexo | Relleno: control = círculo vacío, AD = círculo lleno |
| Diagnóstico, cuando el sexo no está en el panel | Color: control azul, AD rojo |
| Efecto con signo (log2FC, NES) | Divergente azul (baja) – blanco – rojo (sube) |
| Convergencia / divergencia de los sexos | Divergente café (divergen) – blanco – morado (convergen) |
| Similitud entre redes (Jaccard, ρ, NMI) | Escala de grises, con el valor escrito en cada celda |
| Orden de regiones | DLPFC, HCN, PCC, TC, CRB (`config/samples.tsv`) |
| Orden de grupos | F-CT, F-AD, M-CT, M-AD |

Reglas de composición:

- **Una letra por gráfica.** Las facetas por región comparten letra; cada medida o pregunta distinta lleva la suya.
- **Cada leyenda aparece una sola vez** en la figura. Si dos paneles usan los mismos colores, solo uno muestra la leyenda.
- **Ningún color significa dos cosas** en la misma figura. Por eso la similitud va en grises y la convergencia en
  morado/café: rojo, azul, verde y naranja ya son diagnóstico y sexo.
- **El n desigual se muestra, no se esconde.** Las redes y los DEGs dependen del número de muestras, y las mujeres
  con AD son 2–3 veces más que los hombres con AD. Los paneles que comparan sexos escriben el n (AD/CT) junto a cada
  valor.

---

## Figuras principales

### Figura 1. Expresión diferencial de AD en mujeres y hombres (`main/degs_ad_by_sex`)

Script: [plot_degs_ad_by_sex.R](workflow/scripts/plot_degs_ad_by_sex.R). Pregunta: ¿AD cambia la expresión igual en
mujeres y en hombres?

| Panel | Qué muestra | Diseño |
|---|---|---|
| a | Número de DEGs de AD vs CT por sexo y región | Barras por sexo; suben arriba de 0, bajan abajo; n (AD/CT) bajo cada barra |
| b | log2FC de AD en mujeres vs en hombres, cada gen | Dispersión con facetas por región; color = sexo en que el gen es DEG; diagonal = mismo efecto; ρ de Spearman |
| c | NES de AD en mujeres vs en hombres, cada vía (GSEA) | Dispersión con facetas por región; coloreadas las vías significativas en la interacción sexo × diagnóstico |

> **Pendiente de las convenciones nuevas:** en c, "efecto mayor en mujeres / en hombres" está en rojo/azul. Debería
> ir en verde/naranja (sexo).

**Pie propuesto.**
**Figure 1. Effect of Alzheimer's disease on gene expression in women and men across five brain regions.**
(a) Number of differentially expressed genes (DEGs; DESeq2, BH-adjusted p < 0.05, |log2FC| ≥ log2(1.2)) between AD and
control in women and in men. Bars above zero are up-regulated genes, below zero down-regulated. Numbers under each bar
are the AD/control samples of that sex. Women had more AD samples, which alone increases the number of DEGs.
(b) log2 fold change of AD in women versus in men for every tested gene, per region. Colour indicates the sex in which
the gene is a DEG. The dashed line marks equal effects in both sexes; ρ is the Spearman correlation.
(c) Normalized enrichment score (NES) of AD in women versus in men for every pathway (GSEA, genes ranked by the Wald
statistic). Highlighted pathways are significant for the sex × diagnosis interaction. Their names are in
Supplementary Figure S6.
DLPFC, dorsolateral prefrontal cortex; HCN, head of caudate nucleus; PCC, posterior cingulate cortex; TC, temporal
cortex; CRB, cerebellum.

### Figura 2. Recableado de las redes de co-expresión por AD en mujeres y hombres (`main/networks_ad_by_sex`)

Script: [plot_networks_by_sex.R](workflow/scripts/plot_networks_by_sex.R). Pregunta: ¿AD reorganiza igual las redes de
mujeres y hombres? Va de lo global (la red) al gen y luego a la función.

| Panel | Qué muestra | Diseño |
|---|---|---|
| a | Cuánto recablea AD a cada sexo: Jaccard de aristas AD vs CT | Paletas; en cada región mujeres a la izquierda y hombres a la derecha; n (AD/CT) vertical bajo cada paleta. Más bajo = más recableo |
| b | Si los sexos se parecen más o menos con AD: índice de convergencia | Heatmap región × 4 medidas (Jaccard de aristas y nodos, ρ de grado, NMI de comunidades); valor escrito; morado = convergen, café = divergen |
| c | Recableado ajustado de cada gen en mujeres vs en hombres | Dispersión con facetas por región; color = sexo en que el gen está en el top 10% (ninguno, solo mujeres, solo hombres, ambos); diagonal; ρ |
| d | Si los genes recableados son los DEGs del mismo sexo | Paletas de log2 odds ratio por sexo y región; * = padj (BH) < 0.05 de la prueba de permutación |
| e | Función de los genes recableados en un solo sexo | Puntos: términos (Reactome) × región, facetas "solo mujeres" / "solo hombres"; tamaño = −log10(padj). Baja el número de términos si no caben |

Leyendas: Sexo (en a; d usa los mismos colores), convergencia (b), "Rewired by AD in" (c) y −log10(padj) (e).

**Pie propuesto.**
**Figure 2. Rewiring of the gene co-expression networks by Alzheimer's disease in women and men.**
Networks were inferred with ARACNe-AP (p < 10⁻⁸, 100 bootstraps, all genes as candidate regulators), one for each
region, sex and diagnosis.
(a) Edge Jaccard index between the AD and control networks of each sex. Lower values mean more rewiring. Numbers under
each point are the AD/control samples. Networks with more samples have more edges, so differences between sexes partly
reflect sample size (Supplementary Figure S7h).
(b) Convergence index of the two sexes: similarity between the female and male networks with AD minus their similarity
in controls, for four measures (edge and node Jaccard, Spearman ρ of node degree, normalized mutual information of
Leiden communities). Purple, the sexes become more alike with AD; brown, less alike.
(c) Rewiring of each gene by AD in women versus in men. Rewiring is 1 − Jaccard index of the gene's neighbours in the AD
and control networks, adjusted for mean degree and expression (residuals of a linear model). Only genes with degree ≥ 5
in at least one network of both sexes are shown. Colour indicates the sex in which the gene is among the 10% most
rewired. The dashed line marks equal rewiring; ρ is the Spearman correlation.
(d) Overlap between the 10% most rewired genes and the DEGs of the same sex (log2 odds ratio). Significance from 1,000
permutations stratified by expression and degree; *, BH-adjusted p < 0.05. Supplementary Figure S8e shows the 5% and
20% thresholds.
(e) Reactome pathways over-represented among the genes rewired in women only and in men only (universe: genes tested
in both sexes; BH-adjusted p < 0.05). Point size, −log10(adjusted p). The most significant terms per region and sex
are shown.

---

## Figuras suplementarias

### S1. Metadata de las muestras (`supplementary/metadata_panel`)

Script: [0.1plot_metadatas.R](workflow/scripts/0.1plot_metadatas.R).

| Panel | Qué muestra |
|---|---|
| a | Muestras por región, diagnóstico y sexo |
| b | Proporción de mujeres y hombres por región y diagnóstico |
| c | Edad de muerte, PMI y RIN por región, diagnóstico y sexo |
| d | Verificación del sexo: expresión de XIST vs UTY en cada muestra |

> **Pendiente:** c junta tres variables bajo una letra. Con la regla de una letra por medida serían c, d y e (y la
> verificación de sexo pasaría a f).

**Pie propuesto.**
**Supplementary Figure S1. Samples of the five brain regions.**
(a) Number of samples per region, diagnosis and sex. ROSMAP: DLPFC, HCN, PCC. Mayo Clinic: TC, CRB.
(b) Proportion of women and men per region and diagnosis.
(c) Age at death, post-mortem interval (PMI) and RNA integrity number (RIN) per region, diagnosis and sex. Ages above 89
are censored at 90.
(d) Sex check: XIST versus UTY expression (CPM) in every sample.
AD and control follow the harmonized diagnosis of the AMP-AD consortium. In ROSMAP, AD requires a clinical diagnosis of
AD (cogdx 4), Braak ≥ 4 and CERAD ≤ 2. Controls require no cognitive impairment (cogdx 1), Braak ≤ 3 and CERAD ≥ 3.

### S2. Control de calidad de las muestras (`supplementary/sample_qc_panel`)

Script: [plot_sample_qc.R](workflow/scripts/plot_sample_qc.R).

| Panel | Qué muestra |
|---|---|
| a | Tamaño de librería (millones de lecturas en genes) de cada muestra |
| b | Genes detectados (CPM > 1) de cada muestra |
| c | Conectividad estandarizada; bajo la línea punteada (z < −3) la muestra es atípica y se quita, con su ID |
| d | Distribución de la expresión de cada muestra antes y después de quitar las covariables técnicas (línea = 5–95%, línea gruesa = 25–75%, punto = mediana) |
| e | Sesgos de GC y longitud (R², NOISeq) antes y después de la corrección EDASeq |

Las muestras atípicas van en rojo en a–c.

> **Pendiente:** el rojo de las atípicas choca con el rojo de AD. Conviene negro o un contorno.

**Pie propuesto.**
**Supplementary Figure S2. Quality control of the samples.**
(a) Library size (millions of reads assigned to genes) and (b) detected genes (CPM > 1) of every sample.
(c) Standardized connectivity of every sample in its region (Oldham et al., 2012; adjacency ((1 + r)/2)²). Samples below
z = −3 (dashed line) were removed as outliers and are labelled.
(d) Distribution of gene expression in each kept sample, ordered by median, before and after removing technical
covariates (limma removeBatchEffect, protecting the diagnosis × sex design). Thin line, 5th–95th percentile; thick line,
interquartile range; point, median.
(e) Variance of expression explained by GC content and by gene length (R², NOISeq) before and after within-lane
(length, then GC; full quantile) and between-lane (upper quartile) normalization with EDASeq.

### S3 y S4. PCA antes y después del control de calidad (`supplementary/pca_pre_panel`, `supplementary/pca_post_panel`)

Script: [plot_pca.R](workflow/scripts/plot_pca.R). El mismo diseño para las dos figuras: una fila por región.

| Panel | Qué muestra |
|---|---|
| a–e | PC1 vs PC2 de cada región; color = diagnóstico, forma = sexo, una elipse por sexo |
| f–j | Varianza de cada PC explicada por cada covariable (R² de PC ~ covariable), escrita donde padj (BH) < 0.05 |

Leyendas al pie de la figura.

> **Pendiente:** la forma por sexo (círculo/triángulo) es la misma que se cambió en la topología por ser difícil de
> leer. Se podría pasar a color = sexo y relleno = diagnóstico.

**Pie propuesto (S3).**
**Supplementary Figure S3. Principal component analysis before quality control.**
(a–e) First two principal components of the 1,000 most variable genes (variance-stabilized counts) of each region.
Colour, diagnosis; shape, sex; ellipses, 95% normal ellipse of each sex.
(f–j) Variance of each of the first 10 principal components explained by each covariate (R² of a linear model PC ~
covariate). Values are written where the BH-adjusted p < 0.05.

**Pie propuesto (S4).** Igual que S3, con el título *"Principal component analysis of the expression used for the
networks"* y esta frase añadida: *"Expression after removal of outlier samples, GC and length correction and removal of
technical covariates."*

### S5. Expresión diferencial, detalle (`supplementary/de_panel`)

Script: [plot_de_supplementary.R](workflow/scripts/plot_de_supplementary.R).

| Panel | Qué muestra |
|---|---|
| a | Volcanes de AD vs CT, mujeres (fila de arriba) y hombres (fila de abajo), cada región; los 3 DEGs más significativos con nombre |
| b | Volcanes de la interacción sexo × diagnóstico; genes expresados en un solo sexo (chrY, XIST…) en negro, no probados |
| c | Correlación de Spearman del log2FC de AD entre cada región y sexo (heatmap) |
| d | Vías principales de la interacción (GSEA) por región; color = NES |

**Pie propuesto.**
**Supplementary Figure S5. Differential expression by sex and sex × diagnosis interaction.**
(a) Volcano plots of AD versus control in women (top) and men (bottom) in every region. The three most significant DEGs
of each panel are named.
(b) Volcano plots of the sex × diagnosis interaction. Positive log2FC means the effect of AD is larger (more positive)
in women than in men. Genes expressed in only one sex were not tested for the interaction and are shown in black.
(c) Spearman correlation of the AD log2 fold changes between every region and sex.
(d) Most significant pathways of the interaction (GSEA) in each region. Positive NES, the pathway increases more (or
decreases less) with AD in women than in men.
DESeq2 model per region: ~ batch + RIN + % intronic reads + PMI + age at death + group (diagnosis × sex). Significance:
BH-adjusted p < 0.05 and |log2FC| ≥ log2(1.2).

### S6. Vías de AD por sexo (`supplementary/pathways_by_sex`)

Script: [plot_pathways_by_sex.R](workflow/scripts/plot_pathways_by_sex.R).

| Panel | Qué muestra |
|---|---|
| a | Reactome: para cada región y sexo, las vías significativas con mayor \|NES\|, mostradas en todas las regiones y sexos |
| b | GO procesos biológicos, mismo diseño |

Color = NES, tamaño = −log10(padj), transparente = no significativa.

**Pie propuesto.**
**Supplementary Figure S6. Pathways changed by AD in women and in men.**
Gene set enrichment analysis (clusterProfiler; genes ranked by the DESeq2 Wald statistic; gene sets of 15–500 genes) of
AD versus control in each sex and region, for (a) Reactome and (b) GO biological processes. For each region and sex the
most significant pathways with the largest |NES| are selected and shown in every column. Colour, NES (red, up in AD;
blue, down); size, −log10(BH-adjusted p); faded points are not significant.

### S7. Topología de las 20 redes (`supplementary/network_topology_panel`)

Script: [plot_network_topology.R](workflow/scripts/plot_network_topology.R). Color = sexo, vacío = control, lleno = AD
en todos los paneles; una sola leyenda.

| Panel | Qué muestra | Diseño |
|---|---|---|
| a–f | Aristas, densidad, grado medio, clustering, modularidad y número de comunidades | En cada región mujeres a la izquierda y hombres a la derecha; un segmento une control con AD, así que su largo y dirección son el cambio con AD |
| g | Distribución de grado P(k) | Log-log, bins logarítmicos; filas control / AD, facetas por región, color = sexo |
| h | Aristas vs número de muestras | Eje x logarítmico; segmento control→AD por sexo y región; etiqueta de región en el punto de AD |

**Pie propuesto.**
**Supplementary Figure S7. Topology of the 20 co-expression networks.**
(a–f) Number of edges, density, mean degree, global clustering coefficient, modularity and number of Leiden communities
(resolution 1) of each network. In each region, women are on the left and men on the right. Open circles, control;
filled circles, AD; the segment joins the control and AD networks of the same sex.
(g) Degree distribution of each network on log–log axes, with logarithmic bins (probability density per unit degree).
(h) Number of edges versus number of samples of each network (logarithmic x axis). Networks inferred from more samples
have more edges, so comparisons between networks combine biological differences and sample size.

### S8. Comparación de las redes (`supplementary/networks_panel`)

Script: [plot_networks_supplementary.R](workflow/scripts/plot_networks_supplementary.R).

| Panel | Qué muestra | Diseño |
|---|---|---|
| a | Jaccard de aristas de los 6 pares de redes de cada región | Heatmap pares × región, grises, valor escrito |
| b | Jaccard de nodos | Igual |
| c | ρ de Spearman del grado | Igual |
| d | NMI de las comunidades | Igual |
| e | Sensibilidad de recableados × DEG al umbral (top 5, 10 y 20%) | Líneas por sexo, facetas por región; punto lleno = umbral principal (10%) |

Cada heatmap tiene su propia escala de grises, del mínimo al máximo de esa medida, porque los rangos son muy distintos
(el Jaccard de nodos ronda 0.9 y el de aristas 0.1). Por eso no llevan barra de color: el valor está escrito.

**Pie propuesto.**
**Supplementary Figure S8. Similarity between every pair of networks and sensitivity of the rewiring analysis.**
(a–d) Similarity of the six pairs of networks of each region: (a) edge Jaccard index, (b) node Jaccard index,
(c) Spearman correlation of node degree, (d) normalized mutual information of Leiden communities. The shading of each
heatmap spans its own range; values are written in each cell. F, female; M, male; CT, control.
(e) log2 odds ratio of the overlap between rewired genes and DEGs of the same sex when the top 5%, 10% (main analysis,
filled points) or 20% most rewired genes are called rewired.

### S9. Función de las comunidades principales (`supplementary/community_enrichment_panel`)

Script: [plot_community_enrichment.R](workflow/scripts/plot_community_enrichment.R). Un solo panel, sin letra.

- **Filas:** el término más significativo de cada una de las 5 comunidades más grandes de cada red, ordenados por
  similitud de su patrón (clustering jerárquico).
- **Columnas:** las 20 redes (F-CT, F-AD, M-CT, M-AD dentro de cada región).
- **Puntos:** color = sexo, vacío = control, lleno = AD, tamaño = −log10(padj) en la comunidad de esa red donde el
  término está más enriquecido. Sin punto = el término no es significativo en ninguna comunidad de esa red.
- La altura crece con el número de términos.

**Pie propuesto.**
**Supplementary Figure S9. Biological functions of the main network communities.**
Over-representation analysis (Reactome; universe: genes of the network; communities of ≥ 20 genes) of the Leiden
communities of each network. Rows are the most significant term of each of the five largest communities of every
network, ordered by hierarchical clustering of their pattern. A point means the term is significant (BH-adjusted
p < 0.05) in at least one community of that network; its size is −log10(adjusted p) in the most enriched community.
Colour, sex; open circles, control; filled circles, AD.

---

## Pendientes de diseño

1. Figura 1c: rojo/azul para "mayor en mujeres / en hombres" → verde/naranja.
2. S1c: separar edad, PMI y RIN en letras propias.
3. S2: atípicas en rojo → negro o contorno (choca con AD).
4. S3/S4: forma por sexo → color = sexo, relleno = diagnóstico, como en la topología.
5. Revisar con los datos reales que las etiquetas de n (Figura 2a) y los términos de la Figura 2e quepan.
