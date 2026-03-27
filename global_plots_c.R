# ============================================================
# Global network topology analysis
# AD transcriptional landscape - sex differences
# ============================================================

# ---- Libraries ----

pacman::p_load(
  "tidyverse",
  "ggplot2",
  "cowplot",
  "ggrepel",
  "ggsci",
  "vroom",
  "tidytext"    # reorder_within / scale_x_reordered
)

# ---- Paths ----

base_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex"
out_dir  <- file.path(base_dir, "global_plots")
dir.create(out_dir, showWarnings = FALSE)
setwd(out_dir)

# ---- Color palette ----

aaas_cols  <- ggsci::pal_aaas()(10)
col_ctrl   <- aaas_cols[1]    # azul oscuro  → control
col_ad     <- aaas_cols[2]    # rojo         → AD

col_f_ctrl <- "#4393C3"       # azul medio   → female control
col_m_ctrl <- col_ctrl        # azul oscuro  → male control
col_f_ad   <- col_ad          # rojo AAAS    → female AD
col_m_ad   <- "#8B0000"       # rojo oscuro  → male AD

# ============================================================
# 1.  Load & parse data
# ============================================================

data_parsed <- vroom(file.path(base_dir, "topos_results/networks_summary.csv")) %>%
  select(-perc_random_50) %>%
  mutate(file_clean = str_remove(file, "\\.network\\.txt")) %>%
  separate(
    file_clean,
    into = c("cohort", "tech", "sex", "condition", "counts", "region"),
    sep  = "_"
  ) %>%
  select(-cohort, -tech, -counts) %>%
  mutate(
    sex       = factor(sex,       levels = c("female", "male")),
    condition = factor(condition, levels = c("control", "AD")),
    region    = factor(region,    levels = c("CRB", "DLPFC", "HCN", "PCC", "TC"))
  )

# ---- Numeric matrix (raw and scaled) ----

num_data <- data_parsed %>%
  select(where(is.numeric)) %>%
  select(-largest_community_id)

num_scaled <- scale(num_data) %>% as.data.frame()

data_scaled <- bind_cols(
  data_parsed %>% select(region, sex, condition),
  num_scaled
)

# ============================================================
# 2.  PCA
#     shape = sex  → female = círculo (16), male = triángulo (17)
#     color = sex × condition → azules control, rojos AD
#     4 elipses con fill semi-transparente
# ============================================================

pca     <- prcomp(num_data, scale. = TRUE, center = TRUE)
var_exp <- summary(pca)$importance[2, 1:2] * 100

pca_df <- as.data.frame(pca$x) %>%
  mutate(
    sex       = data_parsed$sex,
    condition = data_parsed$condition,
    region    = data_parsed$region,
    sex_cond  = factor(
      paste(sex, condition, sep = "_"),
      levels = c("female_control", "male_control", "female_AD", "male_AD")
    ),
    label = as.character(data_parsed$region)
  )

#Save pca data

#vroom::vroom_write(pca_df, file = "pca_df.csv")

loadings <- as.data.frame(pca$rotation) %>%
  rownames_to_column("variable") %>%
  arrange(desc(abs(PC1)))

#Save loadings

vroom::vroom_write(loadings, file = "pca_loadings.csv")

sex_cond_colors <- c(
  "female_control" = col_f_ctrl,
  "male_control"   = col_m_ctrl,
  "female_AD"      = col_f_ad,
  "male_AD"        = col_m_ad
)

sex_cond_labels <- c(
  "female_control" = "FC",
  "male_control"   = "MC",
  "female_AD"      = "FAD",
  "male_AD"        = "MAD"
)

sex_cond_shapes <- c(
  "female_control" = 16,
  "male_control"   = 17,
  "female_AD"      = 16,
  "male_AD"        = 17
)

pca_global.p <- ggplot(pca_df,
                       aes(PC1, PC2,
                           color = sex_cond,
                           shape = sex_cond)) +
  # FIX 1: alpha en fill de elipses para que no tapen los puntos
  stat_ellipse(aes(group = sex_cond, fill = sex_cond),
               geom      = "polygon",
               alpha     = 0.08,
               linewidth = 0.8,
               level     = 0.85,
               show.legend = FALSE) +
  stat_ellipse(aes(group = sex_cond),
               linewidth = 0.8, level = 0.85,
               show.legend = FALSE) +
  geom_point(size = 5) +
  geom_text_repel(
    aes(label = label),
    color        = "black",
    size         = 3.5,
    show.legend  = FALSE,
    box.padding  = 0.4,
    max.overlaps = 20
  ) +
  scale_color_manual(values = sex_cond_colors, labels = sex_cond_labels, name = NULL) +
  scale_fill_manual(values  = sex_cond_colors, labels = sex_cond_labels, name = NULL) +
  scale_shape_manual(values = sex_cond_shapes, labels = sex_cond_labels, name = NULL) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  ) +
  theme_cowplot(font_size = 13) +
  theme(legend.position = "top") +
  guides(
    color = guide_legend(nrow = 2, override.aes = list(size = 4)),
    shape = guide_legend(nrow = 2)
  )

ggsave("pca_global.jpg", pca_global.p, width = 7, height = 7, dpi = 300)


pca_regions <- pca_global.p +
  facet_wrap(~region, nrow =5)

ggsave("pca_regions.jpg", pca_regions, width = 5, height = 17, dpi = 300)

# ============================================================
# 3.  Contributions to PC1 & PC2
#     Orden: mayor → menor (barra más larga arriba), independiente por PC
# ============================================================

# Calcular contribuciones y ordenar por facet
# fct_reorder ascendente + coord_flip → barra más larga queda arriba
contrib_long <- (sweep(pca$rotation^2, 2, colSums(pca$rotation^2), "/") * 100) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "contribution") %>%
  filter(PC %in% c("PC1", "PC2")) %>%
  mutate(variable = variable %>% str_replace_all("_", " ") %>% str_to_sentence())

# Orden independiente por PC: crear variable con sufijo de PC para free_y
contrib_pc1 <- contrib_long %>%
  filter(PC == "PC1") %>%
  mutate(variable = fct_reorder(variable, contribution))

contrib_pc2 <- contrib_long %>%
  filter(PC == "PC2") %>%
  mutate(variable = fct_reorder(variable, contribution))

contrib_long <- bind_rows(contrib_pc1, contrib_pc2)

contrib.p <- ggplot(contrib_long,
                    aes(x = variable, y = contribution)) +
  geom_col(fill = col_ctrl, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Variable", y = "Contribution (%)") +
  theme_cowplot(font_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave("contrib_pca.jpg", contrib.p, width = 14, height = 10, dpi = 300)

# ============================================================
# 4.  Euclidean distance heatmap
#     FIX 3: etiquetas eje X escalonadas con guide_axis(n.dodge = 2)
# ============================================================

sex_abbr  <- c(female = "F", male = "M")
cond_abbr <- c(control = "C", AD = "AD")

net_labels <- data_parsed %>%
  mutate(net_id = paste0(
    sex_abbr[as.character(sex)],
    cond_abbr[as.character(condition)],
    "-", region
  )) %>%
  pull(net_id)

rownames(num_scaled) <- net_labels

dist_mat <- dist(num_scaled, method = "euclidean") %>% as.matrix()

region_levels <- c("CRB", "DLPFC", "HCN", "PCC", "TC")
group_order   <- c("FAD", "FC", "MAD", "MC")

net_order <- outer(group_order, region_levels,
                   FUN = function(g, r) paste0(g, "-", r)) %>%
  as.vector() %>%
  .[. %in% rownames(dist_mat)]

dist_long <- as.data.frame(dist_mat) %>%
  rownames_to_column("net_i") %>%
  pivot_longer(-net_i, names_to = "net_j", values_to = "distance") %>%
  mutate(
    net_i = factor(net_i, levels = net_order),
    net_j = factor(net_j, levels = net_order)
  )

dist_long <-dist_long %>%
  filter(as.numeric(net_i) < as.numeric(net_j))

region_breaks <- c(4.5, 8.5, 12.5, 16.5)

euclid_heatmap.p <- ggplot(dist_long,
                           aes(x = net_j, y = net_i, fill = distance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = region_breaks, color = "grey30",
             linewidth = 0.9, linetype = "solid") +
  geom_vline(xintercept = region_breaks, color = "grey30",
             linewidth = 0.9, linetype = "solid") +
  scale_fill_gradient(
    low  = "white",
    high = "#631879",
    name = "Euclidean\ndistance"
  ) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +   # etiquetas escalonadas
  labs(x = NULL, y = NULL) +
  theme_cowplot(font_size = 10) +
  theme(
    axis.text.x     = element_text(size = 10),
    axis.text.y     = element_text(size = 13),
    legend.position = "right"
  )

ggsave("euclidean_dist_heatmap.jpg",
       euclid_heatmap.p, width = 11, height = 10, dpi = 300)

# ============================================================
# 5.  Delta heatmap  (AD - Control, en unidades SD)
#     Escala secuencial: morado claro → morado oscuro
# ============================================================

delta_data <- data_scaled %>%
  pivot_longer(where(is.numeric), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  mutate(
    delta  = AD - control,
    metric = metric %>% str_replace_all("_", " ") %>% str_to_sentence(),
    sex    = factor(sex,
                    levels = c("female", "male"),
                    labels = c("Female", "Male"))
  )

delta_matrix <- delta_data %>%
  mutate(group = paste(region, sex, sep = "_")) %>%
  select(metric, group, delta) %>%
  pivot_wider(names_from = group, values_from = delta) %>%
  column_to_rownames("metric") %>%
  as.matrix()

hc_metrics        <- hclust(
  as.dist(1 - cor(t(delta_matrix), use = "pairwise.complete.obs")),
  method = "ward.D2"
)
delta_data$metric <- factor(delta_data$metric,
                            levels = rownames(delta_matrix)[hc_metrics$order])

delta_heatmap.p <- ggplot(delta_data,
                          aes(x = region, y = metric, fill = delta)) +
  geom_tile(color = "white") +
  facet_wrap(~ sex) +
  scale_fill_gradient2(
    low  = "#3B4992",
    mid  = "white",
    high = "#BB0021",
    midpoint = 0,
    name = "AD − Control\n(SD units)"
  ) +
  labs(x = "Region", y = "Network metric") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text  = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid  = element_blank()
  )

ggsave("delta_heatmap.jpg", delta_heatmap.p, width = 8, height = 10, dpi = 300)

# ============================================================
# 6.  Panel final
#     A = PCA  |  B = Contribuciones
#     C = Distancias euclidianas  |  D = Delta heatmap
# ============================================================

top_row    <- plot_grid(pca_global.p,     contrib.p,        labels = c("A", "B"), ncol = 2)
bottom_row <- plot_grid(euclid_heatmap.p, delta_heatmap.p,  labels = c("C", "D"), ncol = 2)
panel_full <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1.1))

ggsave("panel_global.jpg",
       panel_full,
       width  = 22,
       height = 20,
       dpi    = 300)

# ============================================================
# END
# ============================================================