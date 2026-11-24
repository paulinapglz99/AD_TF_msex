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

base_dir <- "/home/tiamat/Documents/AD_TF_sex_data/networks_ad_sex/"
out_dir  <- file.path(base_dir, "global_plots_local")
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

vroom::vroom_write(pca_df, file = "pca_df.csv")

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
  geom_point(size = 8) +
  geom_text_repel(
    aes(label = label),
    color        = "black",
    size         = 6.5,
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
  theme(legend.position = "top", 
        legend.text=element_text(size=20)) +
  guides(
    color = guide_legend(nrow = 2, override.aes = list(size = 4)),
    shape = guide_legend(nrow = 2)
  )

#Vis

pca_global.p

ggsave("pca_global.jpg", pca_global.p, width = 7, height = 7, dpi = 300)

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
  labs(x = element_blank(), y = "Contribution (%)") +
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
  labs(x = element_blank(), y = element_blank()) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text  = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid  = element_blank()
  )

#Vis

delta_heatmap.p

ggsave("delta_heatmap.jpg", delta_heatmap.p, width = 8, height = 10, dpi = 300)

# ============================================================
# ADDITIONS: Sexual symmetry + Delta PC plots
# Mirrors the comparative panel from interactive analysis
# Insert after section 2 (PCA) in your original script
# Requires: pca_df, var_exp, sex_cond_colors, col_f_ctrl,
#           col_m_ctrl, col_f_ad, col_m_ad (already defined)
# ============================================================

# ---- Build trajectory delta data ----
# One row per sex × region, with control and AD positions + deltas

traj_data <- pca_df %>%
  select(PC1, PC2, sex, condition, region) %>%
  pivot_wider(
    id_cols     = c(sex, region),
    names_from  = condition,
    values_from = c(PC1, PC2)
  ) %>%
  rename(
    x_ctrl = PC1_control, y_ctrl = PC2_control,
    x_ad   = PC1_AD,      y_ad   = PC2_AD
  ) %>%
  mutate(
    sex       = factor(sex, levels = c("female", "male")),
    region    = factor(region, levels = c("CRB", "DLPFC", "HCN", "PCC", "TC")),
    delta_PC1 = x_ad - x_ctrl,
    delta_PC2 = y_ad - y_ctrl
  )

# ---- Symmetry data: sex gap in control vs AD ----
# FC − MC and FAD − MAD, for PC1 and PC2

symmetry_data <- pca_df %>%
  select(PC1, PC2, sex, condition, region) %>%
  pivot_wider(
    id_cols     = c(region, condition),
    names_from  = sex,
    values_from = c(PC1, PC2)
  ) %>%
  rename(
    PC1_female = PC1_female, PC1_male = PC1_male,
    PC2_female = PC2_female, PC2_male = PC2_male
  ) %>%
  mutate(
    gap_PC1 = PC1_female - PC1_male,
    gap_PC2 = PC2_female - PC2_male,
    region    = factor(region,    levels = c("CRB", "DLPFC", "HCN", "PCC", "TC")),
    condition = factor(condition, levels = c("control", "AD"),
                       labels = c("Control (FC − MC)", "AD (FAD − MAD)"))
  )

# ============================================================
# Plot A: Δ PC1 by AD — females vs males per region
# ============================================================

delta_PC1_long <- traj_data %>%
  select(sex, region, delta_PC1) %>%
  mutate(
    sex    = factor(sex, labels = c("Female (FC→FAD)", "Male (MC→MAD)")),
    region = factor(region, levels = c("CRB", "DLPFC", "HCN", "PCC", "TC"))
  )

delta_pc1.p <- ggplot(delta_PC1_long,
                      aes(x = region, y = delta_PC1, fill = sex)) +
  geom_col(
    position = position_dodge(width = 0.65),
    width    = 0.55,
    color    = NA
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%+.2f", delta_PC1),
      vjust = ifelse(delta_PC1 >= 0, -0.4, 1.3)
    ),
    position    = position_dodge(width = 0.65),
    size        = 3.2,
    fontface    = "bold",
    color       = "grey20"
  ) +
  scale_fill_manual(
    values = c("Female (FC→FAD)" = col_f_ctrl,
               "Male (MC→MAD)"   = col_m_ctrl),
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  labs(
    x = "Brain region",
    y = expression(Delta * " PC1 (AD − Control)"),
    title = paste0("PC1 — global integration (", round(var_exp[1], 1), "% variance)")
  ) +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = "top",
    plot.title      = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

#Vis 
delta_pc1.p

# ============================================================
# Plot B: Δ PC2 by AD — females vs males per region
# ============================================================

delta_PC2_long <- traj_data %>%
  select(sex, region, delta_PC2) %>%
  mutate(
    sex    = factor(sex, labels = c("Female (FC→FAD)", "Male (MC→MAD)")),
    region = factor(region, levels = c("CRB", "DLPFC", "HCN", "PCC", "TC"))
  )

delta_pc2.p <- ggplot(delta_PC2_long,
                      aes(x = region, y = delta_PC2, fill = sex)) +
  geom_col(
    position = position_dodge(width = 0.65),
    width    = 0.55,
    color    = NA
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%+.2f", delta_PC2),
      vjust = ifelse(delta_PC2 >= 0, -0.4, 1.3)
    ),
    position    = position_dodge(width = 0.65),
    size        = 3.2,
    fontface    = "bold",
    color       = "grey20"
  ) +
  scale_fill_manual(
    values = c("Female (FC→FAD)" = col_f_ctrl,
               "Male (MC→MAD)"   = col_m_ctrl),
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  labs(
    x = "Brain region",
    y = expression(Delta * " PC2 (AD − Control)"),
    title = paste0("PC2 — mesoscale organization (", round(var_exp[2], 1), "% variance)")
  ) +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position = "top",
    plot.title      = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

#Vis 

delta_pc2.p

# ============================================================
# Plot C: Sexual symmetry — sex gap in PC1 (control vs AD)
# ============================================================

sym_pc1.p <- ggplot(symmetry_data,
                    aes(x = region, y = gap_PC1, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.65),
    width    = 0.55,
    color    = NA
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%+.2f", gap_PC1),
      vjust = ifelse(gap_PC1 >= 0, -0.4, 1.3)
    ),
    position    = position_dodge(width = 0.65),
    size        = 6,
    fontface    = "bold",
    color       = "grey20"
  ) +
  scale_fill_manual(
    values = c("Control (FC − MC)" = col_m_ctrl,
               "AD (FAD − MAD)"    = col_f_ad),
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  labs(
    x     = element_blank(),
    y     = "PC1 sex gap (Female − Male)",
    title = "Sexual asymmetry in PC1 — global integration"
  ) +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position    = "top",
    plot.title         = element_text(size = 17, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

sym_pc1.p

# ============================================================
# Plot D: Sexual symmetry — sex gap in PC2 (control vs AD)
# ============================================================

sym_pc2.p <- ggplot(symmetry_data,
                    aes(x = region, y = gap_PC2, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.65),
    width    = 0.55,
    color    = NA
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_text(
    aes(
      label = sprintf("%+.2f", gap_PC2),
      vjust = ifelse(gap_PC2 >= 0, -0.4, 1.3)
    ),
    position    = position_dodge(width = 0.65),
    size        = 6,
    fontface    = "bold",
    color       = "grey20"
  ) +
  scale_fill_manual(
    values = c("Control (FC − MC)" = col_m_ctrl,
               "AD (FAD − MAD)"    = col_f_ad),
    name   = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.18, 0.18))
  ) +
  labs(
    x     = element_blank(),
    y     = "PC2 sex gap (Female − Male)",
    title = "Sexual asymmetry in PC2 — mesoscale organization"
  ) +
  theme_cowplot(font_size = 12) +
  theme(
    legend.position    = "top",
    plot.title         = element_text(size = 17, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

#Vis

sym_pc2.p

# ============================================================
# Save individual plots
# ============================================================

ggsave("delta_pc1.jpg",   delta_pc1.p, width = 7, height = 5, dpi = 300)
ggsave("delta_pc2.jpg",   delta_pc2.p, width = 7, height = 5, dpi = 300)
ggsave("sym_pc1.jpg",     sym_pc1.p,   width = 7, height = 5, dpi = 300)
ggsave("sym_pc2.jpg",     sym_pc2.p,   width = 7, height = 5, dpi = 300)

# ============================================================
# FIGURE 1 — Main results panel
#   A = PCA global (large)
#   B = Delta heatmap (metrics x region x sex)
#   C = Sexual asymmetry PC1
#   D = Sexual asymmetry PC2
# ============================================================

# Fix factor order in symmetry_data if not already done
symmetry_data <- symmetry_data %>%
  mutate(condition = factor(condition,
                            levels = c("Control (FC − MC)", "AD (FAD − MAD)")))

# Rebuild sym plots with consistent x-axis order and cleaned theme
sym_pc1_final.p <- sym_pc1.p +
  scale_x_discrete(limits = c("CRB", "DLPFC", "HCN", "PCC", "TC")) +
  theme(legend.position = "top")

sym_pc2_final.p <- sym_pc2.p +
  scale_x_discrete(limits = c("CRB", "DLPFC", "HCN", "PCC", "TC")) +
  theme(legend.position = "top")

# Rebuild delta heatmap with sex ordered Female → Male
delta_heatmap_final.p <- delta_heatmap.p +
  theme(legend.position = "right")

# ---- Assemble ----
row1 <-  plot_grid(
  pca_global.p,
  delta_heatmap_final.p,
  labels     = c("A", "B"),
  ncol       = 2,
  rel_widths = c(1, 1),
  align      = "hv",
  axis       = "tb"
)

row2 <- plot_grid(
  sym_pc1_final.p,
  sym_pc2_final.p,
  labels     = c("C", "D"),
  ncol       = 2,
  rel_widths = c(1, 1),
  align      = "hv",
  axis       = "tb"
)

fig1 <- plot_grid(
  row1,
  row2,
  labels      = c("A", ""),
  ncol        = 1,
  rel_heights = c(1.2, 1)
)

#Vis
fig1 

#Save
ggsave(
  "topos_PCA_main.jpg",
  fig1,
  width  = 16,
  height = 16,
  dpi    = 300
)

# ============================================================
# SUPPLEMENTARY FIGURE S1
#   Row 1: PCA by region (facet, horizontal — nrow=1)
#   Row 2: PCA contributions (PC1 | PC2)
#   Row 3: ΔPC1 | ΔPC2
# ============================================================

# PCA by region: switch to nrow=1 (horizontal) for supplementary
# so it reads left to right and doesn't compress vertically

pca_regions_horiz <- ggplot(pca_df,
                            aes(PC1, PC2,
                                color = sex_cond,
                                shape = sex_cond)) +
  geom_point(size = 8) +
  scale_color_manual(values = sex_cond_colors, labels = sex_cond_labels, name = NULL) +
  scale_fill_manual(values  = sex_cond_colors, labels = sex_cond_labels, name = NULL) +
  scale_shape_manual(values = sex_cond_shapes, labels = sex_cond_labels, name = NULL) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
    
  ) +
  theme_cowplot(font_size = 13) +
  theme(legend.position = "top", 
        legend.text=element_text(size=20)) +
  guides(
    color = guide_legend(nrow = 2, override.aes = list(size = 4)),
    shape = guide_legend(nrow = 2)
  ) + 
  facet_wrap(~ region, nrow = 1) +
  theme(
    axis.text=element_text(size=12),
    axis.title=element_text(size=18),
    legend.position = "top",
    strip.text      = element_text(face = "bold", size = 19), 
    #strip.background = "gray", 
    #strip.placement = "bottom"
  )

#Vis
pca_regions_horiz

# Contributions: already built as contrib.p (PC1 | PC2 faceted)

# Delta plots: reuse delta_pc1.p and delta_pc2.p

row_s2 <- plot_grid(
  contrib.p,
  labels = c("B"),
  ncol   = 1
)

row_s3 <- plot_grid(
  delta_pc1.p,
  delta_pc2.p,
  labels     = c("C", "D"),
  ncol       = 2,
  align      = "hv",
  axis       = "tb"
)

figS1 <- plot_grid(
  pca_regions_horiz,
  row_s2,
  row_s3,
  labels      = c("A", "", ""),
  ncol        = 1,
  rel_heights = c(0.9, 1.1, 0.9)
)

ggsave(
  "PCA_topos_supplementary.jpg",
  figS1,
  width  = 18,
  height = 20,
  dpi    = 300
)

# ============================================================
# END
# ============================================================