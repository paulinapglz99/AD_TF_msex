#Global plots

#Libraries 

pacman::p_load("tidyverse", 
               "ggplot2", 
               "cowplot", 
               "ggrepel", 
               "ggsci")

out_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/global_plots"

dir.create(out_dir)
setwd(out_dir)

data <- vroom::vroom(file = "/STORAGE/genut/alzheimer/networks_ad_sex/topos_results/networks_summary.csv")
colnames(data)

data_parsed <- data %>%
  mutate(file_clean = str_remove(file, "\\.network\\.txt")) %>%
  separate(
    file_clean,
    into = c("cohort", "tech", "sex", "condition", "counts", "region"),
    sep = "_"
  ) %>%
  select(-cohort, -tech, -counts)%>%
  mutate(
    sex = factor(sex, levels = c("female","male")),
    condition = factor(condition, levels = c("control","AD")),
    region = factor(region, levels = c("CRB","DLPFC","HCN","PCC","TC")))


#Preprocess data

num_data <- data_parsed %>%
  select(where(is.numeric)) %>%
  select(-largest_community_id)


#Global plots

#Libraries 

pacman::p_load("tidyverse", 
               "ggplot2", 
               "cowplot", 
               "ggrepel", 
               "ggsci")

out_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/global_plots"

dir.create(out_dir)
setwd(out_dir)

data <- vroom::vroom(file = "/STORAGE/genut/alzheimer/networks_ad_sex/topos_results/networks_summary.csv") %>% dplyr::select(-perc_random_50)
colnames(data)

data_parsed <- data %>%
  mutate(file_clean = str_remove(file, "\\.network\\.txt")) %>%
  separate(
    file_clean,
    into = c("cohort", "tech", "sex", "condition", "counts", "region"),
    sep = "_"
  ) %>%
  select(-cohort, -tech, -counts)

#Preprocess data

num_data <- data_parsed %>%
  select(where(is.numeric)) %>%
  select(-largest_community_id)

num_scaled <- scale(num_data) %>%
  as.data.frame()

data_scaled <- bind_cols(
  data_parsed %>% select(region, sex, condition),
  num_scaled)

#Get delta data
delta_data <- data_scaled %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "metric",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = value
  ) %>%
  mutate(delta = AD - control) %>%
  mutate(
    metric = metric %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    sex = factor(sex,
                 levels = c("female","male"),
                 labels = c("Female","Male"))
  )

# 4. Ordenar métricas por magnitud promedio de cambio

metric_order <- delta_data %>%
  group_by(metric) %>%
  summarise(mean_abs_delta = mean(abs(delta))) %>%
  arrange(desc(mean_abs_delta)) %>%
  pull(metric)

delta_data$metric <- factor(delta_data$metric,
                            levels = metric_order)

delta_matrix <- delta_data %>%
  mutate(group = paste(region, sex, sep = "_")) %>%
  select(metric, group, delta) %>%
  pivot_wider(names_from = group, values_from = delta) %>%
  column_to_rownames("metric") %>%
  as.matrix()

dist_metrics <- as.dist(1 - cor(t(delta_matrix), use = "pairwise.complete.obs"))

hc_metrics <- hclust(dist_metrics, method = "ward.D2")

metric_order_clustered <- rownames(delta_matrix)[hc_metrics$order]

delta_data$metric <- factor(delta_data$metric,
                            levels = metric_order_clustered)

k <- 4
metric_clusters <- cutree(hc_metrics, k = k)

cluster_df <- data.frame(
  metric = names(metric_clusters),
  cluster = metric_clusters)

#Plot data

max_abs <- max(abs(delta_data$delta), na.rm = TRUE)

heatmap_delta <- ggplot(delta_data,
                        aes(x = region,
                            y = metric,
                            fill = delta)) +
  geom_tile(color = "white") +
  facet_wrap(~ sex) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs)) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()) +
  labs(
    fill = "AD - Control\n(Change in SD units)",
    x = "Region",
    y = "Network metric")

#Rojo = AD mayor
#Azul = AD menor

#Save plot
ggsave("heatmapdelta.jpg",
       heatmap_delta,
       width = 8,
       height = 10,
       dpi = 300)

#PCA

pca <- prcomp(num_data, scale = TRUE, center = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  mutate(
    sex = data_parsed$sex,
    condition = data_parsed$condition,
    region = data_parsed$region
  )

var_exp <- summary(pca)$importance[2, 1:2] * 100

pca.p <- ggplot(pca_df, aes(PC1, PC2,
                           color = sex,
                           shape = condition)) +
  geom_point(size = 6) +
  geom_text(aes(label = region), color = "black",  vjust = -1) +
  theme_minimal(base_size = 13) +
  scale_color_aaas() +
  stat_ellipse(aes(group = sex)) +
  labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_cowplot() +
  theme(legend.position = "top",
    plot.title = element_text(hjust = 0.5))

#Save plot 

ggsave("pca_global.jpg",
       pca.p,
       width = 10,
       height = 10,
       dpi = 300)

#PCA per region

pca_region.p <- ggplot(pca_df, aes(PC1, PC2,
                            color = sex,
                            shape = condition)) +
  geom_point(size = 6) +
  facet_wrap(~ region) +
  theme_minimal(base_size = 13) +
  scale_color_aaas() +
  #stat_ellipse(aes(group = sex)) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_cowplot() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

#Save plot 

ggsave("pca_regions.jpg",
       pca_region.p,
       width = 15,
       height = 12,
       dpi = 300)

#PCA per sex

pca_sex.p <- ggplot(pca_df, aes(PC1, PC2,
                                   color = region, 
                                   shape = condition)) +
  geom_point(size = 6) +
  facet_wrap(~ sex) +
  theme_minimal(base_size = 13) +
  scale_color_aaas() +
  geom_text_repel(aes(label = region),
                  size = 5.5,
                  show.legend = FALSE, 
                  color = "black") +
  #stat_ellipse(aes(group = region)) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_cowplot() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5))

#Save plot 

ggsave("pca_sex.jpg",
       pca_sex.p,
       width = 12,
       height = 9,
       dpi = 300)

#Contribution to PCs

#Global
contrib_global <- sweep(pca$rotation^2, 2, colSums(pca$rotation^2), FUN = "/") * 100

contrib_global.df <- as.data.frame(contrib_global[, 1:2]) %>% 
  rownames_to_column("variable") %>% 
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "contribution"
  ) %>% 
  group_by(PC) %>%
  arrange(desc(contribution), .by_group = TRUE) %>%
  mutate(variable = factor(variable, levels = variable)) %>%
  ungroup() %>%
  mutate(variable = variable %>%
           str_replace_all("_", " ") %>%
           str_to_sentence())


contrib_global.p <- ggplot(contrib_global.df,
                    aes(x = variable, y = contribution, fill = PC)) +
  geom_col(show.legend = TRUE) +
  coord_flip() +
  facet_wrap(~ PC, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "PC1" = "#1446A0",
      "PC2" = "#1446A0"
    )
  ) +
  labs(
    x = "Variable",
    y = "Contribution (%)",
    title = "Variable contributions to PC1 and PC2 - Global"
  ) +
  theme_cowplot()+
  theme(legend.position = "none")


#Save plot

ggsave("contrib_global_pca.jpeg",
       plot = contrib_global.p,
       device = "jpeg",
       width = 15,
       height = 13,
       units = "in",
       dpi = 300)

#Panel PCA

panel1 <- cowplot::plot_grid(pca.p, pca_sex.p, labels = c("A", "B"), ncol = 1)
panel2 <- cowplot::plot_grid(panel1, contrib_global.p, heatmap_delta, labels = c("", "C", "D"), ncol = 3)

ggsave("panel_pca.jpeg",
       plot = panel2,
       device = "jpeg",
       width = 25,
       height = 13,
       units = "in",
       dpi = 300)

#La red está organizada alrededor de un “core” denso que coordina muchos nodos periféricos más débiles?
#“Las redes AD muestran reorganización topológica global diferente por sexo”???

#Nucleus

nucleo_periferia <- ggplot(data, aes(kcore_max, deg_mean)) +
  geom_point(size = 3) +
  theme_bw()

ggsave("nucleo_periferia.jpg",
       nucleo_periferia,
       width = 8,
       height = 10,
       dpi = 300)

#END?