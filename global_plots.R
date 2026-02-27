#Global plots

#Libraries 

pacman::p_load("tidyverse", 
               "ggplot2", 
               "cowplot", 
               "ggrepel")

out_dir <- "/STORAGE/genut/alzheimer/networks_ad_sex/global_plots"

dir.create(out_dir)
setwd(out_dir)

data <- vroom::vroom(file = "/STORAGE/genut/alzheimer/networks_ad_sex/topos_results/networks_summary.csv")

data_parsed <- data %>%
  mutate(file_clean = str_remove(file, "\\.network\\.txt")) %>%
  separate(
    file_clean,
    into = c("cohort", "tech", "sex", "condition", "counts", "region"),
    sep = "_"
  ) %>%
  select(-cohort, -tech, -counts)

#Preprocess data

num_data <- num_data <- data_parsed %>%
  select(where(is.numeric)) %>%
  select(-largest_community_id, -perc_random_50)

data_parsed <- data_parsed %>%
  mutate(
    sex = factor(sex, levels = c("female","male")),
    condition = factor(condition, levels = c("control","AD")),
    region = factor(region, levels = c("CRB","DLPFC","HCN","PCC","TC"))
  )

scaled_data <- scale(num_data) %>%
  as.data.frame() %>%
  mutate(
    region = data_parsed$region,
    sex = data_parsed$sex,
    condition = data_parsed$condition,
    id = paste(region, sex, condition, sep = "_")
  )

#plot heatmap

heatmap_data <- scaled_data %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "metric",
    values_to = "value"
  )

heatmap_data <- heatmap_data %>%
  mutate(metric = metric %>%
           str_replace_all("_", " ") %>%
           str_to_sentence())

heatmap.p <- ggplot(heatmap_data,
            aes(x = condition,
                y = metric,
                fill = value)) +
  geom_tile(color = "white") +
  facet_grid(sex ~ region, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Condition",
    y = "Network metric",
    fill = "Z-score"
  )

#Save heatmap

ggsave("heatmap_global.jpg",
       heatmap.p,
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
  theme_minimal(base_size = 13) +
  scale_color_aaas() +
  stat_ellipse(aes(group = sex)) +
  labs(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
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
                                   color = region)) +
  geom_point(size = 6) +
  facet_wrap(~ sex) +
  theme_minimal(base_size = 13) +
  scale_color_aaas() +
  geom_text_repel(aes(label = region),
                  size = 6,
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
      "PC2" = "#7d1538"
    )
  ) +
  labs(
    x = "Variable",
    y = "Contribution (%)",
    title = "Variable contributions to PC1 and PC2 - Global"
  ) +
  theme_cowplot()

#Save plot

ggsave("contrib_global_pca.jpeg",
       plot = contrib_global.p,
       device = "jpeg",
       width = 15,
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