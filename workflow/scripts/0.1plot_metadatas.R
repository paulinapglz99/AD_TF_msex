# 0.1plot_metadatas.R
# Supplementary figure with the sample metadata of the five brain regions:
#   a  Samples per region, diagnosis and sex
#   b  Proportion of women and men per region and diagnosis
#   c  Age at death, PMI and RIN per region, diagnosis and sex
#   d  Sex check: XIST vs UTY expression in every sample
# Diagnosis (AD / CT) is Sage's harmonized 'diagnosis' column from the
# *_Covariates_ages_censored.tsv files, the same variable used in the rest of the workflow.
# Colors, fonts, group order and sizes come from scripts/R/theme_pub.R.
# Run by Snakemake (rule plot_metadata).

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  source(snakemake@input[["theme"]])
})

regions <- read_tsv(snakemake@input[["regions"]], show_col_types = FALSE)
diagnoses <- unlist(snakemake@params[["diagnoses"]])

# ---- Covariates -------------------------------------------------------------

read_covariates <- function(path, study) {
  read_tsv(path, col_types = cols(.default = col_character())) %>%
    transmute(specimenID,
              individualID,
              study = study,
              covariates_tissue = tissue,
              diagnosis,
              sex,
              age_death = parse_number(age_death), # "90+" (censored) becomes 90
              pmi = as.numeric(pmi),
              RIN = as.numeric(RIN))
}

covariates <- bind_rows(
  read_covariates(snakemake@input[["rosmap_covariates"]], "ROSMAP"),
  read_covariates(snakemake@input[["mayo_covariates"]], "Mayo")
) %>%
  # Sage's tissue labels -> region names of config/samples.tsv (e.g. ACC -> HCN)
  inner_join(select(regions, region, study, covariates_tissue),
             by = c("study", "covariates_tissue")) %>%
  filter(diagnosis %in% diagnoses) %>%
  set_pub_levels()

cat("Samples per region, diagnosis and sex:\n")
print(as.data.frame(count(covariates, region, diagnosis, sex)))

# ---- XIST and UTY counts ----------------------------------------------------

sex_genes <- c(ENSG00000229807 = "XIST", ENSG00000183878 = "UTY")

# Reads only the header and the two genes (grep), not the whole matrix
read_sex_genes <- function(path) {
  pattern <- paste0("^(feature|", paste0(names(sex_genes), "\\.", collapse = "|"), ")")
  read_tsv(pipe(paste("grep -E", shQuote(pattern), shQuote(path))),
           show_col_types = FALSE) %>%
    mutate(gene = sex_genes[sub("\\..*$", "", feature)]) %>%
    select(-feature) %>%
    pivot_longer(-gene, names_to = "counts_id", values_to = "count")
}

sex_expression <- c(snakemake@input[["rosmap_counts"]], snakemake@input[["mayo_counts"]]) %>%
  lapply(read_sex_genes) %>%
  bind_rows() %>%
  pivot_wider(names_from = gene, values_from = count)

# Mayo cerebellum is "_CBE" in Sage's covariates but "_CER" in the count matrix
sex_check <- covariates %>%
  mutate(counts_id = if_else(study == "Mayo",
                             sub("_CBE$", "_CER", specimenID),
                             specimenID)) %>%
  inner_join(sex_expression, by = "counts_id")

cat("\nSamples with XIST/UTY counts:", nrow(sex_check), "of", nrow(covariates), "\n")

# ---- Plots ------------------------------------------------------------------

# Region names on the x axis of panels a and b overlap if horizontal
tilted_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

# The sex legend is drawn once, from panel a; the other panels hide theirs
no_legend <- guides(fill = "none")

n_groups <- covariates %>%
  count(region, diagnosis, sex) %>%
  group_by(region, diagnosis) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# a. Samples per group
p_n <- ggplot(n_groups, aes(region, n, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8),
            vjust = -0.3, size = PUB_LABEL_SIZE) +
  facet_wrap(~diagnosis, labeller = as_labeller(DX_LABELS)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_sex() +
  labs(x = NULL, y = "Samples") +
  theme_pub() +
  tilted_x

# b. Proportion of women and men
p_prop <- ggplot(n_groups, aes(region, prop, fill = sex)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = paste0(round(100 * prop), "%")),
            position = position_stack(vjust = 0.5),
            colour = "white", size = PUB_LABEL_SIZE) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = PUB_LINE_WIDTH) +
  facet_wrap(~diagnosis, labeller = as_labeller(DX_LABELS)) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_sex() +
  no_legend +
  labs(x = NULL, y = "Proportion of samples") +
  theme_pub() +
  tilted_x

# c. Age at death, PMI and RIN
covariates_long <- covariates %>%
  pivot_longer(c(age_death, pmi, RIN), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable,
                           levels = c("age_death", "pmi", "RIN"),
                           labels = c("Age at death (years)", "PMI (hours)", "RIN")))

p_covariates <- ggplot(covariates_long, aes(diagnosis, value, fill = sex)) +
  geom_violin(position = position_dodge(width = 0.85), scale = "width",
              alpha = 0.6, linewidth = PUB_LINE_WIDTH) +
  geom_boxplot(position = position_dodge(width = 0.85), width = 0.15,
               outlier.size = 0.2, linewidth = PUB_LINE_WIDTH, show.legend = FALSE) +
  facet_grid(variable ~ region, scales = "free_y", switch = "y") +
  scale_x_dx() +
  scale_fill_sex() +
  no_legend +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(strip.placement = "outside")

# d. Sex check
p_sex <- ggplot(sex_check, aes(log2(XIST + 1), log2(UTY + 1), fill = sex)) +
  geom_point(shape = 21, colour = "white", stroke = 0.1, size = 0.9, alpha = 0.8) +
  facet_wrap(~study) +
  scale_fill_sex() +
  no_legend +
  labs(x = expression(log[2]("XIST counts" + 1)),
       y = expression(log[2]("UTY counts" + 1))) +
  theme_pub()

panel <- (p_n | p_prop) / p_covariates / p_sex +
  plot_layout(heights = c(1, 2.4, 1.2), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  # patchwork places collected legends on the right unless told otherwise
  theme(legend.position = theme_pub()$legend.position)

save_pub(panel, snakemake@output[["pdf"]], width = "double", height = FIG_MAX_HEIGHT)

cat("\nDone\n")
