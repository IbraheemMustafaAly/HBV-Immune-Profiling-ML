# ============================================================
#  HBV Multi-omics Analysis — R Script
#  Integrative Cytokine + miRNA + SNP Analysis
#  Target: Q1 Journal Publication
# ============================================================
#
#  COHORTS (complete-case, no imputation):
#   Cohort A: SNP only          n=200 (G0=104, G1=96)
#   Cohort B: ELISA only        n=143 (G0=46,  G1=97)
#   Cohort C: miRNA only        n=98  (G0=48,  G1=50)
#   Cohort D: ELISA + miRNA     n=77  (G0=29,  G1=48)
#   Cohort E: ALL features      n=73  (G0=29,  G1=44)
#
#  OUTPUT FILES:
#   - Figures 1-7 (PDF/PNG, publication quality 300dpi)
#   - stats_results.csv
#   - cohort_D_ml_ready.csv  (for Python Colab)
#   - cohort_E_ml_ready.csv  (for Python Colab)
# ============================================================


# ---- 0. PACKAGES -------------------------------------------
# Run once to install:
# install.packages(c("tidyverse","ggpubr","rstatix","corrplot",
#                    "ComplexHeatmap","FactoMineR","factoextra",
#                    "circlize","RColorBrewer","scales","ggrepel",
#                    "HardyWeinberg","epitools","patchwork"))

library(tidyverse)
library(ggpubr)
library(rstatix)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(HardyWeinberg)
library(epitools)
library(patchwork)

# ---- 1. LOAD & CLEAN DATA ----------------------------------

raw <- read.csv("Whole_Data_final_5_cytokine_and_mirs_1_.csv",
                na.strings = c("", " ", "NA", "NaN"),
                stringsAsFactors = FALSE)

# Convert all columns to numeric except Group
raw <- raw %>%
  mutate(across(-Group, ~ suppressWarnings(as.numeric(.))))

# Keep only valid groups (0 and 1)
raw <- raw %>%
  mutate(Group = as.numeric(Group)) %>%
  filter(Group %in% c(0, 1)) %>%
  mutate(Group = factor(Group, levels = c(0, 1),
                        labels = c("Control", "HBV")))

cat("Total valid samples:", nrow(raw), "\n")
cat("Controls:", sum(raw$Group == "Control"), "\n")
cat("HBV:     ", sum(raw$Group == "HBV"), "\n\n")


# ---- 2. DEFINE FEATURE GROUPS ------------------------------

elisa_cols <- c("IL10_ELISA", "TGFb_ELISA", "IL6_level",
                "TNFA_ELISA", "IFN_ELISA")

elisa_labels <- c("IL-10", "TGF-β", "IL-6", "TNF-α", "IFN-γ")

mirna_cols <- c("mir_10","mir_17","mir_21","mir_24","mir_26",
                "mir_122","mir_125","mir_145","mir_146",
                "mir_148","mir_155","mir_221")

mirna_labels <- c("miR-10","miR-17","miR-21","miR-24","miR-26",
                  "miR-122","miR-125","miR-145","miR-146",
                  "miR-148","miR-155","miR-221")

# SNP: use only the _Geno columns (1=homozygous major,
#       2=heterozygous, 3=homozygous minor)
geno_cols <- c("IL10_1082_Geno","IL10_819_Geno",
               "TGFb_800_Geno","TGFb_509_Geno",
               "TGFb_codon10_Geno","TGFb_codon25_Geno",
               "TNFA_863_Geno","TNFA_376_Geno",
               "TNFA_308_Geno","TNFA_857_Geno","TNFA_489_Geno")

geno_labels <- c("IL10 -1082","IL10 -819",
                 "TGFb -800","TGFb -509",
                 "TGFb codon10","TGFb codon25",
                 "TNFa -863","TNFa -376",
                 "TNFa -308","TNFa -857","TNFa -489")


# ---- 3. BUILD COHORTS (complete-case) ----------------------

cohort_A <- raw %>%
  select(Group, all_of(geno_cols)) %>% drop_na()

cohort_B <- raw %>%
  select(Group, all_of(elisa_cols)) %>% drop_na()

cohort_C <- raw %>%
  select(Group, all_of(mirna_cols)) %>% drop_na()

cohort_D <- raw %>%
  select(Group, all_of(elisa_cols), all_of(mirna_cols)) %>% drop_na()

cohort_E <- raw %>%
  select(Group, all_of(geno_cols), all_of(elisa_cols),
         all_of(mirna_cols)) %>% drop_na()

cat("Cohort A (SNP):         n=", nrow(cohort_A),
    " G0=", sum(cohort_A$Group=="Control"),
    " G1=", sum(cohort_A$Group=="HBV"), "\n")
cat("Cohort B (ELISA):       n=", nrow(cohort_B),
    " G0=", sum(cohort_B$Group=="Control"),
    " G1=", sum(cohort_B$Group=="HBV"), "\n")
cat("Cohort C (miRNA):       n=", nrow(cohort_C),
    " G0=", sum(cohort_C$Group=="Control"),
    " G1=", sum(cohort_C$Group=="HBV"), "\n")
cat("Cohort D (ELISA+miRNA): n=", nrow(cohort_D),
    " G0=", sum(cohort_D$Group=="Control"),
    " G1=", sum(cohort_D$Group=="HBV"), "\n")
cat("Cohort E (ALL):         n=", nrow(cohort_E),
    " G0=", sum(cohort_E$Group=="Control"),
    " G1=", sum(cohort_E$Group=="HBV"), "\n\n")


# ---- 4. LOG-TRANSFORM CONTINUOUS FEATURES -----------------
# ELISA values span orders of magnitude → log10 transform
# miRNA → log2 transform (fold-change scale)

log_transform <- function(df, cols, base = "log10") {
  df %>% mutate(across(
    all_of(cols),
    ~ if (base == "log10") log10(. + 0.001) else log2(. + 0.001)
  ))
}

cohort_B_log <- log_transform(cohort_B, elisa_cols, "log10")
cohort_C_log <- log_transform(cohort_C, mirna_cols, "log2")
cohort_D_log <- cohort_D %>%
  log_transform(elisa_cols, "log10") %>%
  log_transform(mirna_cols, "log2")
cohort_E_log <- cohort_E %>%
  log_transform(elisa_cols, "log10") %>%
  log_transform(mirna_cols, "log2")


# ====================================================
#  PART I — STATISTICAL ANALYSIS
# ====================================================

# ---- 5A. MANN-WHITNEY U TEST — ELISA -------------------

cat("===== ELISA Mann-Whitney Results =====\n")
elisa_stats <- map_dfr(seq_along(elisa_cols), function(i) {
  col  <- elisa_cols[i]
  lab  <- elisa_labels[i]
  data <- cohort_B_log %>% select(Group, value = all_of(col))

  test <- wilcox.test(value ~ Group, data = data,
                      exact = FALSE, conf.int = TRUE)

  # Median (IQR) per group
  med <- data %>%
    group_by(Group) %>%
    summarise(med = median(value, na.rm = TRUE),
              q1  = quantile(value, 0.25, na.rm = TRUE),
              q3  = quantile(value, 0.75, na.rm = TRUE),
              .groups = "drop")

  tibble(
    Feature      = lab,
    Control_med  = round(med$med[med$Group == "Control"], 3),
    Control_IQR  = paste0(round(med$q1[med$Group=="Control"],3),
                          "–",
                          round(med$q3[med$Group=="Control"],3)),
    HBV_med      = round(med$med[med$Group == "HBV"], 3),
    HBV_IQR      = paste0(round(med$q1[med$Group=="HBV"],3),
                          "–",
                          round(med$q3[med$Group=="HBV"],3)),
    W_statistic  = test$statistic,
    P_value      = test$p.value
  )
})

# FDR correction (Benjamini-Hochberg)
elisa_stats$P_adj <- p.adjust(elisa_stats$P_value, method = "BH")
elisa_stats$Significance <- case_when(
  elisa_stats$P_adj < 0.001 ~ "***",
  elisa_stats$P_adj < 0.01  ~ "**",
  elisa_stats$P_adj < 0.05  ~ "*",
  TRUE                       ~ "ns"
)

print(elisa_stats)
cat("\n")


# ---- 5B. MANN-WHITNEY U TEST — miRNA -------------------

cat("===== miRNA Mann-Whitney Results =====\n")
mirna_stats <- map_dfr(seq_along(mirna_cols), function(i) {
  col  <- mirna_cols[i]
  lab  <- mirna_labels[i]
  data <- cohort_C_log %>% select(Group, value = all_of(col))

  test <- wilcox.test(value ~ Group, data = data,
                      exact = FALSE, conf.int = TRUE)

  med <- data %>%
    group_by(Group) %>%
    summarise(med = median(value, na.rm = TRUE),
              q1  = quantile(value, 0.25, na.rm = TRUE),
              q3  = quantile(value, 0.75, na.rm = TRUE),
              .groups = "drop")

  tibble(
    Feature      = lab,
    Control_med  = round(med$med[med$Group == "Control"], 3),
    HBV_med      = round(med$med[med$Group == "HBV"], 3),
    W_statistic  = test$statistic,
    P_value      = test$p.value
  )
})

mirna_stats$P_adj <- p.adjust(mirna_stats$P_value, method = "BH")
mirna_stats$Significance <- case_when(
  mirna_stats$P_adj < 0.001 ~ "***",
  mirna_stats$P_adj < 0.01  ~ "**",
  mirna_stats$P_adj < 0.05  ~ "*",
  TRUE                       ~ "ns"
)

print(mirna_stats)
cat("\n")


# ---- 5C. CHI-SQUARE + ODDS RATIO — SNPs ----------------

cat("===== SNP Chi-square + OR Results =====\n")
snp_stats <- map_dfr(seq_along(geno_cols), function(i) {
  col <- geno_cols[i]
  lab <- geno_labels[i]

  data <- cohort_A %>%
    select(Group, geno = all_of(col)) %>%
    drop_na() %>%
    mutate(geno = factor(geno, levels = c(1,2,3),
                         labels = c("Maj/Maj","Het","Min/Min")))

  tbl <- table(data$Group, data$geno)

  # Chi-square
  chi <- tryCatch(
    chisq.test(tbl, simulate.p.value = TRUE, B = 2000),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  # Genotype frequencies
  freq <- prop.table(tbl, margin = 1) * 100

  tibble(
    SNP          = lab,
    Control_GG   = round(freq["Control", "Maj/Maj"], 1),
    Control_GA   = round(freq["Control", "Het"],     1),
    Control_AA   = round(freq["Control", "Min/Min"], 1),
    HBV_GG       = round(freq["HBV", "Maj/Maj"],    1),
    HBV_GA       = round(freq["HBV", "Het"],         1),
    HBV_AA       = round(freq["HBV", "Min/Min"],     1),
    Chi2         = round(as.numeric(chi$statistic), 3),
    P_value      = chi$p.value
  )
})

snp_stats$P_adj <- p.adjust(snp_stats$P_value, method = "BH")
snp_stats$Significance <- case_when(
  snp_stats$P_adj < 0.001 ~ "***",
  snp_stats$P_adj < 0.01  ~ "**",
  snp_stats$P_adj < 0.05  ~ "*",
  TRUE                     ~ "ns"
)

print(snp_stats)
cat("\n")


# ---- 5D. HARDY-WEINBERG EQUILIBRIUM (Controls only) ----

cat("===== Hardy-Weinberg Equilibrium (Controls) =====\n")
hwe_results <- map_dfr(seq_along(geno_cols), function(i) {
  col <- geno_cols[i]
  lab <- geno_labels[i]

  ctrl <- cohort_A %>%
    filter(Group == "Control") %>%
    pull(all_of(col)) %>%
    table()

  # Expects counts: c(AA, AB, BB)
  counts <- c(
    as.numeric(ctrl["1"]),   # homozygous major
    as.numeric(ctrl["2"]),   # heterozygous
    as.numeric(ctrl["3"])    # homozygous minor
  )
  counts[is.na(counts)] <- 0

  hwe <- tryCatch(
    HWChisq(counts, verbose = FALSE),
    error = function(e) list(chisq = NA, pval = NA)
  )

  tibble(
    SNP     = lab,
    n_GG    = counts[1],
    n_GA    = counts[2],
    n_AA    = counts[3],
    HWE_chi = round(hwe$chisq, 3),
    HWE_p   = round(hwe$pval,  4)
  )
})

hwe_results$HWE_status <- ifelse(
  hwe_results$HWE_p > 0.05, "In HWE", "Deviation"
)

print(hwe_results)
cat("\n")


# ====================================================
#  PART II — VISUALIZATION (Publication Figures)
# ====================================================

# Theme for all plots
theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.text        = element_text(color = "black", size = 11),
      axis.title       = element_text(size = 12, face = "bold"),
      legend.position  = "right",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.text      = element_text(size = 10),
      plot.title       = element_text(size = 13, face = "bold",
                                      hjust = 0.5),
      strip.text       = element_text(size = 11, face = "bold"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.4)
    )
}

group_colors <- c("Control" = "#2166AC", "HBV" = "#D6604D")


# ---- FIGURE 1 — ELISA Boxplots -------------------------

# Reshape to long format
elisa_long <- cohort_B_log %>%
  pivot_longer(all_of(elisa_cols),
               names_to  = "Cytokine",
               values_to = "log10_level") %>%
  mutate(Cytokine = factor(Cytokine,
                           levels = elisa_cols,
                           labels = elisa_labels))

# Add significance labels from stats
sig_df <- elisa_stats %>%
  select(Feature, Significance) %>%
  rename(Cytokine = Feature)

# Get y positions for significance brackets
y_pos <- elisa_long %>%
  group_by(Cytokine) %>%
  summarise(y.position = max(log10_level, na.rm = TRUE) + 0.15,
            .groups = "drop")

sig_df <- sig_df %>% left_join(y_pos, by = "Cytokine")

fig1 <- ggplot(elisa_long,
               aes(x = Group, y = log10_level, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
               width = 0.55, alpha = 0.85, linewidth = 0.5) +
  geom_jitter(aes(color = Group), width = 0.15, size = 0.8,
              alpha = 0.4) +
  facet_wrap(~ Cytokine, scales = "free_y", nrow = 1) +
  stat_pvalue_manual(sig_df, label = "Significance",
                     tip.length = 0.02, size = 4.5,
                     bracket.size = 0.5) +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Differential Cytokine Expression: HBV vs Controls",
    x     = NULL,
    y     = "log₁₀ (Cytokine level)",
    fill  = "Group",
    color = "Group"
  ) +
  theme_pub() +
  theme(legend.position = "top")

ggsave("Figure1_ELISA_boxplots.pdf", fig1,
       width = 14, height = 5, dpi = 300)
ggsave("Figure1_ELISA_boxplots.png", fig1,
       width = 14, height = 5, dpi = 300)
cat("Figure 1 saved.\n")


# ---- FIGURE 2 — miRNA Boxplots -------------------------

mirna_long <- cohort_C_log %>%
  pivot_longer(all_of(mirna_cols),
               names_to  = "miRNA",
               values_to = "log2_expr") %>%
  mutate(miRNA = factor(miRNA,
                        levels = mirna_cols,
                        labels = mirna_labels))

sig_mirna <- mirna_stats %>%
  select(Feature, Significance) %>%
  rename(miRNA = Feature)

y_pos_mirna <- mirna_long %>%
  group_by(miRNA) %>%
  summarise(y.position = max(log2_expr, na.rm = TRUE) + 0.2,
            .groups = "drop")

sig_mirna <- sig_mirna %>% left_join(y_pos_mirna, by = "miRNA")

fig2 <- ggplot(mirna_long,
               aes(x = Group, y = log2_expr, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1.5,
               width = 0.55, alpha = 0.85, linewidth = 0.5) +
  geom_jitter(aes(color = Group), width = 0.15, size = 0.8,
              alpha = 0.4) +
  facet_wrap(~ miRNA, scales = "free_y", nrow = 2) +
  stat_pvalue_manual(sig_mirna, label = "Significance",
                     tip.length = 0.02, size = 4) +
  scale_fill_manual(values  = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Differential miRNA Expression: HBV vs Controls",
    x     = NULL,
    y     = "log₂ (relative expression)",
    fill  = "Group",
    color = "Group"
  ) +
  theme_pub() +
  theme(legend.position = "top")

ggsave("Figure2_miRNA_boxplots.pdf", fig2,
       width = 16, height = 8, dpi = 300)
ggsave("Figure2_miRNA_boxplots.png", fig2,
       width = 16, height = 8, dpi = 300)
cat("Figure 2 saved.\n")


# ---- FIGURE 3 — PCA (Cytokines + miRNA) ----------------

# PCA on ELISA (Cohort B, log-transformed)
pca_elisa <- PCA(
  cohort_B_log %>% select(all_of(elisa_cols)) %>%
    scale() %>% as.data.frame(),
  graph = FALSE, ncp = 5
)

pca_elisa_df <- as.data.frame(pca_elisa$ind$coord) %>%
  mutate(Group = cohort_B_log$Group)

var_elisa <- round(pca_elisa$eig[1:2, 2], 1)

p_pca_elisa <- ggplot(pca_elisa_df,
                      aes(x = Dim.1, y = Dim.2,
                          color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values  = group_colors) +
  labs(
    title = "PCA — Cytokine profiles",
    x = paste0("PC1 (", var_elisa[1], "%)"),
    y = paste0("PC2 (", var_elisa[2], "%)")
  ) +
  theme_pub()

# PCA on miRNA (Cohort C, log-transformed)
pca_mirna <- PCA(
  cohort_C_log %>% select(all_of(mirna_cols)) %>%
    scale() %>% as.data.frame(),
  graph = FALSE, ncp = 5
)

pca_mirna_df <- as.data.frame(pca_mirna$ind$coord) %>%
  mutate(Group = cohort_C_log$Group)

var_mirna <- round(pca_mirna$eig[1:2, 2], 1)

p_pca_mirna <- ggplot(pca_mirna_df,
                      aes(x = Dim.1, y = Dim.2,
                          color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values  = group_colors) +
  labs(
    title = "PCA — miRNA profiles",
    x = paste0("PC1 (", var_mirna[1], "%)"),
    y = paste0("PC2 (", var_mirna[2], "%)")
  ) +
  theme_pub()

fig3 <- p_pca_elisa + p_pca_mirna +
  plot_annotation(
    title = "Figure 3 — Multivariate Immune Profiling (PCA)",
    theme = theme(plot.title = element_text(size = 14,
                                            face = "bold",
                                            hjust = 0.5))
  )

ggsave("Figure3_PCA.pdf", fig3, width = 12, height = 5, dpi = 300)
ggsave("Figure3_PCA.png", fig3, width = 12, height = 5, dpi = 300)
cat("Figure 3 saved.\n")


# ---- FIGURE 4 — Cytokine–miRNA Correlation Heatmap -----
# Uses Cohort D (ELISA + miRNA, log-transformed)

cor_matrix <- cor(
  cohort_D_log %>% select(all_of(elisa_cols), all_of(mirna_cols)),
  method = "spearman",
  use    = "complete.obs"
)

# Subset: cytokines (rows) x miRNA (cols)
cor_sub <- cor_matrix[elisa_labels, mirna_labels]
rownames(cor_sub) <- elisa_labels
colnames(cor_sub) <- mirna_labels

pdf("Figure4_CytokineMiRNA_Heatmap.pdf", width = 10, height = 5)
corrplot(cor_sub,
         method     = "color",
         type       = "full",
         col        = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(200),
         addCoef.col = "black",
         number.cex  = 0.7,
         tl.col      = "black",
         tl.srt      = 45,
         cl.lim      = c(-1, 1),
         title       = "Spearman Correlation: Cytokines vs miRNAs",
         mar         = c(0, 0, 2, 0))
dev.off()
cat("Figure 4 saved.\n")


# ---- FIGURE 5 — Global Correlation Heatmap (all features) -

cor_all <- cor(
  cohort_D_log %>%
    select(all_of(elisa_cols), all_of(mirna_cols)) %>%
    rename_with(~ elisa_labels, all_of(elisa_cols)) %>%
    rename_with(~ mirna_labels, all_of(mirna_cols)),
  method = "spearman",
  use    = "complete.obs"
)

pdf("Figure5_Global_Correlation_Heatmap.pdf",
    width = 12, height = 11)
corrplot(cor_all,
         method      = "color",
         type        = "upper",
         order       = "hclust",
         col         = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(200),
         addCoef.col  = "black",
         number.cex   = 0.65,
         tl.col       = "black",
         tl.srt       = 45,
         diag         = FALSE,
         title        = "Global Spearman Correlation: All Features",
         mar          = c(0, 0, 2, 0))
dev.off()
cat("Figure 5 saved.\n")


# ---- FIGURE 6 — Hierarchical Clustering Heatmap --------
# Cohort D, scaled features, annotated by Group

library(pheatmap)

mat <- cohort_D_log %>%
  select(all_of(elisa_cols), all_of(mirna_cols)) %>%
  as.matrix() %>%
  t()

rownames(mat) <- c(elisa_labels, mirna_labels)

# Annotation bar
ann <- data.frame(Group = cohort_D_log$Group)
rownames(ann) <- paste0("S", seq_len(nrow(cohort_D_log)))
colnames(mat) <- rownames(ann)

ann_colors <- list(Group = c("Control" = "#2166AC",
                              "HBV"     = "#D6604D"))

pdf("Figure6_Clustering_Heatmap.pdf", width = 14, height = 6)
pheatmap(mat,
         annotation_col   = ann,
         annotation_colors = ann_colors,
         scale            = "row",
         clustering_method = "ward.D2",
         color            = colorRampPalette(
           rev(brewer.pal(11, "RdBu")))(100),
         show_colnames    = FALSE,
         fontsize_row     = 10,
         main             = "Hierarchical Clustering — Cytokine + miRNA Profiles")
dev.off()
cat("Figure 6 saved.\n")


# ---- FIGURE 7 — SNP Genotype Frequency Barplot ----------

snp_long <- cohort_A %>%
  pivot_longer(all_of(geno_cols),
               names_to  = "SNP",
               values_to = "Genotype") %>%
  mutate(
    SNP      = factor(SNP, levels = geno_cols, labels = geno_labels),
    Genotype = factor(Genotype, levels = c(1,2,3),
                      labels = c("Maj/Maj","Het","Min/Min"))
  ) %>%
  drop_na() %>%
  count(Group, SNP, Genotype) %>%
  group_by(Group, SNP) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

fig7 <- ggplot(snp_long,
               aes(x = SNP, y = pct, fill = Genotype,
                   alpha = Group)) +
  geom_bar(stat = "identity", position = "dodge",
           linewidth = 0.3, color = "white") +
  facet_wrap(~ Group, nrow = 2) +
  scale_fill_manual(values = c("#4DAF4A","#FF7F00","#E41A1C")) +
  scale_alpha_manual(values = c(1, 0.85)) +
  labs(
    title = "SNP Genotype Frequencies: HBV vs Controls",
    x     = "SNP",
    y     = "Frequency (%)",
    fill  = "Genotype"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                   size = 9),
        legend.position = "top")

ggsave("Figure7_SNP_Genotypes.pdf", fig7,
       width = 14, height = 8, dpi = 300)
ggsave("Figure7_SNP_Genotypes.png", fig7,
       width = 14, height = 8, dpi = 300)
cat("Figure 7 saved.\n")


# ====================================================
#  PART III — EXPORT RESULTS & ML-READY CSVs
# ====================================================

# ---- Save statistical results -------------------------

stats_all <- list(
  ELISA = elisa_stats,
  miRNA = mirna_stats,
  SNP   = snp_stats,
  HWE   = hwe_results
)

write.csv(elisa_stats, "stats_ELISA.csv",   row.names = FALSE)
write.csv(mirna_stats, "stats_miRNA.csv",   row.names = FALSE)
write.csv(snp_stats,   "stats_SNP.csv",     row.names = FALSE)
write.csv(hwe_results, "stats_HWE.csv",     row.names = FALSE)
cat("Statistical results saved.\n")


# ---- Export ML-ready CSVs for Python Colab ------------
# Cohort D: ELISA + miRNA (log-transformed, no missing)
# Cohort E: ELISA + miRNA + SNP (log-transformed, no missing)

ml_D <- cohort_D_log %>%
  mutate(Group = as.integer(Group == "HBV"))

ml_E <- cohort_E_log %>%
  mutate(Group = as.integer(Group == "HBV"))

write.csv(ml_D, "cohort_D_ELISA_miRNA_ml_ready.csv",
          row.names = FALSE)
write.csv(ml_E, "cohort_E_ALL_ml_ready.csv",
          row.names = FALSE)

cat("\n=== ML-ready files exported ===\n")
cat("cohort_D_ELISA_miRNA_ml_ready.csv  n=", nrow(ml_D),
    "  features=", ncol(ml_D)-1, "\n")
cat("cohort_E_ALL_ml_ready.csv          n=", nrow(ml_E),
    "  features=", ncol(ml_E)-1, "\n")

cat("\n✓ R analysis complete. All figures and CSVs saved.\n")
cat("  → Upload the two ML CSV files to Google Colab.\n")
