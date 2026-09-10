# Independent non-PAR X signals: combine clumped significant eQTLs across 4 settings (wide), sex-heterogeneity test, second LD clumping, bar plots, supplementary table
rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ieugwasr)
# CONFIG
cis_all_rds <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
rsid_annotation_rds <- "../00.3_rsid_annotation/rsid_gene_annotation.rds"
bfile_prefix      <- "../00.1_X_chr_bfile_QC/chrX_QC_EU.0.G"
afreq_both_path   <- paste0(bfile_prefix, ".afreq")
afreq_female_path <- paste0(bfile_prefix, ".female.afreq")
afreq_male_path   <- paste0(bfile_prefix, ".male.afreq")
clump_dir              <- "../04_significant_eQTL_X/"
clump_both_model_1_rds <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_both_model_1.rds")
clump_both_model_2_rds <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_both_model_2.rds")
clump_female_rds       <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_female_xchr_model_2.rds")
clump_male_rds         <- paste0(clump_dir, "eqtl_sig_clumped_all_X_chr_2stepFDR_after_prune_male_xchr_model_1.rds")
# 1. Annotate cis-eQTL results with rsID, gene, allele frequency
df_eqtl_cis_all <- readRDS(cis_all_rds)
df_eqtl_cis_all$fdr <- NULL
df_rsid_gene_annotation <- readRDS(rsid_annotation_rds)
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,df_rsid_gene_annotation,by='ID')
afreq_both <- fread(afreq_both_path)
afreq_both$OBS_CT <- NULL
colnames(afreq_both)[5] <- "EAF_both"
afreq_female <- fread(afreq_female_path)
afreq_female$OBS_CT <- NULL
colnames(afreq_female)[5] <- "EAF_female"
afreq_male <- fread(afreq_male_path)
afreq_male$OBS_CT <- NULL
colnames(afreq_male)[5] <- "EAF_male"
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_both,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_female,by=c("#CHROM","ID", "REF", "ALT"))
df_eqtl_cis_all <- left_join(df_eqtl_cis_all,afreq_male,by=c("#CHROM","ID", "REF", "ALT"))
# 2. Keep cell-type-gene-variant pairs significant in >= 1 setting; pivot wide across settings
df_clump_sig_both_model_1 <- readRDS(clump_both_model_1_rds)
df_clump_sig_both_model_2 <- readRDS(clump_both_model_2_rds)
df_clump_sig_female_model_2 <- readRDS(clump_female_rds)
df_clump_sig_male_model_1 <- readRDS(clump_male_rds)
df_clump_sig_all_model <- rbind(df_clump_sig_both_model_1,
                                df_clump_sig_both_model_2,
                                df_clump_sig_female_model_2,
                                df_clump_sig_male_model_1)

df_eqtl_cis_clump_sig_all_model <- semi_join(df_eqtl_cis_all,df_clump_sig_all_model,
                                             by=c("ID","cell_type","gene"))
df_eqtl_cis_clump_sig_all_model <- df_eqtl_cis_clump_sig_all_model %>%
  left_join(
    df_clump_sig_all_model %>%
      select(ID, cell_type, gene, setting, significant_by_2step_FDR),
    by = c("ID", "cell_type", "gene", "setting")
  )
df_eqtl_cis_clump_sig_all_model$TEST <- NULL
df_eqtl_cis_clump_sig_all_model$ERRCODE <- NULL
# A1 should equal ALT
identical(df_eqtl_cis_clump_sig_all_model$A1,df_eqtl_cis_clump_sig_all_model$ALT)
df_eqtl_cis_clump_sig_all_model$A1 <- NULL

id_cols <- c(
  "#CHROM", "POS", "ID", "rsid","REF", "ALT",
  "cell_type", "gene", "par_type",
  "EAF_both", "EAF_female", "EAF_male"
)
value_cols <- c(
  "OBS_CT", "BETA", "SE", "T_STAT", "LOG10_P", "pvalue",
  "significant_by_2step_FDR"
)
df_eqtl_wide <- df_eqtl_cis_clump_sig_all_model %>%
  pivot_wider(
    id_cols   = all_of(id_cols),
    names_from  = setting,
    values_from = all_of(value_cols),
    names_glue  = "{.value}_{setting}"
  )
# 3. Min p-value across settings, male/female beta ratio, sex-heterogeneity test
p_cols <- grep("^pvalue_", colnames(df_eqtl_wide), value = TRUE)
df_eqtl_wide$min_pvalue <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) min(x, na.rm = TRUE)
)
df_eqtl_wide$min_pvalue_model <- apply(
  df_eqtl_wide[, p_cols],
  1,
  function(x) {
    p_cols[which.min(x)]
  }
)
df_eqtl_wide$min_pvalue_model <- sub("pvalue_", "", df_eqtl_wide$min_pvalue_model)

df_eqtl_wide$male_female_ratio <- (df_eqtl_wide$BETA_male_xchr_model_1)/(df_eqtl_wide$BETA_female_xchr_model_2)
df_eqtl_wide$abs_male_female_ratio <- abs(df_eqtl_wide$male_female_ratio)

eps <- 1e-6
df_eqtl_wide <- df_eqtl_wide %>%
  mutate(
    # male beta/SE (0/2 coding) rescaled to the female 0/1/2 scale
    male_beta_for_het = BETA_male_xchr_model_1 / 2,
    male_se_for_het   = SE_male_xchr_model_1 / 2,
    
    Z_het = (male_beta_for_het - BETA_female_xchr_model_2) /
      sqrt(male_se_for_het^2 + SE_female_xchr_model_2^2),
    P_het = 2 * pnorm(-abs(Z_het)),
    P_het_FDR = p.adjust(P_het, method = "BH"),
    
    male_sign   = sign(ifelse(abs(male_beta_for_het) < eps, 0, male_beta_for_het)),
    female_sign = sign(ifelse(abs(BETA_female_xchr_model_2) < eps, 0, BETA_female_xchr_model_2)),
    opposite_dir = (male_sign * female_sign) == -1,
    
    delta_abs = abs(male_beta_for_het) - abs(BETA_female_xchr_model_2)
  ) %>%
  mutate(
    sex_het_label_nominal = case_when(
      is.na(P_het) ~ NA_character_,
      P_het >= 0.05 ~ "No sex heterogeneity",
      
      !opposite_dir & delta_abs > 0 ~ "Male-biased (concordant)",
      !opposite_dir & delta_abs < 0 ~ "Female-biased (concordant)",
      
      opposite_dir & delta_abs > 0 ~ "Male-dominant (opposite)",
      opposite_dir & delta_abs < 0 ~ "Female-dominant (opposite)",
      
      TRUE ~ "Unclassified"
    ),
    
    sex_het_label_fdr = case_when(
      is.na(P_het_FDR) ~ NA_character_,
      P_het_FDR >= 0.05 ~ "No sex heterogeneity",
      
      !opposite_dir & delta_abs > 0 ~ "Male-biased (concordant)",
      !opposite_dir & delta_abs < 0 ~ "Female-biased (concordant)",
      
      opposite_dir & delta_abs > 0 ~ "Male-dominant (opposite)",
      opposite_dir & delta_abs < 0 ~ "Female-dominant (opposite)",
      
      TRUE ~ "Unclassified"
    )
  )

saveRDS(df_eqtl_wide,"df_eqtl_wide_before_second_clump.rds")
writexl::write_xlsx(df_eqtl_wide,"df_eqtl_wide_before_second_clump.xlsx")
df_eqtl_wide <- readRDS("df_eqtl_wide_before_second_clump.rds")
# 4. Second LD clumping across settings
column_order <- colnames(df_eqtl_wide)

df_eqtl_wide$snp <- df_eqtl_wide$rsid
df_eqtl_wide$rsid <- df_eqtl_wide$ID
df_eqtl_wide$pval <- df_eqtl_wide$min_pvalue
df_eqtl_wide$id <- paste0(df_eqtl_wide$gene,"_",df_eqtl_wide$cell_type)

df_eqtl_wide <- ld_clump(
  df_eqtl_wide,
  clump_kb  = 1000,
  clump_r2  = 0.1,
  clump_p   = 1,
  bfile     = bfile_prefix,
  plink_bin = plinkbinr::get_plink_exe()
)

df_eqtl_wide$rsid <- df_eqtl_wide$snp
df_eqtl_wide$snp <- NULL
df_eqtl_wide$pval <- NULL
df_eqtl_wide$id <- NULL
df_eqtl_wide <- select(df_eqtl_wide,column_order)

saveRDS(df_eqtl_wide,"df_eqtl_wide_after_second_clump.rds")
writexl::write_xlsx(df_eqtl_wide,"df_eqtl_wide_after_second_clump.xlsx")
# 5. Best-model bar plot per cell type: log2(|male/female beta ratio|), colored by best model
df_eqtl_wide_plot <- readRDS("df_eqtl_wide_after_second_clump.rds")
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  filter(!is.na(sex_het_label_nominal))

df_eqtl_wide_plot <- arrange(df_eqtl_wide_plot,desc(df_eqtl_wide_plot$abs_male_female_ratio))
df_eqtl_wide_plot$log_ratio <- log2(df_eqtl_wide_plot$abs_male_female_ratio)

df_eqtl_wide_plot$idx <- seq_len(nrow(df_eqtl_wide_plot))

df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = c("Exc","Inh","Ast","Oli","OPC","Mic","End")
    )
  ) %>%
  arrange(cell_type) %>%
  mutate(idx = row_number())

writexl::write_xlsx(df_eqtl_wide_plot,"non_PAR_variants_in_barplot_ordered_by_ratio_per_cell_type.xlsx")

cell_breaks <- df_eqtl_wide_plot %>%
  group_by(cell_type) %>%
  summarise(max_idx = max(idx)) %>%
  pull(max_idx)
cell_breaks <- cell_breaks + 0.5
cell_breaks <- c(0.5,cell_breaks)

ggplot(
  df_eqtl_wide_plot,
  aes(
    x = idx,
    y = log_ratio,
    fill = min_pvalue_model
  )
) +
  geom_col(
    width = 1,
    color = "black",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = log2(1),
    linetype = "solid",
    linewidth = 0.15,
    color = "black"
  ) +
  geom_hline(
    yintercept = log2(1.5),
    linetype = "dashed",
    linewidth = 0.15,
    color = "#ff7f00"
  ) +
  scale_fill_manual(
    values = c(
      "male_xchr_model_1"   = "#1f78b4",
      "female_xchr_model_2" = "#e31a1c",
      "both_model_2"        = "#33a02c",
      "both_model_1"        = "#ff7f00"
    ),
    breaks = c("male_xchr_model_1", "female_xchr_model_2", "both_model_1", "both_model_2"),
    labels = c("Male", "Female", "eXCI", "rXCI"),
    name = "Most significant model"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = paste0(" "),
    x = NULL,
    y = "log2 absolute value of\nMale / Female beta ratio"
  )+
  scale_y_continuous(
    breaks = seq(
      -5,
      10,
      by = 1
    )
  )+
  theme(
    legend.position = c(0.02, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.2
    )
  )+
  geom_vline(
    xintercept = cell_breaks,
    linetype = "dashed",
    linewidth = 0.15,
    color = "grey50"
  )

# variants without a male result are not shown
ggsave("Barplot_best_model_7_cell_type_panel_WIDE.pdf",width = 11,height = 3.9,dpi = 300)
# 6. Horizontal bar plot: variants per sex-het category, striped overlay = MAF >= 0.05 in all strata
df_eqtl_wide_plot <- readRDS("df_eqtl_wide_after_second_clump.rds")
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  filter(!is.na(sex_het_label_nominal))

# NA EAF_male counts as not common
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    maf_both    = pmin(EAF_both,   1 - EAF_both),
    maf_female  = pmin(EAF_female, 1 - EAF_female),
    maf_male    = pmin(EAF_male,   1 - EAF_male),
    common_all3 = (maf_both >= 0.05 & maf_female >= 0.05 & maf_male >= 0.05) %in% TRUE
  )

sex_het_colors_5 <- c(
  "Male-biased (concordant)"   = "#1f78b4",
  "Female-biased (concordant)" = "#e31a1c",
  "Male-dominant (opposite)"   = "#08306b",
  "Female-dominant (opposite)" = "#67000d",
  "No sex heterogeneity"       = "grey70"
)

# levels reversed: y axis runs bottom-up
sexhet_levels <- rev(c(
  "Male-biased (concordant)",
  "Male-dominant (opposite)",
  "Female-biased (concordant)",
  "Female-dominant (opposite)",
  "No sex heterogeneity"
))

bar_df2 <- df_eqtl_wide_plot %>%
  group_by(sex_het_label_nominal) %>%
  summarise(
    n_total  = n(),
    n_common = sum(common_all3),
    .groups  = "drop"
  ) %>%
  mutate(
    sex_het_label_nominal = factor(sex_het_label_nominal, levels = sexhet_levels)
  )

p_bar2 <- ggplot(bar_df2, aes(y = sex_het_label_nominal)) +
  geom_col(
    aes(x = n_total, fill = sex_het_label_nominal),
    width = 0.7, color = "black", linewidth = 0.2
  ) +
  geom_col_pattern(
    aes(x = n_common),
    fill            = NA,
    color           = "black",
    linewidth       = 0.2,
    width           = 0.7,
    pattern         = "stripe",
    pattern_fill    = "black",
    pattern_color   = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02
  ) +
  geom_text(
    aes(x = n_total, label = n_total),
    hjust = -0.3, size = 4, fontface = "bold"
  ) +
  scale_fill_manual(
    values = sex_het_colors_5,
    name = "Sex effect pattern"
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of independent lead variants",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p_bar2
ggsave("barplot_sex_het_horizontal_MAF.pdf", p_bar2, width = 5, height = 2.7, dpi = 300)
# 7. Horizontal bar plot: variants per best model, striped overlay = MAF >= 0.05 in all strata
df_eqtl_wide_plot <- readRDS("df_eqtl_wide_after_second_clump.rds")
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  filter(!is.na(sex_het_label_nominal))

# NA EAF_male counts as not common
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    maf_both    = pmin(EAF_both,   1 - EAF_both),
    maf_female  = pmin(EAF_female, 1 - EAF_female),
    maf_male    = pmin(EAF_male,   1 - EAF_male),
    common_all3 = (maf_both >= 0.05 & maf_female >= 0.05 & maf_male >= 0.05) %in% TRUE
  )

# levels reversed: y axis runs bottom-up
model_levels <- c("rXCI", "eXCI", "Female", "Male")

bar_df <- df_eqtl_wide_plot %>%
  group_by(min_pvalue_model) %>%
  summarise(
    n_total  = n(),
    n_common = sum(common_all3),
    .groups  = "drop"
  ) %>%
  mutate(
    model_label = recode(min_pvalue_model,
                         "both_model_1"        = "eXCI",
                         "both_model_2"        = "rXCI",
                         "female_xchr_model_2" = "Female",
                         "male_xchr_model_1"   = "Male"
    ),
    model_label = factor(model_label, levels = model_levels)
  )

p_bar <- ggplot(bar_df, aes(y = model_label)) +
  geom_col(
    aes(x = n_total, fill = model_label),
    width = 0.7, color = "black", linewidth = 0.2
  ) +
  geom_col_pattern(
    aes(x = n_common),
    fill            = NA,
    color           = "black",
    linewidth       = 0.2,
    width           = 0.7,
    pattern         = "stripe",
    pattern_fill    = "black",
    pattern_color   = "black",
    pattern_angle   = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.02
  ) +
  geom_text(
    aes(x = n_total, label = n_total),
    hjust = -0.3, size = 4, fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Male"   = "#1f78b4",
      "Female" = "#e31a1c",
      "rXCI"   = "#33a02c",
      "eXCI"   = "#ff7f00"
    ),
    name = "Most significant model"
  ) +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 14) +
  labs(
    x = "Number of independent lead variants",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p_bar
ggsave("barplot_best_model_variant_horizontal_MAF.pdf", p_bar, width = 5, height = 2.7, dpi = 300)
# 8. Sex-het bar plot per cell type: log2(|male/female beta ratio|), colored by sex-het category, * = P_het FDR < 0.05
df_eqtl_wide_plot <- readRDS("df_eqtl_wide_after_second_clump.rds")
df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  filter(!is.na(sex_het_label_nominal))

df_eqtl_wide_plot <- arrange(df_eqtl_wide_plot,desc(df_eqtl_wide_plot$abs_male_female_ratio))
df_eqtl_wide_plot$log_ratio <- log2(df_eqtl_wide_plot$abs_male_female_ratio)

df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    log_ratio_se = case_when(
      is.na(BETA_male_xchr_model_1) ~ NA_real_,
      is.na(BETA_female_xchr_model_2) ~ NA_real_,
      is.na(SE_male_xchr_model_1) ~ NA_real_,
      is.na(SE_female_xchr_model_2) ~ NA_real_,
      BETA_male_xchr_model_1 == 0 ~ NA_real_,
      BETA_female_xchr_model_2 == 0 ~ NA_real_,
      TRUE ~ sqrt(
        (SE_male_xchr_model_1 / BETA_male_xchr_model_1)^2 +
          (SE_female_xchr_model_2 / BETA_female_xchr_model_2)^2
      ) / log(2)
    ),
    log_ratio_lower_95 = log_ratio - 1.96 * log_ratio_se,
    log_ratio_upper_95 = log_ratio + 1.96 * log_ratio_se
  )

df_eqtl_wide_plot$idx <- seq_len(nrow(df_eqtl_wide_plot))

sex_het_colors_5 <- c(
  "Male-biased (concordant)"   = "#1f78b4",
  "Female-biased (concordant)" = "#e31a1c",
  "Male-dominant (opposite)"   = "#08306b",
  "Female-dominant (opposite)" = "#67000d",
  "No sex heterogeneity"       = "grey70"
)

df_eqtl_wide_plot <- df_eqtl_wide_plot %>%
  mutate(
    cell_type = factor(
      cell_type,
      levels = c("Exc","Inh","Ast","Oli","OPC","Mic","End")
    )
  ) %>%
  arrange(cell_type) %>%
  mutate(idx = row_number())

cell_breaks <- df_eqtl_wide_plot %>%
  group_by(cell_type) %>%
  summarise(max_idx = max(idx)) %>%
  pull(max_idx)
cell_breaks <- cell_breaks + 0.5
cell_breaks <- c(0.5,cell_breaks)

ggplot(
  df_eqtl_wide_plot,
  aes(
    x = idx,
    y = log_ratio,
    fill = sex_het_label_nominal
  )
) +
  geom_col(
    width = 1,
    color = "black",
    linewidth = 0.2
  ) +
  geom_hline(
    yintercept = log2(1.5),
    linetype = "dashed",
    linewidth = 0.15,
    color = "#1f78b4"
  ) +
  geom_hline(
    yintercept = log2(0.5),
    linetype = "dashed",
    linewidth = 0.15,
    color = "#e31a1c"
  ) +
  scale_fill_manual(
    values = sex_het_colors_5,
    name   = "Sex-effect pattern (P-HET<0.05)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = " ",
    x = NULL,
    y = "log2 absolute value of\nMale / Female beta ratio"
  ) +
  scale_y_continuous(
    breaks = seq(
      -5,
      10,
      by = 1
    )
  )+
  theme(
    legend.position = c(0.02, 0),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.2
    )
  )+
  geom_vline(
    xintercept = cell_breaks,
    linetype = "dashed",
    linewidth = 0.15,
    color = "grey50"
  )+
  geom_text(
    data = df_eqtl_wide_plot %>%
      filter(P_het_FDR < 0.05) %>%
      mutate(
        y_star = ifelse(
          log_ratio >= 0,
          log_ratio + 0.05,
          log_ratio - 0.15
        )
      ),
    aes(
      x = idx,
      y = y_star,
      label = "*"
    ),
    size = 2
  )

ggsave("Barplot_sex_het_7_cell_type_panel_WIDE.pdf",width = 11,height = 3.9,dpi = 300)
# 9. Supplementary table
rm(list=ls())
library(tidyverse)
df_eqtl_wide_plot <- readxl::read_excel("non_PAR_variants_in_barplot_ordered_by_ratio_per_cell_type.xlsx")
df_S4 <- df_eqtl_wide_plot %>%
  select(
    Gene        = gene,
    `PAR type`  = par_type,
    `Cell type` = cell_type,
    CHR         = `#CHROM`,
    POS,
    rsid,
    NEA         = REF,
    EA          = ALT,
    EAF_both, EAF_female, EAF_male,
    
    n_male    = OBS_CT_male_xchr_model_1,
    BETA_male = BETA_male_xchr_model_1,
    SE_male   = SE_male_xchr_model_1,
    P_male    = pvalue_male_xchr_model_1,
    
    n_female    = OBS_CT_female_xchr_model_2,
    BETA_female = BETA_female_xchr_model_2,
    SE_female   = SE_female_xchr_model_2,
    P_female    = pvalue_female_xchr_model_2,
    
    n_both_eXCI    = OBS_CT_both_model_1,
    BETA_both_eXCI = BETA_both_model_1,
    SE_both_eXCI   = SE_both_model_1,
    P_both_eXCI    = pvalue_both_model_1,
    
    n_both_rXCI    = OBS_CT_both_model_2,
    BETA_both_rXCI = BETA_both_model_2,
    SE_both_rXCI   = SE_both_model_2,
    P_both_rXCI    = pvalue_both_model_2,
    
    best_model = min_pvalue_model,
    P_best     = min_pvalue,
    
    Z_het,
    P_het,
    P_het_FDR,
    sex_het_label_nominal,
    sex_het_label_fdr
  )

writexl::write_xlsx(df_S4, "supplementary table.xlsx")