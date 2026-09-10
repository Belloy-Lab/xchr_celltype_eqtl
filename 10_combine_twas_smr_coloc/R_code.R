# Merge TWAS + SMR + COLOC; support flags; per-trait heatmaps and a stacked AD/PD/MSA figure with right-side annotations
rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
# CONFIG
twas_rds  <- "../07_twas/7.8_combine_twas_results/twas_res_distinct.rds"
smr_rds   <- "../08_smr/8.4_combine_results/smr_res_distinct.rds"
coloc_rds <- "../09_coloc/9.3_combine_results/coloc_res_distinct.rds"
deg_xlsx  <- "../06.3_sex_DEG/DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx"
eqtl_rds  <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
fdr_smr_cut    <- 0.05   # SMR support
heidi_cut      <- 0.05   # SMR support: p_HEIDI NA or above
coloc_pp4_cut  <- 0.6   # COLOC support
pp4_label_cut  <- 0.5   # PP4 label shown above this
deg_p_cut      <- 0.05
het_p_cut      <- 0.05
zlim_quantile  <- 0.99   # heatmap colour saturation
target_traits  <- c("AD", "PD", "MSA")
# gene groups shown as separate row slices in the stacked figure
box_groups <- list(
  AD  = list(c("KRBOX4", "CHST7", "SLC9A7", "JADE3")),
  PD  = list(c("RAB40A", "MORF4L2")),
  MSA = list(c("SLC6A8", "PLXNA3", "FAM3A", "IKBKG"))
)
cell_types <- c("Ast", "End", "Exc", "Inh", "Mic", "OPC", "Oli")
cell_type_colors <- c(
  Ast = "#1b9e77", End = "#d95f02", Exc = "#7570b3", Inh = "#e7298a",
  Mic = "#66a61e", OPC = "#e6ab02", Oli = "#a6761d"
)
models <- c("female", "male", "eXCI", "rXCI")
model_colors <- c(
  female = "#e31a1c", male = "#1f78b4", eXCI = "#ff7f00", rXCI = "#33a02c"
)
na_color <- "grey85"
twas_low_col  <- "#4B0082"
twas_high_col <- "#B8860B"
deg_levels <- c("Upregulated in female", "Upregulated in male", "No significant difference")
deg_colors <- c(`Upregulated in female`     = "#e31a1c",
                `Upregulated in male`       = "#1f78b4",
                `No significant difference` = na_color)
het_levels <- c("Female-biased (concordant)", "Male-biased (concordant)",
                "Female-dominant (opposite)", "Male-dominant (opposite)",
                "No sex het", "Unclassified")
het_colors <- c(`Female-biased (concordant)` = "#fb9a99",
                `Male-biased (concordant)`    = "#a6cee3",
                `Female-dominant (opposite)`  = "#e31a1c",
                `Male-dominant (opposite)`    = "#1f78b4",
                `No sex het`                  = na_color,
                `Unclassified`                = "grey60")
support_colors <- c(yes = "#333333", no = na_color)
# 1. Merge TWAS + SMR + COLOC
twas_res_distinct  <- readRDS(twas_rds)
smr_res_distinct   <- readRDS(smr_rds)
coloc_res_distinct <- readRDS(coloc_rds)
twas_smr_coloc_merged <- twas_res_distinct %>%
  left_join(smr_res_distinct,
            by = c("Gene", "gwas_model", "cell_type", "trait")) %>%
  left_join(coloc_res_distinct,
            by = c("Gene", "gwas_model", "cell_type", "trait"))
# 2. FDR and support flags per trait x cell type x model
twas_smr_coloc_merged <- twas_smr_coloc_merged %>%
  group_by(trait, cell_type, gwas_model) %>%
  mutate(
    FDR_SMR       = p.adjust(p_SMR, method = "fdr"),
    TWAS_SUPPORT  = !is.na(eqtl_support),
    SMR_SUPPORT   = FDR_SMR < fdr_smr_cut & (is.na(p_HEIDI) | p_HEIDI > heidi_cut),
    COLOC_SUPPORT = PP4 > coloc_pp4_cut
  ) %>%
  ungroup()
twas_smr_coloc_merged <- arrange(twas_smr_coloc_merged, CHR, gene_start)
# final support = TWAS and (SMR or COLOC)
twas_smr_coloc_merged$final_support <-
  twas_smr_coloc_merged$TWAS_SUPPORT &
  (twas_smr_coloc_merged$SMR_SUPPORT | twas_smr_coloc_merged$COLOC_SUPPORT)
# top model per trait x cell type x gene
twas_smr_coloc_merged <- twas_smr_coloc_merged %>%
  group_by(trait, cell_type, Gene) %>%
  mutate(
    top_model_in_group = if (any(final_support, na.rm = TRUE)) {
      final_support &
        abs(TWAS.Z_best) == max(abs(TWAS.Z_best)[final_support], na.rm = TRUE)
    } else {
      FALSE
    }
  ) %>%
  ungroup()
traits <- unique(twas_smr_coloc_merged$trait)
writexl::write_xlsx(twas_smr_coloc_merged, "twas_smr_coloc_merged.xlsx")
# 3. Pad matrices to the full column set
force_all_celltypes <- function(m, all_ct) {
  missing_ct <- setdiff(all_ct, colnames(m))
  if (length(missing_ct) > 0) {
    add_mat <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing_ct),
                      dimnames = list(rownames(m), missing_ct))
    m <- cbind(m, add_mat)
  }
  m[, all_ct, drop = FALSE]
}
force_all_cols <- function(m, all_cols) {
  missing_cols <- setdiff(all_cols, colnames(m))
  if (length(missing_cols) > 0) {
    add_mat <- matrix(NA_real_, nrow = nrow(m), ncol = length(missing_cols),
                      dimnames = list(rownames(m), missing_cols))
    m <- cbind(m, add_mat)
  }
  m[, all_cols, drop = FALSE]
}
# 4. Per-trait heatmaps: all genes / final-support genes / all models
for (i in seq_len(length(traits))) {
  twas_smr_coloc_trait <- filter(twas_smr_coloc_merged, trait == traits[i])
  
  twas_smr_coloc_trait_sig <- filter(twas_smr_coloc_trait, TWAS_SUPPORT == TRUE)
  ctwas_sig_gene <- unique(twas_smr_coloc_trait_sig$Gene)
  
  twas_smr_coloc_trait_multi <- filter(twas_smr_coloc_trait, Gene %in% ctwas_sig_gene)
  
  # strongest model per cell type
  twas_smr_coloc_trait_for_plot <- twas_smr_coloc_trait_multi %>%
    group_by(cell_type, Gene) %>%
    slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  twas_smr_coloc_trait_for_plot <- arrange(twas_smr_coloc_trait_for_plot, CHR, gene_start)
  
  df_loci <- distinct(twas_smr_coloc_trait_for_plot, Gene, gene_start)
  writexl::write_xlsx(df_loci, paste0("loci_", traits[i], ".xlsx"))
  
  if (nrow(twas_smr_coloc_trait_for_plot) == 0) next
  
  df <- twas_smr_coloc_trait_for_plot %>%
    mutate(x_axis = cell_type) %>%
    select(Gene, x_axis, gwas_model, TWAS.Z_best, TWAS_SUPPORT,
           SMR_SUPPORT, COLOC_SUPPORT, PP4, top_model_in_group) %>%
    distinct()
  
  df <- df %>%
    complete(Gene, x_axis,
             fill = list(TWAS.Z_best = NA_real_, TWAS_SUPPORT = FALSE,
                         SMR_SUPPORT = FALSE, COLOC_SUPPORT = FALSE,
                         PP4 = NA_real_, gwas_model = ""))
  
  mat <- df %>%
    select(Gene, x_axis, TWAS.Z_best) %>%
    pivot_wider(names_from = x_axis, values_from = TWAS.Z_best) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  mat <- force_all_celltypes(mat, cell_types)
  
  gene_order <- unique(twas_smr_coloc_trait_for_plot$Gene)
  gene_order <- gene_order[gene_order %in% rownames(mat)]
  mat <- mat[gene_order, , drop = FALSE]
  
  zlim <- as.numeric(quantile(abs(mat), zlim_quantile, na.rm = TRUE))
  if (is.na(zlim) || zlim == 0) zlim <- 1
  
  mat_plot <- mat
  mat_plot[!is.na(mat_plot) & mat_plot >  zlim] <-  zlim
  mat_plot[!is.na(mat_plot) & mat_plot < -zlim] <- -zlim
  
  col_fun <- colorRampPalette(c(twas_low_col, "white", twas_high_col))(100)
  
  df$key <- paste(df$Gene, df$x_axis, sep = "__")
  
  make_cell_fun <- function(mat_ref) {
    function(j, i, x, y, width, height, fill) {
      gene <- rownames(mat_ref)[i]; xax <- colnames(mat_ref)[j]
      key <- paste(gene, xax, sep = "__")
      idx <- which(df$key == key)
      if (length(idx) == 0) return()
      row <- df[idx[1], ]
      
      if (isTRUE(row$top_model_in_group)) {
        grid.rect(x = x, y = y, width = width * 0.95, height = height * 0.95,
                  gp = gpar(col = "black", lwd = 1.5, fill = NA))
      }
      if (row$gwas_model != "") {
        grid.text(row$gwas_model, x, y + unit(7, "pt"),
                  gp = gpar(fontsize = 4.5, col = "black"))
      }
      label_top <- paste0(ifelse(isTRUE(row$TWAS_SUPPORT), "*", ""),
                          ifelse(isTRUE(row$SMR_SUPPORT), "#", ""))
      if (label_top != "") {
        grid.text(label_top, x, y + unit(1, "pt"),
                  gp = gpar(fontsize = 7, col = "black"))
      }
      label_bottom <- ""
      if (isTRUE(row$COLOC_SUPPORT) && !is.na(row$PP4) && row$PP4 > pp4_label_cut) {
        label_bottom <- sprintf("%.2f", row$PP4)
      }
      if (label_bottom != "") {
        grid.text(label_bottom, x, y - unit(6, "pt"),
                  gp = gpar(fontsize = 5.5, col = "black"))
      }
    }
  }
  
  # Figure 1: all genes
  pdf(paste0(traits[i], "_twas_smr_coloc.pdf"), width = 7, height = 9)
  ht <- Heatmap(
    mat_plot, name = "TWAS Z", col = col_fun, na_col = na_color,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), row_names_side = "left",
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
    column_title = traits[i], column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cell_fun = make_cell_fun(mat_plot))
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  # Figure 3: all models, grouped by cell type
  df_multi <- twas_smr_coloc_trait_multi %>%
    mutate(x_axis = paste(gwas_model, cell_type, sep = "_")) %>%
    select(Gene, x_axis, TWAS.Z_best, TWAS_SUPPORT,
           SMR_SUPPORT, COLOC_SUPPORT, PP4, top_model_in_group) %>%
    distinct()
  df_multi$key <- paste(df_multi$Gene, df_multi$x_axis, sep = "__")
  
  x_full <- as.vector(outer(models, cell_types, paste, sep = "_"))
  
  mat_multi <- df_multi %>%
    select(Gene, x_axis, TWAS.Z_best) %>%
    pivot_wider(names_from = x_axis, values_from = TWAS.Z_best) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  mat_multi <- force_all_cols(mat_multi, x_full)
  
  gene_order_multi <- unique(twas_smr_coloc_trait_for_plot$Gene)
  gene_order_multi <- gene_order_multi[gene_order_multi %in% rownames(mat_multi)]
  mat_multi <- mat_multi[gene_order_multi, , drop = FALSE]
  
  zlim_multi <- as.numeric(quantile(abs(mat_multi), zlim_quantile, na.rm = TRUE))
  if (is.na(zlim_multi) || zlim_multi == 0) zlim_multi <- 1
  mat_multi_plot <- mat_multi
  mat_multi_plot[!is.na(mat_multi_plot) & mat_multi_plot >  zlim_multi] <-  zlim_multi
  mat_multi_plot[!is.na(mat_multi_plot) & mat_multi_plot < -zlim_multi] <- -zlim_multi
  
  cell_fun_multi <- function(j, i, x, y, width, height, fill) {
    gene <- rownames(mat_multi_plot)[i]; xax <- colnames(mat_multi_plot)[j]
    key <- paste(gene, xax, sep = "__")
    idx <- which(df_multi$key == key)
    if (length(idx) == 0) return()
    row <- df_multi[idx[1], ]
    
    if (isTRUE(row$top_model_in_group)) {
      grid.rect(x = x, y = y, width = width * 0.95, height = height * 0.95,
                gp = gpar(col = "black", lwd = 1.5, fill = NA))
    }
    label_top <- paste0(ifelse(isTRUE(row$TWAS_SUPPORT), "*", ""),
                        ifelse(isTRUE(row$SMR_SUPPORT), "#", ""))
    if (label_top != "") {
      grid.text(label_top, x, y + unit(2, "pt"),
                gp = gpar(fontsize = 6, col = "black"))
    }
    label_bottom <- ""
    if (isTRUE(row$COLOC_SUPPORT) && !is.na(row$PP4) && row$PP4 > pp4_label_cut) {
      label_bottom <- sprintf("%.2f", row$PP4)
    }
    if (label_bottom != "") {
      grid.text(label_bottom, x, y - unit(4, "pt"),
                gp = gpar(fontsize = 5, col = "black"))
    }
  }
  
  pdf(paste0(traits[i], "_twas_smr_coloc_all_models.pdf"), width = 12, height = 9)
  ht <- Heatmap(
    mat_multi_plot, name = "TWAS Z", col = col_fun, na_col = na_color,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), row_names_side = "left",
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
    column_title = traits[i], column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cell_fun = cell_fun_multi)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  # Figure 2: final-support genes with right annotation
  df_final_support_gene <- filter(twas_smr_coloc_trait_for_plot, final_support == TRUE)
  
  mat_plot_sig <- mat_plot[
    rownames(mat_plot) %in% df_final_support_gene$Gene, , drop = FALSE
  ]
  
  if (nrow(mat_plot_sig) == 0) next
  
  mat_plot_sig <- force_all_celltypes(mat_plot_sig, cell_types)
  
  # right annotation from the best cell type (max |Z|)
  anno_src <- twas_smr_coloc_trait_for_plot %>%
    filter(Gene %in% rownames(mat_plot_sig)) %>%
    group_by(Gene) %>%
    slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Gene, best_cell_type = cell_type, best_model = gwas_model,
           SMR_SUPPORT, COLOC_SUPPORT)
  
  anno_src <- anno_src[match(rownames(mat_plot_sig), anno_src$Gene), ]
  
  best_model <- anno_src$best_model
  smr_vec    <- ifelse(anno_src$SMR_SUPPORT   %in% TRUE, "yes", "no")
  coloc_vec  <- ifelse(anno_src$COLOC_SUPPORT %in% TRUE, "yes", "no")
  
  ha_right <- rowAnnotation(
    `best model` = best_model,
    `SMR`        = smr_vec,
    `COLOC`      = coloc_vec,
    col = list(
      `best model` = model_colors,
      `SMR`        = support_colors,
      `COLOC`      = support_colors
    ),
    gp = gpar(col = "black", lwd = 0.8),
    annotation_name_gp = gpar(fontsize = 7),
    annotation_name_rot = 45,
    simple_anno_size = unit(8, "mm"),
    gap = unit(1, "mm")
  )
  
  pdf(paste0(traits[i], "_twas_smr_coloc_only_significant.pdf"), width = 11, height = 3)
  
  ht <- Heatmap(
    mat_plot_sig,
    name = "TWAS Z",
    col = col_fun,
    na_col = na_color,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    row_names_side = "left",
    
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    
    column_title = traits[i],
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    right_annotation = ha_right,
    cell_fun = make_cell_fun(mat_plot_sig)
  )
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}
# 5. DEG sex bias and eQTL sex-het lookup per trait x gene, taken at the best cell type
gene_best_ct <- twas_smr_coloc_merged %>%
  filter(trait %in% target_traits, final_support == TRUE) %>%
  group_by(trait, Gene) %>%
  slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(trait, Gene, cell_type, eqtl_setting_coloc)
# DEG sex bias (nominal p; logFC > 0 = higher in female)
df_deg_extra <- readxl::read_excel(deg_xlsx) %>%
  mutate(
    DEG_sex_bias = case_when(
      logFC > 0 & P.Value < deg_p_cut ~ "Upregulated in female",
      logFC < 0 & P.Value < deg_p_cut ~ "Upregulated in male",
      TRUE ~ "No significant difference"
    )
  ) %>%
  select(cell_type, Gene = updated_Gene_symbol, DEG_sex_bias) %>%
  distinct(cell_type, Gene, .keep_all = TRUE)
# eQTL sex heterogeneity: top variant of each of the 4 models, female/male beta paired on the same variant, male /2, P_het
df_eqtl_results <- read_rds(eqtl_rds)
df_eqtl_results_sub <- semi_join(
  df_eqtl_results, gene_best_ct, by = c("gene" = "Gene", "cell_type")
)
het_top <- gene_best_ct %>% select(cell_type, Gene) %>% distinct()
het_top$top_rxci_variant   <- NA_character_
het_top$top_exci_variant   <- NA_character_
het_top$top_female_variant <- NA_character_
het_top$top_male_variant   <- NA_character_
for (j in seq_len(nrow(het_top))) {
  sub_j <- filter(df_eqtl_results_sub,
                  cell_type == het_top$cell_type[j],
                  gene      == het_top$Gene[j])
  sub_j <- arrange(sub_j, pvalue)
  het_top$top_rxci_variant[j]   <- filter(sub_j, setting == 'both_model_2')$ID[1]
  het_top$top_exci_variant[j]   <- filter(sub_j, setting == 'both_model_1')$ID[1]
  het_top$top_female_variant[j] <- filter(sub_j, setting == 'female_xchr_model_2')$ID[1]
  het_top$top_male_variant[j]   <- filter(sub_j, setting == 'male_xchr_model_1')$ID[1]
}
variant_long <- het_top %>%
  select(Gene, cell_type,
         top_rxci_variant, top_exci_variant,
         top_female_variant, top_male_variant) %>%
  pivot_longer(cols = starts_with("top_"), names_to = "model", values_to = "ID") %>%
  filter(!is.na(ID))
variant_long2 <- variant_long %>%
  group_by(Gene, cell_type, ID) %>%
  summarise(model = paste(unique(model), collapse = ";"), .groups = "drop")
df_eqtl_results_sub2 <- left_join(
  variant_long2, df_eqtl_results_sub,
  by = c("cell_type", "Gene" = "gene", "ID" = "ID")
)
eqtl_female <- df_eqtl_results_sub2 %>%
  filter(setting == "female_xchr_model_2") %>%
  select(Gene, cell_type, ID, beta_f = BETA, se_f = SE)
eqtl_male <- df_eqtl_results_sub2 %>%
  filter(setting == "male_xchr_model_1") %>%
  select(Gene, cell_type, ID, beta_m = BETA, se_m = SE)
eqtl_het <- inner_join(eqtl_female, eqtl_male,
                       by = c("Gene", "cell_type", "ID"))
eqtl_het$beta_m <- eqtl_het$beta_m / 2
eqtl_het$se_m   <- eqtl_het$se_m / 2
eqtl_het <- eqtl_het %>%
  mutate(
    Z_het = (beta_f - beta_m) / sqrt(se_f^2 + se_m^2),
    P_het = 2 * pnorm(-abs(Z_het)),
    opposite_dir = sign(beta_f) != sign(beta_m),
    delta_abs    = abs(beta_m) - abs(beta_f),
    sex_het = case_when(
      is.na(P_het)                     ~ NA_character_,
      P_het >= het_p_cut               ~ "No sex het",
      !opposite_dir & delta_abs > 0    ~ "Male-biased (concordant)",
      !opposite_dir & delta_abs < 0    ~ "Female-biased (concordant)",
      opposite_dir  & delta_abs > 0    ~ "Male-dominant (opposite)",
      opposite_dir  & delta_abs < 0    ~ "Female-dominant (opposite)",
      TRUE                             ~ "Unclassified"
    )
  )
# most significant het variant per cell type x gene
eqtl_het <- eqtl_het %>%
  arrange(P_het) %>%
  distinct(cell_type, Gene, .keep_all = TRUE) %>%
  select(cell_type, Gene, sex_het)
gene_extra <- gene_best_ct %>%
  left_join(df_deg_extra, by = c("cell_type", "Gene")) %>%
  left_join(eqtl_het,     by = c("cell_type", "Gene")) %>%
  mutate(
    DEG_sex_bias = ifelse(is.na(DEG_sex_bias), "No significant difference", DEG_sex_bias),
    sex_het      = ifelse(is.na(sex_het), "No sex het", sex_het)
  ) %>%
  select(trait, Gene, DEG_sex_bias, sex_het)
# 6. Stacked AD / PD / MSA figure, final-support genes only
target_traits_present <- target_traits[target_traits %in% traits]
make_cf <- function(mat_ref, df_ref) {
  function(j, i, x, y, width, height, fill) {
    gene <- rownames(mat_ref)[i]; xax <- colnames(mat_ref)[j]
    key  <- paste(gene, xax, sep = "__")
    idx  <- which(df_ref$key == key)
    if (length(idx) == 0) return()
    row <- df_ref[idx[1], ]
    if (isTRUE(row$top_model_in_group)) {
      grid.rect(x = x, y = y, width = width * 0.95, height = height * 0.95,
                gp = gpar(col = "black", lwd = 1.5, fill = NA))
    }
    label_top <- paste0(ifelse(isTRUE(row$TWAS_SUPPORT), "*", ""),
                        ifelse(isTRUE(row$SMR_SUPPORT), "#", ""))
    if (label_top != "") grid.text(label_top, x, y + unit(1, "pt"),
                                   gp = gpar(fontsize = 7, col = "black"))
    if (isTRUE(row$COLOC_SUPPORT) && !is.na(row$PP4) && row$PP4 > pp4_label_cut) {
      grid.text(sprintf("%.2f", row$PP4), x, y - unit(6, "pt"),
                gp = gpar(fontsize = 5.5, col = "black"))
    }
  }
}
build_trait_block <- function(trait_name) {
  d <- filter(twas_smr_coloc_merged, trait == trait_name)
  sig_gene <- unique(filter(d, TWAS_SUPPORT == TRUE)$Gene)
  d_for_plot <- filter(d, Gene %in% sig_gene) %>%
    group_by(cell_type, Gene) %>%
    slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(CHR, gene_start)
  if (nrow(d_for_plot) == 0) return(NULL)
  
  df_local <- d_for_plot %>%
    mutate(x_axis = cell_type) %>%
    select(Gene, x_axis, gwas_model, TWAS.Z_best, TWAS_SUPPORT,
           SMR_SUPPORT, COLOC_SUPPORT, PP4, top_model_in_group) %>%
    distinct() %>%
    complete(Gene, x_axis,
             fill = list(TWAS.Z_best = NA_real_, TWAS_SUPPORT = FALSE,
                         SMR_SUPPORT = FALSE, COLOC_SUPPORT = FALSE,
                         PP4 = NA_real_, gwas_model = ""))
  df_local$key <- paste(df_local$Gene, df_local$x_axis, sep = "__")
  
  m <- df_local %>%
    select(Gene, x_axis, TWAS.Z_best) %>%
    pivot_wider(names_from = x_axis, values_from = TWAS.Z_best) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  m <- force_all_celltypes(m, cell_types)
  
  gene_order <- unique(d_for_plot$Gene)
  gene_order <- gene_order[gene_order %in% rownames(m)]
  m <- m[gene_order, , drop = FALSE]
  
  sig_final <- unique(filter(d_for_plot, final_support == TRUE)$Gene)
  m <- m[rownames(m) %in% sig_final, , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  m <- force_all_celltypes(m, cell_types)
  
  anno_src <- d_for_plot %>%
    filter(Gene %in% rownames(m)) %>%
    group_by(Gene) %>%
    slice_max(abs(TWAS.Z_best), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Gene, best_cell_type = cell_type, best_model = gwas_model,
           SMR_SUPPORT, COLOC_SUPPORT)
  
  extra_this <- filter(gene_extra, trait == trait_name) %>%
    select(Gene, DEG_sex_bias, sex_het)
  anno_src <- left_join(anno_src, extra_this, by = "Gene")
  
  anno_src <- anno_src[match(rownames(m), anno_src$Gene), ]
  
  list(mat = m, df = df_local, anno = anno_src)
}
blocks <- lapply(target_traits_present, build_trait_block)
names(blocks) <- target_traits_present
present <- target_traits_present[!sapply(blocks, is.null)]
# row_split: each box_groups gene group in its own slice
make_box_split <- function(gene_names, groups) {
  n <- length(gene_names)
  if (n == 0) return(factor(integer(0)))
  if (is.null(groups)) groups <- list()
  membership <- rep(NA_integer_, n)
  for (gi in seq_along(groups)) {
    membership[gene_names %in% groups[[gi]]] <- gi
  }
  grp <- integer(n)
  g <- 1L; grp[1] <- 1L
  if (n > 1) {
    for (idx in 2:n) {
      prev <- membership[idx - 1]; cur <- membership[idx]
      same <- (is.na(prev) && is.na(cur)) ||
        (!is.na(prev) && !is.na(cur) && prev == cur)
      if (!same) g <- g + 1L
      grp[idx] <- g
    }
  }
  factor(grp, levels = unique(grp))
}
if (length(present) > 0) {
  # symmetric colour range from max |Z| across traits
  all_abs <- unlist(lapply(present, function(t) abs(as.vector(blocks[[t]]$mat))))
  zmax <- max(all_abs, na.rm = TRUE)
  if (is.na(zmax) || zmax == 0) zmax <- 1
  col_fun_combined <- colorRamp2(c(-zmax, 0, zmax), c(twas_low_col, "white", twas_high_col))
  
  ht_list <- NULL
  for (k in seq_along(present)) {
    t  <- present[k]
    bk <- blocks[[t]]
    m_plot <- bk$mat
    
    split_fac <- make_box_split(rownames(m_plot), box_groups[[t]])
    nsl <- nlevels(split_fac)
    rt  <- rep("", nsl)
    rt[ceiling(nsl / 2)] <- t   # trait name on the middle slice only
    
    is_first <- (k == 1)
    is_last  <- (k == length(present))
    
    ha_right <- rowAnnotation(
      `Best locus model`    = factor(bk$anno$best_model,     levels = models),
      `DEG sex bias`        = factor(bk$anno$DEG_sex_bias,   levels = deg_levels),
      `eQTL sex heterogeneity` = factor(bk$anno$sex_het,     levels = het_levels),
      `SMR`                 = factor(ifelse(bk$anno$SMR_SUPPORT   %in% TRUE, "yes", "no"), levels = c("yes","no")),
      `COLOC`               = factor(ifelse(bk$anno$COLOC_SUPPORT %in% TRUE, "yes", "no"), levels = c("yes","no")),
      col = list(`Best locus model`       = model_colors,
                 `DEG sex bias`           = deg_colors,
                 `eQTL sex heterogeneity` = het_colors,
                 `SMR`                    = support_colors,
                 `COLOC`                  = support_colors),
      gp = gpar(col = "black", lwd = 0.8),
      annotation_name_gp = gpar(fontsize = 7), annotation_name_rot = 45,
      simple_anno_size = unit(8, "mm"), gap = unit(1, "mm"),
      show_annotation_name = is_last,
      show_legend = is_first,
      annotation_legend_param = list(
        `Best locus model`       = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `DEG sex bias`           = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `eQTL sex heterogeneity` = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `SMR`                    = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8)),
        `COLOC`                  = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8))
      )
    )
    
    ht <- Heatmap(
      m_plot,
      name = if (is_first) "TWAS Z" else paste0("TWAS Z_", t),
      col = col_fun_combined, na_col = na_color,
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_split = split_fac,
      row_gap   = unit(0.8, "mm"),
      show_row_names = TRUE, row_names_gp = gpar(fontsize = 6), row_names_side = "left",
      show_column_names = is_last, column_names_gp = gpar(fontsize = 8),
      column_names_rot = 45,
      row_title = rt, row_title_gp = gpar(fontsize = 11, fontface = "bold"),
      right_annotation = ha_right,
      cell_fun = make_cf(m_plot, bk$df),
      show_heatmap_legend = is_first,
      heatmap_legend_param = list(
        direction = "horizontal",
        title_gp  = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8)
      )
    )
    
    ht_list <- if (is.null(ht_list)) ht else ht_list %v% ht
  }
  
  # height scales with gene count
  total_rows <- sum(sapply(present, function(t) nrow(blocks[[t]]$mat)))
  n_gaps_total <- sum(sapply(present, function(t)
    nlevels(make_box_split(rownames(blocks[[t]]$mat), box_groups[[t]])) - 1L))
  pdf("AD_PD_MSA_twas_smr_coloc_stacked.pdf",
      width = 11, height = max(6, 1.5 + total_rows * 0.3 + n_gaps_total * 0.08))
  draw(ht_list,
       heatmap_legend_side    = "top",
       annotation_legend_side = "right",
       column_title = NULL,
       column_title_gp = gpar(fontsize = 13, fontface = "bold"))
  dev.off()
}

# After running the code, the color of the best locus model was manually adjusted.
