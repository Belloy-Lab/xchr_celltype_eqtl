# Sex DEG (female vs male) on chrX eQTL input genes per cell type, plus figures; case/control-adjusted model used only for the logFC correlation
rm(list = ls())
library(tidyverse)
library(edgeR)
library(limma)
library(data.table)
library(HGNChelper)
library(EnhancedVolcano)
library(patchwork)
library(ggsignif)
# CONFIG
metadata_xlsx <- "path/to/metadata.xlsx"
pseudobulk_dir <- "path/to/pseudobulk/"
cov_pheno_dir  <- "../01.2_cov_pheno/"
# order of the pseudobulk rds files
cell_types_analysis <- c("Ast", "End", "Exc", "Inh", "Mic", "Oli", "OPC")
# display order
cell_types <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "End")
deg_xlsx <- "DEG_female_vs_male_X_chr_gene_non_adj_disease.xlsx"
# 1. Metadata
meta_data <- readxl::read_excel(metadata_xlsx)
meta_data$msex <- factor(meta_data$msex, levels = c("1", "0"), labels = c("Male", "Female"))
meta_data$case_control <- ifelse(meta_data$cogdx == 1, "control",
                                 ifelse(meta_data$cogdx %in% 2:6, "case", NA))
meta_data <- select(meta_data, specimenID, batch, Study, age_death, msex, pmi, case_control)

count_files   <- dir(pseudobulk_dir, pattern = "rds")
analysis_list <- tibble(
  cell_types = cell_types_analysis,
  count_path = paste0(pseudobulk_dir, count_files)
)
analysis_list$cov_path <- paste0(cov_pheno_dir, "both_", analysis_list$cell_types, ".cov.txt")
# 2. DEG function: limma-voom, coef = msex; optional case/control adjustment
run_deg <- function(adjust_case_control = FALSE) {
  res_list <- vector("list", nrow(analysis_list))
  for (i in seq_len(nrow(analysis_list))) {
    cell_type <- analysis_list$cell_types[i]
    cat("Run DEG analysis (sex) in", cell_type,
        if (adjust_case_control) "[case/control-adjusted]" else "", "\n")
    df <- as.data.frame(read_rds(analysis_list$count_path[i]))
    
    # samples with < 10 cells are NA in the cov files
    valid_samples <- fread(analysis_list$cov_path[i])$IID
    meta <- meta_data %>%
      filter(specimenID %in% valid_samples,
             !is.na(Study), !is.na(msex), !is.na(age_death), !is.na(batch))
    meta$pmi[is.na(meta$pmi)] <- median(meta$pmi, na.rm = TRUE)
    df <- select(df, meta$specimenID)
    
    design <- if (adjust_case_control) {
      model.matrix(~ msex + Study + pmi + age_death + batch + case_control, data = meta)
    } else {
      model.matrix(~ msex + Study + pmi + age_death + batch, data = meta)
    }
    
    dge <- DGEList(counts = df)
    dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    v <- voom(dge, design)
    
    # drop mean log2CPM < 2, re-voom
    logcpm <- v$E[apply(v$E, 1, mean) > 2.0, ]
    dge$counts <- dge$counts[rownames(logcpm), , drop = FALSE]
    v <- voom(dge, design)
    
    fit <- eBayes(lmFit(v, design))
    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    res <- rownames_to_column(res, var = "Gene")
    res$gene_cell <- paste0(res$Gene, "_", cell_type)
    tbl <- summary(meta$msex)
    res$compare_summary <- paste(names(tbl), tbl, sep = ": ", collapse = ", ")
    
    # update symbols; duplicates keep the HGNC-approved one
    gsu <- checkGeneSymbols(res$Gene)
    res$updated_Gene_symbol <- gsu$Suggested.Symbol
    res$Approved <- gsu$Approved
    if (any(duplicated(res$updated_Gene_symbol))) {
      res <- res %>%
        group_by(updated_Gene_symbol) %>%
        arrange(desc(Approved == TRUE)) %>%
        slice(1) %>%
        ungroup()
    }
    res <- res[!duplicated(res$updated_Gene_symbol), ]
    res$gene_cell <- paste0(res$updated_Gene_symbol, "_", cell_type)
    res$cell_type <- cell_type
    colnames(res)[1] <- "Old_gene_symbol"
    
    # restrict to eQTL input genes; FDR within them
    pheno <- fread(paste0(cov_pheno_dir, "both_", cell_type, ".pheno_chrX.tsv"))
    pheno$FID <- NULL; pheno$IID <- NULL
    res <- filter(res, updated_Gene_symbol %in% colnames(pheno))
    
    res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
    res_list[[i]] <- res
  }
  rbindlist(res_list, fill = TRUE)
}
# 3. Non-adjusted DEG
DEG_all <- run_deg(adjust_case_control = FALSE)
writexl::write_xlsx(DEG_all, deg_xlsx)
DEG_all <- readxl::read_excel(deg_xlsx)   # normalize column types
# 4. Volcano plots; non-End panels use compressed x/y axes beyond a breakpoint
brk      <- 1
cmp      <- 0.05
x_breaks <- c(-2, -1, 0, 1, 2, 3, 4, 5, 6)
y_brk    <- 100
y_cmp    <- 0.05
y_breaks <- c(0, 50, 100, 200, 300)
y_max    <- 300
# End: real coordinates
end_x_breaks <- c(0.5)
end_x_right  <- 0.6
# margins in pt
panel_margin_t <- 2.5
panel_margin_r <- 3
panel_margin_b <- 2.5
panel_margin_l <- 2.5
title_margin_b <- 1.5
xtitle_margin_t <-  0
ytitle_margin_r <- -8
axis_text_size <- 15
squish_x <- scales::trans_new(
  name      = "squish_x",
  transform = function(x) ifelse(x <= brk, x, brk + (x - brk) * cmp),
  inverse   = function(x) ifelse(x <= brk, x, brk + (x - brk) / cmp)
)
squish_y <- scales::trans_new(
  name      = "squish_y",
  transform = function(y) ifelse(y <= y_brk, y, y_brk + (y - y_brk) * y_cmp),
  inverse   = function(y) ifelse(y <= y_brk, y, y_brk + (y - y_brk) / y_cmp)
)
plot_list <- list()
for (ct in cell_types) {
  cat("generating volcano plot for:", ct, "\n")
  
  df_plot <- dplyr::filter(DEG_all, cell_type == ct)
  
  df_plot$color_group <- "Not sig"
  df_plot$color_group[df_plot$adj.P.Val < 0.05 & df_plot$logFC > 0] <- "Female high expression"
  df_plot$color_group[df_plot$adj.P.Val < 0.05 & df_plot$logFC < 0] <- "Male high expression"
  
  keyvals <- ifelse(df_plot$adj.P.Val >= 0.05, "grey80",
                    ifelse(df_plot$logFC > 0, "red", "blue"))
  names(keyvals)[keyvals == "red"]    <- "Female-upregulated DEGs"
  names(keyvals)[keyvals == "blue"]   <- "Male-upregulated DEGs"
  names(keyvals)[keyvals == "grey80"] <- "No differential expression"
  
  p <- EnhancedVolcano(
    df_plot,
    lab = df_plot$updated_Gene_symbol,
    x = "logFC",
    y = "adj.P.Val",
    pCutoff = 0.05,
    FCcutoff = 0,
    colCustom = keyvals,
    pointSize = 1.5,
    labSize = 3,
    axisLabSize = axis_text_size,
    title = ct,
    subtitle = NULL,
    drawConnectors = TRUE,
    caption = NULL,
    max.overlaps = 20,
    ylab = expression(-Log[10]("adj.P.Val"))
  )
  
  if (ct != "End") {
    p <- p + ggplot2::scale_x_continuous(
      trans  = squish_x,
      breaks = x_breaks,
      expand = ggplot2::expansion(mult = 0.05)
    )
    p <- p + ggplot2::scale_y_continuous(
      trans  = squish_y,
      breaks = y_breaks,
      limits = c(0, y_max),
      expand = ggplot2::expansion(mult = 0.02)
    )
  } else {
    end_p   <- -log10(df_plot$adj.P.Val); end_p <- end_p[is.finite(end_p)]
    x_left  <- min(df_plot$logFC, na.rm = TRUE) - 0.1
    p <- p +
      ggplot2::scale_x_continuous(
        breaks = end_x_breaks,
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::coord_cartesian(
        xlim = c(x_left, end_x_right),
        ylim = c(0, max(end_p, na.rm = TRUE) + 5)
      )
  }
  
  plot_list[[ct]] <- p
}
combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(margin = ggplot2::margin(t = xtitle_margin_t)),
    axis.title.y = element_text(margin = ggplot2::margin(r = ytitle_margin_r)),
    plot.margin  = ggplot2::margin(panel_margin_t, panel_margin_r,
                                   panel_margin_b, panel_margin_l),
    plot.title   = element_text(margin = ggplot2::margin(b = title_margin_b))
  )
ggsave("volcanoplot_DEG_female_vs_male_X_chr_gene_adj_disease_fdr_significant_version_2.pdf",
       combined_plot, width = 17, height = 15.8, dpi = 300)
# 5. Stacked bar: Female-up / Male-up / Non-DEG per cell type
n_gene_total <- as.data.frame(table(DEG_all$cell_type))
colnames(n_gene_total) <- c("cell_type", "n_total")
DEG_sig_dir <- DEG_all %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(
    type = dplyr::case_when(
      logFC > 0 ~ "Female-upregulated",
      logFC < 0 ~ "Male-upregulated",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(type)) %>%
  dplyr::group_by(cell_type, type) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop")
n_sig_per_ct <- DEG_sig_dir %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(n_sig = sum(count), .groups = "drop")
df_nonsig <- n_gene_total %>%
  dplyr::left_join(n_sig_per_ct, by = "cell_type") %>%
  dplyr::mutate(
    n_sig = ifelse(is.na(n_sig), 0, n_sig),
    count = n_total - n_sig,
    type = "Non-DEG"
  ) %>%
  dplyr::select(cell_type, type, count)
df_long2 <- dplyr::bind_rows(DEG_sig_dir, df_nonsig)
df_long2$cell_type <- factor(df_long2$cell_type, levels = cell_types)
df_long2$type <- factor(
  df_long2$type,
  levels = c("Non-DEG", "Female-upregulated", "Male-upregulated")
)
ggplot(df_long2, aes(x = cell_type, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +
  geom_text(
    aes(label = ifelse(count > 0, count, "")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Female-upregulated" = "red",
      "Male-upregulated" = "blue",
      "Non-DEG" = "grey80"
    ),
    breaks = c("Female-upregulated", "Male-upregulated", "Non-DEG")
  ) +
  labs(x = NULL, y = "Number of genes", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12)
  )
ggsave("barplot_DEG_across_model_split_direction.pdf", width = 6.2, height = 3.5, dpi = 300)
# 6. Female- vs male-biased DEG percentage per cell type; chi-square vs reference cell type
DEG_direction <- DEG_all %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::mutate(
    DEG_direction = dplyr::case_when(
      logFC > 0 ~ "Female-upregulated",
      logFC < 0 ~ "Male-upregulated",
      TRUE ~ "No direction"
    )
  ) %>%
  dplyr::filter(DEG_direction %in% c("Female-upregulated", "Male-upregulated"))
DEG_direction_summary <- DEG_direction %>%
  dplyr::group_by(cell_type, DEG_direction) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(cell_type) %>%
  dplyr::mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent, 1), "%\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()
lv <- cell_types
DEG_direction_summary$cell_type <- factor(DEG_direction_summary$cell_type, levels = lv)

ref_cell     <- "Mic"
use_adjusted <- T   # stars use FDR-adjusted p
get_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}
# chi-square on counts
DEG_wide <- DEG_direction_summary %>%
  dplyr::select(cell_type, DEG_direction, n) %>%
  tidyr::pivot_wider(names_from = DEG_direction, values_from = n, values_fill = 0) %>%
  as.data.frame()
for (col in c("Female-upregulated", "Male-upregulated")) {
  if (!col %in% colnames(DEG_wide)) DEG_wide[[col]] <- 0
}
miss_ct <- setdiff(lv, as.character(DEG_wide$cell_type))
if (length(miss_ct) > 0) {
  DEG_wide <- dplyr::bind_rows(
    DEG_wide,
    data.frame(cell_type = miss_ct,
               `Female-upregulated` = 0,
               `Male-upregulated`   = 0,
               check.names = FALSE)
  )
}
rownames(DEG_wide) <- as.character(DEG_wide$cell_type)
exc_F <- DEG_wide[ref_cell, "Female-upregulated"]
exc_M <- DEG_wide[ref_cell, "Male-upregulated"]
other_cells <- lv[lv != ref_cell]
chisq_results <- lapply(other_cells, function(ct) {
  ct_F <- DEG_wide[ct, "Female-upregulated"]
  ct_M <- DEG_wide[ct, "Male-upregulated"]
  m <- matrix(
    c(exc_F, exc_M,
      ct_F,  ct_M),
    nrow = 2, byrow = TRUE,
    dimnames = list(c(ref_cell, ct), c("Female", "Male"))
  )
  test <- chisq.test(m)
  data.frame(
    comparison     = paste0(ref_cell, "_vs_", ct),
    cell_type      = ct,
    exc_F          = exc_F,
    exc_M          = exc_M,
    exc_female_pct = exc_F / (exc_F + exc_M) * 100,
    ct_F           = ct_F,
    ct_M           = ct_M,
    ct_female_pct  = ct_F / (ct_F + ct_M) * 100,
    chisq          = unname(test$statistic),
    df             = unname(test$parameter),
    p_value        = test$p.value
  )
}) %>% bind_rows()
chisq_results$p_adj <- p.adjust(chisq_results$p_value, method = "fdr")
print(chisq_results)
writexl::write_xlsx(chisq_results, "chisq_Exc_vs_others_DEG_direction.xlsx")

# significance brackets
pos <- setNames(seq_along(lv), lv)
p_for_star <- if (use_adjusted) chisq_results$p_adj else chisq_results$p_value
star_labels <- vapply(p_for_star, get_stars, character(1))
max_h  <- 100
base_y <- max_h * 1.02
step_y <- max_h * 0.05
bracket_df <- data.frame(
  cell_type  = other_cells,
  xmin       = pos[ref_cell],
  xmax       = pos[other_cells],
  annotation = star_labels,
  stringsAsFactors = FALSE
)
bracket_df <- bracket_df[order(bracket_df$xmax), ]
bracket_df$y_position <- base_y + step_y * (seq_len(nrow(bracket_df)) - 1)

p <- ggplot(DEG_direction_summary, aes(x = cell_type, y = percent, fill = DEG_direction)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  geom_signif(
    annotations = bracket_df$annotation,
    xmin        = bracket_df$xmin,
    xmax        = bracket_df$xmax,
    y_position  = bracket_df$y_position,
    tip_length  = 0.005,
    size        = 0.3,
    textsize    = 2.5,
    vjust       = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Female-upregulated" = "red",
      "Male-upregulated" = "blue"
    )
  ) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(
    ylim = c(0, max(bracket_df$y_position) + step_y),
    clip = "off"
  ) +
  labs(x = NULL, y = "Percentage of significant DEGs", fill = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
  )
print(p)
ggsave(
  "barplot_DEG_female_male_percentage_by_cell_type.pdf",
  plot = p,
  width = 4.5,
  height = 5.4,
  dpi = 300
)
# 7. Cell-type-specific vs shared DEGs: A) sharing histogram, B) heatmap of shared genes, C) heatmap of specific genes
gene_col <- "updated_Gene_symbol"
stopifnot(gene_col %in% colnames(DEG_all))

DEG_sig <- DEG_all %>%
  dplyr::mutate(gene = .data[[gene_col]]) %>%
  dplyr::filter(adj.P.Val < 0.05, !is.na(logFC), logFC != 0) %>%
  dplyr::mutate(direction = ifelse(logFC > 0, "Female", "Male")) %>%
  dplyr::group_by(gene, cell_type) %>%
  dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

gene_class <- DEG_sig %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_celltype = dplyr::n_distinct(cell_type),
    celltypes  = paste(cell_type[order(match(cell_type, lv))], collapse = ", "),
    n_female   = sum(direction == "Female"),
    n_male     = sum(direction == "Male"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    category = dplyr::case_when(
      n_celltype == 1             ~ "Cell-type-specific",
      n_female == 0 | n_male == 0 ~ "Shared (concordant)",
      TRUE                        ~ "Shared (discordant)"
    )
  ) %>%
  dplyr::arrange(dplyr::desc(n_celltype), category, gene)
cat("\n===== gene classification summary =====\n"); print(table(gene_class$category))
writexl::write_xlsx(gene_class, "DEG_gene_classification_specific_vs_shared.xlsx")

prop_specific_shared <- gene_class %>%
  dplyr::mutate(group = ifelse(n_celltype == 1, "Cell-type-specific", "Shared")) %>%
  dplyr::count(group) %>%
  dplyr::mutate(percent = round(n / sum(n) * 100, 1))
prop_shared_dir <- gene_class %>%
  dplyr::filter(n_celltype >= 2) %>%
  dplyr::count(category) %>%
  dplyr::mutate(percent_within_shared = round(n / sum(n) * 100, 1))
cat("\n===== specific vs shared =====\n"); print(prop_specific_shared)
cat("\n===== shared: concordant vs discordant =====\n"); print(prop_shared_dir)

col_map <- c("Female"    = "red",
             "Male"      = "blue",
             "NS"        = "#9E9E9E",
             "NotTested" = "#ECECEC")
legend_param <- list(
  at     = c("Female", "Male", "NS", "NotTested"),
  labels = c("Female-upregulated DEGs", "Male-upregulated DEGs",
             "No differential expression", "Not tested")
)
legend_param_top <- c(legend_param, list(direction = "horizontal", nrow = 1))
# gene x cell type status matrix (Female / Male / NS / NotTested)
build_status_mat <- function(genes) {
  if (length(genes) == 0) return(NULL)
  long <- DEG_all %>%
    dplyr::mutate(gene = .data[[gene_col]]) %>%
    dplyr::filter(gene %in% genes) %>%
    dplyr::group_by(gene, cell_type) %>%
    dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(status = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Female",
      adj.P.Val < 0.05 & logFC < 0 ~ "Male",
      TRUE                         ~ "NS"
    )) %>%
    dplyr::select(gene, cell_type, status) %>%
    tidyr::complete(gene = genes, cell_type = lv) %>%
    dplyr::mutate(status = ifelse(is.na(status), "NotTested", status),
                  cell_type = factor(cell_type, levels = lv))
  m <- long %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = status) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  m[, lv, drop = FALSE]
}

# A) sharing histogram, shared gene names on bar tops
col_purple <- "#6A3D9A"   # concordant
col_gold   <- "#C8920A"   # discordant
col_grey   <- "grey75"    # specific
gene_fontsize <- 9
y_axis_max    <- 60
share_summary <- gene_class %>%
  dplyr::mutate(concord = dplyr::case_when(
    n_celltype == 1             ~ "Specific",
    n_female == 0 | n_male == 0 ~ "Concordant",
    TRUE                        ~ "Discordant"
  )) %>%
  dplyr::count(n_celltype, concord)
share_summary$concord <- factor(share_summary$concord,
                                levels = c("Specific", "Concordant", "Discordant"))
bar_top <- share_summary %>%
  dplyr::group_by(n_celltype) %>%
  dplyr::summarise(bar_h = sum(n), .groups = "drop")
off <- max(bar_top$bar_h) * 0.02
label_rich <- gene_class %>%
  dplyr::filter(n_celltype >= 2) %>%
  dplyr::mutate(
    concord = ifelse(grepl("discordant", category), "Discordant", "Concordant"),
    col     = ifelse(concord == "Discordant", col_gold, col_purple),
    tag     = paste0("<span style='color:", col, "'><i>", gene, "</i></span>")
  ) %>%
  dplyr::arrange(n_celltype, concord, gene) %>%
  dplyr::group_by(n_celltype) %>%
  dplyr::summarise(label = paste(tag, collapse = "<br>"), .groups = "drop") %>%
  dplyr::left_join(bar_top, by = "n_celltype") %>%
  dplyr::mutate(y_anchor = bar_h + off)
y_max <- y_axis_max
p_A <- ggplot(share_summary, aes(x = factor(n_celltype), y = n, fill = concord)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  ggtext::geom_richtext(
    data = label_rich,
    aes(x = factor(n_celltype), y = y_anchor, label = label),
    inherit.aes = FALSE, vjust = 0, hjust = 0.5,
    size = gene_fontsize / ggplot2::.pt, lineheight = 1.0,
    label.color = NA, fill = NA
  ) +
  scale_fill_manual(
    values = c("Specific"   = col_grey,
               "Concordant" = col_purple,
               "Discordant" = col_gold),
    breaks = c("Specific", "Concordant", "Discordant"),
    labels = c("Specific", "Shared (same direction)", "Shared (opposite direction)")
  ) +
  scale_y_continuous(breaks = seq(0, y_axis_max, 20)) +
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  labs(x = "Number of cell types in which a gene is DEG",
       y = "Number of genes", fill = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        plot.margin = margin(t = 15, r = 12, b = 10, l = 12))
n_label_max <- gene_class %>% dplyr::filter(n_celltype >= 2) %>%
  dplyr::count(n_celltype) %>% dplyr::pull(n)
n_label_max <- if (length(n_label_max) == 0) 1 else max(n_label_max)
fig_h_A <- max(5, n_label_max * 0.16 + 3.5)

# B) heatmap of shared genes (>= 2 cell types)
shared_info  <- gene_class %>% dplyr::filter(n_celltype >= 2)
shared_genes <- shared_info$gene
if (length(shared_genes) >= 1) {
  mat_shared <- build_status_mat(shared_genes)
  ord <- order(-shared_info$n_celltype[match(rownames(mat_shared), shared_info$gene)])
  mat_shared   <- mat_shared[ord, , drop = FALSE]
  mat_shared_t <- t(mat_shared)
  
  nB    <- ncol(mat_shared_t)
  showB <- nB <= 120
  wB <- if (showB) max(6, nB * 0.16 + 2.5) else min(16, max(8, nB * 0.035 + 3))
  hB <- if (showB) 4.2 else 3.2
  
  htB <- ComplexHeatmap::Heatmap(
    mat_shared_t, name = "Status", col = col_map,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_names = showB,
    column_names_gp  = grid::gpar(fontsize = 8, fontface = "italic"),
    column_names_rot = 45,
    row_names_side   = "left", row_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "white", lwd = 1),
    heatmap_legend_param = legend_param_top
  )
}

# C) heatmap of cell-type-specific genes
specific_info <- gene_class %>%
  dplyr::filter(n_celltype == 1) %>%
  dplyr::mutate(sig_ct  = celltypes,
                sig_dir = ifelse(n_female > 0, "Female", "Male"))
specific_genes <- specific_info$gene
if (length(specific_genes) >= 1) {
  mat_spec <- build_status_mat(specific_genes)
  idx       <- match(rownames(mat_spec), specific_info$gene)
  sig_ct_v  <- specific_info$sig_ct[idx]
  sig_dir_v <- specific_info$sig_dir[idx]
  ordc      <- order(match(sig_ct_v, lv), sig_dir_v, rownames(mat_spec))
  mat_spec  <- mat_spec[ordc, , drop = FALSE]
  
  idx2        <- match(rownames(mat_spec), specific_info$gene)
  col_split_c <- droplevels(factor(specific_info$sig_ct[idx2], levels = lv))
  mat_spec_t  <- t(mat_spec)
  
  nC    <- ncol(mat_spec_t)
  showC <- nC <= 120
  wC <- if (showC) max(6, nC * 0.16 + 2.5) else min(18, max(8, nC * 0.03 + 3))
  hC <- if (showC) 4.5 else 3.5
  
  htC <- ComplexHeatmap::Heatmap(
    mat_spec_t, name = "Status", col = col_map,
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_split    = col_split_c,
    column_title_gp = grid::gpar(fontsize = 11),
    column_gap      = grid::unit(1.5, "mm"),
    show_column_names = showC,
    column_names_gp  = grid::gpar(fontsize = 7, fontface = "italic"),
    column_names_rot = 45,
    row_names_side   = "left", row_names_gp = grid::gpar(fontsize = 10),
    rect_gp = grid::gpar(col = "white", lwd = 0.5),
    show_heatmap_legend = FALSE
  )
}

# combined: A | B over C
if (exists("htB") && exists("htC")) {
  wA <- 6
  gB <- grid::grid.grabExpr(ComplexHeatmap::draw(htB, heatmap_legend_side = "top"))
  gC <- grid::grid.grabExpr(ComplexHeatmap::draw(htC))
  
  row1 <- cowplot::plot_grid(p_A, gB, ncol = 2, rel_widths = c(wA, wB),
                             labels = c("A", "B"), label_size = 16)
  combo <- cowplot::plot_grid(row1, gC, ncol = 1, rel_heights = c(fig_h_A, hC),
                              labels = c("", "C"), label_size = 16)
  
  ggsave("DEG_gene_sharing_combined.pdf", combo,
         width = wA + wB, height = fig_h_A + hC, dpi = 300, limitsize = FALSE)
}
# 8. logFC correlation: case/control-adjusted vs non-adjusted
DEG_all_non_adj <- DEG_all
DEG_all_adj     <- run_deg(adjust_case_control = TRUE)

identical(DEG_all_adj$gene_cell, DEG_all_non_adj$gene_cell)

x <- DEG_all_adj$logFC
y <- DEG_all_non_adj$logFC
# Pearson if normal, else Spearman
norm_ok <- if (length(x) <= 5000) {
  shapiro.test(x)$p.value > 0.05 && shapiro.test(y)$p.value > 0.05
} else {
  # Shapiro not available for n > 5000; use skewness
  sk <- function(v) mean((v - mean(v))^3) / sd(v)^3
  abs(sk(x)) < 1 && abs(sk(y)) < 1
}
method <- if (norm_ok) "pearson" else "spearman"
ct_cor <- cor.test(x, y, method = method, exact = FALSE)
r_val <- as.numeric(ct_cor$estimate)
label_method <- if (method == "pearson") "Pearson r" else "Spearman rho"
cat(sprintf("Normality passed: %s -> using %s\n", norm_ok, method))
cat(sprintf("%s = %.3f, p = %.3e\n", label_method, r_val, ct_cor$p.value))

lim <- range(c(x, y))
df_cor <- data.frame(x = x, y = y)
p <- ggplot(df_cor, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.5, size = 1.6, color = "#2c7fb8") +
  geom_smooth(method = "lm", se = TRUE, color = "#d95f02", fill = "#d95f02", alpha = 0.15) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 3.8,
           label = sprintf("%s = %.6f\nn = %d", label_method, r_val, length(x))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "logFC (adjusted case/control)",
       y = "logFC (non-adjusted case/control)",
       title = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p)
ggsave("DEG_logFC_correlation_adj_vs_nonadj.pdf", p, width = 5.5, height = 5.5)