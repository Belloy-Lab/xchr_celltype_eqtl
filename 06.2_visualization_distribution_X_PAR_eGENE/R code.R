# eGene distribution across cell types: UpSet plot, pie chart, eGene proportion bar plot (chi-square vs Exc), eGenes per cell type
rm(list = ls())
library(data.table)
library(tidyverse)
library(ComplexUpset)
library(ggsci)
library(ggsignif)
# CONFIG
indep_x_rds   <- "../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds"
indep_par_rds <- "../05_independent_signal_PAR/df_eqtl_wide_after_second_clump.rds"
cov_pheno_dir <- "../01.2_cov_pheno/"
n_cis_x_xlsx   <- "../04_significant_eQTL_X/N_of_cis_eQTL_X.xlsx"
n_cis_par_xlsx <- "../04_significant_eQTL_PAR/N_of_cis_eQTL_PAR.xlsx"
# 1. Independent eGene-cell type pairs; expression pattern = number of cell types the gene was tested in
df_X_final_independent_variants <- readRDS(indep_x_rds)
df_PAR_final_independent_variants <- readRDS(indep_par_rds)
df_X_final_independent_variants <- select(df_X_final_independent_variants,cell_type,gene,par_type)
df_PAR_final_independent_variants <- select(df_PAR_final_independent_variants,cell_type,gene,par_type)
df_X_PAR_final_independent_variants <- rbind(df_X_final_independent_variants,df_PAR_final_independent_variants)
df_X_PAR_final_independent_variants$gene_cell <- paste0(df_X_PAR_final_independent_variants$gene,"_",
                                                        df_X_PAR_final_independent_variants$cell_type)
df_X_PAR_final_independent_variants <-  distinct(df_X_PAR_final_independent_variants,gene_cell,.keep_all = T)

df_pheno_ast <- fread(paste0(cov_pheno_dir, "both_Ast.pheno_chrX.tsv"))
input_gene_ast <- colnames(df_pheno_ast)[-c(1,2)]
df_pheno_end <- fread(paste0(cov_pheno_dir, "both_End.pheno_chrX.tsv"))
input_gene_end <- colnames(df_pheno_end)[-c(1,2)]
df_pheno_exc <- fread(paste0(cov_pheno_dir, "both_Exc.pheno_chrX.tsv"))
input_gene_exc <- colnames(df_pheno_exc)[-c(1,2)]
df_pheno_inh <- fread(paste0(cov_pheno_dir, "both_Inh.pheno_chrX.tsv"))
input_gene_inh <- colnames(df_pheno_inh)[-c(1,2)]
df_pheno_mic <- fread(paste0(cov_pheno_dir, "both_Mic.pheno_chrX.tsv"))
input_gene_mic <- colnames(df_pheno_mic)[-c(1,2)]
df_pheno_oli <- fread(paste0(cov_pheno_dir, "both_Oli.pheno_chrX.tsv"))
input_gene_oli <- colnames(df_pheno_oli)[-c(1,2)]
df_pheno_opc <- fread(paste0(cov_pheno_dir, "both_OPC.pheno_chrX.tsv"))
input_gene_opc <- colnames(df_pheno_opc)[-c(1,2)]
all_input_gene <- c(input_gene_ast,
                    input_gene_end,
                    input_gene_exc,
                    input_gene_inh,
                    input_gene_mic,
                    input_gene_oli,
                    input_gene_opc)
all_input_gene <- as.data.frame(table(all_input_gene))
all_input_gene$gene_expression_pattern <- ifelse(all_input_gene$Freq == 1,"Cell type-specific","Multiple cell types")

rm(df_X_final_independent_variants,
   df_PAR_final_independent_variants,
   input_gene_ast,
   input_gene_end,
   input_gene_exc,
   input_gene_inh,
   input_gene_mic,
   input_gene_oli,
   input_gene_opc,
   df_pheno_ast,
   df_pheno_end,
   df_pheno_exc,
   df_pheno_inh,
   df_pheno_mic,
   df_pheno_oli,
   df_pheno_opc)
colnames(all_input_gene)[1] <- "gene"
all_input_gene$Freq <- NULL
df_X_PAR_final_independent_variants <- left_join(df_X_PAR_final_independent_variants,all_input_gene)
writexl::write_xlsx(df_X_PAR_final_independent_variants,"Upsetplot_input_df.xlsx")
# 2. UpSet plot
df_gene <- df_X_PAR_final_independent_variants %>%
  distinct(cell_type, gene, gene_expression_pattern)
df_binary <- df_gene %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = cell_type,
    values_from = value,
    values_fill = 0
  )
df_binary <- df_binary %>%
  mutate(
    n_celltype = rowSums(across(where(is.numeric)))
  )
df_binary_filtered <- df_binary
celltype_cols <- setdiff(
  colnames(df_binary_filtered),
  c("gene", "gene_expression_pattern", "n_celltype")
)

df_plot <- df_binary_filtered %>%
  mutate(
    intersection_name = apply(
      select(., all_of(celltype_cols)),
      1,
      function(x) {
        paste(celltype_cols[which(x == 1)], collapse = "&")
      }
    ),
    degree = rowSums(select(., all_of(celltype_cols)))
  )

# order intersections by degree, then size
intersection_size_df <- df_plot %>%
  group_by(intersection_name) %>%
  summarise(
    size = n(),
    degree = first(degree),
    .groups = "drop"
  ) %>%
  arrange(degree, desc(size))
intersection_list <- lapply(
  intersection_size_df$intersection_name,
  function(x) strsplit(x, "&")[[1]]
)

p <- upset(
  df_plot,
  intersect = celltype_cols,
  intersections = intersection_list,
  sort_intersections = FALSE,
  width_ratio = 0.2,
  base_annotations = list(
    "# of eGenes \n in intersections" =
      intersection_size(
        mapping = aes(fill = gene_expression_pattern),
        text = list(vjust = -0.5, colour = "black")
      ) +
      scale_fill_manual(
        name   = "Expression pattern",
        values = c(
          "Cell type-specific"  = "grey70",
          "Multiple cell types" = "#2C6BA0"
        ),
        breaks = c("Cell type-specific", "Multiple cell types"),
        labels = c(
          "Cell type-specific"  = "Gene available in 1 cell-type only",
          "Multiple cell types" = "Gene available in ≥2 cell-types"
        )
      ) +
      theme(
        legend.position = c(0.78, 0.85),
        legend.background = element_blank(),
        legend.key = element_blank()
      )
  ),
  set_sizes = upset_set_size()
) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p
ggsave("UpSetplot.pdf", p, width = 10.5, height = 6, dpi = 300)
# 3. Pie chart: genes by number of cell types
df_count <- df_gene %>%
  distinct(cell_type, gene) %>%
  group_by(gene) %>%
  summarise(n_celltype = n(), .groups = "drop") %>%
  count(n_celltype) %>%
  arrange(n_celltype)
df_count$n_celltype <- factor(df_count$n_celltype,
                              levels = df_count$n_celltype)
df_count$n_celltype <- factor(
  df_count$n_celltype,
  levels = rev(levels(df_count$n_celltype))
)
p <- ggplot(df_count, aes(x = "", y = n, fill = n_celltype)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  geom_text(aes(label = ""),
            position = position_stack(vjust = 0.5),
            size = 6,
            color = "white") +
  scale_fill_brewer(palette = "Blues") +
  theme(legend.position = "none")
print(p)
df_count
ggsave("pie_chart.pdf",width = 6,height = 6,dpi = 300)
# 4. eGene proportion per cell type; chi-square Exc vs each other cell type
df_X_PAR_final_independent_variants <- readxl::read_excel("Upsetplot_input_df.xlsx")
n_independent_variant_count_per_cell_type <- as.data.frame(table(df_X_PAR_final_independent_variants$cell_type))
colnames(n_independent_variant_count_per_cell_type) <- c("cell_type","n_sig_independent")
df_1 <- readxl::read_excel(n_cis_x_xlsx)
df_2 <- readxl::read_excel(n_cis_par_xlsx)
df_1 <- distinct(df_1,cell_type,n_gene_included,.keep_all = T)
df_2 <- distinct(df_2,cell_type,n_gene_included,.keep_all = T)
if (identical(df_1$cell_type,df_2$cell_type)) {
  df_12 <- left_join(df_1,df_2,by="cell_type")
}
df_12$n_input_gene_X_PAR <- df_12$n_gene_included.x + df_12$n_gene_included.y
df_12 <- select(df_12,cell_type,n_input_gene_X_PAR)
df_12 <- left_join(df_12,n_independent_variant_count_per_cell_type,by="cell_type")
df_12$n_non_sig <- df_12$n_input_gene_X_PAR-df_12$n_sig_independent

ref_cell     <- "Exc"
use_adjusted <- T   # stars use FDR-adjusted p
get_stars <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}
exc_row    <- df_12[df_12$cell_type == ref_cell, ]
exc_sig    <- exc_row$n_sig_independent
exc_nonsig <- exc_row$n_non_sig
lv <- c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "End")
other_cells <- lv[lv != ref_cell]
chisq_results <- lapply(other_cells, function(ct) {
  row <- df_12[df_12$cell_type == ct, ]
  m <- matrix(
    c(exc_sig,               exc_nonsig,
      row$n_sig_independent, row$n_non_sig),
    nrow = 2, byrow = TRUE,
    dimnames = list(c(ref_cell, ct), c("eGene", "non_eGene"))
  )
  test <- chisq.test(m)
  data.frame(
    comparison = paste0(ref_cell, "_vs_", ct),
    cell_type  = ct,
    exc_n_sig    = exc_sig,
    exc_n_total  = exc_row$n_input_gene_X_PAR,
    ct_n_sig     = row$n_sig_independent,
    ct_n_total   = row$n_input_gene_X_PAR,
    exc_prop     = exc_sig / exc_row$n_input_gene_X_PAR,
    ct_prop      = row$n_sig_independent / row$n_input_gene_X_PAR,
    chisq        = unname(test$statistic),
    df           = unname(test$parameter),
    p_value      = test$p.value
  )
}) %>% bind_rows()
chisq_results$p_adj <- p.adjust(chisq_results$p_value, method = "fdr")
print(chisq_results)
writexl::write_xlsx(chisq_results, "chisq_Exc_vs_others_eGene_proportion.xlsx")

# significance brackets
pos <- setNames(seq_along(lv), lv)
p_for_star <- if (use_adjusted) chisq_results$p_adj else chisq_results$p_value
star_labels <- vapply(p_for_star, get_stars, character(1))
max_h  <- max(df_12$n_input_gene_X_PAR)
base_y <- max_h * 1.02
step_y <- max_h * 0.04
bracket_df <- data.frame(
  cell_type  = other_cells,
  xmin       = pos[ref_cell],
  xmax       = pos[other_cells],
  annotation = star_labels,
  stringsAsFactors = FALSE
)
bracket_df <- bracket_df[order(bracket_df$xmax), ]
bracket_df$y_position <- base_y + step_y * (seq_len(nrow(bracket_df)) - 1)

df_long <- df_12 %>%
  pivot_longer(
    cols = c(n_sig_independent, n_non_sig),
    names_to = "type",
    values_to = "count"
  )
df_long$cell_type <- factor(df_long$cell_type, levels = lv)
p <- ggplot(df_long, aes(x = cell_type, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = count),
    position = position_stack(vjust = 0.5),
    size = 4,
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
      "n_sig_independent" = "#E64B35",
      "n_non_sig" = "grey80"
    ),
    labels = c(
      "n_sig_independent" = "eGene",
      "n_non_sig" = "Non-eGene"
    )
  ) +
  scale_y_continuous(breaks = c(0,100, 200, 300, 400, 500)) +
  coord_cartesian(clip = "off") +
  labs(
    x = NULL,
    y = "Number of genes",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    plot.margin = margin(t = 60, r = 10, b = 10, l = 10)
  )
print(p)
ggsave("barplot_eGene_X_PAR_independent_across_model.pdf", plot = p, width = 6.5, height = 5.5, dpi = 600)
# 5. eGenes per cell type, grouped by number of independent eQTLs (starts fresh; literal paths)
rm(list=ls())
df_X_final_independent_variants <- readRDS("../05_independent_signal_X/df_eqtl_wide_after_second_clump.rds")
df_PAR_final_independent_variants <- readRDS("../05_independent_signal_PAR/df_eqtl_wide_after_second_clump.rds")
df_X_final_independent_variants$abs_BETA <- abs(df_X_final_independent_variants$BETA_both_model_2)
df_X_final_independent_variants$MAF <- pmin(df_X_final_independent_variants$EAF_both,
                                            1 - df_X_final_independent_variants$EAF_both)
df_PAR_final_independent_variants$abs_BETA <- abs(df_PAR_final_independent_variants$BETA_both_par)
df_PAR_final_independent_variants$MAF <- pmin(df_PAR_final_independent_variants$EAF_both,
                                              1 - df_PAR_final_independent_variants$EAF_both)
df_X_final_independent_variants <- select(df_X_final_independent_variants,
                                          gene,
                                          cell_type)
df_PAR_final_independent_variants <- select(df_PAR_final_independent_variants,
                                            gene,
                                            cell_type)
df_X_PAR_final_independent_variants <- rbind(df_X_final_independent_variants,df_PAR_final_independent_variants)
df_count <- df_X_PAR_final_independent_variants %>%
  count(cell_type, gene)
df_count <- df_count %>%
  mutate(freq_group = as.character(n))
df_plot <- df_count %>%
  count(cell_type, freq_group)
df_plot$cell_type <- factor(
  df_plot$cell_type,
  levels = c("Exc", "Inh", "Ast", "Oli", "OPC", "Mic", "End")
)
pd <- position_dodge2(width = 0.65, preserve = "single")
ggplot(
  df_plot,
  aes(
    x = cell_type,
    y = n,
    fill = freq_group
  )
) +
  geom_col(
    position = pd,
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  geom_text(
    aes(label = n),
    position = pd,
    vjust = -0.35,
    size = 3
  ) +
  scale_fill_npg(name = "# of independent eQTL") +
  guides(fill = guide_legend(ncol = 3, byrow = FALSE)) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.2
    )
  ) +
  labs(
    x = "Cell type",
    y = "Number of eGenes"
  ) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)))
ggsave("N of eGenes.pdf",width = 6.5,height = 4,dpi = 300)