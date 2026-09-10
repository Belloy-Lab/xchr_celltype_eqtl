# LocusZoom (locuszoomr) female vs male eQTL panels for example sex-heterogeneous genes in Exc
rm(list = ls())
library(tidyverse)
library(locuszoomr)
library(data.table)
library(EnsDb.Hsapiens.v86)
# CONFIG
pl <- "plink1.9"
bfile_dir  <- "../00.1_X_chr_bfile_QC/"
bfile_name <- "chrX_QC_EU.0.G"
afreq_dir <- "../00.1_X_chr_bfile_QC/"
eqtl_rds <- "../03_merge_eQTL_res_X/df_eqtl_cis_all.rds"
rsid_rds <- "../00.3_rsid_annotation/rsid_gene_annotation.rds"
sex_het_xlsx <- "../05_independent_signal_X/non_PAR_variants_in_barplot_ordered_by_ratio_per_cell_type.xlsx"
target_genes  <- c("INTS6L", "GYG2", "MID2")
display_cell_type <- "Exc"
flank_bp      <- 150000   # half-window around the index variant
ens_db        <- "EnsDb.Hsapiens.v86"
highlight_col <- "red"
# 1. Annotate cis-eQTL results with rsID and allele frequency
df_eqtl_all <- readRDS(eqtl_rds)
df_rsid <- readRDS(rsid_rds)
df_eqtl_all <- left_join(df_eqtl_all, dplyr::select(df_rsid, ID, rsid), by = 'ID')
afreq_both <- fread(paste0(afreq_dir, bfile_name, ".afreq"))
afreq_both$OBS_CT <- NULL
colnames(afreq_both)[5] <- "EAF_both"
afreq_female <- fread(paste0(afreq_dir, bfile_name, ".female.afreq"))
afreq_female$OBS_CT <- NULL
colnames(afreq_female)[5] <- "EAF_female"
afreq_male <- fread(paste0(afreq_dir, bfile_name, ".male.afreq"))
afreq_male$OBS_CT <- NULL
colnames(afreq_male)[5] <- "EAF_male"
df_eqtl_all <- left_join(df_eqtl_all, afreq_both,   by = c("#CHROM","ID", "REF", "ALT"))
df_eqtl_all <- left_join(df_eqtl_all, afreq_female, by = c("#CHROM","ID", "REF", "ALT"))
df_eqtl_all <- left_join(df_eqtl_all, afreq_male,   by = c("#CHROM","ID", "REF", "ALT"))
# 2. Index variants of the target genes
df_sex_het <- readxl::read_excel(sex_het_xlsx)
df_sex_het <- dplyr::arrange(df_sex_het, min_pvalue)
df_sex_het <- dplyr::filter(df_sex_het,
                            gene %in% target_genes,
                            cell_type == display_cell_type)
df_sex_het <- dplyr::distinct(df_sex_het, cell_type, gene, .keep_all = T)
analysis_list <- distinct(df_sex_het, gene, cell_type, ID, rsid, POS,
                          EAF_both, EAF_female, EAF_male)
writexl::write_xlsx(analysis_list, "index_variants.xlsx")
# 3. Per index variant: LD, locus objects, female/male panels
dir.create("temp", showWarnings = FALSE)
for (i in 1:nrow(analysis_list)) {
  eqtl_female <- dplyr::filter(df_eqtl_all, gene == analysis_list$gene[i] & cell_type == analysis_list$cell_type[i],
                               setting == "female_xchr_model_2")
  eqtl_male <- dplyr::filter(df_eqtl_all, gene == analysis_list$gene[i] & cell_type == analysis_list$cell_type[i],
                             setting == "male_xchr_model_1")
  eqtl_female <- dplyr::semi_join(eqtl_female, eqtl_male, by = 'ID')
  eqtl_male <- dplyr::semi_join(eqtl_male, eqtl_female, by = 'ID')
  
  eqtl_female <- dplyr::select(eqtl_female,
                               chrom = `#CHROM`,
                               pos = POS,
                               rsid = ID,
                               other_allele = REF,
                               effect_allele = ALT,
                               p = pvalue,
                               beta = BETA,
                               se = SE)
  eqtl_male <- dplyr::select(eqtl_male,
                             chrom = `#CHROM`,
                             pos = POS,
                             rsid = ID,
                             other_allele = REF,
                             effect_allele = ALT,
                             p = pvalue,
                             beta = BETA,
                             se = SE)
  
  eqtl_female <- arrange(eqtl_female, p)
  eqtl_male   <- arrange(eqtl_male, p)
  
  top_snp <- analysis_list$ID[i]
  
  # r2 to the index variant
  write.table(
    unique(c(eqtl_female$rsid, eqtl_male$rsid)),
    paste0("temp/snp_list", i, ".txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  command <- paste(
    pl,
    "--bfile", paste0(bfile_dir, bfile_name),
    "--allow-no-sex",
    "--r2",
    "--ld-snp", top_snp,
    "--ld-window 99999",
    "--ld-window-kb 99999",
    "--ld-window-r2 0",
    "--extract", paste0(getwd(), "/temp/snp_list", i, ".txt"),
    "--keep-allele-order",
    "--out", paste0(getwd(), "/temp/ld_res", i),
    sep = " "
  )
  
  system(command)
  
  ld_df <- fread(paste0("temp/ld_res", i, ".ld"))
  ld_df <- dplyr::select(ld_df, rsid = SNP_B, R2)
  
  eqtl_female <- left_join(eqtl_female, ld_df)
  eqtl_male   <- left_join(eqtl_male, ld_df)
  
  xrange_min <- analysis_list$POS[i] - flank_bp
  xrange_max <- analysis_list$POS[i] + flank_bp
  
  loc_female <- locus(
    eqtl_female,
    gene = NULL,
    seqname = eqtl_female$chrom[1],
    xrange = c(xrange_min, xrange_max),
    index_snp = top_snp,
    LD = "R2",
    ens_db = ens_db
  )
  
  loc_male <- locus(
    eqtl_male,
    gene = NULL,
    seqname = eqtl_female$chrom[1],
    xrange = c(xrange_min, xrange_max),
    index_snp = top_snp,
    LD = "R2",
    ens_db = ens_db
  )
  
  loc_female <- link_recomb(loc_female, genome = "hg38")
  loc_male   <- link_recomb(loc_male, genome = "hg38")
  
  loc_female$TX <- dplyr::filter(loc_female$TX, gene_biotype == "protein_coding")
  loc_male$TX   <- dplyr::filter(loc_male$TX,   gene_biotype == "protein_coding")
  
  pdf(
    paste0(
      "LCzoom_",
      analysis_list$gene[i],
      "_",
      analysis_list$cell_type[i],
      "_female_male.pdf"
    ),
    width = 3,
    height = 4.875
  )
  
  layout(matrix(1:3, ncol = 1), heights = c(1.8, 1.8, 0.6))
  
  ylim_eqtl <- c(
    0,
    max(
      -log10(min(eqtl_female$p)),
      -log10(min(eqtl_male$p))
    ) + 0.1
  )
  
  scatter_plot(
    loc_female,
    ylim = ylim_eqtl,
    labels = NULL,
    legend_pos = "topleft",
    xaxt = "n",
    xlab = ""
  )
  title(paste0("Female eQTL ", analysis_list$gene[i], " (", analysis_list$cell_type[i], ")",
               ", index variant: ", analysis_list$rsid[i],
               ", EAF female: ", analysis_list$EAF_female[i]), line = 2.7)
  
  scatter_plot(
    loc_male,
    ylim = ylim_eqtl,
    labels = NULL,
    legend_pos = "topleft",
    xaxt = "n",
    xlab = ""
  )
  title(paste0("Male eQTL ", analysis_list$gene[i], " (", analysis_list$cell_type[i], ")",
               ", index variant: ", analysis_list$rsid[i],
               ", EAF male: ", analysis_list$EAF_male[i]), line = 2.7)
  
  # gene track, target gene highlighted
  genetracks(
    loc_female,
    highlight     = analysis_list$gene[i],
    highlight_col = highlight_col
  )
  
  dev.off()
}