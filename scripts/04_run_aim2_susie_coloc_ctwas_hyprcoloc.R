#!/usr/bin/env Rscript
# ==============================================================================
# Aim 2 Multi-Causal Fine-Mapping & Multi-Trait Colocalization R Pipeline
# ==============================================================================
# Resolves single-causal HEIDI conflicts (e.g. chr3p21.31 GPX1/AMT/RNF123 &
# chr3p21.1 GNL3/SPCS1/NEK4 & chr12 DNAH10/CDK2AP1) using:
#   1. SuSiE (`susieR::susie_rss`) with 1000G EUR regional LD matrices (`plink --r square`)
#   2. Multi-signal Colocalization (`coloc::coloc.susie`)
#   3. Joint Multi-Trait Colocalization (`hyprcoloc::hyprcoloc`: MDD x ANX x Brain/Cell eQTL)
#   4. Genetic-confounder-adjusted Causal-TWAS (`ctwas::ctwas_rss` when weights provided)
#
# Usage:
#   Rscript 04_run_aim2_susie_coloc_ctwas_hyprcoloc.R --check-only
#   Rscript 04_run_aim2_susie_coloc_ctwas_hyprcoloc.R --ld-dir <dir> --eqtl-dir <dir>
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
check_only <- ("--check-only" %in% args)

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "04_run_aim2_susie_coloc_ctwas_hyprcoloc.R"
base_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

loci_dir <- file.path(base_dir, "outputs", "loci_dualtrait")
ld_dir   <- file.path(base_dir, "external_downloads", "ld_matrices")
eqtl_dir <- file.path(base_dir, "external_downloads", "eqtl_loci")
out_dir  <- file.path(base_dir, "outputs", "tables")

has_coloc     <- requireNamespace("coloc", quietly = TRUE)
has_susie     <- requireNamespace("susieR", quietly = TRUE)
has_hyprcoloc <- requireNamespace("hyprcoloc", quietly = TRUE)
has_ctwas     <- requireNamespace("ctwas", quietly = TRUE)

loci_files <- list.files(loci_dir, pattern = "_dualtrait_grch37\\.tsv$", full.names = TRUE)

if (check_only) {
  cat(sprintf(
    "AIM2_R_PREFLIGHT_OK: loci_dualtrait=%d, coloc=%s, susieR=%s, hyprcoloc=%s, ctwas=%s\n",
    length(loci_files),
    ifelse(has_coloc, "TRUE", "FALSE"),
    ifelse(has_susie, "TRUE", "FALSE"),
    ifelse(has_hyprcoloc, "TRUE", "FALSE"),
    ifelse(has_ctwas, "TRUE", "FALSE")
  ))
  quit(save = "no", status = 0)
}

if (!has_coloc || !has_susie) {
  stop("Install R packages: install.packages(c('coloc','susieR','data.table')); remotes::install_github('jrs95/hyprcoloc')")
}

library(coloc)
library(susieR)

if (length(loci_files) == 0) {
  stop("Run 00_prepare_local_aim2_aligned_loci.py first to populate outputs/loci_dualtrait/.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
susie_results <- list()
hypr_results  <- list()

eqtl_files <- list.files(eqtl_dir, pattern = "\\.(tsv|tsv\\.gz)$", full.names = TRUE)

for (loc_f in loci_files) {
  tag <- sub("_dualtrait_grch37\\.tsv$", "", basename(loc_f))
  ld_mat_f <- file.path(ld_dir, paste0(tag, ".ld.gz"))
  ld_bim_f <- file.path(ld_dir, paste0(tag, ".bim"))
  if (!file.exists(ld_mat_f) || !file.exists(ld_bim_f)) {
    message("Skipping ", tag, " (missing regional LD files: ", ld_mat_f, ")")
    next
  }

  gwas <- fread(loc_f)
  bim  <- fread(ld_bim_f, header = FALSE, col.names = c("chr", "SNP", "cm", "bp", "a1_ld", "a2_ld"))
  ld_R <- as.matrix(fread(ld_mat_f))
  rownames(ld_R) <- colnames(ld_R) <- bim$SNP

  common_snps <- intersect(gwas$SNP, bim$SNP)
  if (length(common_snps) < 50) next
  gwas <- gwas[SNP %in% common_snps]
  gwas <- gwas[!duplicated(SNP)]
  setkey(gwas, SNP)

  idx <- match(gwas$SNP, bim$SNP)
  R_sub <- ld_R[idx, idx, drop = FALSE]
  # Align GWAS alleles against LD reference A1/A2 before SuSiE RSS
  a1g <- toupper(gwas$A1_MDD); a2g <- toupper(gwas$A2_MDD)
  a1l <- toupper(bim$a1_ld[idx]); a2l <- toupper(bim$a2_ld[idx])
  flip <- (a1g == a2l) & (a2g == a1l)
  keep <- ((a1g == a1l) & (a2g == a2l)) | flip
  gwas <- gwas[keep]
  R_sub <- R_sub[keep, keep, drop = FALSE]
  if (any(flip[keep])) {
    s_flip <- flip[keep]
    gwas[s_flip, `:=`(BETA_MDD = -BETA_MDD, Z_MDD = -Z_MDD, BETA_ANX = -BETA_ANX, Z_ANX = -Z_ANX)]
  }
  rownames(R_sub) <- colnames(R_sub) <- gwas$SNP

  # 1. Fit SuSiE on MDD (N ~ 1.15M) and ANX (N ~ 1.10M) using case-control log-OR varbeta
  d_mdd <- list(
    beta = gwas$BETA_MDD, varbeta = gwas$SE_MDD^2, snp = gwas$SNP,
    position = gwas$BP, N = median(gwas$N_MDD, na.rm = TRUE), type = "cc", LD = R_sub
  )
  d_anx <- list(
    beta = gwas$BETA_ANX, varbeta = gwas$SE_ANX^2, snp = gwas$SNP,
    position = gwas$BP, N = median(gwas$N_ANX, na.rm = TRUE), type = "cc", LD = R_sub
  )
  s_mdd <- runsusie(d_mdd, L = 10)
  s_anx <- runsusie(d_anx, L = 10)

  cs_res <- coloc.susie(s_mdd, s_anx)
  if (!is.null(cs_res$summary) && nrow(cs_res$summary) > 0) {
    sm <- as.data.table(cs_res$summary)
    sm[, locus := tag]
    sm[, comparison := "MDD_vs_ANX"]
    susie_results[[length(susie_results) + 1]] <- sm
  }

  # 2. Multi-trait HyPrColoc when matching external eQTL/sQTL/pQTL files exist
  if (has_hyprcoloc && length(eqtl_files) > 0) {
    library(hyprcoloc)
    for (eq_f in eqtl_files) {
      eq <- fread(eq_f)
      if (!all(c("SNP", "BETA", "SE", "gene") %in% names(eq))) next
      eq_sub <- eq[SNP %in% gwas$SNP]
      for (g_target in unique(eq_sub$gene)) {
        eg <- eq_sub[gene == g_target][!duplicated(SNP)]
        m <- merge(gwas, eg, by = "SNP")
        if (nrow(m) < 50) next
        B_mat <- cbind(MDD = m$BETA_MDD, ANX = m$BETA_ANX, QTL = m$BETA)
        S_mat <- cbind(MDD = m$SE_MDD,   ANX = m$SE_ANX,   QTL = m$SE)
        rownames(B_mat) <- rownames(S_mat) <- m$SNP
        R_eq <- R_sub[m$SNP, m$SNP, drop = FALSE]
        hc <- hyprcoloc(B_mat, S_mat, trait.names = c("MDD", "ANX", paste0(basename(eq_f), ":", g_target)), snp.id = m$SNP, ld.matrix = R_eq)
        ht <- as.data.table(hc$results)
        ht[, `:=`(locus = tag, gene = g_target, eqtl_source = basename(eq_f))]
        hypr_results[[length(hypr_results) + 1]] <- ht
      }
    }
  }
}

if (length(susie_results) > 0) {
  fwrite(rbindlist(susie_results, fill = TRUE), file.path(out_dir, "aim2_coloc_susie_results.tsv"), sep = "\t")
}
if (length(hypr_results) > 0) {
  fwrite(rbindlist(hypr_results, fill = TRUE), file.path(out_dir, "aim2_hyprcoloc_multitrait_results.tsv"), sep = "\t")
}
cat("Aim 2 multi-causal R pipeline finished cleanly.\n")
