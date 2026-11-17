#!/usr/bin/env Rscript
# 12a_run_p0_coloc_abf_sensitivity.R
# Evaluates coloc.abf prior sensitivity (p12 = 1e-6, 1e-5, 1e-4) across all 7 dual-trait loci
# and renders native coloc::sensitivity() diagnostic plot for SLC12A5.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

has_pkgs <- all(sapply(c("data.table", "coloc"), requireNamespace, quietly = TRUE))
if (check_only) {
  cat(sprintf("P0_COLOC_PREFLIGHT: coloc=%s\n", ifelse(has_pkgs, "TRUE", "FALSE")))
  quit(save = "no", status = 0)
}
if (!has_pkgs) stop("Requires data.table and coloc.")

suppressPackageStartupMessages({ library(data.table); library(coloc) })
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
loci_dir <- file.path(script_dir, "outputs", "loci_dualtrait")
out_tab  <- file.path(script_dir, "outputs", "tables", "p0_aim2_all7loci_coloc_abf_p12_sensitivity.tsv")
fig_dir  <- file.path(script_dir, "outputs", "figures")

loci_files <- sort(list.files(loci_dir, pattern = "_dualtrait_grch37[.]tsv$", full.names = TRUE))
rows <- list()
for (lf in loci_files) {
  tag <- sub("_dualtrait_grch37[.]tsv$", "", basename(lf))
  gw <- fread(lf)[!is.na(BETA_MDD) & SE_MDD > 0 & !is.na(BETA_ANX) & SE_ANX > 0]
  d_mdd <- list(beta = gw$BETA_MDD, varbeta = gw$SE_MDD^2, snp = gw$SNP, position = gw$BP, N = median(gw$N_MDD, na.rm = TRUE), s = 0.34, type = "cc")
  d_anx <- list(beta = gw$BETA_ANX, varbeta = gw$SE_ANX^2, snp = gw$SNP, position = gw$BP, N = median(gw$N_ANX, na.rm = TRUE), s = 0.25, type = "cc")
  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    suppressWarnings({ res <- coloc.abf(dataset1 = d_mdd, dataset2 = d_anx, p12 = p12) })
    sm <- as.list(res$summary)
    rows[[length(rows) + 1]] <- data.table(
      locus = tag, nsnps = sm$nsnps, p12_prior = p12,
      PP.H0 = sm$PP.H0.abf, PP.H1 = sm$PP.H1.abf, PP.H2 = sm$PP.H2.abf,
      PP.H3 = sm$PP.H3.abf, PP.H4 = sm$PP.H4.abf,
      mdd_lead = gw$SNP[which.min(gw$P_MDD)], mdd_min_p = min(gw$P_MDD),
      anx_lead = gw$SNP[which.min(gw$P_ANX)], anx_min_p = min(gw$P_ANX)
    )
    # Render native coloc::sensitivity plot for SLC12A5 at standard prior p12 = 1e-5
    if (tag == "SLC12A5_chr20" && abs(p12 - 1e-5) < 1e-12) {
      png(file.path(fig_dir, "Figure_Aim2B_Coloc_Sensitivity_SLC12A5.png"), width = 2600, height = 1400, res = 300)
      suppressWarnings(sensitivity(res, rule = "H4 > 0.8", plot.manhattans = FALSE))
      dev.off()
      pdf(file.path(fig_dir, "Figure_Aim2B_Coloc_Sensitivity_SLC12A5.pdf"), width = 8.6, height = 4.6)
      suppressWarnings(sensitivity(res, rule = "H4 > 0.8", plot.manhattans = FALSE))
      dev.off()
    }
  }
}
fwrite(rbindlist(rows), out_tab, sep = "\t")
cat("[SUCCESS] Saved coloc p12 sensitivity table & native SLC12A5 sensitivity plot.\n")
