#!/usr/bin/env Rscript
# Real-Data Bayesian Colocalization (coloc.abf) Runner for Scheme 2
# Reads real fine-mapped GRCh37 GWAS locus files from ../outputs/loci_extracted/
# and harmonizes against user-supplied brain/blood eQTL summary statistics placed in
# ../inputs/external_eqtl/. Never fabricates synthetic eQTL statistics.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[1]) else "scripts/03_run_coloc_abf_when_eqtl_available.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
base_dir <- dirname(script_dir)
in_dir <- file.path(base_dir, "inputs")
eqtl_dir <- file.path(in_dir, "external_eqtl")
loci_dir <- file.path(base_dir, "outputs", "loci_extracted")
out_tab <- file.path(base_dir, "outputs", "tables", "coloc_abf_real_eqtl_results.tsv")
check_only <- ("--check-only" %in% args)

eqtl_files <- list.files(eqtl_dir, pattern = "\\.(tsv|tsv\\.gz|txt|txt\\.gz)$", full.names = TRUE)
loci_files <- list.files(loci_dir, pattern = "_gwas_grch37\\.tsv$", full.names = TRUE)
has_coloc <- requireNamespace("coloc", quietly = TRUE)

if (check_only) {
  cat(sprintf(
    "PREFLIGHT_OK: coloc_pkg=%s, loci_extracted=%d, external_eqtl=%d\n",
    ifelse(has_coloc, "TRUE", "FALSE"), length(loci_files), length(eqtl_files)
  ))
  quit(save = "no", status = 0)
}

if (!has_coloc) {
  stop(paste(
    "R package 'coloc' is not installed in this R environment.",
    "Install via: Rscript -e \"install.packages('coloc', repos='https://cloud.r-project.org')\"",
    sep = "\n"
  ))
}
library(coloc)

if (length(eqtl_files) == 0) {
  stop(paste(
    "No external eQTL summary files found in:", eqtl_dir,
    "See ../inputs/external_eqtl/README_DOWNLOAD_EQTL.txt for required columns",
    "(rs, chr, pos, ea, nea, beta, se, maf, n, gene) matching GRCh37 coordinates.",
    sep = "\n"
  ))
}

if (length(loci_files) == 0) {
  stop("Run 01_extract_gwas_loci_and_finemap.py first to populate outputs/loci_extracted/.")
}

cat("Found", length(loci_files), "GWAS loci and", length(eqtl_files), "eQTL datasets.\n")
N_GWAS <- 500000L
S_PROP <- 0.34

is_palindromic <- function(a1, a2) {
  p <- paste0(toupper(a1), toupper(a2))
  p %in% c("AT", "TA", "CG", "GC")
}

results <- list()
for (eq_f in eqtl_files) {
  eq <- fread(eq_f)
  req_cols <- c("rs", "beta", "se", "maf", "n", "gene")
  if (!all(req_cols %in% names(eq))) {
    warning("Skipping ", basename(eq_f), " (missing required columns: ", paste(setdiff(req_cols, names(eq)), collapse = ", "), ")")
    next
  }
  for (loc_f in loci_files) {
    gw <- fread(loc_f)
    m <- merge(gw, eq, by = "rs", suffixes = c("_gwas", "_eqtl"), allow.cartesian = TRUE)
    m <- m[
      !is.na(beta_gwas) & !is.na(beta_eqtl) &
      se_gwas > 0 & se_eqtl > 0 &
      af > 0 & af < 1 &
      maf > 0 & maf < 1
    ]
    if (nrow(m) < 20) next

    # Harmonize alleles when supplied; exclude ambiguous palindromic SNPs near MAF 0.5
    if (all(c("allele1", "allele0", "ea", "nea") %in% names(m))) {
      pal <- is_palindromic(m$allele1, m$allele0) & (abs(m$af - 0.5) < 0.08 | abs(m$maf - 0.5) < 0.08)
      m <- m[!pal]
      flip <- (toupper(m$allele1) == toupper(m$nea)) & (toupper(m$allele0) == toupper(m$ea))
      same <- (toupper(m$allele1) == toupper(m$ea)) & (toupper(m$allele0) == toupper(m$nea))
      m <- m[same | flip]
      m[flip, beta_eqtl := -beta_eqtl]
    }
    if (nrow(m) < 20) next

    for (gname in unique(m$gene)) {
      sub <- m[gene == gname]
      sub <- sub[!duplicated(rs)]
      if (nrow(sub) < 20) next
      maf_gwas <- pmin(sub$af, 1 - sub$af)
      maf_eqtl <- pmin(sub$maf, 1 - sub$maf)
      d1 <- list(beta = sub$beta_gwas, varbeta = sub$se_gwas^2, snp = sub$rs, MAF = maf_gwas, N = N_GWAS, s = S_PROP, type = "cc")
      d2 <- list(beta = sub$beta_eqtl, varbeta = sub$se_eqtl^2, snp = sub$rs, MAF = maf_eqtl, N = median(sub$n, na.rm = TRUE), type = "quant")
      res <- coloc.abf(dataset1 = d1, dataset2 = d2)
      s <- res$summary
      results[[length(results) + 1]] <- data.table(
        locus_file = basename(loc_f),
        eqtl_file = basename(eq_f),
        gene = gname,
        n_shared_snps = as.integer(s["nsnps"]),
        PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"],
        PP.H2 = s["PP.H2.abf"], PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"]
      )
    }
  }
}

if (length(results) > 0) {
  out_dt <- rbindlist(results)[order(-PP.H4)]
  fwrite(out_dt, out_tab, sep = "\t")
  cat("Saved real coloc.abf results to:", out_tab, "\n")
} else {
  cat("No overlapping loci/genes met >=20 shared SNP threshold.\n")
}
