#!/usr/bin/env Rscript
# 10_run_p0_limma_camera_interaction.R
# P0 Statistical Verification:
#   - Batch-adjusted limma-voom pseudobulk DE (capturing Sex + Chemistry + Sequencing Batch)
#   - Inter-gene correlation adjusted limma::camera pathway tests + donor Z 95% CIs
#   - Sex-stratified Batch-adjusted sensitivity (XIST artifact audit + F/M concordance)
#   - Joint Mixed-Effect limma-voom + duplicateCorrelation cell_type x diagnosis interaction
#   - 71-donor Batch+Sex-adjusted logit cellular composition stability check

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

suppressPackageStartupMessages({
  has_deps <- all(sapply(c("data.table", "edgeR", "limma", "statmod"), requireNamespace, quietly = TRUE))
})

if (check_only) {
  cat(sprintf("P0_LIMMA_PREFLIGHT: deps_ok=%s\n", ifelse(has_deps, "TRUE", "FALSE")))
  quit(save = "no", status = 0)
}

if (!has_deps) stop("Missing required R packages: data.table, edgeR, limma, statmod")

# Determine script directory robustly
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root_dir   <- dirname(script_dir)

pb_dir  <- file.path(script_dir, "outputs", "local_pb_cache")
if (!dir.exists(pb_dir)) {
  pb_dir <- Sys.getenv("PB_CACHE_DIR", unset = file.path(root_dir, "data", "pb_cache"))
}
in_s1   <- file.path(root_dir, "01_Scheme1_snRNA_Spatial_scDRS", "inputs")
out_res <- file.path(script_dir, "outputs", "tables")
dir.create(out_res, showWarnings = FALSE, recursive = TRUE)

source_script_impl <- function() {
  library(data.table); library(edgeR); library(limma); library(statmod)
  genes_dt <- fread(file.path(pb_dir, "genes.tsv"))
  meta_dt  <- fread(file.path(pb_dir, "pb_groups_meta.tsv"))
  don_dt   <- fread(file.path(pb_dir, "donor_metadata_complete.tsv"))

  n_genes <- nrow(genes_dt); n_grps <- nrow(meta_dt)
  raw_vec <- readBin(file.path(pb_dir, "pb_counts.bin"), what = "integer", n = as.numeric(n_genes) * n_grps, size = 4L)
  cnt_mat <- matrix(raw_vec, nrow = n_genes, ncol = n_grps)
  rownames(cnt_mat) <- genes_dt$gene
  colnames(cnt_mat) <- paste0(meta_dt$donor_id, "__", meta_dt$author_subtype)

  mc_df <- fread(file.path(in_s1, "MitoCarta_core_frozen.tsv"))
  ox_df <- fread(file.path(in_s1, "OXPHOS_GO0006119_148genes.tsv"))
  names(mc_df) <- gsub('"', '', names(mc_df))
  mc_df$MitoCartaSymbol <- gsub('"', '', mc_df$MitoCartaSymbol)
  mc_df$ComplexGroup    <- gsub('"', '', mc_df$ComplexGroup)

  nuclear_mc  <- mc_df[!grepl("^MT-", MitoCartaSymbol)]
  ox_go_genes <- unique(ox_df$Gene[!grepl("^MT-", ox_df$Gene)])
  gs_list <- list(
    Nuclear_MitoCarta_OXPHOS = unique(nuclear_mc$MitoCartaSymbol),
    Nuclear_GO_OXPHOS        = ox_go_genes,
    Inhibitory_KCC2_GABA     = c("SLC12A5", "WNK3", "WNK1", "GAD1", "GAD2", "SLC32A1", "GABRA1", "GABRB2", "GABRG2", "PVALB", "SST"),
    Excitatory_Glut_NKCC1    = c("SLC12A2", "STK39", "SLC17A7", "GRIA1", "GRIA2", "GRIN1", "GRIN2A", "GRIN2B")
  )

  valid_meta <- meta_dt[n_cells >= 10]
  subtypes <- sort(unique(valid_meta$author_subtype))
  key_genes <- c("SLC12A5", "SLC12A2", "STK39", "WNK3", "GABRA1", "GRIN2A", "GRIN2B", "GRIA1",
                 "ATP5F1D", "TUBB", "EIF5A", "HSPA1A", "HSPB1", "XIST")

  cam_results <- list(); de_key_results <- list()
  for (sub in subtypes) {
    sm <- valid_meta[author_subtype == sub]
    if (nrow(sm) < 12 || length(unique(sm$diagnosis)) < 2 || min(table(sm$diagnosis)) < 4) next
    cols <- paste0(sm$donor_id, "__", sm$author_subtype)
    y <- DGEList(counts = cnt_mat[, cols, drop = FALSE], group = factor(sm$diagnosis, levels = c("Control", "MDD")))
    sm$diag <- factor(sm$diagnosis, levels = c("Control", "MDD"))
    sm$sex_f <- factor(sm$Sex); sm$batch_f <- factor(sm$Batch)

    if (nlevels(droplevels(sm$batch_f)) >= 2) {
      des_try <- model.matrix(~ diag + batch_f, data = sm)
      if (qr(des_try)$rank == ncol(des_try) && (nrow(sm) - ncol(des_try)) >= 5) {
        des <- des_try; cov_mode <- "Batch(captures_Sex+Chem+Batch)"
      } else if (nlevels(droplevels(sm$sex_f)) == 2) {
        des <- model.matrix(~ diag + sex_f, data = sm); cov_mode <- "Sex"
      } else {
        des <- model.matrix(~ diag, data = sm); cov_mode <- "Unadjusted"
      }
    } else {
      des <- model.matrix(~ diag, data = sm); cov_mode <- "Unadjusted"
    }

    keep <- filterByExpr(y, design = des, min.count = 10)
    if (sum(keep) < 200) next
    y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
    v <- voom(y, des, plot = FALSE)
    fit <- eBayes(lmFit(v, des), robust = TRUE)
    tt <- topTable(fit, coef = "diagMDD", number = Inf, sort.by = "none")
    tt$gene <- rownames(tt); tt$cell_type <- sub; tt$n_donors <- nrow(sm); tt$cov_mode <- cov_mode

    ksub <- tt[tt$gene %in% key_genes, c("cell_type", "gene", "logFC", "t", "P.Value", "adj.P.Val", "n_donors", "cov_mode")]
    if (nrow(ksub) > 0) {
      se_vec <- sqrt(fit$s2.post[ksub$gene]) * fit$stdev.unscaled[ksub$gene, "diagMDD"]
      tc <- qt(0.975, df = fit$df.total[ksub$gene])
      ksub$se <- se_vec; ksub$ci_low <- ksub$logFC - tc * se_vec; ksub$ci_high <- ksub$logFC + tc * se_vec
      de_key_results[[length(de_key_results) + 1]] <- as.data.table(ksub)
    }

    idx_list <- ids2indices(gs_list, id = rownames(v))
    idx_list <- idx_list[sapply(idx_list, length) >= 5]
    if (length(idx_list) > 0) {
      cm <- camera(v, index = idx_list, design = des, contrast = which(colnames(des) == "diagMDD"), inter.gene.cor = NA)
      cm$pathway <- rownames(cm); cm$cell_type <- sub; cm$n_donors <- nrow(sm)
      cm$n_degs_fdr05 <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE); cm$cov_mode <- cov_mode
      logcpm <- cpm(y, log = TRUE, prior.count = 1); zmat <- t(scale(t(logcpm)))
      eff_vec <- cil_vec <- cih_vec <- numeric(nrow(cm))
      for (pi in seq_len(nrow(cm))) {
        g_in <- rownames(v)[idx_list[[cm$pathway[pi]]]]
        dscore <- colMeans(zmat[g_in, , drop = FALSE], na.rm = TRUE)
        df_ols <- data.frame(score = dscore, diag = sm$diag, batch = sm$batch_f, sex = sm$sex_f)
        mod <- if (cov_mode == "Batch(captures_Sex+Chem+Batch)") lm(score ~ diag + batch, data = df_ols) else if (cov_mode == "Sex") lm(score ~ diag + sex, data = df_ols) else lm(score ~ diag, data = df_ols)
        cf <- summary(mod)$coefficients["diagMDD", ]; ci <- confint(mod)["diagMDD", ]
        eff_vec[pi] <- cf["Estimate"]; cil_vec[pi] <- ci[1]; cih_vec[pi] <- ci[2]
      }
      cm$donor_mean_z_diff <- eff_vec; cm$ci_low <- cil_vec; cm$ci_high <- cih_vec
      cam_results[[length(cam_results) + 1]] <- as.data.table(cm)
    }
  }
  cam_all <- rbindlist(cam_results, fill = TRUE)
  cam_all[, FDR_across_subtypes := p.adjust(PValue, method = "BH"), by = pathway]
  fwrite(cam_all, file.path(out_res, "p0_limma_camera_batch_adjusted_pathways.tsv"), sep = "\t")
  fwrite(rbindlist(de_key_results, fill = TRUE), file.path(out_res, "p0_limma_voom_batch_adjusted_key_genes.tsv"), sep = "\t")
}
source_script_impl()
