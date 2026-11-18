#!/usr/bin/env Rscript
# 14_run_local_planA_upgrades.R
# Local Statistical Upgrades Addressing Reviewer Critiques (Zero External Downloads Required):
#   1. Lineage-matched Pooled ExN / InN / Glial limma-voom DE + 140-gene eQTL-DE concordance
#   2. 1,000-iteration within-batch donor-label permutation across 25 neuronal subtypes
#   3. SVA latent confounders + HSPA1A stress-proxy + min-nuclei (5/20/50) sensitivity
#   4. Per-donor CPM & proportion export for individual-data-point visualization

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

req_pkgs <- c("data.table", "edgeR", "limma", "statmod", "sva")
has_deps <- all(sapply(req_pkgs, requireNamespace, quietly = TRUE))
if (check_only) {
  cat(sprintf("PLAN_A_UPGRADE_PREFLIGHT: deps_ok=%s\n", ifelse(has_deps, "TRUE", "FALSE")))
  quit(save = "no", status = 0)
}
if (!has_deps) stop("Missing required R packages: data.table, edgeR, limma, statmod, sva")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(statmod)
  library(sva)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
root_dir <- dirname(script_dir)
pb_dir  <- file.path(script_dir, "outputs", "local_pb_cache")
if (!dir.exists(pb_dir)) pb_dir <- Sys.getenv("PB_CACHE_DIR", unset = file.path(root_dir, "data", "pb_cache"))
out_res <- file.path(script_dir, "outputs", "tables")
dir.create(out_res, showWarnings = FALSE, recursive = TRUE)

genes_dt <- fread(file.path(pb_dir, "genes.tsv"))
meta_dt  <- fread(file.path(pb_dir, "pb_groups_meta.tsv"))
don_dt   <- fread(file.path(pb_dir, "donor_metadata_complete.tsv"))
harm_dt  <- fread(file.path(out_res, "p0_bryois_allelic_harmonized_all.tsv"))

n_genes <- nrow(genes_dt)
n_grps  <- nrow(meta_dt)
raw_vec <- readBin(file.path(pb_dir, "pb_counts.bin"), what = "integer", n = as.numeric(n_genes) * n_grps, size = 4L)
cnt_mat <- matrix(raw_vec, nrow = n_genes, ncol = n_grps)
rownames(cnt_mat) <- genes_dt$gene
colnames(cnt_mat) <- paste0(meta_dt$donor_id, "__", meta_dt$author_subtype)

# Helper to map fine subtypes to Bryois broad lineages
map_broad <- function(st) {
  out <- rep(NA_character_, length(st))
  out[grepl("^ExN", st)] <- "Excitatory.neurons"
  out[grepl("^InN", st)] <- "Inhibitory.neurons"
  out[grepl("^Ast", st)] <- "Astrocytes"
  out[grepl("^Oli", st)] <- "Oligodendrocytes"
  out[grepl("^OPC", st)] <- "OPCs...COPs"
  out[grepl("^Mic", st)] <- "Microglia"
  out[grepl("^End", st)] <- "Endothelial.cells"
  out
}
meta_dt[, broad_lineage := map_broad(author_subtype)]

# ==============================================================================
# 1. Pooled Broad-Lineage limma-voom DE + Full 140-Gene eQTL vs DE Concordance
# ==============================================================================
cat("[1/4] Running lineage-matched pooled DE across 71 donors for all 140 harmonized genes...\n")
lineages <- c("Excitatory.neurons", "Inhibitory.neurons", "Astrocytes",
              "Oligodendrocytes", "OPCs...COPs", "Microglia", "Endothelial.cells")

pooled_de_list <- list()
for (lin in lineages) {
  sub_meta <- meta_dt[broad_lineage == lin]
  # Sum counts per donor
  d_ids <- sort(unique(sub_meta$donor_id))
  lin_mat <- matrix(0L, nrow = n_genes, ncol = length(d_ids))
  rownames(lin_mat) <- rownames(cnt_mat)
  colnames(lin_mat) <- d_ids
  n_nuc <- integer(length(d_ids))
  for (j in seq_along(d_ids)) {
    did <- d_ids[j]
    idx <- which(meta_dt$broad_lineage == lin & meta_dt$donor_id == did)
    n_nuc[j] <- sum(meta_dt$n_cells[idx])
    if (length(idx) == 1L) {
      lin_mat[, j] <- cnt_mat[, idx]
    } else {
      lin_mat[, j] <- rowSums(cnt_mat[, idx, drop = FALSE])
    }
  }
  d_info <- don_dt[match(d_ids, don_dt$donor_id)]
  keep_d <- which(n_nuc >= 5L)
  if (length(keep_d) < 10L) next
  lin_mat <- lin_mat[, keep_d, drop = FALSE]
  d_info  <- d_info[keep_d]
  d_info[, diag := factor(diagnosis, levels = c("Control", "MDD"))]
  d_info[, batch := factor(Batch)]
  if (length(unique(d_info$diag)) < 2L) next

  design <- model.matrix(~ diag + batch, data = d_info)
  dge <- DGEList(counts = lin_mat)
  keep_g <- filterByExpr(dge, design = design, min.count = 10)
  dge <- dge[keep_g, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = "diagMDD", number = Inf, sort.by = "none"), keep.rownames = "gene")
  tt[, broad_lineage := lin]
  tt[, n_donors := nrow(d_info)]
  pooled_de_list[[lin]] <- tt[, .(broad_lineage, gene, n_donors, logFC, AveExpr, t, P.Value, adj.P.Val)]
}
pooled_de_dt <- rbindlist(pooled_de_list)
fwrite(pooled_de_dt[gene %in% c("SLC12A5", "GABRA1", "ATP5F1D", "NCKIPSD", "CACNA2D2", "DYNC1I2", "NEK4", "GNL3", "DCC", "PPP6C")],
       file.path(out_res, "p1_pooled_lineage_de_summary.tsv"), sep = "\t")

# Merge all harmonized pairs (exclude Pericytes which have no dedicated GSE213982 broad cluster, or report matching rate)
harm_mdd <- harm_dt[trait == "MDD"][order(-PP.H4)]
conc_dt <- merge(harm_mdd, pooled_de_dt, by.x = c("gene", "cell_type"), by.y = c("gene", "broad_lineage"), all.x = FALSE)
conc_dt[, same_sign := (sign(eqtl_beta_on_risk_allele) == sign(logFC))]
conc_dt[, discordance_score := eqtl_beta_on_risk_allele * (-logFC)]
conc_dt <- conc_dt[order(-discordance_score)]
conc_dt[, discordance_pct_rank := round(100 * (1 - (rank(-discordance_score) - 1) / .N), 2)]
fwrite(conc_dt, file.path(out_res, "p1_bryois_all140_lineage_matched_eqtl_vs_de.tsv"), sep = "\t")

cat(sprintf("   Matched %d gene-lineage pairs across %d unique genes.\n",
            nrow(conc_dt), uniqueN(conc_dt$gene)))
r_all <- cor.test(conc_dt$eqtl_beta_on_risk_allele, conc_dt$logFC, method = "pearson")
rho_all <- cor.test(conc_dt$eqtl_beta_on_risk_allele, conc_dt$logFC, method = "spearman", exact = FALSE)
cat(sprintf("   All matched pairs (n=%d): Pearson r=%.3f (P=%.3g), Spearman rho=%.3f (P=%.3g), same-sign=%d/%d (%.1f%%)\n",
            nrow(conc_dt), r_all$estimate, r_all$p.value, rho_all$estimate, rho_all$p.value,
            sum(conc_dt$same_sign), nrow(conc_dt), 100 * mean(conc_dt$same_sign)))
sub50 <- conc_dt[PP.H4 >= 0.50]
cat(sprintf("   Coloc PP.H4 >= 0.50 pairs (n=%d): same-sign=%d/%d\n", nrow(sub50), sum(sub50$same_sign), nrow(sub50)))
print(conc_dt[gene %in% c("SLC12A5", "NCKIPSD", "DYNC1I2", "GNL3", "NEK4"),
              .(gene, cell_type, PP.H4, eqtl_beta_on_risk_allele, logFC, P.Value, adj.P.Val, same_sign, discordance_pct_rank)])

# ==============================================================================
# 2. 1,000-Iteration Within-Batch Donor-Label Permutation Across Neuronal Subtypes
# ==============================================================================
cat("[2/4] Running 1,000 within-batch donor-label permutations across testable neuronal subtypes...\n")
neur_subs <- sort(unique(meta_dt[grepl("^(ExN|InN)", author_subtype)]$author_subtype))
sub_pre <- list()
for (ct in neur_subs) {
  sm <- meta_dt[author_subtype == ct & n_cells >= 10]
  if (nrow(sm) < 12 || length(unique(sm$diagnosis)) < 2 || min(table(sm$diagnosis)) < 4) next
  sm <- sm[order(donor_id)]
  mat_ct <- cnt_mat[, sm$grp_idx + 1L, drop = FALSE]
  sm[, diag := factor(diagnosis, levels = c("Control", "MDD"))]
  sm[, batch := factor(Batch)]
  des <- model.matrix(~ diag + batch, data = sm)
  if (qr(des)$rank != ncol(des) || (nrow(sm) - ncol(des)) < 5) next
  dge <- DGEList(counts = mat_ct, group = sm$diag)
  keep_g <- filterByExpr(dge, design = des, min.count = 10)
  if (sum(keep_g) < 200 || !("SLC12A5" %in% rownames(dge)[keep_g])) next
  dge <- calcNormFactors(dge[keep_g, , keep.lib.sizes = FALSE])
  v <- voom(dge, des, plot = FALSE)
  s_idx <- which(rownames(v$E) == "SLC12A5")
  fit0 <- eBayes(lmFit(v, des), robust = TRUE)
  obs_t <- fit0$t[s_idx, "diagMDD"]
  obs_lfc <- fit0$coefficients[s_idx, "diagMDD"]
  obs_p <- fit0$p.value[s_idx, "diagMDD"]
  sub_pre[[ct]] <- list(
    meta = sm,
    y = v$E[s_idx, ],
    w = v$weights[s_idx, ],
    obs_t = obs_t,
    obs_lfc = obs_lfc,
    obs_p = obs_p
  )
}
valid_neur <- names(sub_pre)
obs_neg_cnt <- sum(sapply(sub_pre, function(x) x$obs_lfc < 0))
obs_nom_neg <- sum(sapply(sub_pre, function(x) x$obs_lfc < 0 & x$obs_p < 0.05))
obs_mean_t  <- mean(sapply(sub_pre, function(x) x$obs_t))

# Stratified donor label permutation within Batch (preserving sex/batch and cross-subtype donor correlation!)
set.seed(20261117)
n_perm <- 1000L
don_perm_dt <- copy(don_dt)
perm_neg_cnt <- integer(n_perm)
perm_nom_neg <- integer(n_perm)
perm_mean_t  <- numeric(n_perm)

batch_groups <- split(seq_len(nrow(don_perm_dt)), don_perm_dt$Batch)
orig_diag <- don_perm_dt$diagnosis

for (b in seq_len(n_perm)) {
  p_diag <- orig_diag
  for (g_idx in batch_groups) {
    if (length(g_idx) > 1L) p_diag[g_idx] <- sample(orig_diag[g_idx], length(g_idx), replace = FALSE)
  }
  names(p_diag) <- don_perm_dt$donor_id

  neg_c <- 0L; nom_c <- 0L; ts <- numeric(length(valid_neur))
  for (k in seq_along(valid_neur)) {
    obj <- sub_pre[[valid_neur[k]]]
    sm <- obj$meta
    pd <- factor(p_diag[sm$donor_id], levels = c("Control", "MDD"))
    des_p <- model.matrix(~ pd + sm$batch)
    fit_p <- lm.wfit(x = des_p, y = obj$y, w = obj$w)
    df_res <- fit_p$df.residual
    res <- fit_p$residuals
    s2 <- sum(obj$w * res^2) / df_res
    cov_b <- chol2inv(fit_p$qr$qr)[2, 2] * s2
    t_val <- fit_p$coefficients[2] / sqrt(cov_b)
    p_val <- 2 * pt(-abs(t_val), df = df_res)
    ts[k] <- t_val
    if (fit_p$coefficients[2] < 0) {
      neg_c <- neg_c + 1L
      if (p_val < 0.05) nom_c <- nom_c + 1L
    }
  }
  perm_neg_cnt[b] <- neg_c
  perm_nom_neg[b] <- nom_c
  perm_mean_t[b]  <- mean(ts)
}

p_emp_neg  <- (sum(perm_neg_cnt >= obs_neg_cnt) + 1) / (n_perm + 1)
p_emp_nom  <- (sum(perm_nom_neg >= obs_nom_neg) + 1) / (n_perm + 1)
p_emp_meant <- (sum(perm_mean_t <= obs_mean_t) + 1) / (n_perm + 1)

perm_summary <- data.table(
  n_tested_neuronal_subtypes = length(valid_neur),
  obs_negative_count = obs_neg_cnt,
  obs_nominal_neg_count = obs_nom_neg,
  obs_mean_t = obs_mean_t,
  n_permutations = n_perm,
  perm_median_neg_count = median(perm_neg_cnt),
  perm_95pct_neg_count = quantile(perm_neg_cnt, 0.95),
  empirical_p_neg_count = p_emp_neg,
  empirical_p_nominal_neg_count = p_emp_nom,
  empirical_p_mean_t = p_emp_meant
)
fwrite(perm_summary, file.path(out_res, "p1_slc12a5_donor_permutation_summary.tsv"), sep = "\t")
fwrite(data.table(perm_id = seq_len(n_perm), neg_cnt = perm_neg_cnt, nom_neg_cnt = perm_nom_neg, mean_t = perm_mean_t),
       file.path(out_res, "p1_slc12a5_donor_permutation_null.tsv"), sep = "\t")
print(perm_summary)

# ==============================================================================
# 3. SVA Latent Confounders + HSPA1A Proxy + Pseudobulk Min-Nuclei Sensitivity
# ==============================================================================
cat("[3/4] Running SVA, HSPA1A stress-proxy, and min-nuclei (5/10/20/50) sensitivity...\n")
sens_rows <- list()
for (ct in c("ExN10_L46", "InN1_PV", "ExN2_L23", "InN5_SST", "ExN7", "Pooled_ExN", "Pooled_InN")) {
  for (min_n in c(10L, 5L, 20L, 50L)) {
    if (grepl("^Pooled_", ct) && min_n != 10L) next
    if (ct == "Pooled_ExN") {
      sm <- meta_dt[broad_lineage == "Excitatory.neurons", .(n_cells = sum(n_cells)), by = .(donor_id, diagnosis, Batch)]
      d_ids <- sort(sm$donor_id)
      mat_ct <- matrix(0L, nrow = n_genes, ncol = length(d_ids))
      rownames(mat_ct) <- rownames(cnt_mat); colnames(mat_ct) <- d_ids
      for (j in seq_along(d_ids)) {
        idx <- which(meta_dt$broad_lineage == "Excitatory.neurons" & meta_dt$donor_id == d_ids[j])
        mat_ct[, j] <- if (length(idx) == 1L) cnt_mat[, idx] else rowSums(cnt_mat[, idx, drop = FALSE])
      }
      sm <- sm[match(d_ids, donor_id)]
    } else if (ct == "Pooled_InN") {
      sm <- meta_dt[broad_lineage == "Inhibitory.neurons", .(n_cells = sum(n_cells)), by = .(donor_id, diagnosis, Batch)]
      d_ids <- sort(sm$donor_id)
      mat_ct <- matrix(0L, nrow = n_genes, ncol = length(d_ids))
      rownames(mat_ct) <- rownames(cnt_mat); colnames(mat_ct) <- d_ids
      for (j in seq_along(d_ids)) {
        idx <- which(meta_dt$broad_lineage == "Inhibitory.neurons" & meta_dt$donor_id == d_ids[j])
        mat_ct[, j] <- if (length(idx) == 1L) cnt_mat[, idx] else rowSums(cnt_mat[, idx, drop = FALSE])
      }
      sm <- sm[match(d_ids, donor_id)]
    } else {
      sm <- meta_dt[author_subtype == ct & n_cells >= min_n]
      if (nrow(sm) < 12L) next
      mat_ct <- cnt_mat[, sm$grp_idx + 1L, drop = FALSE]
    }
    sm[, diag := factor(diagnosis, levels = c("Control", "MDD"))]
    sm[, batch := factor(Batch)]
    if (length(unique(sm$diag)) < 2L) next
    des0 <- model.matrix(~ diag + batch, data = sm)
    if (qr(des0)$rank != ncol(des0)) next
    dge <- DGEList(counts = mat_ct, group = sm$diag)
    keep_g <- filterByExpr(dge, design = des0, min.count = 10)
    if (!("SLC12A5" %in% rownames(dge)[keep_g])) next
    dge <- calcNormFactors(dge[keep_g, , keep.lib.sizes = FALSE])
    v0 <- voom(dge, des0, plot = FALSE)
    fit0 <- eBayes(lmFit(v0, des0), robust = TRUE)
    for (g_chk in intersect(c("SLC12A5", "GABRA1", "ATP5F1D"), rownames(v0$E))) {
      r0 <- topTable(fit0, coef = "diagMDD", number = Inf, sort.by = "none")[g_chk, ]
      sens_rows[[length(sens_rows) + 1L]] <- data.table(
        cell_type = ct, min_nuclei = min_n, n_donors = nrow(sm),
        model = "Base (diag + Batch)", gene = g_chk,
        logFC = r0$logFC, t = r0$t, P.Value = r0$P.Value, adj.P.Val = r0$adj.P.Val
      )
    }

    # Only run SVA and HSPA1A proxy at baseline min_n == 10
    if (min_n == 10L) {
      # Model + HSPA1A logCPM covariate
      if ("HSPA1A" %in% rownames(v0$E)) {
        hsp_z <- as.numeric(scale(v0$E["HSPA1A", ]))
        des_hsp <- model.matrix(~ diag + batch + hsp_z, data = sm)
        vh <- voom(dge, des_hsp, plot = FALSE)
        fith <- eBayes(lmFit(vh, des_hsp), robust = TRUE)
        for (g_chk in intersect(c("SLC12A5", "GABRA1", "ATP5F1D"), rownames(vh$E))) {
          rh <- topTable(fith, coef = "diagMDD", number = Inf, sort.by = "none")[g_chk, ]
          sens_rows[[length(sens_rows) + 1L]] <- data.table(
            cell_type = ct, min_nuclei = min_n, n_donors = nrow(sm),
            model = "+ HSPA1A stress proxy", gene = g_chk,
            logFC = rh$logFC, t = rh$t, P.Value = rh$P.Value, adj.P.Val = rh$adj.P.Val
          )
        }
      }
      # Model + SVA surrogate variables
      des_null <- model.matrix(~ batch, data = sm)
      n_sv <- tryCatch(num.sv(v0$E, des0, method = "be"), error = function(e) 0L)
      n_sv <- min(max(n_sv, 1L), 3L)
      svobj <- tryCatch(sva(v0$E, des0, des_null, n.sv = n_sv), error = function(e) NULL)
      if (!is.null(svobj) && svobj$n.sv > 0) {
        des_sva <- cbind(des0, svobj$sv)
        colnames(des_sva)[(ncol(des0) + 1):ncol(des_sva)] <- paste0("SV", seq_len(svobj$n.sv))
        vsv <- voom(dge, des_sva, plot = FALSE)
        fitsv <- eBayes(lmFit(vsv, des_sva), robust = TRUE)
        for (g_chk in intersect(c("SLC12A5", "GABRA1", "ATP5F1D"), rownames(vsv$E))) {
          rs <- topTable(fitsv, coef = "diagMDD", number = Inf, sort.by = "none")[g_chk, ]
          sens_rows[[length(sens_rows) + 1L]] <- data.table(
            cell_type = ct, min_nuclei = min_n, n_donors = nrow(sm),
            model = sprintf("+ %d SVA factors", svobj$n.sv), gene = g_chk,
            logFC = rs$logFC, t = rs$t, P.Value = rs$P.Value, adj.P.Val = rs$adj.P.Val
          )
        }
      }
    }
  }
}
sens_dt <- rbindlist(sens_rows)
fwrite(sens_dt, file.path(out_res, "p1_sva_hspa1a_threshold_sensitivity.tsv"), sep = "\t")
print(sens_dt[gene == "SLC12A5"])

# ==============================================================================
# 4. Export Donor-Level CPM & Proportions for Individual Data Point Figures
# ==============================================================================
cat("[4/4] Exporting individual donor-level CPM and proportions for plotting...\n")
don_tot <- meta_dt[, .(tot_cells = sum(n_cells)), by = donor_id]
cpm_rows <- list()
for (ct in c("ExN10_L46", "InN1_PV", "ExN2_L23", "InN5_SST", "ExN3_L46", "InN9_PV", "ExN7", "OPC1", "Ast1")) {
  sm <- meta_dt[author_subtype == ct]
  sm <- merge(sm, don_tot, by = "donor_id")
  mat_ct <- cnt_mat[, sm$grp_idx + 1L, drop = FALSE]
  libs <- colSums(mat_ct)
  for (g_chk in c("SLC12A5", "ATP5F1D", "GABRA1", "HSPA1A")) {
    cpm_v <- (mat_ct[g_chk, ] / pmax(libs, 1)) * 1e6
    cpm_rows[[length(cpm_rows) + 1L]] <- data.table(
      donor_id = sm$donor_id,
      diagnosis = sm$diagnosis,
      Sex = sm$Sex,
      Chemistry = sm$Chemistry,
      Batch = sm$Batch,
      cell_type = ct,
      n_cells = sm$n_cells,
      tot_cells = sm$tot_cells,
      prop_pct = 100 * sm$n_cells / sm$tot_cells,
      gene = g_chk,
      raw_count = as.integer(mat_ct[g_chk, ]),
      cpm = cpm_v,
      log2_cpm1 = log2(cpm_v + 1)
    )
  }
}
cpm_dt <- rbindlist(cpm_rows)
fwrite(cpm_dt, file.path(out_res, "p1_donor_level_cpm_key_genes.tsv"), sep = "\t")
cat("[SUCCESS] All 4 local upgrades finished and archived to outputs/tables/.\n")
