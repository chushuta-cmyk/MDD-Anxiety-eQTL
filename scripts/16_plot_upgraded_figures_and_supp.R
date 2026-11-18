#!/usr/bin/env Rscript
# 16_plot_upgraded_figures_and_supp.R
# Renders Molecular Psychiatry-compliant Main Figures 2-4 (with individual donor points)
# and Supplementary Figures S1-S3 (GWAS overlap sensitivity, 30-gene 20q13.12 map,
# donor permutation null + ATP5F1D/HSPA1A audit).

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

req_pkgs <- c("data.table", "ggplot2", "patchwork", "ggrepel", "scales")
has_pkgs <- all(sapply(req_pkgs, requireNamespace, quietly = TRUE))
if (check_only) {
  cat(sprintf("UPGRADED_PLOT_PREFLIGHT: all_pkgs=%s\n", ifelse(has_pkgs, "TRUE", "FALSE")))
  quit(save = "no", status = 0)
}
if (!has_pkgs) stop("Missing required R plotting packages.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
tab_dir <- file.path(script_dir, "outputs", "tables")
fig_dir <- file.path(script_dir, "outputs", "figures")
art_dir <- Sys.getenv("ARTIFACT_DIR", unset = "")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_mp <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.subtitle = element_text(size = 8.4, color = "#4B5563", hjust = 0),
    axis.title = element_text(face = "bold", size = 9.2),
    axis.text = element_text(color = "#111827", size = 8.2),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.0),
    strip.background = element_rect(fill = "#F3F4F6", color = "#D1D5DB"),
    strip.text = element_text(face = "bold", size = 8.5)
  )

# ==============================================================================
# FIGURE 2: Bryois 8-Cell-Type Matrix + Full 177-Pair eQTL vs DE Concordance
# ==============================================================================
bdf <- fread(file.path(tab_dir, "aim3_bryois_8celltype_coloc_MDD_and_ANX.tsv"))
top_genes <- unique(bdf[PP.H4 >= 0.40]$gene)
sub_mdd <- bdf[gene %in% top_genes & trait == "MDD", .(PP.H4 = max(PP.H4, na.rm = TRUE)), by = .(gene, cell_type)]

pref_cols <- c("Excitatory.neurons", "Inhibitory.neurons", "Astrocytes", "Oligodendrocytes",
               "OPCs...COPs", "Microglia", "Endothelial.cells", "Pericytes")
col_lbls <- c("ExN", "InN", "Ast", "Oli", "OPC", "Mic", "End", "Per")
col_map <- setNames(col_lbls, pref_cols)

max_g <- sub_mdd[, .(mx = max(PP.H4, na.rm = TRUE)), by = gene][order(-mx)]
grid_dt <- CJ(gene = max_g$gene, cell_type = pref_cols)
grid_dt <- merge(grid_dt, sub_mdd, by = c("gene", "cell_type"), all.x = TRUE)
grid_dt[, cell_short := factor(col_map[cell_type], levels = col_lbls)]
grid_dt[, gene_f := factor(gene, levels = rev(max_g$gene))]
grid_dt[, lbl := ifelse(is.na(PP.H4), "NT", sprintf("%.2f", PP.H4))]
grid_dt[, txt_col := ifelse(is.na(PP.H4), "#6B7280", ifelse(PP.H4 > 0.60, "white", "#111827"))]

p2A <- ggplot(grid_dt, aes(x = cell_short, y = gene_f, fill = PP.H4)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = lbl, color = txt_col), size = 2.6, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradientn(colors = c("#F7FCF5", "#C7E9C0", "#74C476", "#238B45", "#00441B"),
                       na.value = "#E5E7EB", limits = c(0, 1), name = "PP.H4") +
  labs(title = "Cell-Type cis-eQTL Colocalization (MDD)",
       subtitle = "Bryois et al. (N = 192); NT = absent from P < 1e-3 summary file",
       x = NULL, y = NULL) +
  theme_mp +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")

conc_dt <- fread(file.path(tab_dir, "p1_bryois_all140_lineage_matched_eqtl_vs_de.tsv"))
conc_dt[, cat := ifelse(gene == "SLC12A5", "SLC12A5 (Top Discordance, PP.H4 = 0.991)",
                        ifelse(PP.H4 >= 0.50, "Other PP.H4 >= 0.50", "All Other Pairs (P_eQTL < 1e-3)"))]
conc_dt[, pt_lbl := ifelse(gene == "SLC12A5", sprintf("SLC12A5 (%s)", gsub("\\.neurons", "", cell_type)),
                           ifelse(gene %in% c("NCKIPSD", "NEK4", "GNL3") & PP.H4 >= 0.50,
                                  sprintf("%s (%s)", gene, substr(cell_type, 1, 3)),
                                  ifelse(gene == "DYNC1I2" & cell_type == "Excitatory.neurons", "DYNC1I2 (Exc)", NA)))]

r_all <- cor.test(conc_dt$eqtl_beta_on_risk_allele, conc_dt$logFC)
sub_txt <- sprintf("177 matched pairs (109 genes): r = %.3f (P = %.2f); 45.8%% concordant sign",
                   r_all$estimate, r_all$p.value)

p2B <- ggplot(conc_dt, aes(x = eqtl_beta_on_risk_allele, y = logFC)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FEF3C7", alpha = 0.35) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#FEF3C7", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.5) +
  geom_point(data = conc_dt[cat == "All Other Pairs (P_eQTL < 1e-3)"], aes(color = cat), alpha = 0.45, size = 1.8) +
  geom_point(data = conc_dt[cat != "All Other Pairs (P_eQTL < 1e-3)"], aes(color = cat), size = 3.2) +
  geom_label_repel(data = conc_dt[!is.na(pt_lbl)], aes(label = pt_lbl, color = cat),
                   size = 2.7, fontface = "bold", box.padding = 0.4, max.overlaps = 15, show.legend = FALSE) +
  scale_color_manual(values = c("SLC12A5 (Top Discordance, PP.H4 = 0.991)" = "#D95F02",
                                "Other PP.H4 >= 0.50" = "#3C5488",
                                "All Other Pairs (P_eQTL < 1e-3)" = "#9CA3AF")) +
  labs(title = "Lineage-Matched eQTL vs. Postmortem DE (109 Genes)",
       subtitle = sub_txt,
       x = "Risk-Allele cis-eQTL Beta (Bryois 2022)",
       y = "Postmortem Pooled Lineage log2FC (GSE213982)") +
  guides(color = guide_legend(nrow = 2)) +
  theme_mp + theme(legend.position = "bottom")

fig2 <- (p2A | p2B) + plot_layout(widths = c(1.05, 1.15)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Figure2_Bryois_CellType_and_177Pair_Polarity.png"), fig2, width = 13.6, height = 6.4, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "Figure2_Bryois_CellType_and_177Pair_Polarity.png"), fig2, width = 13.6, height = 6.4, dpi = 300)

# ==============================================================================
# FIGURE 3: Interaction, CAMERA Audit, Individual Donor SLC12A5 CPM, Individual Donor scDRS
# ==============================================================================
df_int <- fread(file.path(tab_dir, "p0_celltype_x_diagnosis_interaction.tsv"))
df_cam <- fread(file.path(tab_dir, "p0_limma_camera_batch_adjusted_pathways.tsv"))
cpm_dt <- fread(file.path(tab_dir, "p1_donor_level_cpm_key_genes.tsv"))
sc_don <- fread(file.path(tab_dir, "p1_scdrs_donor_subtype_means.tsv"))
sc_sum <- fread(file.path(tab_dir, "aim3_scdrs_38subtypes_summary.tsv"))

keep_c <- c("PV_eff", "ExN10_eff", "Int_PV_vs_ExN10", "Int_PV_vs_ExN7")
keep_g <- c("SLC12A5", "GABRA1", "ATP5F1D")
pa_dt <- df_int[contrast %in% keep_c & gene %in% keep_g]
pa_dt[, se := abs(logFC / t)]
pa_dt[t == 0, se := 0]
pa_dt[, ci_lo := logFC - 1.96 * se]
pa_dt[, ci_hi := logFC + 1.96 * se]
c_lbls <- c("PV_eff" = "InN1_PV (MDD vs Ctrl)", "ExN10_eff" = "ExN10_L46 (MDD vs Ctrl)",
            "Int_PV_vs_ExN10" = "Int: InN1_PV x ExN10_L46", "Int_PV_vs_ExN7" = "Int: InN1_PV x ExN7")
pa_dt[, contrast_f := factor(c_lbls[contrast], levels = rev(unname(c_lbls)))]
pa_dt[, gene := factor(gene, levels = keep_g)]
pa_dt[, sig_lbl := ifelse(adj.P.Val < 0.05, "*", "")]

p3A <- ggplot(pa_dt, aes(x = logFC, y = contrast_f, color = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.6) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.55), size = 0.55, linewidth = 0.6) +
  geom_text(aes(x = ci_lo - 0.04, label = sig_lbl), position = position_dodge(width = 0.55), hjust = 1, size = 3.8, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("SLC12A5" = "#00A087", "GABRA1" = "#3C5488", "ATP5F1D" = "#E64B35")) +
  coord_cartesian(xlim = c(-1.75, 0.72)) +
  labs(title = "Subtype x Diagnosis Mixed-Effect Model",
       subtitle = "duplicateCorrelation (rho = 0.211; *FDR < 0.05)",
       x = "log2 Fold Change (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 1)) +
  theme_mp + theme(legend.position = "bottom")

show_cts <- c("InN4_VIP", "ExN10_L46", "InN1_PV", "End1", "InN9_PV", "ExN15_L56", "InN5_SST", "ExN2_L23", "ExN7", "ExN4_L35")
pb_dt <- df_cam[pathway == "Nuclear_MitoCarta_OXPHOS" & cell_type %in% show_cts]
pb_dt[, cell_f := factor(cell_type, levels = rev(show_cts))]
pb_dt[, sig_cat := ifelse(PValue < 0.05, "P < 0.05 (uncorr)", ifelse(PValue < 0.10, "0.05 <= P < 0.10", "P >= 0.10"))]
pb_dt[, lbl_p := sprintf("P=%.3f (r=%.2f)", PValue, Correlation)]

p3B <- ggplot(pb_dt, aes(x = donor_mean_z_diff, y = cell_f, color = sig_cat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.6) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.55, linewidth = 0.6) +
  geom_text(aes(x = 0.28, label = lbl_p), hjust = 0, size = 2.75, color = "#1F2937", show.legend = FALSE) +
  scale_color_manual(values = c("P < 0.05 (uncorr)" = "#E64B35", "0.05 <= P < 0.10" = "#F39B7F", "P >= 0.10" = "#6B7280")) +
  coord_cartesian(xlim = c(-0.70, 0.75)) +
  labs(title = "Correlation-Adjusted limma::camera OXPHOS",
       subtitle = "0/34 subtypes pass FDR < 0.05; ExN7 resets to P = 0.676",
       x = "Donor-Level Nuclear OXPHOS Z Diff (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 1)) +
  theme_mp + theme(legend.position = "bottom")

# Panel C: Individual donor SLC12A5 log2(CPM+1) across ExN10_L46, InN1_PV, ExN7 by Chemistry wave
slc_don <- cpm_dt[gene == "SLC12A5" & cell_type %in% c("ExN10_L46", "InN1_PV", "ExN7") & n_cells >= 10]
slc_don[, wave := ifelse(Chemistry == "v3", "Female (v3)", "Male (v2)")]
slc_don[, diagnosis := factor(diagnosis, levels = c("Control", "MDD"))]
slc_don[, cell_f := factor(cell_type, levels = c("ExN10_L46", "InN1_PV", "ExN7"))]

p3C <- ggplot(slc_don, aes(x = cell_f, y = log2_cpm1, fill = diagnosis, color = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, position = position_dodge(width = 0.75), width = 0.58) +
  geom_point(position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75), size = 1.25, alpha = 0.75) +
  facet_wrap(~wave, ncol = 2) +
  scale_fill_manual(values = c("Control" = "#4B5563", "MDD" = "#E64B35")) +
  scale_color_manual(values = c("Control" = "#374151", "MDD" = "#C0392B")) +
  labs(title = "Individual Donor SLC12A5 Expression by Wave",
       subtitle = "Points = individual donors (n_nuclei >= 10); ExN10_L46 FDR = 0.019",
       x = NULL, y = "SLC12A5 log2(CPM + 1)") +
  theme_mp + theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1))

# Panel D: Donor-level scDRS with individual donor dots overlaid
show_sc <- c("InN1_PV", "InN10_ADARB2", "InN4_VIP", "InN9_PV", "OPC1", "ExN15_L56", "ExN10_L46", "ExN7")
sc_pts <- sc_don[author_subtype %in% show_sc]
sc_pts[, sub_f := factor(author_subtype, levels = rev(show_sc))]
sc_sum_sub <- sc_sum[author_subtype %in% show_sc]
sc_sum_sub[, sub_f := factor(author_subtype, levels = rev(show_sc))]
sc_sum_sub[, class := ifelse(grepl("^InN", author_subtype), "InN", ifelse(grepl("^OPC", author_subtype), "OPC", "ExN"))]

p3D <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.5) +
  geom_jitter(data = sc_pts, aes(x = norm_score, y = sub_f), height = 0.18, width = 0, size = 1.0, alpha = 0.28, color = "#6B7280") +
  geom_errorbar(data = sc_sum_sub, aes(y = sub_f, xmin = donor_mean_norm_score - donor_se_norm_score,
                                       xmax = donor_mean_norm_score + donor_se_norm_score, color = class),
                orientation = "y", width = 0.25, linewidth = 0.85) +
  geom_point(data = sc_sum_sub, aes(x = donor_mean_norm_score, y = sub_f, color = class), size = 3.0) +
  scale_color_manual(values = c("InN" = "#00A087", "OPC" = "#E7298A", "ExN" = "#3C5488")) +
  labs(title = "Donor-Level Polygenic Risk (scDRS)",
       subtitle = "Gray dots = individual donors; colored point/bar = mean +/- SE",
       x = "Mean Donor scDRS Normalized Score", y = NULL) +
  theme_mp + theme(legend.position = "bottom")

fig3 <- (p3A | p3B) / (p3C | p3D) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Figure3_Interaction_Camera_Sex_scDRS_Points.png"), fig3, width = 13.8, height = 10.4, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "Figure3_Interaction_Camera_Sex_scDRS_Points.png"), fig3, width = 13.8, height = 10.4, dpi = 300)

# ==============================================================================
# FIGURE 4: Individual Donor Proportions + SLC12A5 SVA/Threshold Sensitivity
# ==============================================================================
prop_dt <- cpm_dt[gene == "SLC12A5" & cell_type %in% c("OPC1", "Ast1", "ExN10_L46", "InN1_PV", "ExN7")]
prop_dt[, diagnosis := factor(diagnosis, levels = c("Control", "MDD"))]
prop_dt[, cell_f := factor(cell_type, levels = c("OPC1", "Ast1", "ExN10_L46", "InN1_PV", "ExN7"))]
fdr_lbls <- data.table(
  cell_f = factor(c("OPC1", "Ast1", "ExN10_L46", "InN1_PV", "ExN7"), levels = c("OPC1", "Ast1", "ExN10_L46", "InN1_PV", "ExN7")),
  lbl = c("FDR=0.002", "FDR=0.004", "P=0.019\nFDR=0.102", "FDR=0.487", "P=0.026\nFDR=0.108"),
  y_pos = c(12.8, 28.6, 8.2, 17.2, 13.4)
)

p4A <- ggplot(prop_dt, aes(x = cell_f, y = prop_pct, fill = diagnosis, color = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, position = position_dodge(width = 0.75), width = 0.58) +
  geom_point(position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75), size = 1.2, alpha = 0.7) +
  geom_text(data = fdr_lbls, aes(x = cell_f, y = y_pos, label = lbl), inherit.aes = FALSE, size = 2.6, fontface = "bold", color = "#111827") +
  scale_fill_manual(values = c("Control" = "#4B5563", "MDD" = "#E64B35")) +
  scale_color_manual(values = c("Control" = "#374151", "MDD" = "#C0392B")) +
  coord_cartesian(ylim = c(0, 31)) +
  labs(title = "Individual Donor Nuclear Proportions (n = 71)",
       subtitle = "Batch- and sex-adjusted logit regression; nominal depletion in ExN10_L46 & ExN7",
       x = NULL, y = "Subtype Proportion (% of QC-Pass Nuclei)") +
  theme_mp + theme(legend.position = "bottom")

sens_dt <- fread(file.path(tab_dir, "p1_sva_hspa1a_threshold_sensitivity.tsv"))
sens_slc <- sens_dt[gene == "SLC12A5" & cell_type %in% c("ExN10_L46", "InN1_PV", "Pooled_ExN", "Pooled_InN")]
sens_slc[, se := abs(logFC / t)]
sens_slc[, mod_lbl := ifelse(min_nuclei == 10L, model, sprintf("Min >= %d nuclei", min_nuclei))]
mod_levels <- c("Base (diag + Batch)", "+ HSPA1A stress proxy", "+ 3 SVA factors", "Min >= 20 nuclei", "Min >= 50 nuclei", "Min >= 5 nuclei")
sens_slc <- sens_slc[mod_lbl %in% mod_levels]
sens_slc[, mod_f := factor(mod_lbl, levels = rev(mod_levels))]
sens_slc[, cell_f := factor(cell_type, levels = c("ExN10_L46", "InN1_PV", "Pooled_InN", "Pooled_ExN"))]

p4B <- ggplot(sens_slc, aes(x = logFC, y = mod_f, color = cell_f)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.55) +
  geom_pointrange(aes(xmin = logFC - 1.96 * se, xmax = logFC + 1.96 * se),
                  position = position_dodge(width = 0.6), size = 0.5, linewidth = 0.55) +
  scale_color_manual(values = c("ExN10_L46" = "#D95F02", "InN1_PV" = "#00A087",
                                "Pooled_InN" = "#3C5488", "Pooled_ExN" = "#7570B3")) +
  labs(title = "SLC12A5 Sensitivity Across Thresholds, Stress Proxy & SVA",
       subtitle = "Persists across min-nuclei & HSPA1A proxy; attenuates under 3 unconstrained SVs",
       x = "log2 Fold Change (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 2)) +
  theme_mp + theme(legend.position = "bottom")

fig4 <- (p4A | p4B) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Figure4_Composition_Points_and_Sensitivity.png"), fig4, width = 13.6, height = 6.2, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "Figure4_Composition_Points_and_Sensitivity.png"), fig4, width = 13.6, height = 6.2, dpi = 300)

# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Full 30-Gene x 8-Cell-Type Heatmap at 20q13.12
# ==============================================================================
w20 <- fread(file.path(tab_dir, "p1_20q13_all30genes_bryois_coloc.tsv"))
w20_mdd <- w20[trait == "MDD", .(PP.H4 = max(PP.H4, na.rm = TRUE)), by = .(gene, cell_type)]
all_g20 <- w20_mdd[, .(mx = max(PP.H4, na.rm = TRUE)), by = gene][order(-mx)]
grid20 <- CJ(gene = all_g20$gene, cell_type = pref_cols)
grid20 <- merge(grid20, w20_mdd, by = c("gene", "cell_type"), all.x = TRUE)
grid20[, cell_short := factor(col_map[cell_type], levels = col_lbls)]
grid20[, gene_f := factor(gene, levels = rev(all_g20$gene))]
grid20[, lbl := ifelse(is.na(PP.H4), "NT", sprintf("%.2f", PP.H4))]
grid20[, txt_col := ifelse(is.na(PP.H4), "#6B7280", ifelse(PP.H4 > 0.60, "white", "#111827"))]

fig_s2 <- ggplot(grid20, aes(x = cell_short, y = gene_f, fill = PP.H4)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = lbl, color = txt_col), size = 2.4, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradientn(colors = c("#F7FCF5", "#C7E9C0", "#74C476", "#238B45", "#00441B"),
                       na.value = "#E5E7EB", limits = c(0, 1), name = "MDD PP.H4") +
  labs(title = "Complete 20q13.12 (chr20:44.15-45.20 Mb) Cell-Type cis-eQTL Colocalization (30 Genes)",
       subtitle = "CD40 passes P_eQTL < 1e-3 only in Pericytes (PP.H4 = 0.054); 29 neighbors show PP.H4 <= 0.198",
       x = NULL, y = NULL) +
  theme_mp + theme(axis.line = element_blank(), axis.ticks = element_blank())

ggsave(file.path(fig_dir, "FigureS2_20q13_All30Genes.png"), fig_s2, width = 8.8, height = 8.4, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "FigureS2_20q13_All30Genes.png"), fig_s2, width = 8.8, height = 8.4, dpi = 300)

# ==============================================================================
# SUPPLEMENTARY FIGURE S3: Donor Permutation Null + Donor-Level ATP5F1D & HSPA1A
# ==============================================================================
perm_null <- fread(file.path(tab_dir, "p1_slc12a5_donor_permutation_null.tsv"))
perm_sum  <- fread(file.path(tab_dir, "p1_slc12a5_donor_permutation_summary.tsv"))

pS3A <- ggplot(perm_null, aes(x = mean_t)) +
  geom_histogram(bins = 35, fill = "#9CA3AF", color = "white") +
  geom_vline(xintercept = perm_sum$obs_mean_t[1], color = "#E64B35", linewidth = 1.1, linetype = "dashed") +
  annotate("text", x = perm_sum$obs_mean_t[1] + 0.08, y = 85, hjust = 0,
           label = sprintf("Observed mean t = %.3f\nEmp P = %.3f (1,000 perms)\n7/25 P < 0.05 (Emp P = %.3f)",
                           perm_sum$obs_mean_t[1], perm_sum$empirical_p_mean_t[1], perm_sum$empirical_p_nominal_neg_count[1]),
           color = "#B91C1C", fontface = "bold", size = 3.0) +
  labs(title = "Within-Batch Donor Permutation Null (25 Subtypes)",
       subtitle = "Preserves cross-subtype donor correlation structure across 71 donors",
       x = "Mean Neuronal SLC12A5 t-Statistic", y = "Permutations") +
  theme_mp

atp_hsp <- cpm_dt[gene %in% c("ATP5F1D", "HSPA1A") & cell_type %in% c("InN1_PV", "ExN10_L46") & n_cells >= 10]
atp_hsp[, wave := ifelse(Chemistry == "v3", "Female (v3)", "Male (v2)")]
atp_hsp[, diagnosis := factor(diagnosis, levels = c("Control", "MDD"))]
atp_hsp[, grp_lbl := paste0(cell_type, "\n", wave)]

pS3B <- ggplot(atp_hsp, aes(x = grp_lbl, y = log2_cpm1, fill = diagnosis, color = diagnosis)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, position = position_dodge(width = 0.75), width = 0.58) +
  geom_point(position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75), size = 1.15, alpha = 0.7) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#4B5563", "MDD" = "#E64B35")) +
  scale_color_manual(values = c("Control" = "#374151", "MDD" = "#C0392B")) +
  labs(title = "Donor-Level ATP5F1D & HSPA1A Across Waves",
       subtitle = "Zero zero-count donors in InN1_PV (69/69); ATP5F1D tracks HSPA1A drop in both waves",
       x = NULL, y = "log2(CPM + 1)") +
  theme_mp + theme(legend.position = "bottom")

fig_s3 <- (pS3A | pS3B) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "FigureS3_Permutation_and_ATP5F1D_Audit.png"), fig_s3, width = 13.2, height = 6.8, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "FigureS3_Permutation_and_ATP5F1D_Audit.png"), fig_s3, width = 13.2, height = 6.8, dpi = 300)

# ==============================================================================
# SUPPLEMENTARY FIGURE S1: GWAS Overlap Sensitivity Forest Across 7 Loci
# ==============================================================================
ov_dt <- fread(file.path(tab_dir, "p1_gwas_overlap_coloc_sensitivity.tsv"))
ov_dt[, loc_short := gsub("_chr.*", "", gsub("chr3p21_31_.*", "chr3p21.31", gsub("chr3p21_1_.*", "chr3p21.1", locus)))]
ov_dt[, loc_short := gsub("_CDK2AP1", "", loc_short)]
ov_dt[, model_lbl := factor(model,
                            levels = rev(c("MDD_all_vs_ANX", "MDD_noUKBB_vs_ANX", "MDD_all_vs_ANX_se_inflated_LDSCint", "MDD_Clin_vs_ANX")),
                            labels = rev(c("MDD all vs ANX", "MDD no-UKBB vs ANX", "MDD all + LDSC SE penalty", "MDD Clinical-only vs ANX")))]

fig_s1_ov <- ggplot(ov_dt, aes(x = PP.H4, y = model_lbl, color = model_lbl)) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "#6B7280", linewidth = 0.5) +
  geom_point(size = 3.0) +
  geom_segment(aes(x = 0, xend = PP.H4, y = model_lbl, yend = model_lbl), linewidth = 0.9) +
  facet_wrap(~loc_short, ncol = 4) +
  scale_color_manual(values = c("MDD all vs ANX" = "#1E3A5F", "MDD no-UKBB vs ANX" = "#00A087",
                                "MDD all + LDSC SE penalty" = "#D95F02", "MDD Clinical-only vs ANX" = "#7570B3")) +
  coord_cartesian(xlim = c(0, 1.05)) +
  labs(title = "Sample-Overlap & Cohort-Subsetting Sensitivity (Single-Causal coloc.abf)",
       subtitle = "SLC12A5 maintains PP.H4 >= 0.982 without UKBB or under LDSC-intercept (0.323) SE inflation",
       x = "coloc.abf PP.H4", y = NULL) +
  theme_mp + theme(legend.position = "none")

ggsave(file.path(fig_dir, "FigureS1_GWAS_Overlap_Sensitivity.png"), fig_s1_ov, width = 13.0, height = 5.4, dpi = 300)
if (nzchar(art_dir)) ggsave(file.path(art_dir, "FigureS1_GWAS_Overlap_Sensitivity.png"), fig_s1_ov, width = 13.0, height = 5.4, dpi = 300)

cat("[SUCCESS] Rendered Main Figures 2-4 (with individual donor points) & Supp Figures S1-S3!\n")
