#!/usr/bin/env Rscript
# 13_plot_publication_figures_ggplot2_pheatmap.R
# Renders authentic publication-grade figures using native R pheatmap, ggplot2, ggrepel, and patchwork.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check-only" %in% args

req_pkgs <- c("data.table", "ggplot2", "patchwork", "pheatmap", "ggrepel", "RColorBrewer", "scales")
has_pkgs <- all(sapply(req_pkgs, requireNamespace, quietly = TRUE))
if (check_only) {
  cat(sprintf("P0_PLOT_R_PREFLIGHT: all_pkgs=%s\n", ifelse(has_pkgs, "TRUE", "FALSE")))
  quit(save = "no", status = 0)
}
if (!has_pkgs) stop("Missing required R plotting packages.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(RColorBrewer)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
tab_dir <- file.path(script_dir, "outputs", "tables")
fig_dir <- file.path(script_dir, "outputs", "figures")
art_dir <- fig_dir
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Native pheatmap: Bryois 2022 8-Brain-Cell-Type cis-eQTL Colocalization
# ==============================================================================
bdf <- fread(file.path(tab_dir, "aim3_bryois_8celltype_coloc_MDD_and_ANX.tsv"))
top_genes <- unique(bdf[PP.H4 >= 0.40]$gene)
sub <- bdf[gene %in% top_genes & trait == "MDD"]
# Aggregate max PP.H4 per gene x cell_type
sub_agg <- sub[, .(PP.H4 = max(PP.H4, na.rm = TRUE), locus = locus[1], chr = chr[1]), by = .(gene, cell_type)]
wide <- dcast(sub_agg, gene ~ cell_type, value.var = "PP.H4")
mat <- as.matrix(wide[, -1, with = FALSE])
rownames(mat) <- wide$gene

pref_cols <- c("Excitatory.neurons", "Inhibitory.neurons", "Astrocytes", "Oligodendrocytes",
               "OPCs...COPs", "Microglia", "Endothelial.cells", "Pericytes")
mat <- mat[, intersect(pref_cols, colnames(mat)), drop = FALSE]
colnames(mat) <- gsub("\\.\\.\\.", "/", gsub("\\.", " ", colnames(mat)))

# Sort rows cleanly by max PP.H4 descending
ord <- order(apply(mat, 1, function(x) max(x, na.rm = TRUE)), decreasing = TRUE)
mat <- mat[ord, , drop = FALSE]

# Build locus chromosome row annotation
gene_chr <- sub_agg[, .(Chr = paste0("chr", chr[1])), by = gene]
ann_row <- data.frame(Locus = factor(gene_chr[match(rownames(mat), gene)]$Chr))
rownames(ann_row) <- rownames(mat)

num_mat <- ifelse(is.na(mat), "NT", sprintf("%.2f", mat))
pal_ylgn <- colorRampPalette(brewer.pal(9, "YlGn"))(100)
breaks_vec <- seq(0, 1, length.out = 101)

render_bryois_pheatmap <- function(out_path, is_pdf = FALSE) {
  if (is_pdf) {
    pdf(out_path, width = 9.0, height = 6.9)
  } else {
    png(out_path, width = 2700, height = 2070, res = 300)
  }
  pheatmap(
    mat,
    color = pal_ylgn,
    breaks = breaks_vec,
    na_col = "#ECEFF4",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = num_mat,
    number_color = ifelse(is.na(mat), "#7A869A", ifelse(mat > 0.60, "white", "#1F2937")),
    fontsize = 9.0,
    fontsize_number = 7.8,
    border_color = "white",
    cellwidth = 36,
    cellheight = 17,
    angle_col = 45,
    annotation_row = ann_row,
    main = "Bryois 2022 Cell-Type cis-eQTL Coloc (MDD PP.H4; NT = Not Tested)"
  )
  dev.off()
}

render_bryois_pheatmap(file.path(fig_dir, "Figure_Aim3_Bryois_CellType_eQTL_Coloc.png"), FALSE)
render_bryois_pheatmap(file.path(fig_dir, "Figure_Aim3_Bryois_CellType_eQTL_Coloc.pdf"), TRUE)
render_bryois_pheatmap(file.path(art_dir, "Figure_Aim3_Bryois_CellType_eQTL_Coloc.png"), FALSE)

# ==============================================================================
# 2. Native ggplot2 + patchwork: Joint Interaction + CAMERA OXPHOS + Sex + scDRS
# ==============================================================================
df_int <- fread(file.path(tab_dir, "p0_celltype_x_diagnosis_interaction.tsv"))
df_cam <- fread(file.path(tab_dir, "p0_limma_camera_batch_adjusted_pathways.tsv"))
df_sex <- fread(file.path(tab_dir, "p0_sex_stratified_batch_adjusted_audit.tsv"))
df_sc  <- fread(file.path(tab_dir, "aim3_scdrs_38subtypes_summary.tsv"))

theme_nat <- theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.subtitle = element_text(size = 8.5, color = "#4B5563", hjust = 0),
    axis.title = element_text(face = "bold", size = 9.5),
    axis.text = element_text(color = "#111827", size = 8.5),
    legend.background = element_rect(fill = "white", color = "#D1D5DB", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 8.2),
    strip.background = element_blank()
  )

# Panel A: Mixed-Model Interaction Forest Plot (Pointrange with approximate 95% CI = logFC +/- 1.96*se)
c_keep <- c("PV_eff", "ExN10_eff", "Int_PV_vs_ExN10", "Int_PV_vs_ExN7")
c_labs <- c("PV_eff" = "InN1_PV\n(MDD vs Ctrl)", "ExN10_eff" = "ExN10_L46\n(MDD vs Ctrl)",
            "Int_PV_vs_ExN10" = "Interaction:\nInN1_PV vs ExN10", "Int_PV_vs_ExN7" = "Interaction:\nInN1_PV vs ExN7")
pa_dt <- df_int[gene %in% c("SLC12A5", "GABRA1", "ATP5F1D") & contrast %in% c_keep]
pa_dt[, se := abs(logFC) / pmax(abs(t), 1e-4)]
pa_dt[, `:=`(ci_lo = logFC - 1.96 * se, ci_hi = logFC + 1.96 * se)]
pa_dt[, contrast_f := factor(contrast, levels = rev(c_keep), labels = rev(c_labs[c_keep]))]
pa_dt[, sig_lbl := ifelse(adj.P.Val < 0.001, "***", ifelse(adj.P.Val < 0.01, "**", ifelse(adj.P.Val < 0.05, "*", "")))]
pa_dt[, gene := factor(gene, levels = c("SLC12A5", "GABRA1", "ATP5F1D"))]

pA <- ggplot(pa_dt, aes(x = logFC, y = contrast_f, color = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.6) +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), position = position_dodge(width = 0.55), size = 0.55, fatten = 2.6) +
  geom_text(aes(x = ci_lo - 0.04, label = sig_lbl), position = position_dodge(width = 0.55), hjust = 1, size = 3.8, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = c("SLC12A5" = "#00A087", "GABRA1" = "#3C5488", "ATP5F1D" = "#E64B35")) +
  coord_cartesian(xlim = c(-1.75, 0.72)) +
  labs(title = "Subtype x Diagnosis Mixed-Effect Model",
       subtitle = "Batch-adjusted voom + duplicateCorrelation (rho = 0.211; *FDR < 0.05)",
       x = "log2 Fold Change (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 1)) +
  theme_nat + theme(legend.position = "bottom")

# Panel B: limma::camera Nuclear-OXPHOS Forest Plot
show_cts <- c("InN4_VIP", "ExN10_L46", "InN1_PV", "End1", "InN9_PV", "ExN15_L56", "InN5_SST", "ExN2_L23", "ExN7", "ExN4_L35")
pb_dt <- df_cam[pathway == "Nuclear_MitoCarta_OXPHOS" & cell_type %in% show_cts]
pb_dt[, cell_f := factor(cell_type, levels = rev(show_cts))]
pb_dt[, sig_cat := ifelse(PValue < 0.05, "P < 0.05", ifelse(PValue < 0.10, "0.05 <= P < 0.10", "P >= 0.10"))]
pb_dt[, lbl_p := sprintf("P=%.3f (r=%.2f)", PValue, Correlation)]

pB <- ggplot(pb_dt, aes(x = donor_mean_z_diff, y = cell_f, color = sig_cat)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280", linewidth = 0.6) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.55, fatten = 2.6) +
  geom_text(aes(x = 0.28, label = lbl_p), hjust = 0, size = 2.85, color = "#1F2937", show.legend = FALSE) +
  scale_color_manual(values = c("P < 0.05" = "#E64B35", "0.05 <= P < 0.10" = "#F39B7F", "P >= 0.10" = "#6B7280")) +
  coord_cartesian(xlim = c(-0.70, 0.72)) +
  labs(title = "Correlation-Adjusted limma::camera OXPHOS Audit",
       subtitle = "ExN7 naive artifact resets to P = 0.676; modest reduction in InN/ExN10",
       x = "Donor-Level Nuclear OXPHOS Z Diff (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 1)) +
  theme_nat + theme(legend.position = "bottom")

# Panel C: Sex-Stratified Replication of SLC12A5 & ATP5F1D
sex_subs <- c("ExN10_L46", "InN1_PV", "ExN2_L23", "InN5_SST", "ExN3_L46", "InN4_VIP")
pc_dt <- df_sex[gene == "SLC12A5" & cell_type %in% sex_subs]
pc_dt[, cell_f := factor(cell_type, levels = rev(sex_subs))]
pc_dt[, stratum := ifelse(sex_stratum == "Female", "Female (n = 38)", "Male (n = 33)")]

pC <- ggplot(pc_dt, aes(x = logFC, y = cell_f, fill = stratum)) +
  geom_vline(xintercept = 0, color = "#111827", linewidth = 0.5) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Female (n = 38)" = "#E7298A", "Male (n = 33)" = "#3C5488")) +
  labs(title = "Sex-Stratified SLC12A5 Reproducibility",
       subtitle = "Batch-adjusted limma-voom within independent sex strata",
       x = "Batch-Adjusted log2 Fold Change (SLC12A5)", y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_nat + theme(legend.position = "bottom")

# Panel D: Donor-Level scDRS (n >= 200 cells, N >= 20 donors)
pd_dt <- df_sc[reliable_n200_d20 == TRUE][order(-donor_mean_norm_score)][1:12]
pd_dt[, class := ifelse(grepl("^InN", author_subtype), "Inhibitory (InN)",
                        ifelse(grepl("^OPC", author_subtype), "OPC", "Excitatory (ExN)"))]
pd_dt[, sub_f := factor(author_subtype, levels = rev(author_subtype))]

pD <- ggplot(pd_dt, aes(x = donor_mean_norm_score, y = sub_f, fill = class)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = donor_mean_norm_score - donor_se_norm_score,
                     xmax = donor_mean_norm_score + donor_se_norm_score), height = 0.28, linewidth = 0.45) +
  scale_fill_manual(values = c("Inhibitory (InN)" = "#00A087", "OPC" = "#E7298A", "Excitatory (ExN)" = "#3C5488")) +
  coord_cartesian(xlim = c(0, 1.45)) +
  labs(title = "Donor-Level Polygenic Disease Score (scDRS)",
       subtitle = "Subtypes with >= 200 cells across 71 donors (excluding sparse ExN17 n = 12)",
       x = "Mean scDRS Score Across Donors (+/- SE)", y = NULL) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_nat + theme(legend.position = "bottom")

fig_mech <- (pA | pB) / (pC | pD) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Figure_Aim3_EI_Chloride_Mitochondria_Mechanism.png"), fig_mech, width = 13.6, height = 10.2, dpi = 300)
ggsave(file.path(fig_dir, "Figure_Aim3_EI_Chloride_Mitochondria_Mechanism.pdf"), fig_mech, width = 13.6, height = 10.2)
ggsave(file.path(art_dir, "Figure_Aim3_EI_Chloride_Mitochondria_Mechanism.png"), fig_mech, width = 13.6, height = 10.2, dpi = 300)

# ==============================================================================
# 3. Native ggplot2 + ggrepel + patchwork: Allelic Reconciliation & Composition
# ==============================================================================
df_lim  <- fread(file.path(tab_dir, "p0_limma_voom_batch_adjusted_key_genes.tsv"))
df_comp <- fread(file.path(tab_dir, "p0_celltype_composition_71donors.tsv"))

s_ex10 <- df_lim[gene == "SLC12A5" & cell_type == "ExN10_L46"]$logFC[1]
s_pv   <- df_lim[gene == "SLC12A5" & cell_type == "InN1_PV"]$logFC[1]

pts_dt <- data.table(
  label  = c("SLC12A5 (ExN10_L46)", "SLC12A5 (InN1_PV)", "NCKIPSD (ExN3_L46)",
             "CACNA2D2 (ANX InN -> PV)", "DYNC1I2 (ExN2_L23)", "GNL3 (InN -> ExN3)"),
  eqtl_b = c(0.230, 0.303, -0.222, -0.200, 0.288, 0.232),
  de_lfc = c(s_ex10, s_pv, -0.424, -0.189, 0.190, -0.458),
  status = c("Discordant (+eQTL / -DE)", "Discordant (+eQTL / -DE)", "Concordant (eQTL & DE sign match)",
             "Concordant (eQTL & DE sign match)", "Concordant (eQTL & DE sign match)", "Lineage Mismatch (InN vs ExN)")
)

pE <- ggplot(pts_dt, aes(x = eqtl_b, y = de_lfc, color = status)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#FEF3C7", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#D1FAE5", alpha = 0.35) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#D1FAE5", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#6B7280") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280") +
  geom_point(size = 3.6, stroke = 0.8) +
  geom_label_repel(aes(label = label), size = 3.0, fontface = "bold", box.padding = 0.45, point.padding = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("Discordant (+eQTL / -DE)" = "#D95F02",
                                "Concordant (eQTL & DE sign match)" = "#00A087",
                                "Lineage Mismatch (InN vs ExN)" = "#3C5488")) +
  coord_cartesian(xlim = c(-0.34, 0.42), ylim = c(-0.55, 0.32)) +
  labs(title = "Allelic Direction Reconciliation (eQTL vs Postmortem DE)",
       subtitle = "Green: concordant genetic & case-control sign; Amber: SLC12A5 paradox",
       x = "Risk-Allele cis-eQTL Beta (Bryois 2022)", y = "Postmortem DLPFC Case-Control log2FC") +
  guides(color = guide_legend(nrow = 2)) +
  theme_nat + theme(legend.position = "bottom")

comp_sub <- rbind(df_comp[order(P.Value)][1:7], df_comp[cell_type %in% c("InN1_PV", "InN5_SST", "ExN3_L46", "ExN2_L23")])
comp_sub <- unique(comp_sub, by = "cell_type")[order(P.Value)]
comp_sub[, sig := ifelse(FDR < 0.05, "FDR < 0.05 (Depleted)", "Stable (FDR >= 0.05)")]
comp_sub[, lbl := sprintf("%s (FDR = %.3f)", cell_type, FDR)]
comp_sub[, cell_f := factor(lbl, levels = rev(lbl))]

pF <- ggplot(comp_sub, aes(x = logit_beta, y = cell_f, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#6B7280") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.55, fatten = 2.6) +
  scale_color_manual(values = c("FDR < 0.05 (Depleted)" = "#E64B35", "Stable (FDR >= 0.05)" = "#4B5563")) +
  labs(title = "71-Donor Cellular Composition Stability",
       subtitle = "InN1_PV & pyramidal proportions unchanged; per-donor LIANA+ had 0/30,235 FDR < 0.05",
       x = "Batch & Sex-Adjusted Logit Proportion Beta (95% CI)", y = NULL) +
  guides(color = guide_legend(nrow = 1)) +
  theme_nat + theme(legend.position = "bottom")

fig_rec <- (pE | pF) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "Figure_Aim3_Synaptic_LR_and_Translational_Targets.png"), fig_rec, width = 13.0, height = 6.1, dpi = 300)
ggsave(file.path(fig_dir, "Figure_Aim3_Synaptic_LR_and_Translational_Targets.pdf"), fig_rec, width = 13.0, height = 6.1)
ggsave(file.path(art_dir, "Figure_Aim3_Synaptic_LR_and_Translational_Targets.png"), fig_rec, width = 13.0, height = 6.1, dpi = 300)

cat("[SUCCESS] All native R pheatmap & ggplot2/patchwork/ggrepel publication figures saved!\n")
