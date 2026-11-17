#!/usr/bin/env python3
"""
Publication-Quality Figure for Aim 2:
Multi-Causal Fine-Mapping & Cross-Trait Colocalization (MDD vs ANX)

Produces:
  - Figure_Aim2_SuSiE_Colocalization_MDD_vs_ANX.png (300 DPI)
  - Figure_Aim2_SuSiE_Colocalization_MDD_vs_ANX.pdf (Vector graphics)
"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Paths
script_dir = os.path.dirname(os.path.abspath(__file__))
tables_dir = os.path.join(script_dir, 'outputs', 'tables')
loci_dir   = os.path.join(script_dir, 'outputs', 'loci_dualtrait')
out_dir    = os.path.join(script_dir, 'outputs', 'figures')
os.makedirs(out_dir, exist_ok=True)

susie_tsv = os.path.join(tables_dir, 'aim2_coloc_susie_MDD_vs_ANX.tsv')
slc_tsv   = os.path.join(loci_dir, 'SLC12A5_chr20_dualtrait_grch37.tsv')
dcc_tsv   = os.path.join(loci_dir, 'DCC_chr18_dualtrait_grch37.tsv')

df_susie = pd.read_csv(susie_tsv, sep='\t')
df_slc   = pd.read_csv(slc_tsv, sep='\t')
df_dcc   = pd.read_csv(dcc_tsv, sep='\t')

# High-end journal aesthetic style (Nature/Cell standard)
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 9
plt.rcParams['axes.linewidth'] = 0.8

fig = plt.figure(figsize=(13, 10.5), dpi=300)
gs = gridspec.GridSpec(2, 2, height_ratios=[1.1, 1.3], width_ratios=[1.1, 1.0], hspace=0.38, wspace=0.28)

# ==================== Panel A: Coloc SuSiE Posterior Probability Breakdown ====================
axA = fig.add_subplot(gs[0, :])
top_loci = df_susie.sort_values('PP.H4.abf', ascending=True).copy()
labels = [f"{r['locus']}\n({r['hit1']} vs {r['hit2']})" for _, r in top_loci.iterrows()]
y_pos = np.arange(len(top_loci))

bar_h4 = axA.barh(y_pos, top_loci['PP.H4.abf'], height=0.55, color='#2ca02c', alpha=0.9, label='PP.H4 (Shared Causal Signal)')
bar_h3 = axA.barh(y_pos, top_loci['PP.H3.abf'], left=top_loci['PP.H4.abf'], height=0.55, color='#1f77b4', alpha=0.75, label='PP.H3 (Distinct Causal Signals)')
bar_other = axA.barh(y_pos, 1.0 - (top_loci['PP.H4.abf'] + top_loci['PP.H3.abf']), left=top_loci['PP.H4.abf'] + top_loci['PP.H3.abf'], height=0.55, color='#d62728', alpha=0.3, label='PP.H0/H1/H2 (Single-trait / No assoc)')

axA.set_yticks(y_pos)
axA.set_yticklabels(labels, fontsize=8)
axA.set_xlabel('Posterior Probability (coloc.susie multi-causal model)', fontweight='bold')
axA.set_xlim(0, 1.02)
axA.axvline(0.8, color='black', linestyle='--', linewidth=0.8, label='Shared Causal Threshold (PP.H4 >= 0.80)')
axA.set_title('A. Multi-Causal Fine-Mapping & Cross-Trait Colocalization between MDD and ANX', loc='left', fontweight='bold', fontsize=11, pad=12)
axA.legend(loc='lower right', frameon=True, fontsize=8, edgecolor='none')
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)

# Annotate SLC12A5 standout result
axA.text(0.985, len(top_loci)-1, ' PP.H4 = 98.1% (KCC2 shared causal variant)', va='center', ha='right', color='white', fontweight='bold', fontsize=8)

# ==================== Panel B: SLC12A5 (chr20) Mirror Manhattan Plot ====================
axB = fig.add_subplot(gs[1, 0])
bp_mb = df_slc['BP'] / 1e6
logp_mdd = -np.log10(np.clip(df_slc['P_MDD'], 1e-300, 1.0))
logp_anx = -np.log10(np.clip(df_slc['P_ANX'], 1e-300, 1.0))

# Upper mirror: MDD
axB.scatter(bp_mb, logp_mdd, s=12, color='#1b9e77', alpha=0.7, edgecolors='none', label='MDD (N=1.15M)')
# Lower mirror: ANX (inverted)
axB.scatter(bp_mb, -logp_anx, s=12, color='#d95f02', alpha=0.7, edgecolors='none', label='ANX (N=1.10M)')

axB.axhline(0, color='gray', linewidth=0.8)
axB.axhline(-np.log10(5e-8), color='#1b9e77', linestyle=':', linewidth=0.7)
axB.axhline(np.log10(5e-8), color='#d95f02', linestyle=':', linewidth=0.7)

# Set comfortable y-limits with ample headroom for title and labels
axB.set_ylim(-12, 23)

# Highlight shared lead hits
lead_mdd = df_slc.loc[df_slc['SNP'] == 'rs9074']
lead_anx = df_slc.loc[df_slc['SNP'] == 'rs6032660']

if not lead_mdd.empty:
    x_val = lead_mdd['BP'].values[0] / 1e6
    y_val = -np.log10(lead_mdd['P_MDD'].values[0])
    axB.scatter([x_val], [y_val], s=70, color='black', edgecolors='gold', linewidth=1.5, zorder=5)
    # Position MDD label clearly in the open upper-left quadrant below peak with zero title collision
    axB.annotate(
        'rs9074 (MDD lead)\nPP.H4 = 0.981',
        xy=(x_val, y_val),
        xytext=(x_val - 0.25, y_val - 2.8),
        arrowprops=dict(arrowstyle='->', lw=1.2, color='#1b9e77', shrinkA=2, shrinkB=4),
        fontsize=8,
        fontweight='bold',
        color='#1b9e77',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#1b9e77', alpha=0.95, lw=0.9)
    )

if not lead_anx.empty:
    x_val = lead_anx['BP'].values[0] / 1e6
    y_val = -(-np.log10(lead_anx['P_ANX'].values[0]))
    axB.scatter([x_val], [y_val], s=70, color='black', edgecolors='cyan', linewidth=1.5, zorder=5)
    # Position ANX label symmetrically
    axB.annotate(
        'rs6032660 (ANX lead)',
        xy=(x_val, y_val),
        xytext=(x_val - 0.25, y_val - 2.2),
        arrowprops=dict(arrowstyle='->', lw=1.2, color='#d95f02', shrinkA=2, shrinkB=4),
        fontsize=8,
        fontweight='bold',
        color='#d95f02',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='#d95f02', alpha=0.95, lw=0.9)
    )

axB.set_title('B. SLC12A5 (KCC2, chr20) Mirror Association Plot', loc='left', fontweight='bold', fontsize=11, pad=14)
axB.set_xlabel('Genomic Position on Chromosome 20 (Mb, GRCh37)', fontweight='bold')
axB.set_ylabel(r'ANX $\leftarrow$   $-\log_{10}(P)$   $\rightarrow$ MDD', fontweight='bold')
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
axB.legend(loc='upper right', frameon=True, fontsize=8)

# ==================== Panel C: DCC (chr18) Multi-Signal Dissection Plot ====================
axC = fig.add_subplot(gs[1, 1])
bp_mb_dcc = df_dcc['BP'] / 1e6
logp_mdd_dcc = -np.log10(np.clip(df_dcc['P_MDD'], 1e-300, 1.0))
logp_anx_dcc = -np.log10(np.clip(df_dcc['P_ANX'], 1e-300, 1.0))

axC.scatter(bp_mb_dcc, logp_mdd_dcc, s=10, color='#7570b3', alpha=0.6, label='MDD (CS=2)')
axC.scatter(bp_mb_dcc, -logp_anx_dcc, s=10, color='#e7298a', alpha=0.6, label='ANX (CS=2)')
axC.axhline(0, color='gray', linewidth=0.8)
axC.axhline(-np.log10(5e-8), color='#7570b3', linestyle=':', linewidth=0.7)
axC.axhline(np.log10(5e-8), color='#e7298a', linestyle=':', linewidth=0.7)
axC.set_ylim(-13, 26)

# Highlight the 4 distinct CS lead SNPs in DCC
cs_dcc_snps = [
    ('rs17411339', 'MDD CS1', '#7570b3', 1),
    ('rs11083029', 'MDD CS2', '#7570b3', 1),
    ('rs7229097', 'ANX CS1', '#e7298a', -1),
    ('rs4245258', 'ANX CS2', '#e7298a', -1)
]
for snp, label, col, direction in cs_dcc_snps:
    row = df_dcc[df_dcc['SNP'] == snp]
    if not row.empty:
        xv = row['BP'].values[0] / 1e6
        pcol = 'P_MDD' if direction == 1 else 'P_ANX'
        yv = direction * (-np.log10(row[pcol].values[0]))
        axC.scatter([xv], [yv], s=55, color=col, edgecolors='black', linewidth=1.2, zorder=6)

axC.set_title('C. DCC (Netrin-1 Receptor, chr18) Allelic Heterogeneity', loc='left', fontweight='bold', fontsize=11, pad=14)
axC.set_xlabel('Genomic Position on Chromosome 18 (Mb, GRCh37)', fontweight='bold')
axC.set_ylabel(r'ANX $\leftarrow$   $-\log_{10}(P)$   $\rightarrow$ MDD', fontweight='bold')
axC.spines['top'].set_visible(False)
axC.spines['right'].set_visible(False)
axC.legend(loc='upper right', frameon=True, fontsize=8)

# Explanatory text box
axC.text(0.04, 0.08, 'SuSiE Result:\nMDD: 2 Distinct Credible Sets\nANX: 2 Distinct Credible Sets\nPP.H3 > 99.4% (Distinct Causal Variants)',
         transform=axC.transAxes, fontsize=7.5, bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0', edgecolor='gray', alpha=0.9))

png_path = os.path.join(out_dir, 'Figure_Aim2_SuSiE_Colocalization_MDD_vs_ANX.png')
pdf_path = os.path.join(out_dir, 'Figure_Aim2_SuSiE_Colocalization_MDD_vs_ANX.pdf')
plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
print(f'Successfully generated {png_path} and {pdf_path} in 300 DPI!')
