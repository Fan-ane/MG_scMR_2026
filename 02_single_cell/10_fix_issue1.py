"""
Fix #1: Split panel D into CD4_NC (left) + CD8_NC (right) control distributions.
Loads all results from saved CSVs — no recomputation.
"""
import pandas as pd
import numpy as np
import os, warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats

data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
out_dir = os.path.join(data_dir, "final_results")
fig_dir = os.path.join(data_dir, "figures_meeting")

C_CD4 = "#1F78B4"
C_CD8 = "#E31A1C"
C_SESN3 = "#D62728"
C_CONTROL = "#AAAAAA"

np.random.seed(440)

# ---- Load all results ----
print("Loading saved results...")
ko_cd4 = pd.read_csv(os.path.join(out_dir, "KO_CD4_NC.csv"))
ko_cd8 = pd.read_csv(os.path.join(out_dir, "KO_CD8_NC.csv"))
gsea_cd4 = pd.read_csv(os.path.join(out_dir, "GSEA_CD4_NC.csv"))
gsea_cd8 = pd.read_csv(os.path.join(out_dir, "GSEA_CD8_NC.csv"))
overlap_df = pd.read_csv(os.path.join(out_dir, "overlap_scDEG.csv"), index_col=0)
controls_df = pd.read_csv(os.path.join(out_dir, "negative_controls_100.csv"))

# Find CD8_NC overlap row (relaxed)
ov_cd8 = overlap_df.loc[[i for i in overlap_df.index if "CD8_NC" in i and "relaxed" in i][0]]
ov_cd4 = overlap_df.loc[[i for i in overlap_df.index if "CD4_NC" in i and "relaxed" in i][0]]

# Extract control metrics for both cell types
for ct in ["CD4_NC", "CD8_NC"]:
    ct_ctrl = controls_df[controls_df["cell_type"] == ct]
    sesn3_row = ct_ctrl[ct_ctrl["gene"] == "SESN3"]
    other_rows = ct_ctrl[ct_ctrl["gene"] != "SESN3"]
    if ct == "CD4_NC":
        sesn3_cd4_top10 = float(sesn3_row["top10_mean_abs_z"].values[0])
        control_cd4_top10 = other_rows["top10_mean_abs_z"].dropna().values.tolist()
    else:
        sesn3_cd8_top10 = float(sesn3_row["top10_mean_abs_z"].values[0])
        control_cd8_top10 = other_rows["top10_mean_abs_z"].dropna().values.tolist()

# Empirical P recalc
n_high_cd4 = sum(1 for v in control_cd4_top10 if v >= sesn3_cd4_top10)
emp_p_cd4 = (1 + n_high_cd4) / (1 + len(control_cd4_top10))
n_high_cd8 = sum(1 for v in control_cd8_top10 if v >= sesn3_cd8_top10)
emp_p_cd8 = (1 + n_high_cd8) / (1 + len(control_cd8_top10))

# KO stats
FREQ_HIGH, FREQ_LOW = 0.50, 0.35
n_cd4_all = (ko_cd4["frequency"] >= FREQ_LOW).sum()
n_cd4_hi = (ko_cd4["frequency"] >= FREQ_HIGH).sum()
n_cd8_all = (ko_cd8["frequency"] >= FREQ_LOW).sum()
n_cd8_hi = (ko_cd8["frequency"] >= FREQ_HIGH).sum()

# ---- FIGURE ----
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7,
    'figure.dpi': 300, 'savefig.dpi': 600,
})

fig = plt.figure(figsize=(21, 13.5))

# ======= A: Workflow (unchanged) =======
ax = fig.add_subplot(2, 3, 1)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
steps = [
    (1.0, 8.2, 8.0, 1.6, "HC CD4_NC / CD8_NC T cells\n(3 donors, 12 refined cell types)", "#E3F2FD"),
    (1.0, 6.3, 8.0, 1.4, "Donor-balanced sampling\n(500 cells/donor × 20 iterations)", "#E3F2FD"),
    (1.0, 4.5, 8.0, 1.3, "scTenifoldKnk SESN3 KO\nZero SESN3 edges in gene-gene\ncorrelation network", "#FFEBEE"),
    (1.0, 2.7, 8.0, 1.3, "Perturbation Z-scores\n(median ± IQR, recurrence frequency)", "#FFEBEE"),
    (1.0, 0.9, 8.0, 1.3, "Validation: GSEA · 100 controls\nOverlap with scDEGs", "#E8F5E9"),
]
for x, y, w, h, text, color in steps:
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='#AAAAAA', linewidth=0.7))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=7, family='monospace')
for i in range(len(steps) - 1):
    _, y1, _, h1, _, _ = steps[i]; _, y2, _, h2, _, _ = steps[i+1]
    ax.annotate('', xy=(5.0, y2 + h2), xytext=(5.0, y1),
                arrowprops=dict(arrowstyle='->', color='#666666', lw=1.5))
ax.set_title("A. Analysis Workflow", fontweight='bold', fontsize=11, loc='left')

# ======= B: CD4_NC KO (unchanged) =======
ax = fig.add_subplot(2, 3, 2)
top = ko_cd4[ko_cd4["gene"] != "SESN3"].head(15)
z_vals, z_err = top["median_z"].values, top["iqr_z"].values / 2
freqs = top["frequency"].values
colors_b = [C_CD4 if f >= FREQ_HIGH else (C_CD4 if f >= FREQ_LOW else '#BBDEFB') for f in freqs]
edge_b = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in freqs]
y_pos = range(len(top))
ax.barh(y_pos, z_vals, xerr=z_err, color=colors_b, alpha=0.9, capsize=2,
        height=0.7, edgecolor=edge_b, linewidth=0.8)
ax.set_yticks(y_pos); ax.set_yticklabels(top["gene"].values, fontsize=7.5)
ax.set_xlabel("Network Perturbation Z-score (± IQR/2)")
ax.axvline(x=0, color='black', lw=0.5); ax.invert_yaxis()
for i, (z, f) in enumerate(zip(z_vals, freqs)):
    if f >= FREQ_HIGH:
        ax.text(z + z_err[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
    elif f >= FREQ_LOW:
        ax.text(z + z_err[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
leg_b = [Line2D([0],[0], marker='s', color='w', markerfacecolor=C_CD4, markersize=10, label='Freq ≥ 0.50'),
         Line2D([0],[0], marker='s', color='w', markerfacecolor=C_CD4, markersize=10, alpha=0.55, label='Freq 0.35–0.50')]
ax.legend(handles=leg_b, fontsize=6.5, loc='lower right', framealpha=0.9)
ax.set_title(f"B. CD4_NC SESN3-KO Perturbed Genes\n({n_cd4_all} recurrent, {n_cd4_hi} highly stable)",
             fontweight='bold', fontsize=10)

# ======= C: CD8_NC KO (unchanged) =======
ax = fig.add_subplot(2, 3, 3)
top8 = ko_cd8[ko_cd8["gene"] != "SESN3"].head(15)
z8, z8e = top8["median_z"].values, top8["iqr_z"].values / 2
f8 = top8["frequency"].values
cb8 = [C_CD8 if f >= FREQ_HIGH else (C_CD8 if f >= FREQ_LOW else '#FFCDD2') for f in f8]
eb8 = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in f8]
yp8 = range(len(top8))
ax.barh(yp8, z8, xerr=z8e, color=cb8, alpha=0.9, capsize=2,
        height=0.7, edgecolor=eb8, linewidth=0.8)
ax.set_yticks(yp8); ax.set_yticklabels(top8["gene"].values, fontsize=7.5)
ax.set_xlabel("Network Perturbation Z-score (± IQR/2)")
ax.axvline(x=0, color='black', lw=0.5); ax.invert_yaxis()
for i, (z, f) in enumerate(zip(z8, f8)):
    if f >= FREQ_HIGH:
        ax.text(z + z8e[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
    elif f >= FREQ_LOW:
        ax.text(z + z8e[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
ax.set_title(f"C. CD8_NC SESN3-KO Perturbed Genes\n({n_cd8_all} recurrent, {n_cd8_hi} highly stable)",
             fontweight='bold', fontsize=10)

# ======= D: Negative Controls — SPLIT CD4_NC (left) + CD8_NC (right) =======
ax = fig.add_subplot(2, 3, 4)
ax.axis('off')

# Create two inset axes within panel D
# Left: CD4_NC
ax_d1 = ax.inset_axes([0.0, 0.15, 0.48, 0.80])
ax_d1.hist(control_cd4_top10, bins=22, color=C_CD4, edgecolor='#333333',
           alpha=0.25, linewidth=0.5)
ax_d1.hist(control_cd4_top10, bins=22, color=C_CD4, edgecolor='#333333',
           alpha=0.6, histtype='step', linewidth=1.5)
ax_d1.axvline(x=sesn3_cd4_top10, color=C_SESN3, linewidth=2.5, linestyle='--',
              label=f'SESN3 = {sesn3_cd4_top10:.2f}')
p_str_d1 = f"P = {emp_p_cd4:.4f}"
ax_d1.text(0.98, 0.95, f"CD4_NC\nEmpirical {p_str_d1}\n"
           f"{n_high_cd4}/{len(control_cd4_top10)} controls ≥ SESN3",
           transform=ax_d1.transAxes, ha='right', va='top', fontsize=7.5,
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=C_CD4, alpha=0.9))
ax_d1.set_xlabel("Mean |Z| Top-10 Genes", fontsize=7.5)
ax_d1.set_ylabel("Frequency", fontsize=7.5)
ax_d1.set_title("CD4_NC", fontweight='bold', fontsize=9, color=C_CD4)
ax_d1.legend(fontsize=6.5, framealpha=0.9)

# Right: CD8_NC
ax_d2 = ax.inset_axes([0.52, 0.15, 0.48, 0.80])
ax_d2.hist(control_cd8_top10, bins=22, color=C_CD8, edgecolor='#333333',
           alpha=0.25, linewidth=0.5)
ax_d2.hist(control_cd8_top10, bins=22, color=C_CD8, edgecolor='#333333',
           alpha=0.6, histtype='step', linewidth=1.5)
ax_d2.axvline(x=sesn3_cd8_top10, color=C_SESN3, linewidth=2.5, linestyle='--',
              label=f'SESN3 = {sesn3_cd8_top10:.2f}')
p_str_d2 = f"P = {emp_p_cd8:.4f}"
ax_d2.text(0.98, 0.95, f"CD8_NC\nEmpirical {p_str_d2}\n"
           f"{n_high_cd8}/{len(control_cd8_top10)} controls ≥ SESN3",
           transform=ax_d2.transAxes, ha='right', va='top', fontsize=7.5,
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=C_CD8, alpha=0.9))
ax_d2.set_xlabel("Mean |Z| Top-10 Genes", fontsize=7.5)
ax_d2.set_ylabel("Frequency", fontsize=7.5)
ax_d2.set_title("CD8_NC", fontweight='bold', fontsize=9, color=C_CD8)
ax_d2.legend(fontsize=6.5, framealpha=0.9)

# Panel title
ax.text(0.5, 1.02, "D. Negative Control Distribution\nTop-Target Perturbation Magnitude (n = 100 matched controls)",
        transform=ax.transAxes, ha='center', va='bottom', fontweight='bold', fontsize=10)

# ======= E: GSEA (unchanged) =======
ax = fig.add_subplot(2, 3, 5)
gsea_comb = []
for df, ct in [(gsea_cd4, "CD4_NC"), (gsea_cd8, "CD8_NC")]:
    if len(df) > 0:
        dc = df.copy(); dc["cell_type"] = ct; gsea_comb.append(dc)
gsea_all = pd.concat(gsea_comb).sort_values("p_value")
gsea_plot = gsea_all.drop_duplicates(subset="pathway").head(16)
neg_log_p = -np.log10(gsea_plot["p_value"].values.clip(min=1e-300))
colors_g = [C_CD4 if ct == "CD4_NC" else C_CD8 for ct in gsea_plot["cell_type"]]
ypg = range(len(gsea_plot))
ax.barh(ypg, neg_log_p, color=colors_g, alpha=0.85, height=0.7,
        edgecolor='white', linewidth=0.3)
ax.set_yticks(ypg)
labels_g = [f"{r['pathway']} [{r['cell_type'].replace('_',' ')}]" for _, r in gsea_plot.iterrows()]
ax.set_yticklabels(labels_g, fontsize=6.5); ax.invert_yaxis()
ax.set_xlabel("-log10(P-value)")
ax.axvline(x=-np.log10(0.05), color='grey', linestyle='--', alpha=0.6, lw=0.8)
ax.axvline(x=-np.log10(0.01), color='grey', linestyle=':', alpha=0.4, lw=0.5)
for i, (_, r) in enumerate(gsea_plot.iterrows()):
    if r["FDR"] < 0.10:
        ax.text(neg_log_p[i] + 0.1, i, '★', va='center', fontsize=9, color='#D4A017')
leg_g = [Patch(facecolor=C_CD4, alpha=0.85, label='CD4_NC'),
         Patch(facecolor=C_CD8, alpha=0.85, label='CD8_NC'),
         Line2D([0],[0], marker='*', color='w', markerfacecolor='#D4A017', markersize=10, label='FDR < 0.10')]
ax.legend(handles=leg_g, fontsize=7, loc='lower right', framealpha=0.9)
ax.set_title("E. Ranked Gene-Set Enrichment Analysis\nof SESN3 Perturbation Scores",
             fontweight='bold', fontsize=10)

# ======= F: Overlap Forest Plot (unchanged) =======
ax = fig.add_subplot(2, 3, 6)
overs = [
    {"label": "CD4_NC", "OR": float(ov_cd4["OR"]), "OR_CI_low": float(ov_cd4["OR_CI_low"]),
     "OR_CI_high": float(ov_cd4["OR_CI_high"]), "p_fisher": float(ov_cd4["p_fisher"]),
     "n_overlap": int(ov_cd4["n_overlap"]), "n_ko": int(ov_cd4["n_ko"]), "color": C_CD4},
    {"label": "CD8_NC", "OR": float(ov_cd8["OR"]), "OR_CI_low": float(ov_cd8["OR_CI_low"]),
     "OR_CI_high": float(ov_cd8["OR_CI_high"]), "p_fisher": float(ov_cd8["p_fisher"]),
     "n_overlap": int(ov_cd8["n_overlap"]), "n_ko": int(ov_cd8["n_ko"]), "color": C_CD8},
]
for i, ov in enumerate(overs):
    or_val, ci_l, ci_h = ov["OR"], ov["OR_CI_low"], ov["OR_CI_high"]
    if np.isnan(or_val) or or_val <= 0:
        or_disp, ci_valid, or_label = 0.5, False, "N/A"
    elif np.isinf(or_val) or or_val > 100:
        or_disp, ci_valid, or_label = 50, not np.isnan(ci_l), "∞"
    else:
        or_disp, ci_valid, or_label = or_val, not np.isnan(ci_l), f"{or_val:.1f}"
    y = 1 - i; mk = 'D' if i == 0 else 's'
    ax.plot(or_disp, y, marker=mk, color=ov["color"], markersize=12,
            markeredgewidth=1.5, markeredgecolor='white', zorder=5)
    if ci_valid and not np.isinf(ci_l) and not np.isinf(ci_h):
        ci_lc = max(0.05, min(ci_l, 50)); ci_hc = min(ci_h, 50)
        ax.plot([ci_lc, ci_hc], [y, y], '-', color=ov["color"], linewidth=3, zorder=4)
    p_str = f"{ov['p_fisher']:.1e}" if ov['p_fisher'] < 0.01 else f"{ov['p_fisher']:.3f}"
    ax.text(55, y, f"OR={or_label}  P={p_str}\nk={ov['n_overlap']}/{ov['n_ko']} (KO∩DEG)",
            va='center', fontsize=7, family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=ov["color"], alpha=0.85))
ax.set_ylim(-0.8, 1.8); ax.set_yticks([0, 1])
ax.set_yticklabels(["CD4_NC", "CD8_NC"], fontsize=10, fontweight='bold')
ax.axvline(x=1, color='black', linestyle='-', lw=0.8, zorder=1)
ax.set_xlabel("Odds Ratio (95% CI)")
ax.set_xscale('symlog', linthresh=1, linscale=0.5); ax.set_xlim(0.03, 75)
ax.set_title("F. SESN3-KO Perturbed ∩ scDEGs\n(Fisher exact, |log2FC|>0.5, P<0.01)",
             fontweight='bold', fontsize=10)

# ---- Final ----
plt.suptitle("In-silico SESN3 Perturbation Predicts Altered Regulatory Networks\n"
             "in MG-Associated Naïve/Central Memory T Cells",
             fontweight='bold', fontsize=14, y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])

fig_path = os.path.join(fig_dir, "Fig8_SESN3_virtual_KO_final")
for fmt, dpi_val, kw in [("png", 600, {}), ("pdf", None, {}),
                          ("tiff", 600, {"pil_kwargs": {"compression": "tiff_lzw"}})]:
    plt.savefig(f"{fig_path}.{fmt}", dpi=dpi_val, bbox_inches='tight',
                facecolor='white', edgecolor='none', **kw)
plt.close()
print(f"[OK] Figure saved: {fig_path}.{{png, pdf, tiff}}")
print("\nChanges in this version:")
print("  #1: Panel D → split CD4_NC (left, blue) + CD8_NC (right, red)")
print("  #3: Title → 'In-silico ... Predicts' (removed 'Reveals')")
print("  #4: X-axis → 'Network Perturbation Z-score' (panels B, C)")
print("  #8: E title → 'Ranked Gene-Set Enrichment Analysis of SESN3 Perturbation Scores'")
