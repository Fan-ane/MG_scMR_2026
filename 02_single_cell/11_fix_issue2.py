"""
Fix #2: Replace F panel overlap forest plot →
        Rank-based disease signature enrichment (GSEA framework)
        Tests: Are MG scDEGs enriched at top of SESN3 perturbation ranking?
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

# ---- Load all saved data ----
print("Loading...")
ko_cd4 = pd.read_csv(os.path.join(out_dir, "KO_CD4_NC.csv"))
ko_cd8 = pd.read_csv(os.path.join(out_dir, "KO_CD8_NC.csv"))
gsea_cd4 = pd.read_csv(os.path.join(out_dir, "GSEA_CD4_NC.csv"))
gsea_cd8 = pd.read_csv(os.path.join(out_dir, "GSEA_CD8_NC.csv"))
controls_df = pd.read_csv(os.path.join(out_dir, "negative_controls_100.csv"))
scdeg_cd4 = pd.read_csv(os.path.join(out_dir, "scDEG_CD4_NC.csv"))
scdeg_cd8 = pd.read_csv(os.path.join(out_dir, "scDEG_CD8_NC.csv"))

# ---- Extract control data ----
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

n_high_cd4 = sum(1 for v in control_cd4_top10 if v >= sesn3_cd4_top10)
emp_p_cd4 = (1 + n_high_cd4) / (1 + len(control_cd4_top10))
n_high_cd8 = sum(1 for v in control_cd8_top10 if v >= sesn3_cd8_top10)
emp_p_cd8 = (1 + n_high_cd8) / (1 + len(control_cd8_top10))

# ============================================================
# NEW: Disease Signature GSEA (replaces overlap forest plot)
# ============================================================
print("\nComputing disease signature enrichment...")

def disease_gsea(ko_summary, scdeg, label, n_perm=2000):
    """
    GSEA: rank genes by KO perturbation |Z| (high→low).
    Test if MG DEGs (|log2FC|>0.5, P<0.01) are enriched at top of KO ranking.

    Returns NES, P-value, leading edge genes, and running enrichment data for plotting.
    """
    # Build ranked list by |median_z|
    z_ranked = ko_summary.set_index("gene")["median_z"].dropna()
    z_ranked = z_ranked.abs().sort_values(ascending=False)

    # Define disease gene set: scDEGs passing threshold
    deg_genes = set(scdeg[scdeg["pass_relaxed"]]["names"].values)

    # Filter to genes in both
    common = [g for g in z_ranked.index if g in deg_genes or True]  # all ranked genes
    bg = set(z_ranked.index)
    pw_in_bg = [g for g in deg_genes if g in bg]

    if len(pw_in_bg) < 5:
        return {"label": label, "n_deg": len(pw_in_bg), "NES": 0, "p_value": 1.0,
                "leading_edge": "", "ES": 0, "running": None, "n_total": 0}

    n_total = len(z_ranked)
    in_pw = np.array([g in deg_genes for g in z_ranked.index], dtype=int)
    n_pw = in_pw.sum()

    # Weighted KS (weight = |Z|)
    hit_w = np.abs(z_ranked.values)
    hit_w[in_pw == 0] = 0
    h_sum = hit_w.sum()
    miss_w = np.ones(n_total)
    miss_w[in_pw == 1] = 0
    m_sum = miss_w.sum()

    running = np.cumsum(hit_w / h_sum - miss_w / m_sum)
    es = running.max() if abs(running.max()) > abs(running.min()) else running.min()

    # Permutation
    perm_es = []
    for _ in range(n_perm):
        pi = np.random.permutation(n_total)
        p_in = in_pw[pi]
        p_hit = hit_w[pi]; p_hit[p_in == 0] = 0
        ph_sum = p_hit.sum()
        if ph_sum == 0: continue
        p_miss = miss_w[pi]; p_miss[p_in == 1] = 0
        pm_sum = p_miss.sum()
        if pm_sum == 0: continue
        p_run = np.cumsum(p_hit / ph_sum - p_miss / pm_sum)
        perm_es.append(p_run.max() if abs(p_run.max()) > abs(p_run.min()) else p_run.min())

    ss = [e for e in perm_es if e * es > 0]
    mean_perm = np.mean(np.abs(ss)) if ss else 1.0
    nes = es / mean_perm if mean_perm > 0 else es

    n_valid = len(perm_es)
    p_val = ((sum(1 for e in perm_es if e >= es) + 1) / (n_valid + 1) if es > 0
             else (sum(1 for e in perm_es if e <= es) + 1) / (n_valid + 1))

    # Leading edge: genes before ES peak
    es_peak_idx = np.argmax(running) if es > 0 else np.argmin(running)
    le_genes = [z_ranked.index[i] for i in range(es_peak_idx + 1) if in_pw[i] == 1][:10]

    return {
        "label": label, "n_deg": int(n_pw), "n_total": n_total,
        "ES": es, "NES": nes, "p_value": p_val,
        "running": running, "es_peak_idx": int(es_peak_idx),
        "leading_edge": ", ".join(le_genes),
    }

dgsea_cd4 = disease_gsea(ko_cd4, scdeg_cd4, "CD4_NC")
dgsea_cd8 = disease_gsea(ko_cd8, scdeg_cd8, "CD8_NC")

for dg in [dgsea_cd4, dgsea_cd8]:
    sig = "***" if dg["p_value"] < 0.01 else ("**" if dg["p_value"] < 0.05 else "")
    print(f"  {dg['label']}: NES={dg['NES']:+.2f}, P={dg['p_value']:.4f} {sig}, "
          f"DEGs tested={dg['n_deg']}, LE={dg['leading_edge'][:60]}")

# ---- KO stats ----
FREQ_HIGH, FREQ_LOW = 0.50, 0.35
n_cd4_all = (ko_cd4["frequency"] >= FREQ_LOW).sum()
n_cd4_hi = (ko_cd4["frequency"] >= FREQ_HIGH).sum()
n_cd8_all = (ko_cd8["frequency"] >= FREQ_LOW).sum()
n_cd8_hi = (ko_cd8["frequency"] >= FREQ_HIGH).sum()

# ============================================================
# FIGURE
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'legend.fontsize': 7,
    'figure.dpi': 300, 'savefig.dpi': 600,
})

fig = plt.figure(figsize=(21, 13.5))

# ======= A: Workflow =======
ax = fig.add_subplot(2, 3, 1)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
steps = [
    (1.0, 8.2, 8.0, 1.6, "HC CD4_NC / CD8_NC T cells\n(3 donors, 12 refined cell types)", "#E3F2FD"),
    (1.0, 6.3, 8.0, 1.4, "Donor-balanced sampling\n(500 cells/donor × 20 iterations)", "#E3F2FD"),
    (1.0, 4.5, 8.0, 1.3, "scTenifoldKnk SESN3 KO\nZero SESN3 edges in gene-gene\ncorrelation network", "#FFEBEE"),
    (1.0, 2.7, 8.0, 1.3, "Network perturbation Z-scores\n(median ± IQR, recurrence frequency)", "#FFEBEE"),
    (1.0, 0.9, 8.0, 1.3, "Validation: GSEA · 100 controls\nDisease signature enrichment", "#E8F5E9"),
]
for x, y, w, h, text, color in steps:
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='#AAAAAA', linewidth=0.7))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=7, family='monospace')
for i in range(len(steps) - 1):
    _, y1, _, h1, _, _ = steps[i]; _, y2, _, h2, _, _ = steps[i+1]
    ax.annotate('', xy=(5.0, y2 + h2), xytext=(5.0, y1),
                arrowprops=dict(arrowstyle='->', color='#666666', lw=1.5))
ax.set_title("A. Analysis Workflow", fontweight='bold', fontsize=11, loc='left')
# Footnote: why HC only
ax.text(5.0, -0.6, "HC T-cell networks used as baseline to predict SESN3 perturbation consequences.",
        ha='center', va='top', fontsize=7, fontstyle='italic', color='#666666',
        transform=ax.transData)

# ======= B: CD4_NC KO =======
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

# ======= C: CD8_NC KO =======
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

# ======= D: Negative Controls =======
ax = fig.add_subplot(2, 3, 4)
ax.axis('off')

ax_d1 = ax.inset_axes([0.0, 0.15, 0.48, 0.80])
ax_d1.hist(control_cd4_top10, bins=22, color=C_CD4, edgecolor='#333333',
           alpha=0.25, linewidth=0.5)
ax_d1.hist(control_cd4_top10, bins=22, color=C_CD4, edgecolor='#333333',
           alpha=0.6, histtype='step', linewidth=1.5)
ax_d1.axvline(x=sesn3_cd4_top10, color=C_SESN3, linewidth=2.5, linestyle='--',
              label=f'SESN3 = {sesn3_cd4_top10:.2f}')
ax_d1.text(0.98, 0.95, f"CD4_NC\nEmpirical P = {emp_p_cd4:.4f}\n"
           f"{n_high_cd4}/{len(control_cd4_top10)} controls ≥ SESN3",
           transform=ax_d1.transAxes, ha='right', va='top', fontsize=7.5,
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=C_CD4, alpha=0.9))
ax_d1.set_xlabel("Mean |Z| Top-10 Genes", fontsize=7.5)
ax_d1.set_ylabel("Frequency", fontsize=7.5)
ax_d1.set_title("CD4_NC", fontweight='bold', fontsize=9, color=C_CD4)
ax_d1.legend(fontsize=6.5, framealpha=0.9)

ax_d2 = ax.inset_axes([0.52, 0.15, 0.48, 0.80])
ax_d2.hist(control_cd8_top10, bins=22, color=C_CD8, edgecolor='#333333',
           alpha=0.25, linewidth=0.5)
ax_d2.hist(control_cd8_top10, bins=22, color=C_CD8, edgecolor='#333333',
           alpha=0.6, histtype='step', linewidth=1.5)
ax_d2.axvline(x=sesn3_cd8_top10, color=C_SESN3, linewidth=2.5, linestyle='--',
              label=f'SESN3 = {sesn3_cd8_top10:.2f}')
ax_d2.text(0.98, 0.95, f"CD8_NC\nEmpirical P = {emp_p_cd8:.4f}\n"
           f"{n_high_cd8}/{len(control_cd8_top10)} controls ≥ SESN3",
           transform=ax_d2.transAxes, ha='right', va='top', fontsize=7.5,
           bbox=dict(boxstyle='round', facecolor='white', edgecolor=C_CD8, alpha=0.9))
ax_d2.set_xlabel("Mean |Z| Top-10 Genes", fontsize=7.5)
ax_d2.set_title("CD8_NC", fontweight='bold', fontsize=9, color=C_CD8)
ax_d2.legend(fontsize=6.5, framealpha=0.9)

ax.text(0.5, 1.02, "D. Negative Control Distribution\nTop-Target Perturbation Magnitude (n = 100 matched controls)",
        transform=ax.transAxes, ha='center', va='bottom', fontweight='bold', fontsize=10)

# ======= E: GSEA =======
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
ax.set_title("E. GSEA Pathway Enrichment\n(Ranked by SESN3 Perturbation Scores)",
             fontweight='bold', fontsize=10)

# ======= F: Disease Signature Enrichment (NEW) =======
ax = fig.add_subplot(2, 3, 6)

# Create two sub-panels: top = CD4_NC running enrichment, bottom = CD8_NC
# Use inset_axes for the running ES plots
for idx, (dg, ct_label, clr) in enumerate([
    (dgsea_cd4, "CD4_NC", C_CD4),
    (dgsea_cd8, "CD8_NC", C_CD8),
]):
    # Position: top 55% for CD4, bottom 45% for CD8
    if idx == 0:
        sub_ax = ax.inset_axes([0.08, 0.52, 0.90, 0.43])
    else:
        sub_ax = ax.inset_axes([0.08, 0.04, 0.90, 0.43])

    if dg["running"] is not None:
        running = dg["running"]
        n_total = len(running)
        x_vals = np.arange(n_total) / n_total * 100  # percentile

        # Running ES line
        sub_ax.plot(x_vals, running, color=clr, linewidth=1.2)
        sub_ax.fill_between(x_vals, 0, running, color=clr, alpha=0.15)

        # Mark ES peak
        peak = dg["es_peak_idx"]
        peak_x = peak / n_total * 100
        sub_ax.axvline(x=peak_x, color=clr, linestyle=':', linewidth=0.8, alpha=0.7)

        # Annotate
        p_str = f"P = {dg['p_value']:.4f}" if dg['p_value'] >= 0.001 else f"P = {dg['p_value']:.2e}"
        sub_ax.text(0.02, 0.92, f"{ct_label}\nNES = {dg['NES']:+.2f}  {p_str}\n"
                    f"DEGs = {dg['n_deg']}  Leading: {dg['leading_edge'][:40]}...",
                    transform=sub_ax.transAxes, fontsize=6.5, va='top',
                    bbox=dict(boxstyle='round', facecolor='white', edgecolor=clr, alpha=0.85))

        # Horizontal line at 0
        sub_ax.axhline(y=0, color='grey', linestyle='-', linewidth=0.3)

        # Labels
        sub_ax.set_ylabel("Running ES", fontsize=7)
        if idx == 1:
            sub_ax.set_xlabel("Genes Ranked by SESN3 Perturbation |Z| (percentile)", fontsize=7)
        else:
            sub_ax.set_xticklabels([])
    else:
        sub_ax.text(0.5, 0.5, f"{ct_label}: No enrichment", ha='center', va='center',
                    fontsize=8, transform=sub_ax.transAxes)

# Panel title
ax.axis('off')
ax.text(0.5, 1.02,
        "F. Disease Signature Enrichment\n"
        "MG scDEGs Enriched Against SESN3 Perturbation Rank\n"
        "(Weighted KS test, 2000 permutations)",
        transform=ax.transAxes, ha='center', va='bottom', fontweight='bold', fontsize=10)

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

print(f"\n[OK] Figure saved: {fig_path}.{{png, pdf, tiff}}")
print("\nPanel F changes:")
print("  → Replaced overlap forest plot with GSEA running enrichment curves")
print("  → CD4_NC and CD8_NC disease signature enrichment shown as ES traces")
print("  → Weighted KS test, genes ranked by SESN3 perturbation |Z|")
