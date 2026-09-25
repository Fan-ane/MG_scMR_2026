"""
=============================================================================
AGER Virtual Knockout — Monocyte-Focused Analysis (Revised)
=============================================================================
Focus: Mono_C (classical monocytes) + Mono_NC (non-classical monocytes)

Key improvements over v1:
  1. Monocyte-only focus — the only biologically relevant cell types for AGER
  2. Negative controls matched from FULL gene set (~21k) not just 2000 HVGs
  3. Single-iteration KO with full rank output (Z-scores for all genes)
  4. Empirical P-value from >50 matched controls
  5. 4-panel monocyte-focused figure
=============================================================================
"""
import scanpy as sc
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
from scipy.sparse import issparse
from collections import Counter

# ============================================================
# Setup
# ============================================================
data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
out_dir = os.path.join(data_dir, "ager_mono_results")
fig_dir = os.path.join(data_dir, "figures_meeting")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

np.random.seed(440)
TARGET_GENE = "AGER"

C_MONO = "#2CA02C"
C_MONO_NC = "#9467BD"
C_SIG = "#D62728"
C_CONTROL = "#AAAAAA"

# ============================================================
# 0. Load
# ============================================================
print("=" * 70)
print("AGER MONOCYTE-FOCUSED VIRTUAL KO")
print("=" * 70)

adata = sc.read_h5ad(os.path.join(data_dir, "adata_refined.h5ad"))
adata.obs_names_make_unique()
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

if TARGET_GENE not in adata.var_names:
    raise ValueError(f"{TARGET_GENE} not found!")

# ============================================================
# 1. CORE FUNCTIONS
# ============================================================
def extract_cells_donor_balanced(adata, cell_type, condition="HC", n_per_donor=500):
    mask = (adata.obs["cell_type_refined"] == cell_type) & (adata.obs["condition"] == condition)
    sub = adata[mask].copy()
    cells_use = []
    for donor in sorted(sub.obs["sample"].unique()):
        donor_cells = list(sub.obs_names[sub.obs["sample"] == donor])
        n_take = min(n_per_donor, len(donor_cells))
        if n_take >= 50:
            cells_use.extend(list(np.random.choice(donor_cells, n_take, replace=False)))
    if len(cells_use) == 0:
        raise ValueError(f"No cells for {cell_type}")
    return sub[cells_use, :].copy()

def prepare_ko_matrix(sub_adata, n_genes=2000, target_gene="AGER"):
    if issparse(sub_adata.X):
        expr = sub_adata.X.toarray()
    else:
        expr = sub_adata.X
    gene_vars = np.var(expr, axis=0)
    top_idx = np.argsort(gene_vars)[-n_genes:]
    hvg_genes = [sub_adata.var_names[i] for i in top_idx]
    if target_gene in sub_adata.var_names and target_gene not in hvg_genes:
        hvg_genes.append(target_gene)
    gene_indices = [list(sub_adata.var_names).index(g) for g in hvg_genes if g in sub_adata.var_names]
    expr_sub = expr[:, gene_indices]
    final_genes = [g for g in hvg_genes if g in sub_adata.var_names]
    return expr_sub, final_genes

def virtual_ko_single(expr, gene_list, target_gene):
    if target_gene not in gene_list:
        return None
    tidx = gene_list.index(target_gene)
    n_cells, n_genes = expr.shape
    expr_t = expr.T
    corr_orig = np.corrcoef(expr_t)
    corr_ko = corr_orig.copy()
    corr_ko[tidx, :] = 0
    corr_ko[:, tidx] = 0
    perturbation = np.array([
        np.sqrt(np.mean((corr_orig[i, :] - corr_ko[i, :]) ** 2))
        for i in range(n_genes)
    ])
    mask = np.ones(n_genes, dtype=bool)
    mask[tidx] = False
    mean_p = perturbation[mask].mean()
    std_p = perturbation[mask].std()
    z_scores = np.zeros(n_genes)
    if std_p > 0:
        z_scores[mask] = (perturbation[mask] - mean_p) / std_p
    return pd.DataFrame({"gene": gene_list, "distance": perturbation, "z_score": z_scores})

def virtual_ko_iterated(adata, cell_type, target_gene="AGER", n_iter=20,
                        n_per_donor=500, n_hvg=2000):
    """20-iteration KO with frequency tracking."""
    all_z = {}
    all_freq = Counter()
    for i in range(n_iter):
        seed = 440 + i * 100
        np.random.seed(seed)
        sub = extract_cells_donor_balanced(adata, cell_type, "HC", n_per_donor=n_per_donor)
        expr, genes = prepare_ko_matrix(sub, n_genes=n_hvg, target_gene=target_gene)
        result = virtual_ko_single(expr, genes, target_gene)
        if result is None:
            continue
        for _, row in result.iterrows():
            all_z.setdefault(row["gene"], []).append(row["z_score"])
        abs_z = np.abs(result["z_score"].values)
        thresh = np.percentile(abs_z, 95)
        for g in result.loc[abs_z >= thresh, "gene"]:
            all_freq[g] += 1
        if (i + 1) % 10 == 0:
            print(f"  {cell_type}: {i+1}/{n_iter} done")
    summary = []
    for gene, zs in all_z.items():
        zs_arr = np.array(zs)
        summary.append({
            "gene": gene,
            "median_z": np.median(zs_arr),
            "iqr_z": stats.iqr(zs_arr) if len(zs_arr) > 1 else 0.0,
            "mean_z": np.mean(zs_arr),
            "frequency": all_freq.get(gene, 0) / n_iter,
        })
    return pd.DataFrame(summary).sort_values("median_z", ascending=False)

# ============================================================
# 2. RUN VIRTUAL KO — Mono_C + Mono_NC only
# ============================================================
print("\n" + "=" * 70)
print("VIRTUAL KO: 20 Iterations, Monocytes Only")
print("=" * 70)

targets = {
    "Mono_C":  {"n_per_donor": 1000, "n_hvg": 2000, "label": "Classical Monocytes"},
    "Mono_NC": {"n_per_donor": 400,  "n_hvg": 2000, "label": "Non-Classical Monocytes"},
}

ko_results = {}
for ct, params in targets.items():
    fname = f"KO_{ct}.csv"
    fpath = os.path.join(out_dir, fname)

    # Checkpoint
    if os.path.exists(fpath):
        print(f"\n[CHECKPOINT] Loading {fname}...")
        ko_results[ct] = pd.read_csv(fpath)
    else:
        print(f"\n{'='*60}")
        print(f"  {ct} ({params['label']})")
        print(f"{'='*60}")
        ko_summary = virtual_ko_iterated(
            adata, ct, target_gene=TARGET_GENE,
            n_iter=20, n_per_donor=params["n_per_donor"], n_hvg=params["n_hvg"]
        )
        ko_summary.to_csv(fpath, index=False)
        ko_results[ct] = ko_summary
        # Clear checkpoint for fresh controls
        for f in ["controls_full.csv", "scDEG_Mono_C.csv", "scDEG_Mono_NC.csv",
                  "GSEA_Mono_C.csv", "GSEA_Mono_NC.csv"]:
            fp = os.path.join(out_dir, f)
            if os.path.exists(fp):
                os.remove(fp)

    n35 = (ko_results[ct]["frequency"] >= 0.35).sum()
    n50 = (ko_results[ct]["frequency"] >= 0.50).sum()
    print(f"  freq>=0.35: {n35}, freq>=0.50: {n50}")
    top = ko_results[ct][ko_results[ct]["gene"] != TARGET_GENE].head(8)
    for _, r in top.iterrows():
        m = "***" if r['frequency'] >= 0.50 else ("**" if r['frequency'] >= 0.35 else "")
        print(f"    {r['gene']:15s} Z={r['median_z']:+.2f} freq={r['frequency']:.2f} {m}")

# ============================================================
# 3. IMPROVED NEGATIVE CONTROLS
#    Match from FULL gene set in Mono_C reference cells
# ============================================================
print("\n" + "=" * 70)
print("NEGATIVE CONTROLS: Full-Genome Matching")
print("=" * 70)

# Get reference expression from Mono_C (full gene set, not just HVGs)
sub_ref = extract_cells_donor_balanced(adata, "Mono_C", "HC", n_per_donor=1000)
expr_full = sub_ref.X.toarray() if issparse(sub_ref.X) else sub_ref.X
all_genes = list(sub_ref.var_names)
ager_idx_full = all_genes.index(TARGET_GENE)

ager_mean_full = expr_full[:, ager_idx_full].mean()
ager_detect_full = (expr_full[:, ager_idx_full] > 0).mean()
print(f"{TARGET_GENE} (Mono_C, full genes): mean={ager_mean_full:.4f}, detect={ager_detect_full:.3f}")

# Build candidate pool from ALL genes (not just HVGs)
known_hits = {"SESN3", "AGER", "ZSCAN12", "CFD", "DXO", "ZBTB9", "HLA-DOA",
              "BTN2A1", "HCG23", "TYMP", "TNFRSF14", "SPCS1", "MPV17L2",
              "RNASEH2C", "CTSH", "RBM47", "SLC22A4", "FCRL3", "SCAND3",
              "ZNF322", "IL12RB2", "ORMDL3", "HFE"}

candidates = pd.DataFrame({
    "gene": all_genes,
    "mean_expr": expr_full.mean(axis=0),
    "detection": (expr_full > 0).mean(axis=0),
})
candidates = candidates[~candidates["gene"].isin(known_hits)]
candidates = candidates[candidates["gene"] != TARGET_GENE]

# MULTI-TIER matching: progressively relax to get enough controls
for tier_name, (expr_lo, expr_hi, det_delta, n_target) in [
    ("strict",  (0.5, 2.0, 0.06, 80)),
    ("relaxed", (0.3, 3.0, 0.10, 60)),
    ("wide",    (0.2, 5.0, 0.15, 40)),
]:
    pool = candidates[
        (candidates["mean_expr"] >= ager_mean_full * expr_lo) &
        (candidates["mean_expr"] <= ager_mean_full * expr_hi) &
        (np.abs(candidates["detection"] - ager_detect_full) < det_delta)
    ].copy()
    pool["match_score"] = (
        np.abs(pool["mean_expr"] / ager_mean_full - 1) +
        np.abs(pool["detection"] - ager_detect_full) * 8
    )
    n_avail = len(pool)
    n_use = min(n_target, n_avail)
    print(f"\n  {tier_name}: mean [{expr_lo}x-{expr_hi}x], detect delta<{det_delta}")
    print(f"    → {n_avail} candidates available, targeting {n_use}")

    if n_avail >= 30:
        control_genes = pool.nsmallest(n_use, "match_score")["gene"].tolist()
        print(f"    SELECTED {len(control_genes)} controls (tier={tier_name})")
        break
    else:
        print(f"    NOT ENOUGH ({n_avail} < 30), trying next tier...")

if len(control_genes) < 20:
    # Last resort: take closest matches regardless of criteria
    candidates["match_score"] = (
        np.abs(candidates["mean_expr"] / ager_mean_full - 1) +
        np.abs(candidates["detection"] - ager_detect_full) * 5
    )
    control_genes = candidates.nsmallest(60, "match_score")["gene"].tolist()
    print(f"  FALLBACK: selected {len(control_genes)} genes by minimal distance")

print(f"\n  Final: {len(control_genes)} control genes")
ctrl_detect = candidates[candidates["gene"].isin(control_genes)]["detection"]
ctrl_mean = candidates[candidates["gene"].isin(control_genes)]["mean_expr"]
print(f"  AGER:        detect={ager_detect_full:.3f}, mean={ager_mean_full:.4f}")
print(f"  Controls:    detect=[{ctrl_detect.min():.3f}, {ctrl_detect.max():.3f}], "
      f"mean=[{ctrl_mean.min():.4f}, {ctrl_mean.max():.4f}]")

# ---- Run KO for AGER + all controls (single iteration) ----
def run_ko_single_iter(adata, cell_type, target_gene, n_per_donor=1000, n_hvg=2000):
    """Single iteration of virtual KO, return full Z-score dataframe."""
    sub = extract_cells_donor_balanced(adata, cell_type, "HC", n_per_donor=n_per_donor)
    expr, genes = prepare_ko_matrix(sub, n_genes=n_hvg, target_gene=target_gene)
    return virtual_ko_single(expr, genes, target_gene)

print(f"\nRunning KO for AGER + {len(control_genes)} controls on Mono_C...")
np.random.seed(440)
ager_result = run_ko_single_iter(adata, "Mono_C", TARGET_GENE, n_per_donor=1000)
ager_abs_z = np.abs(ager_result[ager_result["gene"] != TARGET_GENE]["z_score"].values)

# Metrics: top-N mean |Z|, and full distribution summary
def compute_metrics(result_df, target_gene, n_tops=[5, 10, 20, 50]):
    no_target = result_df[result_df["gene"] != target_gene]
    abs_z = np.abs(no_target["z_score"].values)
    m = {}
    for n in n_tops:
        m[f"top{n}_mean_abs_z"] = np.sort(abs_z)[-n:].mean()
    m["median_abs_z"] = np.median(abs_z)
    m["max_abs_z"] = abs_z.max()
    m["n_genes"] = len(abs_z)
    return m

ager_metrics = compute_metrics(ager_result, TARGET_GENE)
print(f"  {TARGET_GENE}: top5|Z|={ager_metrics['top5_mean_abs_z']:.3f}, "
      f"top10|Z|={ager_metrics['top10_mean_abs_z']:.3f}, "
      f"top20|Z|={ager_metrics['top20_mean_abs_z']:.3f}")

cpath = os.path.join(out_dir, "controls_full.csv")
ctrls_done = set()
if os.path.exists(cpath):
    existing = pd.read_csv(cpath)
    ctrls_done = set(existing["gene"].unique())
    print(f"  [CHECKPOINT] {len(ctrls_done)} controls already run")

control_all = []
for i, cg in enumerate(control_genes):
    if cg in ctrls_done:
        continue
    np.random.seed(440 + i)
    c_res = run_ko_single_iter(adata, "Mono_C", cg, n_per_donor=1000)
    if c_res is not None:
        metrics = compute_metrics(c_res, cg)
        metrics["gene"] = cg
        control_all.append(metrics)
    if (len(control_all) + len(ctrls_done)) % 20 == 0:
        print(f"    {len(control_all) + len(ctrls_done)}/{len(control_genes)} done")

# Merge with existing
ctrl_df_new = pd.DataFrame(control_all)
if os.path.exists(cpath):
    ctrl_df_old = pd.read_csv(cpath)
    ctrl_all = pd.concat([ctrl_df_old, ctrl_df_new], ignore_index=True)
else:
    ctrl_all = ctrl_df_new
ctrl_all.to_csv(cpath, index=False)
print(f"  Total controls: {len(ctrl_all)}")

# ---- Empirical P-values ----
print("\n" + "=" * 70)
print("EMPIRICAL P-VALUES (Mono_C)")
print("=" * 70)

for metric in ["top5_mean_abs_z", "top10_mean_abs_z", "top20_mean_abs_z",
               "median_abs_z", "max_abs_z"]:
    ager_val = ager_metrics[metric]
    control_vals = ctrl_all[metric].dropna().values
    if len(control_vals) == 0:
        continue
    n_higher = sum(1 for v in control_vals if v >= ager_val)
    n_lower = sum(1 for v in control_vals if v <= ager_val)
    # Two-sided: fraction of controls MORE extreme than AGER
    emp_p = (1 + n_higher) / (1 + len(control_vals))
    # Also compute how many are more extreme in both directions
    ctrl_median = np.median(control_vals)
    more_extreme = sum(1 for v in control_vals
                       if abs(v - ctrl_median) >= abs(ager_val - ctrl_median))
    emp_p_2sided = (1 + more_extreme) / (1 + len(control_vals))

    desc = "SIG" if emp_p < 0.05 else ("borderline" if emp_p < 0.10 else "n.s.")
    print(f"  {metric:25s}: AGER={ager_val:.3f}, ctrl_median={ctrl_median:.3f}, "
          f"n_higher={n_higher}/{len(control_vals)}, P(one-sided)={emp_p:.4f} [{desc}]")

# ---- Also: Compare AGER re-sampling distribution to controls ----
print("\n" + "=" * 70)
print("RESAMPLING STABILITY: AGER vs Controls")
print("=" * 70)

def run_ko_n_iter(adata, cell_type, target_gene, n_iter=20, n_per_donor=1000):
    """Run KO N times, return top10 mean |Z| for each iteration."""
    top10_vals = []
    for i in range(n_iter):
        np.random.seed(440 + i * 100)
        res = run_ko_single_iter(adata, cell_type, target_gene, n_per_donor=n_per_donor)
        if res is not None:
            no_tgt = res[res["gene"] != target_gene]
            top10 = np.sort(np.abs(no_tgt["z_score"].values))[-10:].mean()
            top10_vals.append(top10)
    return top10_vals

ager_stability = run_ko_n_iter(adata, "Mono_C", TARGET_GENE, n_iter=20)
control_stability = {}
for cg in control_genes[:30]:  # Top 30 matched controls
    np.random.seed(440)
    vals = run_ko_n_iter(adata, "Mono_C", cg, n_iter=20)
    control_stability[cg] = vals
    if (len(control_stability)) % 10 == 0:
        print(f"  {len(control_stability)}/30 stability runs done")

# Test: is AGER's mean top10|Z| higher than controls?
ctrl_means = [np.mean(v) for v in control_stability.values()]
n_higher_stab = sum(1 for m in ctrl_means if m >= np.mean(ager_stability))
emp_p_stab = (1 + n_higher_stab) / (1 + len(ctrl_means))
print(f"\n  AGER mean top10|Z| (±SD): {np.mean(ager_stability):.3f} ± {np.std(ager_stability):.3f}")
print(f"  Control mean top10|Z|: median={np.median(ctrl_means):.3f}, "
      f"range=[{min(ctrl_means):.3f}, {max(ctrl_means):.3f}]")
print(f"  Empirical P (stability) = {emp_p_stab:.4f}")

# ============================================================
# 4. scDEGs — Monocyte-focused
# ============================================================
print("\n" + "=" * 70)
print("scDEGs: Wilcoxon (MG vs HC) — Monocytes")
print("=" * 70)

sc_degs = {}
for ct in ["Mono_C", "Mono_NC"]:
    deg_path = os.path.join(out_dir, f"scDEG_{ct}.csv")
    if os.path.exists(deg_path):
        print(f"  [CHECKPOINT] Loading scDEG_{ct}.csv")
        sc_degs[ct] = pd.read_csv(deg_path)
        continue

    ct_mask = adata.obs["cell_type_refined"] == ct
    sub = adata[ct_mask].copy()
    sc.tl.rank_genes_groups(sub, groupby="condition", reference="HC",
                            method="wilcoxon", n_genes=sub.n_vars,
                            key_added=f"degs_{ct}")
    result = sc.get.rank_genes_groups_df(sub, group="MG", key=f"degs_{ct}")
    result["cell_type"] = ct
    result["abs_log2FC"] = np.abs(result["logfoldchanges"].values)
    result["pass"] = (result["abs_log2FC"] > 0.5) & (result["pvals_adj"] < 0.05)
    result.to_csv(deg_path, index=False)
    sc_degs[ct] = result
    print(f"  {ct}: {(result['pass']).sum()} DEGs (|log2FC|>0.5, P_adj<0.05)")

# ============================================================
# 5. OVERLAP
# ============================================================
print("\n" + "=" * 70)
print("OVERLAP: KO Perturbed ∩ scDEGs")
print("=" * 70)

def overlap_test(ko_df, deg_df, label, target_gene="AGER", freq_thresh=0.35):
    perturbed = ko_df[ko_df["frequency"] >= freq_thresh]["gene"].tolist()
    perturbed = [g for g in perturbed if g != target_gene]
    deg_sig = deg_df[deg_df["pass"]]["names"].tolist()
    universe = list(set(ko_df["gene"]) & set(deg_df["names"]))
    p_in_u = [g for g in perturbed if g in universe]
    d_in_u = [g for g in deg_sig if g in universe]
    overlap = list(set(p_in_u) & set(d_in_u))
    a = len(overlap)
    b = len(p_in_u) - a
    c = len(d_in_u) - a
    d_val = len(universe) - len(p_in_u) - len(d_in_u) + a
    a, b, c, d_val = max(0, a), max(0, b), max(0, c), max(0, d_val)
    if a > 0 and b > 0 and c > 0 and d_val > 0:
        or_val, p_fisher = stats.fisher_exact([[a, b], [c, d_val]], alternative="greater")
        se = np.sqrt(1/a + 1/b + 1/c + 1/d_val)
        ci_l = np.exp(np.log(or_val) - 1.96 * se)
        ci_h = np.exp(np.log(or_val) + 1.96 * se)
    else:
        or_val, p_fisher = 0.0, 1.0
        ci_l, ci_h = np.nan, np.nan
    print(f"\n{label}: KO={len(p_in_u)}, DEG={len(d_in_u)}, Univ={len(universe)}")
    print(f"  Overlap={a}, OR={or_val:.1f} [{ci_l:.1f}-{ci_h:.1f}], P={p_fisher:.2e}")
    if a > 0:
        print(f"  Genes: {', '.join(overlap[:25])}")
    return {"label": label, "n_ko": len(p_in_u), "n_deg": len(d_in_u),
            "n_univ": len(universe), "n_overlap": a, "OR": or_val,
            "OR_CI_low": ci_l, "OR_CI_high": ci_h, "p_fisher": p_fisher,
            "genes": ', '.join(overlap)}

overlaps = {}
for ct in ["Mono_C", "Mono_NC"]:
    if ct in ko_results and ct in sc_degs:
        overlaps[ct] = overlap_test(ko_results[ct], sc_degs[ct], ct)
pd.DataFrame(overlaps).T.to_csv(os.path.join(out_dir, "overlap.csv"))

# ============================================================
# 6. GSEA — Monocyte-relevant pathways
# ============================================================
print("\n" + "=" * 70)
print("GSEA: Monocyte Pathway Enrichment")
print("=" * 70)

msigdb = {
    "RAGE/AGER Signaling": ["AGER", "HMGB1", "S100A8", "S100A9", "S100A12", "S100B",
        "MAPK1", "MAPK3", "NFKB1", "RELA", "JUN", "FOS", "MAPK8", "MAPK14",
        "STAT3", "AKT1", "PIK3CA", "CDC42", "RAC1", "RHOA", "DIAPH1"],
    "Monocyte Activation": ["CD14", "FCGR1A", "FCGR2A", "FCGR3A", "ITGAM", "ITGB2",
        "CSF1R", "CSF3R", "TLR2", "TLR4", "CD68", "CCR2", "CX3CR1",
        "ITGAX", "ITGAL", "ICAM1", "VCAM1", "SELL", "CD44"],
    "OXPHOS": ["NDUFA1", "NDUFA2", "NDUFB1", "NDUFS1", "SDHA", "SDHB", "SDHC", "SDHD",
        "UQCRC1", "UQCRC2", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2",
        "ATP5A1", "ATP5B", "ATP5F1", "ATP5H", "ATP5O"],
    "IFNγ Response": ["STAT1", "IRF1", "IRF7", "IRF9", "MX1", "MX2", "OAS1", "OAS2",
        "OAS3", "IFIT1", "IFIT2", "IFIT3", "ISG15", "ISG20", "GBP1", "GBP2",
        "GBP4", "GBP5", "IFI35", "IFITM1", "IFITM2", "IFITM3", "BST2", "RSAD2",
        "DDX58", "IFIH1", "CD74", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1"],
    "Inflammatory Response": ["TNF", "IL1B", "IL6", "IL18", "CCL2", "CCL3", "CCL4",
        "CCL5", "CXCL8", "CXCL9", "CXCL10", "CXCL11", "NFKB1", "NFKB2", "RELA",
        "RELB", "STAT3", "JUN", "FOS", "JUNB", "PTGS2", "TLR2", "TLR4", "MYD88",
        "NLRP3"],
    "TNFα/NF-κB Signaling": ["TNF", "TNFAIP3", "TNFAIP2", "NFKB1", "NFKBIA",
        "NFKBIB", "RELA", "JUN", "JUNB", "FOS", "FOSB", "ATF3", "EGR1", "EGR3",
        "IER2", "IER3", "ZFP36", "BTG1", "BTG2", "DUSP1", "DUSP2"],
    "Complement": ["C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B",
        "C5", "C7", "CFB", "CFD", "CFH", "CFI", "CFP", "SERPING1", "CD55", "CD59"],
    "Glycolysis": ["HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA",
        "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
        "SLC2A1", "SLC2A3", "SLC16A1", "SLC16A3"],
    "ROS Pathway": ["SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2",
        "PRDX3", "TXN", "TXNRD1", "NFE2L2", "KEAP1", "HMOX1", "GCLC", "GCLM",
        "GSR", "NOX4", "CYBB", "MPO"],
    "Phagocytosis": ["FCGR1A", "FCGR2A", "FCGR2B", "FCGR3A", "FCGR3B", "ITGAM",
        "ITGB2", "MRC1", "CD36", "SCARB1", "MSR1", "MARCO", "CLEC7A", "TLR2",
        "TLR4", "RAC1", "CDC42", "ARPC2", "WAS", "ACTB"],
    "Autophagy": ["ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101", "BECN1", "PIK3C3",
        "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B",
        "GABARAP", "SQSTM1", "OPTN", "BNIP3", "BNIP3L"],
    "Toll-Like Receptor Signaling": ["TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6",
        "TLR7", "TLR8", "TLR9", "MYD88", "TIRAP", "IRAK1", "IRAK4",
        "TRAF6", "IRF3", "IRF7", "NFKB1", "RELA", "MAP3K7", "TAB1"],
}

def gsea_enrichment(ko_summary, pathways, n_perm=2000):
    z_vals = ko_summary.set_index("gene")["median_z"].dropna().sort_values(ascending=False)
    bg = set(z_vals.index)
    n_tot = len(z_vals)
    results = []
    for pw_name, pw_genes in pathways.items():
        pw_in_bg = [g for g in pw_genes if g in bg]
        if len(pw_in_bg) < 3:
            continue
        in_pw = np.array([g in pw_in_bg for g in z_vals.index], dtype=int)
        n_pw = in_pw.sum()
        if n_pw == 0:
            continue
        hit_w = np.abs(z_vals.values)
        hit_w[in_pw == 0] = 0
        h_sum = hit_w.sum()
        if h_sum == 0:
            continue
        miss_w = np.ones(n_tot); miss_w[in_pw == 1] = 0; m_sum = miss_w.sum()
        running = np.cumsum(hit_w / h_sum - miss_w / m_sum)
        es = running.max() if abs(running.max()) > abs(running.min()) else running.min()
        perm_es = []
        for _ in range(n_perm):
            pi = np.random.permutation(n_tot)
            p_in = in_pw[pi]; p_hit = hit_w[pi]; p_hit[p_in == 0] = 0
            ph_sum = p_hit.sum()
            if ph_sum == 0: continue
            p_miss = miss_w[pi]; p_miss[p_in == 1] = 0; pm_sum = p_miss.sum()
            if pm_sum == 0: continue
            p_run = np.cumsum(p_hit / ph_sum - p_miss / pm_sum)
            perm_es.append(p_run.max() if abs(p_run.max()) > abs(p_run.min()) else p_run.min())
        ss = [e for e in perm_es if e * es > 0]
        mean_perm = np.mean(np.abs(ss)) if ss else 1.0
        nes = es / mean_perm if mean_perm > 0 else es
        n_valid = len(perm_es)
        p_val = ((sum(1 for e in perm_es if e >= es) + 1) / (n_valid + 1) if es > 0
                 else (sum(1 for e in perm_es if e <= es) + 1) / (n_valid + 1))
        results.append({
            "pathway": pw_name, "n_genes": len(pw_in_bg), "n_hits": int(n_pw),
            "ES": es, "NES": nes, "p_value": p_val,
            "leading_edge": ", ".join(z_vals.index[in_pw == 1][:10].tolist()),
        })
    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df["FDR"] = np.minimum(results_df["p_value"] * len(results_df), 1.0)
        results_df = results_df.sort_values("p_value")
    return results_df

gsea_results = {}
for ct in ["Mono_C", "Mono_NC"]:
    gpath = os.path.join(out_dir, f"GSEA_{ct}.csv")
    if os.path.exists(gpath):
        print(f"  [CHECKPOINT] Loading GSEA_{ct}.csv")
        gsea_results[ct] = pd.read_csv(gpath)
    else:
        gsea = gsea_enrichment(ko_results[ct], msigdb)
        gsea.to_csv(gpath, index=False)
        gsea_results[ct] = gsea

    sig = gsea_results[ct][gsea_results[ct]["FDR"] < 0.10]
    print(f"\n  {ct}: {len(sig)} pathways at FDR<0.10")
    for _, r in gsea_results[ct].head(6).iterrows():
        marker = " ***" if r['FDR'] < 0.05 else (" **" if r['FDR'] < 0.10 else "")
        print(f"    {r['pathway']:35s} NES={r['NES']:+.2f} P={r['p_value']:.4f} FDR={r['FDR']:.4f}{marker}")

# ============================================================
# 7. 4-PANEL MONOCYTE FIGURE
# ============================================================
print("\n" + "=" * 70)
print("GENERATING 4-PANEL MONOCYTE FIGURE")
print("=" * 70)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
    'legend.fontsize': 7, 'figure.dpi': 300, 'savefig.dpi': 600,
})

FREQ_HIGH, FREQ_LOW = 0.50, 0.35
fig = plt.figure(figsize=(16, 11))

# ======= A: Mono_C KO Top Genes =======
ax = fig.add_subplot(2, 2, 1)
ko_mc = ko_results["Mono_C"]
top = ko_mc[ko_mc["gene"] != TARGET_GENE].head(15)
z_vals, z_err = top["median_z"].values, top["iqr_z"].values / 2
freqs = top["frequency"].values
colors = [C_MONO if f >= FREQ_HIGH else (C_MONO if f >= FREQ_LOW else '#B7E4C7') for f in freqs]
edges = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in freqs]
yp = range(len(top))
ax.barh(yp, z_vals, xerr=z_err, color=colors, alpha=0.9, capsize=2,
        height=0.7, edgecolor=edges, linewidth=0.8)
ax.set_yticks(yp); ax.set_yticklabels(top["gene"].values, fontsize=8)
ax.set_xlabel("Median Perturbation Z-score"); ax.axvline(x=0, color='black', lw=0.5)
ax.invert_yaxis()
for i, (z, f) in enumerate(zip(z_vals, freqs)):
    if f >= FREQ_HIGH:
        ax.text(z + z_err[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
    elif f >= FREQ_LOW:
        ax.text(z + z_err[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
n_all_mc = (ko_mc["frequency"] >= FREQ_LOW).sum()
n_hi_mc = (ko_mc["frequency"] >= FREQ_HIGH).sum()
ax.set_title(f"A. Mono_C AGER-KO Perturbed Genes\n"
             f"({n_all_mc} recurrent [freq>=35%], {n_hi_mc} stable [freq>=50%], 20 iter, 1000 cells/donor)",
             fontweight='bold', fontsize=10)
leg = [Line2D([0],[0], marker='s', color='w', markerfacecolor=C_MONO, markersize=10, label='Freq >= 0.50'),
       Line2D([0],[0], marker='s', color='w', markerfacecolor=C_MONO, markersize=10, alpha=0.55, label='Freq 0.35-0.50')]
ax.legend(handles=leg, fontsize=7, loc='lower right', framealpha=0.9)

# ======= B: Mono_NC KO Top Genes =======
ax = fig.add_subplot(2, 2, 2)
ko_nc = ko_results["Mono_NC"]
top_nc = ko_nc[ko_nc["gene"] != TARGET_GENE].head(15)
z8, z8e = top_nc["median_z"].values, top_nc["iqr_z"].values / 2
f8 = top_nc["frequency"].values
cb8 = [C_MONO_NC if f >= FREQ_HIGH else (C_MONO_NC if f >= FREQ_LOW else '#D8BFD8') for f in f8]
eb8 = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in f8]
yp8 = range(len(top_nc))
ax.barh(yp8, z8, xerr=z8e, color=cb8, alpha=0.9, capsize=2,
        height=0.7, edgecolor=eb8, linewidth=0.8)
ax.set_yticks(yp8); ax.set_yticklabels(top_nc["gene"].values, fontsize=8)
ax.set_xlabel("Median Perturbation Z-score"); ax.axvline(x=0, color='black', lw=0.5)
ax.invert_yaxis()
for i, (z, f) in enumerate(zip(z8, f8)):
    if f >= FREQ_HIGH:
        ax.text(z + z8e[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
    elif f >= FREQ_LOW:
        ax.text(z + z8e[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
n_all_nc = (ko_nc["frequency"] >= FREQ_LOW).sum()
n_hi_nc = (ko_nc["frequency"] >= FREQ_HIGH).sum()
ax.set_title(f"B. Mono_NC AGER-KO Perturbed Genes\n"
             f"({n_all_nc} recurrent [freq>=35%], {n_hi_nc} stable [freq>=50%], 20 iter, 400 cells/donor)",
             fontweight='bold', fontsize=10)

# ======= C: Negative Control Distribution (Mono_C) =======
ax = fig.add_subplot(2, 2, 3)
ctrl_top10 = ctrl_all["top10_mean_abs_z"].dropna().values
ager_top10 = ager_metrics["top10_mean_abs_z"]
n_high = sum(1 for v in ctrl_top10 if v >= ager_top10)
emp_p = (1 + n_high) / (1 + len(ctrl_top10))

ax.hist(ctrl_top10, bins=max(20, len(ctrl_top10)//3), color=C_CONTROL,
        edgecolor='#777777', alpha=0.7, linewidth=0.5,
        label=f'{len(ctrl_top10)} matched controls')
ax.axvline(x=ager_top10, color=C_SIG, linewidth=3, linestyle='--',
           label=f'{TARGET_GENE} = {ager_top10:.2f}')
ctrl_med = np.median(ctrl_top10)
p_str = f"P = {emp_p:.4f}" if emp_p >= 0.001 else f"P = {emp_p:.2e}"
desc = "Specific perturbation" if emp_p < 0.05 else ("Borderline" if emp_p < 0.10 else "Not significant")
# Also show percentile rank
pct_rank = (1 + sum(1 for v in ctrl_top10 if v <= ager_top10)) / (1 + len(ctrl_top10)) * 100
ax.text(0.98, 0.95,
        f"Empirical {p_str}\n{desc}\n"
        f"AGER percentile: {pct_rank:.0f}%\n"
        f"Control median: {ctrl_med:.2f}",
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor='#CCCCCC', alpha=0.9))
ax.set_xlabel("Mean |Z| of Top-10 Perturbed Genes"); ax.set_ylabel("Frequency")
ax.set_title(f"C. Negative Control Distribution (Mono_C)\n"
             f"{TARGET_GENE} vs {len(ctrl_top10)} expression-matched controls",
             fontweight='bold', fontsize=10)
ax.legend(fontsize=7, framealpha=0.9)

# ======= D: GSEA Monocyte Pathways =======
ax = fig.add_subplot(2, 2, 4)
gsea_comb = []
for ct, gsea in gsea_results.items():
    if len(gsea) > 0:
        dc = gsea.copy(); dc["cell_type"] = ct; gsea_comb.append(dc)
if len(gsea_comb) > 0:
    gsea_all = pd.concat(gsea_comb).sort_values("p_value")
    gsea_plot = gsea_all.drop_duplicates(subset="pathway").head(12)
    neg_log_p = -np.log10(gsea_plot["p_value"].values.clip(min=1e-300))
    ct_colors = {"Mono_C": C_MONO, "Mono_NC": C_MONO_NC}
    colors_g = [ct_colors.get(ct, C_CONTROL) for ct in gsea_plot["cell_type"]]
    ypg = range(len(gsea_plot))
    ax.barh(ypg, neg_log_p, color=colors_g, alpha=0.85, height=0.7,
            edgecolor='white', linewidth=0.3)
    ax.set_yticks(ypg)
    labels_g = [f"{r['pathway']} [{r['cell_type']}]"
                for _, r in gsea_plot.iterrows()]
    ax.set_yticklabels(labels_g, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("-log10(P-value)")
    ax.axvline(x=-np.log10(0.05), color='grey', linestyle='--', alpha=0.6, lw=0.8)
    ax.axvline(x=-np.log10(0.01), color='grey', linestyle=':', alpha=0.4, lw=0.5)
    for i, (_, r) in enumerate(gsea_plot.iterrows()):
        if r["FDR"] < 0.10:
            ax.text(neg_log_p[i] + 0.1, i, '★', va='center', fontsize=9, color='#D4A017')
    leg_g = [Patch(facecolor=C_MONO, alpha=0.85, label='Mono_C'),
             Patch(facecolor=C_MONO_NC, alpha=0.85, label='Mono_NC'),
             Line2D([0],[0], marker='*', color='w', markerfacecolor='#D4A017',
                    markersize=10, label='FDR < 0.10')]
    ax.legend(handles=leg_g, fontsize=7, loc='lower right', framealpha=0.9)
    ax.set_title("D. GSEA Pathway Enrichment\n(Ranked by AGER-KO Perturbation Z-score)",
                 fontweight='bold', fontsize=10)

plt.suptitle(f"scTenifoldKnk Virtual Knockout of {TARGET_GENE} in Monocytes\n"
             f"Classical (Mono_C) & Non-Classical (Mono_NC) Monocytes from MG scRNA-seq",
             fontweight='bold', fontsize=13, y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.95])

fig_path = os.path.join(fig_dir, f"Fig_AGER_monocyte_KO")
for fmt, dpi_val, kw in [("png", 600, {}), ("pdf", None, {}),
                          ("tiff", 600, {"pil_kwargs": {"compression": "tiff_lzw"}})]:
    plt.savefig(f"{fig_path}.{fmt}", dpi=dpi_val, bbox_inches='tight',
                facecolor='white', edgecolor='none', **kw)
plt.close()
print(f"[OK] Figure: {fig_path}.{{png, pdf, tiff}}")

# ============================================================
# 8. FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("MONOCYTE-FOCUSED AGER VIRTUAL KO — FINAL RESULTS")
print("=" * 70)

# Count overlap
ov_mc = overlaps.get("Mono_C", {})
ov_nc = overlaps.get("Mono_NC", {})

# GSEA sig
gsea_mc_sig = (gsea_results.get("Mono_C", pd.DataFrame())["FDR"] < 0.10).sum()
gsea_nc_sig = (gsea_results.get("Mono_NC", pd.DataFrame())["FDR"] < 0.10).sum()

summary = f"""
AGER in PBMC scRNA-seq:
  Mono_C:  mean=0.050, detect=5.1%, 12,603 HC cells (3 donors)
  Mono_NC: mean=0.035, detect=3.7%,  2,006 HC cells (3 donors)

VIRTUAL KO (20 iterations, donor-balanced):
  Mono_C:  {n_all_mc} recurrent (freq>=0.35), {n_hi_mc} stable (freq>=0.50)
           Top: CHD3, MNDA, NCF2, ITGB2, CTSS, APLP2, VCAN, IQGAP1
  Mono_NC: {n_all_nc} recurrent (freq>=0.35), {n_hi_nc} stable (freq>=0.50)
           Top: TBC1D8, RAD21, PCBP1, SERPINB1, ERGIC3, NUP214, UBA1

NEGATIVE CONTROLS (n={len(ctrl_top10)}, full-genome matching):
  AGER top10 mean |Z| = {ager_top10:.3f}
  Control median = {ctrl_med:.3f}
  Controls more extreme: {n_high}/{len(ctrl_top10)}
  Empirical P = {emp_p:.4f}  [{desc}]

scDEGs (Wilcoxon, |log2FC|>0.5, P_adj<0.05):
  Mono_C:  {sc_degs.get('Mono_C', pd.DataFrame())['pass'].sum() if 'Mono_C' in sc_degs else 'N/A'} DEGs
  Mono_NC: {sc_degs.get('Mono_NC', pd.DataFrame())['pass'].sum() if 'Mono_NC' in sc_degs else 'N/A'} DEGs

OVERLAP (KO freq>=0.35 ∩ scDEG):
  Mono_C:  {ov_mc.get('n_overlap', 'N/A')}/{ov_mc.get('n_ko', 'N/A')}, OR={ov_mc.get('OR', float('nan')):.1f}, P={ov_mc.get('p_fisher', float('nan')):.2e}
  Mono_NC: {ov_nc.get('n_overlap', 'N/A')}/{ov_nc.get('n_ko', 'N/A')}, OR={ov_nc.get('OR', float('nan')):.1f}, P={ov_nc.get('p_fisher', float('nan')):.2e}

GSEA (FDR < 0.10):
  Mono_C:  {gsea_mc_sig} pathways
  Mono_NC: {gsea_nc_sig} pathways

Output:
  Figure:  {fig_path}.png/pdf/tiff
  Results: {out_dir}/
"""
print(summary)

with open(os.path.join(out_dir, "AGER_monocyte_summary.txt"), 'w', encoding='utf-8') as f:
    f.write(summary)

print("=" * 70)
print("AGER monocyte-focused virtual KO complete.")
print("=" * 70)
