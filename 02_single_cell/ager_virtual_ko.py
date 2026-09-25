"""
=============================================================================
AGER Virtual Knockout — scTenifoldKnk Analysis
=============================================================================
Analogous to SESN3 virtual KO validation work.

Background:
  AGER (Advanced Glycosylation End-product Specific Receptor, RAGE) was
  identified as a risk-increasing MR/colocalization hit in MG. It encodes
  the receptor for advanced glycation endproducts (RAGE), involved in
  inflammatory signaling.

AGER Expression in this dataset:
  Mono_C:   mean=0.046, detect=4.0%, 12,603 HC cells → PRIMARY TARGET
  Mono_NC:  mean=0.039, detect=3.6%,  2,006 HC cells → secondary
  CD4_NC:   mean=0.035, detect=2.1%, 18,536 HC cells → comparison
  CD8_NC:   mean=0.033, detect=2.0%,  1,606 HC cells → comparison

Note: AGER detection is lower than SESN3 (~2-4% vs ~16-19% in T cells),
so virtual KO results may be noisier. We use more cells (1000/donor for
Mono_C) to compensate.

Methodology (same as SESN3):
  1. Extract HC cells of target type, donor-balanced sampling
  2. Build gene-gene Pearson correlation network (2000 HVGs + AGER)
  3. Simulate AGER KO by zeroing all edges connected to AGER
  4. Compute perturbation Z-scores per gene
  5. 20 iterations for stability assessment (freq >= 0.35)
  6. Rank-based GSEA
  7. 100 matched negative controls
  8. Overlap with MG vs HC scDEGs
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
out_dir = os.path.join(data_dir, "ager_ko_results")
fig_dir = os.path.join(data_dir, "figures_meeting")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

np.random.seed(440)

TARGET_GENE = "AGER"

# Color palette
C_MONO = "#2CA02C"    # green for monocytes
C_MONO_NC = "#9467BD" # purple for non-classical monocytes
C_CD4 = "#1F78B4"
C_CD8 = "#E31A1C"
C_AGER = "#D62728"
C_CONTROL = "#AAAAAA"

# ============================================================
# 0. Load data & Check AGER
# ============================================================
print("=" * 70)
print("LOADING DATA & CHECKING AGER EXPRESSION")
print("=" * 70)

adata = sc.read_h5ad(os.path.join(data_dir, "adata_refined.h5ad"))
adata.obs_names_make_unique()
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")

# Verify AGER exists
if TARGET_GENE not in adata.var_names:
    raise ValueError(f"{TARGET_GENE} not found in dataset!")

# AGER expression summary
if issparse(adata.X):
    X = adata.X.toarray()
else:
    X = adata.X
ager_idx = list(adata.var_names).index(TARGET_GENE)
ager_expr = X[:, ager_idx]

print(f"\n{TARGET_GENE} overall: mean={ager_expr.mean():.4f}, "
      f"detection={(ager_expr>0).mean():.3f}")

# Per refined cell type (HC only)
print(f"\n{TARGET_GENE} expression in HC cells per cell type:")
cell_type_stats = {}
for ct in sorted(adata.obs["cell_type_refined"].unique()):
    mask = (adata.obs["cell_type_refined"] == ct) & (adata.obs["condition"] == "HC")
    n_hc = mask.sum()
    ct_expr = ager_expr[mask.values]
    ct_mean = ct_expr.mean()
    ct_detect = (ct_expr > 0).mean()
    donors = adata[mask].obs["sample"].unique()
    print(f"  {ct:15s}: HC={n_hc:6d}, mean={ct_mean:.4f}, detect={ct_detect:.3f}, donors={list(donors)}")
    cell_type_stats[ct] = {"n_hc": n_hc, "mean": ct_mean, "detect": ct_detect,
                           "n_donors": len(donors), "donors": list(donors)}

# ============================================================
# 1. CORE FUNCTIONS
# ============================================================
def extract_cells_donor_balanced(adata, cell_type, condition="HC", n_per_donor=500):
    """Extract donor-balanced cells of given refined cell type."""
    mask = (adata.obs["cell_type_refined"] == cell_type) & (adata.obs["condition"] == condition)
    sub = adata[mask].copy()
    cells_use = []
    for donor in sorted(sub.obs["sample"].unique()):
        donor_cells = list(sub.obs_names[sub.obs["sample"] == donor])
        n_take = min(n_per_donor, len(donor_cells))
        if n_take >= 50:
            cells_use.extend(list(np.random.choice(donor_cells, n_take, replace=False)))
    if len(cells_use) == 0:
        raise ValueError(f"No cells for {cell_type} in {condition}")
    return sub[cells_use, :].copy()

def prepare_ko_matrix(sub_adata, n_genes=2000, target_gene="AGER"):
    """Select HVGs + ensure target gene included."""
    if issparse(sub_adata.X):
        expr = sub_adata.X.toarray()
    else:
        expr = sub_adata.X

    gene_vars = np.var(expr, axis=0)
    top_idx = np.argsort(gene_vars)[-n_genes:]
    hvg_genes = [sub_adata.var_names[i] for i in top_idx]

    if target_gene in sub_adata.var_names and target_gene not in hvg_genes:
        hvg_genes.append(target_gene)

    gene_indices = [list(sub_adata.var_names).index(g) for g in hvg_genes
                    if g in sub_adata.var_names]
    expr_sub = expr[:, gene_indices]
    final_genes = [g for g in hvg_genes if g in sub_adata.var_names]
    return expr_sub, final_genes

def virtual_ko_single(expr, gene_list, target_gene):
    """Virtual KO: zero out target edges, compute perturbation per gene."""
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

    return pd.DataFrame({
        "gene": gene_list, "distance": perturbation, "z_score": z_scores,
    })

def virtual_ko_iterated(adata, cell_type, target_gene="AGER", n_iter=20,
                        n_per_donor=500, n_hvg=2000):
    """Run virtual KO with n_iter iterations."""
    all_z = {}
    all_freq = Counter()

    for i in range(n_iter):
        seed = 440 + i * 100
        np.random.seed(seed)
        sub = extract_cells_donor_balanced(adata, cell_type, "HC", n_per_donor=n_per_donor)
        expr, genes = prepare_ko_matrix(sub, n_genes=n_hvg, target_gene=target_gene)
        result = virtual_ko_single(expr, genes, target_gene)
        if result is None:
            print(f"  WARNING: {target_gene} not in gene list for iter {i}")
            continue
        for _, row in result.iterrows():
            g = row["gene"]
            all_z.setdefault(g, []).append(row["z_score"])
        # Mark genes in top 5% by |Z| as "detected" this iteration
        abs_z = np.abs(result["z_score"].values)
        thresh = np.percentile(abs_z, 95)
        sig_genes = result.loc[abs_z >= thresh, "gene"].tolist()
        for g in sig_genes:
            all_freq[g] += 1
        if (i + 1) % 10 == 0:
            print(f"  {cell_type}: {i+1}/{n_iter} done ({len(sig_genes)} top-5% genes)")

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

    summary_df = pd.DataFrame(summary).sort_values("median_z", ascending=False)
    return summary_df, all_z

# ============================================================
# 2. RUN VIRTUAL KO — Primary targets
# ============================================================
print("\n" + "=" * 70)
print(f"VIRTUAL KO: 20 Iterations for {TARGET_GENE}")
print("=" * 70)

# Determine target cell types with sufficient AGER expression and cells
# Mono_C: best monocyte population (12,603 HC cells, detection=4.0%)
# Mono_NC: non-classical monocytes (2,006 HC cells, detection=3.6%)
# CD4_NC: for comparison with SESN3 (18,536 HC cells, detection=2.1%)
# CD8_NC: for comparison with SESN3 (1,606 HC cells, detection=2.0%)

target_cell_types = {
    "Mono_C":   {"n_per_donor": 1000, "label": "Classical Monocytes"},
    "Mono_NC":  {"n_per_donor": 400,  "label": "Non-Classical Monocytes"},
    "CD4_NC":   {"n_per_donor": 500,  "label": "CD4 Naive/Central Memory T"},
    "CD8_NC":   {"n_per_donor": 400,  "label": "CD8 Naive/Central Memory T"},
}

ko_results = {}
for ct, params in target_cell_types.items():
    n_hc = cell_type_stats.get(ct, {}).get("n_hc", 0)
    n_donors = cell_type_stats.get(ct, {}).get("n_donors", 0)
    ct_detect = cell_type_stats.get(ct, {}).get("detect", 0)

    if n_hc < 100 or n_donors < 3:
        print(f"\n  SKIP {ct}: insufficient HC cells ({n_hc}) or donors ({n_donors})")
        continue
    if ct_detect < 0.01:
        print(f"\n  WARNING {ct}: very low detection ({ct_detect:.3f}) — results may be noisy")

    # CHECKPOINT: load existing results if available
    fname = f"KO_{ct}.csv"
    fpath = os.path.join(out_dir, fname)
    if os.path.exists(fpath):
        print(f"\n  [CHECKPOINT] Loading existing {fname}...")
        ko_summary = pd.read_csv(fpath)
        ko_results[ct] = ko_summary
        n_freq35 = (ko_summary["frequency"] >= 0.35).sum()
        n_freq50 = (ko_summary["frequency"] >= 0.50).sum()
        print(f"  freq>=0.35: {n_freq35}, freq>=0.50: {n_freq50}")
        continue

    print(f"\n{'='*60}")
    print(f"  {ct} ({params['label']}): {n_hc} HC cells, "
          f"{n_donors} donors, detect={ct_detect:.3f}")
    print(f"  n_per_donor={params['n_per_donor']}, n_iter=20, n_hvg=2000")
    print(f"{'='*60}")

    ko_summary, ko_z = virtual_ko_iterated(
        adata, ct, target_gene=TARGET_GENE,
        n_iter=20, n_per_donor=params["n_per_donor"], n_hvg=2000
    )

    # Save
    ko_summary.to_csv(fpath, index=False)
    ko_results[ct] = ko_summary

    # Summary stats
    n_freq35 = (ko_summary["frequency"] >= 0.35).sum()
    n_freq50 = (ko_summary["frequency"] >= 0.50).sum()
    n_freq70 = (ko_summary["frequency"] >= 0.70).sum()
    print(f"  freq>=0.35: {n_freq35}, freq>=0.50: {n_freq50}, freq>=0.70: {n_freq70}")

    # Top hits
    top_genes = ko_summary[ko_summary["gene"] != TARGET_GENE].head(15)
    print(f"  Top perturbed genes:")
    for _, r in top_genes.iterrows():
        marker = " ***" if r['frequency'] >= 0.50 else (" **" if r['frequency'] >= 0.35 else "")
        print(f"    {r['gene']:15s} Z={r['median_z']:+.2f} IQR={r['iqr_z']:.2f} freq={r['frequency']:.2f}{marker}")

print(f"\nKO results saved to: {out_dir}/")

# ============================================================
# 3. NEGATIVE CONTROLS — 100 matched genes
# ============================================================
print("\n" + "=" * 70)
print(f"NEGATIVE CONTROLS: 100 genes matched to {TARGET_GENE}")
print("=" * 70)

# Characterize AGER from Mono_C (best expressing)
primary_ct = "Mono_C"
sub_ref = extract_cells_donor_balanced(adata, primary_ct, "HC", n_per_donor=1000)
expr_ref, genes_ref = prepare_ko_matrix(sub_ref, target_gene=TARGET_GENE)
ager_idx_ref = genes_ref.index(TARGET_GENE)
ager_mean = expr_ref[:, ager_idx_ref].mean()
ager_detect = (expr_ref[:, ager_idx_ref] > 0).mean()
print(f"{TARGET_GENE} ({primary_ct}): mean_expr={ager_mean:.4f}, detection={ager_detect:.3f}")

# Known hits to exclude (all MR/colocalization hits)
known_hits = {"SESN3", "AGER", "ZSCAN12", "CFD", "DXO", "ZBTB9", "HLA-DOA",
              "BTN2A1", "HCG23", "TYMP", "TNFRSF14", "SPCS1", "MPV17L2",
              "RNASEH2C", "CTSH", "RBM47", "SLC22A4", "FCRL3", "SCAND3",
              "ZNF322", "IL12RB2", "ORMDL3", "HFE"}

gene_props = pd.DataFrame({
    "gene": genes_ref,
    "mean_expr": expr_ref.mean(axis=0),
    "detection": (expr_ref > 0).mean(axis=0),
})
gene_props = gene_props[~gene_props["gene"].isin(known_hits)]
gene_props = gene_props[gene_props["gene"] != TARGET_GENE]

# Tighter matching for AGER (wider range since expression is lower)
matched = gene_props[
    (gene_props["mean_expr"] >= ager_mean * 0.5) &
    (gene_props["mean_expr"] <= ager_mean * 2.0) &
    (np.abs(gene_props["detection"] - ager_detect) < 0.10)
].copy()
matched["match_score"] = (
    np.abs(matched["mean_expr"] / ager_mean - 1) +
    np.abs(matched["detection"] - ager_detect) * 10
)

n_controls = min(100, len(matched))
if n_controls < 30:
    print(f"WARNING: Only {n_controls} matched controls — relaxing criteria")
    matched = gene_props[
        (gene_props["mean_expr"] >= ager_mean * 0.3) &
        (gene_props["mean_expr"] <= ager_mean * 3.0) &
        (np.abs(gene_props["detection"] - ager_detect) < 0.15)
    ].copy()
    matched["match_score"] = (
        np.abs(matched["mean_expr"] / ager_mean - 1) +
        np.abs(matched["detection"] - ager_detect) * 10
    )
    n_controls = min(100, len(matched))

control_genes = matched.nsmallest(n_controls, "match_score")["gene"].tolist()
print(f"Selected {len(control_genes)} matched controls")
print(f"  Detection range: {matched.nsmallest(n_controls, 'match_score')['detection'].min():.3f} - "
      f"{matched.nsmallest(n_controls, 'match_score')['detection'].max():.3f}")

# Run controls: compute top-N mean |Z| for AGER and each control
def compute_top_z_metrics(adata, cell_type, target_gene, n_tops=[5, 10, 20],
                          n_per_donor=500):
    """Run virtual KO and return dict of mean |Z| for top-N genes."""
    sub = extract_cells_donor_balanced(adata, cell_type, "HC", n_per_donor=n_per_donor)
    expr, genes = prepare_ko_matrix(sub, target_gene=target_gene)
    result = virtual_ko_single(expr, genes, target_gene)
    if result is None:
        return None
    result_no_target = result[result["gene"] != target_gene]
    abs_z = np.abs(result_no_target["z_score"].values)
    metrics = {}
    for n in n_tops:
        metrics[f"top{n}_mean_abs_z"] = np.sort(abs_z)[-n:].mean()
    return metrics

# Run for primary cell type (Mono_C)
print(f"\nRunning {len(control_genes)} controls on {primary_ct} (est. ~15 min)...")
control_metrics = {}

# AGER itself
ager_m = compute_top_z_metrics(adata, primary_ct, TARGET_GENE, n_per_donor=1000)
if ager_m:
    print(f"\n  {TARGET_GENE} {primary_ct}: top5={ager_m['top5_mean_abs_z']:.3f}, "
          f"top10={ager_m['top10_mean_abs_z']:.3f}, top20={ager_m['top20_mean_abs_z']:.3f}")
    control_metrics[f"{TARGET_GENE}_{primary_ct}"] = ager_m

for i, cg in enumerate(control_genes):
    np.random.seed(440 + i)
    cm = compute_top_z_metrics(adata, primary_ct, cg, n_per_donor=1000)
    if cm is not None:
        control_metrics[f"{cg}_{primary_ct}"] = cm
    if (i + 1) % 50 == 0:
        print(f"    {i+1}/{len(control_genes)} controls done")

# Also run for CD4_NC (comparison with SESN3)
if "CD4_NC" in ko_results:
    print(f"\nRunning controls for CD4_NC (comparison)...")
    ager_cd4_m = compute_top_z_metrics(adata, "CD4_NC", TARGET_GENE, n_per_donor=500)
    if ager_cd4_m:
        control_metrics[f"{TARGET_GENE}_CD4_NC"] = ager_cd4_m
    for i, cg in enumerate(control_genes):
        np.random.seed(440 + i)
        cm = compute_top_z_metrics(adata, "CD4_NC", cg, n_per_donor=500)
        if cm is not None:
            control_metrics[f"{cg}_CD4_NC"] = cm
        if (i + 1) % 50 == 0:
            print(f"    {i+1}/{len(control_genes)} controls done for CD4_NC")

# Empirical P-values
print("\n" + "=" * 70)
print("EMPIRICAL P-VALUES")
print("=" * 70)

for ct_label in [primary_ct, "CD4_NC"]:
    ager_key = f"{TARGET_GENE}_{ct_label}"
    if ager_key not in control_metrics:
        continue
    control_keys = [k for k in control_metrics if k.endswith(f"_{ct_label}")
                    and not k.startswith(TARGET_GENE)]

    if len(control_keys) == 0:
        continue

    for metric in ["top5_mean_abs_z", "top10_mean_abs_z", "top20_mean_abs_z"]:
        ager_val = control_metrics[ager_key][metric]
        control_vals = [control_metrics[k][metric] for k in control_keys]
        n_higher = sum(1 for v in control_vals if v >= ager_val)
        emp_p = (1 + n_higher) / (1 + len(control_vals))
        desc = "SIGNIFICANT" if emp_p < 0.05 else ("borderline" if emp_p < 0.10 else "n.s.")
        print(f"  {ct_label} {metric}: {TARGET_GENE}={ager_val:.3f}, "
              f"control_median={np.median(control_vals):.3f}, "
              f"n_higher={n_higher}/{len(control_vals)}, "
              f"Empirical P={emp_p:.4f} [{desc}]")

# Save controls
ctrl_df = pd.DataFrame([
    {"gene": k.rsplit("_", 1)[0], "cell_type": k.rsplit("_", 1)[1], **v}
    for k, v in control_metrics.items()
])
ctrl_df.to_csv(os.path.join(out_dir, "negative_controls.csv"), index=False)
print(f"\nControls saved: {out_dir}/negative_controls.csv")

# ============================================================
# 4. scDEGs — Wilcoxon rank-sum per refined cell type
# ============================================================
print("\n" + "=" * 70)
print("scDEGs: Wilcoxon Rank-Sum Test (MG vs HC)")
print("=" * 70)

deg_targets = [ct for ct in ["Mono_C", "Mono_NC", "CD4_NC", "CD8_NC"]
               if ct in ko_results]

sc_degs = {}
for ct in deg_targets:
    print(f"\nComputing DEGs for {ct}...")
    ct_mask = adata.obs["cell_type_refined"] == ct
    sub = adata[ct_mask].copy()

    n_hc = (sub.obs["condition"] == "HC").sum()
    n_mg = (sub.obs["condition"] == "MG").sum()
    print(f"  HC={n_hc}, MG={n_mg}")

    if n_hc < 10 or n_mg < 10:
        print(f"  SKIP: insufficient cells")
        continue

    sc.tl.rank_genes_groups(
        sub, groupby="condition", reference="HC",
        method="wilcoxon", n_genes=sub.n_vars,
        key_added=f"degs_{ct}"
    )

    result = sc.get.rank_genes_groups_df(sub, group="MG", key=f"degs_{ct}")
    result["cell_type"] = ct
    result["abs_log2FC"] = np.abs(result["logfoldchanges"].values)
    result["pass_sc"] = (result["abs_log2FC"] > 0.5) & (result["pvals_adj"] < 0.05)
    result["pass_relaxed"] = (result["abs_log2FC"] > 0.5) & (result["pvals"] < 0.01)

    sc_degs[ct] = result
    print(f"  DEGs (|log2FC|>0.5 + P_adj<0.05): {result['pass_sc'].sum()}")
    print(f"  DEGs (|log2FC|>0.5 + P<0.01): {result['pass_relaxed'].sum()}")
    result.to_csv(os.path.join(out_dir, f"scDEG_{ct}.csv"), index=False)

# ============================================================
# 5. OVERLAP ANALYSIS: KO Perturbed ∩ scDEGs
# ============================================================
print("\n" + "=" * 70)
print("OVERLAP: KO-Perturbed ∩ scDEGs")
print("=" * 70)

def overlap_test(ko_summary, deg_result, label, target_gene=TARGET_GENE, freq_thresh=0.35):
    """Test overlap with scDEGs."""
    results = {}
    for deg_label, deg_col in [
        ("scDEG strict (|log2FC|>0.5, P_adj<0.05)", "pass_sc"),
        ("scDEG relaxed (|log2FC|>0.5, P<0.01)", "pass_relaxed"),
    ]:
        perturbed = ko_summary[ko_summary["frequency"] >= freq_thresh]["gene"].tolist()
        perturbed = [g for g in perturbed if g != target_gene]
        deg_sig = deg_result[deg_result[deg_col]]["names"].tolist()
        universe = list(set(ko_summary["gene"]) & set(deg_result["names"]))

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

        key = f"{label} | {deg_label}"
        results[key] = {
            "label": label, "deg_type": deg_label,
            "n_ko": len(p_in_u), "n_deg": len(d_in_u),
            "n_univ": len(universe), "n_overlap": a,
            "OR": or_val, "OR_CI_low": ci_l, "OR_CI_high": ci_h,
            "p_fisher": p_fisher, "genes": ', '.join(overlap),
        }

        sig = "***" if p_fisher < 0.01 else ("**" if p_fisher < 0.05 else "")
        print(f"\n{key}:")
        print(f"  KO={len(p_in_u)}, DEG={len(d_in_u)}, Universe={len(universe)}")
        print(f"  Overlap={a}, OR={or_val:.1f} [{ci_l:.1f}-{ci_h:.1f}], P={p_fisher:.2e} {sig}")
        if a > 0:
            print(f"  Genes: {', '.join(overlap[:20])}")

    return results

all_overlaps = {}
for ct in deg_targets:
    if ct in ko_results and ct in sc_degs:
        all_overlaps.update(overlap_test(ko_results[ct], sc_degs[ct], ct))

ov_df = pd.DataFrame(all_overlaps).T
ov_df.to_csv(os.path.join(out_dir, "overlap_scDEG.csv"))
print(f"\nOverlap results saved: {out_dir}/overlap_scDEG.csv")

# ============================================================
# 6. GSEA PATHWAY ENRICHMENT
# ============================================================
print("\n" + "=" * 70)
print("GSEA: Rank-Based Pathway Enrichment")
print("=" * 70)

# MSigDB hallmark + immune pathways
msigdb = {
    "IFNγ Response": ["STAT1", "IRF1", "IRF7", "IRF9", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
        "IFIT1", "IFIT2", "IFIT3", "ISG15", "ISG20", "GBP1", "GBP2", "GBP4", "GBP5",
        "IFI35", "IFITM1", "IFITM2", "IFITM3", "BST2", "RSAD2", "DDX58", "IFIH1",
        "CD74", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "PSMB9"],
    "Inflammatory Response": ["TNF", "IL1B", "IL6", "IL18", "CCL2", "CCL3", "CCL4", "CCL5",
        "CXCL8", "CXCL9", "CXCL10", "CXCL11", "NFKB1", "NFKB2", "RELA", "RELB",
        "STAT3", "JUN", "FOS", "JUNB", "PTGS2", "TLR2", "TLR4", "MYD88", "NLRP3"],
    "TNFα/NF-κB Signaling": ["TNF", "TNFAIP3", "TNFAIP2", "NFKB1", "NFKBIA", "NFKBIB",
        "RELA", "JUN", "JUNB", "FOS", "FOSB", "ATF3", "EGR1", "EGR3",
        "IER2", "IER3", "ZFP36", "BTG1", "BTG2", "DUSP1", "DUSP2"],
    "OXPHOS": ["NDUFA1", "NDUFA2", "NDUFB1", "NDUFS1", "SDHA", "SDHB", "SDHC", "SDHD",
        "UQCRC1", "UQCRC2", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2",
        "ATP5A1", "ATP5B", "ATP5F1", "ATP5H", "ATP5O"],
    "ROS Pathway": ["SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2", "PRDX3",
        "TXN", "TXNRD1", "NFE2L2", "KEAP1", "HMOX1", "GCLC", "GCLM", "GSR",
        "NOX4", "CYBB", "MPO"],
    "mTORC1 Signaling": ["MTOR", "RPTOR", "RICTOR", "AKT1", "AKT1S1", "DEPTOR", "MLST8",
        "TSC1", "TSC2", "RHEB", "EIF4EBP1", "RPS6KB1", "RPS6", "EIF4E",
        "ULK1", "ATG13", "RB1CC1", "SESN1", "SESN2"],
    "Autophagy": ["ULK1", "ULK2", "ATG13", "RB1CC1", "ATG101", "BECN1", "PIK3C3",
        "ATG3", "ATG5", "ATG7", "ATG12", "ATG16L1", "MAP1LC3A", "MAP1LC3B",
        "GABARAP", "SQSTM1", "OPTN", "BNIP3", "BNIP3L"],
    "Apoptosis": ["BCL2", "BCL2L1", "BCL2L11", "BAX", "BAK1", "BAD", "BID",
        "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CYCS", "APAF1",
        "XIAP", "BIRC2", "BIRC3", "TP53", "PMAIP1", "BBC3"],
    "IL-2/STAT5 Signaling": ["IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B",
        "LCK", "FYN", "CISH", "SOCS1", "SOCS2", "SOCS3", "BCL2", "BCL2L1",
        "CCND2", "CCND3", "MYC", "PIM1", "IRF4", "GZMB", "IFNG"],
    "Complement": ["C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B",
        "C5", "C7", "CFB", "CFD", "CFH", "CFI", "CFP", "SERPING1", "CD55", "CD59"],
    "TCR Signaling": ["CD3D", "CD3E", "CD3G", "CD247", "LCK", "FYN", "ZAP70",
        "LAT", "ITK", "PLCG1", "PRKCQ", "CARD11", "MALT1", "BCL10",
        "NFATC1", "NFATC2", "NFKB1", "RELA", "JUN", "FOS"],
    "IL-6/JAK/STAT3": ["IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3",
        "CRP", "SAA1", "SAA2", "HP", "FGB", "FGA", "FGG"],
    "Glycolysis": ["HK1", "HK2", "HK3", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA",
        "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
        "SLC2A1", "SLC2A3", "SLC16A1", "SLC16A3"],
    "Hypoxia": ["HIF1A", "ARNT", "EPAS1", "HIF3A", "VEGFA", "VEGFB", "PGK1",
        "LDHA", "SLC2A1", "BNIP3", "BNIP3L", "EGLN1", "EGLN2", "EGLN3"],
    "p53 Pathway": ["TP53", "CDKN1A", "BBC3", "PMAIP1", "BAX", "FAS", "MDM2",
        "RRM2B", "SESN1", "SESN2", "SESN3", "GADD45A", "GADD45B",
        "FDXR", "TRIAP1", "ZMAT3"],
    "Toll-Like Receptor Signaling": ["TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6",
        "TLR7", "TLR8", "TLR9", "MYD88", "TIRAP", "IRAK1", "IRAK4",
        "TRAF6", "IRF3", "IRF7", "NFKB1", "RELA", "MAP3K7", "TAB1"],
    "RAGE/AGER Signaling": ["AGER", "HMGB1", "S100A8", "S100A9", "S100A12", "S100B",
        "MAPK1", "MAPK3", "NFKB1", "RELA", "JUN", "FOS", "MAPK8", "MAPK14",
        "STAT3", "AKT1", "PIK3CA", "CDC42", "RAC1", "RHOA", "DIAPH1"],
    "Monocyte Activation": ["CD14", "FCGR1A", "FCGR2A", "FCGR3A", "ITGAM", "ITGB2",
        "CSF1R", "CSF3R", "TLR2", "TLR4", "CD68", "CCR2", "CX3CR1",
        "ITGAX", "ITGAL", "ICAM1", "VCAM1", "SELL", "CD44"],
}

def gsea_enrichment(ko_summary, pathways, n_perm=2000):
    """GSEA with weighted KS statistic and permutation-based P-value."""
    z_vals = ko_summary.set_index("gene")["median_z"].dropna()
    z_vals = z_vals.sort_values(ascending=False)
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
        miss_w = np.ones(n_tot)
        miss_w[in_pw == 1] = 0
        m_sum = miss_w.sum()

        running = np.cumsum(hit_w / h_sum - miss_w / m_sum)
        es = running.max() if abs(running.max()) > abs(running.min()) else running.min()

        perm_es = []
        for _ in range(n_perm):
            pi = np.random.permutation(n_tot)
            p_in = in_pw[pi]
            p_hit = hit_w[pi]
            p_hit[p_in == 0] = 0
            ph_sum = p_hit.sum()
            if ph_sum == 0:
                continue
            p_miss = miss_w[pi]
            p_miss[p_in == 1] = 0
            pm_sum = p_miss.sum()
            if pm_sum == 0:
                continue
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
for ct in target_cell_types:
    if ct not in ko_results:
        continue
    print(f"\nGSEA for {ct}...")
    gsea = gsea_enrichment(ko_results[ct], msigdb)
    gsea.to_csv(os.path.join(out_dir, f"GSEA_{ct}.csv"), index=False)
    gsea_results[ct] = gsea

    n_sig = (gsea["FDR"] < 0.10).sum()
    print(f"  Pathways at FDR<0.10: {n_sig}")
    for _, r in gsea.head(6).iterrows():
        sig = "***" if r['FDR'] < 0.05 else ("**" if r['FDR'] < 0.10 else "")
        print(f"    {r['pathway']:35s} NES={r['NES']:+.2f} P={r['p_value']:.3f} FDR={r['FDR']:.3f} {sig}")

# ============================================================
# 7. PUBLICATION-QUALITY 6-PANEL FIGURE
# ============================================================
print("\n" + "=" * 70)
print("GENERATING PUBLICATION FIGURE")
print("=" * 70)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

FREQ_HIGH, FREQ_LOW = 0.50, 0.35
fig = plt.figure(figsize=(21, 13.5))

# ======= A: Workflow =======
ax = fig.add_subplot(2, 3, 1)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
steps = [
    (1.0, 8.2, 8.0, 1.6, "HC Monocytes / T cells\n(3 donors, 12 refined cell types)", "#E8F5E9"),
    (1.0, 6.3, 8.0, 1.4, "Donor-balanced sampling\n(400-1000 cells/donor × 20 iterations)", "#E8F5E9"),
    (1.0, 4.5, 8.0, 1.3, "scTenifoldKnk AGER KO\nZero AGER edges in gene-gene\ncorrelation network", "#FFEBEE"),
    (1.0, 2.7, 8.0, 1.3, "Perturbation Z-scores\n(median ± IQR, recurrence frequency)", "#FFEBEE"),
    (1.0, 0.9, 8.0, 1.3, "Validation: GSEA · 100 controls\nOverlap with scDEGs", "#E3F2FD"),
]
for x, y, w, h, text, color in steps:
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='#AAAAAA', linewidth=0.7))
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=7, family='monospace')
for i in range(len(steps) - 1):
    _, y1, _, h1, _, _ = steps[i]; _, y2, _, h2, _, _ = steps[i+1]
    ax.annotate('', xy=(5.0, y2 + h2), xytext=(5.0, y1),
                arrowprops=dict(arrowstyle='->', color='#666666', lw=1.5))
ax.set_title(f"A. AGER Virtual KO Analysis Workflow", fontweight='bold', fontsize=11, loc='left')

# ======= B: Mono_C KO Results (Primary) =======
ax = fig.add_subplot(2, 3, 2)
if "Mono_C" in ko_results:
    ko_mc = ko_results["Mono_C"]
    top = ko_mc[ko_mc["gene"] != TARGET_GENE].head(15)
    z_vals, z_err = top["median_z"].values, top["iqr_z"].values / 2
    freqs = top["frequency"].values
    colors_b = [C_MONO if f >= FREQ_HIGH else (C_MONO if f >= FREQ_LOW else '#B7E4C7') for f in freqs]
    edge_b = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in freqs]
    y_pos = range(len(top))
    ax.barh(y_pos, z_vals, xerr=z_err, color=colors_b, alpha=0.9, capsize=2,
            height=0.7, edgecolor=edge_b, linewidth=0.8)
    ax.set_yticks(y_pos); ax.set_yticklabels(top["gene"].values, fontsize=7.5)
    ax.set_xlabel("Median Perturbation Z-score (± IQR/2)"); ax.axvline(x=0, color='black', lw=0.5)
    ax.invert_yaxis()
    for i, (z, f) in enumerate(zip(z_vals, freqs)):
        if f >= FREQ_HIGH:
            ax.text(z + z_err[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
        elif f >= FREQ_LOW:
            ax.text(z + z_err[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
    n_all = (ko_mc["frequency"] >= FREQ_LOW).sum()
    n_hi = (ko_mc["frequency"] >= FREQ_HIGH).sum()
    ax.set_title(f"B. Mono_C AGER-KO Perturbed Genes\n({n_all} recurrent [freq≥35%], {n_hi} stable [freq≥50%])",
                 fontweight='bold', fontsize=10)
    leg_b = [Line2D([0],[0], marker='s', color='w', markerfacecolor=C_MONO, markersize=10, label='Freq ≥ 0.50'),
             Line2D([0],[0], marker='s', color='w', markerfacecolor=C_MONO, markersize=10, alpha=0.55, label='Freq 0.35–0.50')]
    ax.legend(handles=leg_b, fontsize=6.5, loc='lower right', framealpha=0.9)
else:
    ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', transform=ax.transAxes)
    ax.set_title("B. Mono_C AGER-KO", fontweight='bold', fontsize=10)

# ======= C: CD4_NC KO Results (Comparison with SESN3) =======
ax = fig.add_subplot(2, 3, 3)
if "CD4_NC" in ko_results:
    ko_cd4 = ko_results["CD4_NC"]
    top_cd4 = ko_cd4[ko_cd4["gene"] != TARGET_GENE].head(15)
    z4, z4e = top_cd4["median_z"].values, top_cd4["iqr_z"].values / 2
    f4 = top_cd4["frequency"].values
    cb4 = [C_CD4 if f >= FREQ_HIGH else (C_CD4 if f >= FREQ_LOW else '#BBDEFB') for f in f4]
    eb4 = ['#222222' if f >= FREQ_HIGH else ('#666666' if f >= FREQ_LOW else '#CCCCCC') for f in f4]
    yp4 = range(len(top_cd4))
    ax.barh(yp4, z4, xerr=z4e, color=cb4, alpha=0.9, capsize=2,
            height=0.7, edgecolor=eb4, linewidth=0.8)
    ax.set_yticks(yp4); ax.set_yticklabels(top_cd4["gene"].values, fontsize=7.5)
    ax.set_xlabel("Median Perturbation Z-score (± IQR/2)"); ax.axvline(x=0, color='black', lw=0.5)
    ax.invert_yaxis()
    for i, (z, f) in enumerate(zip(z4, f4)):
        if f >= FREQ_HIGH:
            ax.text(z + z4e[i] + 0.15, i, '★', va='center', fontsize=8, color='#D4A017')
        elif f >= FREQ_LOW:
            ax.text(z + z4e[i] + 0.15, i, '◆', va='center', fontsize=7, color='#999999')
    n4a = (ko_cd4["frequency"] >= FREQ_LOW).sum()
    n4h = (ko_cd4["frequency"] >= FREQ_HIGH).sum()
    ax.set_title(f"C. CD4_NC AGER-KO Perturbed Genes\n({n4a} recurrent [freq≥35%], {n4h} stable [freq≥50%])",
                 fontweight='bold', fontsize=10)
else:
    ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', transform=ax.transAxes)
    ax.set_title("C. CD4_NC AGER-KO", fontweight='bold', fontsize=10)

# ======= D: Negative Control Distribution (Mono_C) =======
ax = fig.add_subplot(2, 3, 4)
ager_key = f"{TARGET_GENE}_{primary_ct}"
control_keys_mc = [k for k in control_metrics if k.endswith(f"_{primary_ct}") and not k.startswith(TARGET_GENE)]
if ager_key in control_metrics and len(control_keys_mc) > 0:
    ager_top10 = control_metrics[ager_key]["top10_mean_abs_z"]
    control_top10 = [control_metrics[k]["top10_mean_abs_z"] for k in control_keys_mc]
    n_high = sum(1 for v in control_top10 if v >= ager_top10)
    emp_p = (1 + n_high) / (1 + len(control_top10))

    ax.hist(control_top10, bins=25, color=C_CONTROL, edgecolor='#777777',
            alpha=0.7, linewidth=0.5, label=f'{len(control_top10)} matched controls')
    ax.axvline(x=ager_top10, color=C_AGER, linewidth=3, linestyle='--',
               label=f'{TARGET_GENE} = {ager_top10:.2f}')
    p_str = f"P = {emp_p:.4f}" if emp_p >= 0.001 else f"P = {emp_p:.2e}"
    desc = "Specific perturbation" if emp_p < 0.05 else ("Borderline" if emp_p < 0.10 else "Not significant")
    ax.text(0.98, 0.95, f"{primary_ct} Empirical {p_str}\n{desc}",
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='#CCCCCC', alpha=0.9))
    ax.set_xlabel(f"Mean |Z| of Top-10 Perturbed Genes"); ax.set_ylabel("Frequency")
    ax.set_title(f"D. Negative Control Distribution ({primary_ct})\n"
                 f"{TARGET_GENE} vs {len(control_top10)} matched controls",
                 fontweight='bold', fontsize=10)
    ax.legend(fontsize=8, framealpha=0.9)
else:
    ax.text(0.5, 0.5, "No control data", ha='center', va='center', transform=ax.transAxes)
    ax.set_title("D. Negative Controls", fontweight='bold', fontsize=10)

# ======= E: GSEA Pathway Enrichment =======
ax = fig.add_subplot(2, 3, 5)
gsea_comb = []
for ct, gsea in gsea_results.items():
    if len(gsea) > 0:
        dc = gsea.copy(); dc["cell_type"] = ct; gsea_comb.append(dc)
if len(gsea_comb) > 0:
    gsea_all = pd.concat(gsea_comb).sort_values("p_value")
    gsea_plot = gsea_all.drop_duplicates(subset="pathway").head(16)
    neg_log_p = -np.log10(gsea_plot["p_value"].values.clip(min=1e-300))

    # Map cell types to colors
    ct_color_map = {"Mono_C": C_MONO, "Mono_NC": C_MONO_NC,
                    "CD4_NC": C_CD4, "CD8_NC": C_CD8}
    colors_g = [ct_color_map.get(ct, C_CONTROL) for ct in gsea_plot["cell_type"]]

    ypg = range(len(gsea_plot))
    ax.barh(ypg, neg_log_p, color=colors_g, alpha=0.85, height=0.7,
            edgecolor='white', linewidth=0.3)
    ax.set_yticks(ypg)
    labels_g = [f"{r['pathway']} [{r['cell_type']}]" for _, r in gsea_plot.iterrows()]
    ax.set_yticklabels(labels_g, fontsize=6.5); ax.invert_yaxis()
    ax.set_xlabel("-log10(P-value)")
    ax.axvline(x=-np.log10(0.05), color='grey', linestyle='--', alpha=0.6, lw=0.8)
    ax.axvline(x=-np.log10(0.01), color='grey', linestyle=':', alpha=0.4, lw=0.5)
    for i, (_, r) in enumerate(gsea_plot.iterrows()):
        if r["FDR"] < 0.10:
            ax.text(neg_log_p[i] + 0.1, i, '★', va='center', fontsize=9, color='#D4A017')
    leg_g = [Patch(facecolor=C_MONO, alpha=0.85, label='Mono_C'),
             Patch(facecolor=C_CD4, alpha=0.85, label='CD4_NC'),
             Line2D([0],[0], marker='*', color='w', markerfacecolor='#D4A017', markersize=10, label='FDR < 0.10')]
    ax.legend(handles=leg_g, fontsize=7, loc='lower right', framealpha=0.9)
    ax.set_title(f"E. GSEA Pathway Enrichment\n(Ranked by {TARGET_GENE}-KO Perturbation Z-score)",
                 fontweight='bold', fontsize=10)
else:
    ax.text(0.5, 0.5, "No GSEA data", ha='center', va='center', transform=ax.transAxes)
    ax.set_title("E. GSEA", fontweight='bold', fontsize=10)

# ======= F: Overlap Forest Plot =======
ax = fig.add_subplot(2, 3, 6)
# Show overlap for Mono_C and CD4_NC
plot_overlaps = []
plot_labels = []
plot_colors = []
for ct, clr in [("Mono_C", C_MONO), ("CD4_NC", C_CD4)]:
    key = f"{ct} | scDEG relaxed (|log2FC|>0.5, P<0.01)"
    if key in all_overlaps:
        plot_overlaps.append(all_overlaps[key])
        plot_labels.append(ct)
        plot_colors.append(clr)

if len(plot_overlaps) > 0:
    for i, (ov, lbl, clr) in enumerate(zip(plot_overlaps, plot_labels, plot_colors)):
        or_val = ov["OR"]
        ci_l, ci_h = ov["OR_CI_low"], ov["OR_CI_high"]
        if np.isnan(or_val) or or_val <= 0:
            or_disp, ci_valid, or_label = 0.5, False, "N/A"
        elif np.isinf(or_val) or or_val > 100:
            or_disp, ci_valid, or_label = 50, not np.isnan(ci_l), "∞"
        else:
            or_disp, ci_valid, or_label = or_val, not np.isnan(ci_l), f"{or_val:.1f}"

        y = len(plot_overlaps) - 1 - i
        mk = 'D' if i == 0 else 's'
        ax.plot(or_disp, y, marker=mk, color=clr, markersize=12,
                markeredgewidth=1.5, markeredgecolor='white', zorder=5)
        if ci_valid and not np.isinf(ci_l) and not np.isinf(ci_h):
            ci_lc = max(0.05, min(ci_l, 50)); ci_hc = min(ci_h, 50)
            ax.plot([ci_lc, ci_hc], [y, y], '-', color=clr, linewidth=3, zorder=4)

        p_str = f"{ov['p_fisher']:.1e}" if ov['p_fisher'] < 0.01 else f"{ov['p_fisher']:.3f}"
        ax.text(55, y, f"OR={or_label}  P={p_str}\nk={ov['n_overlap']}/{ov['n_ko']} (KO∩DEG)",
                va='center', fontsize=7, family='monospace',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=clr, alpha=0.85))

    ax.set_ylim(-0.8, len(plot_overlaps) - 0.2)
    ax.set_yticks(range(len(plot_overlaps)))
    ax.set_yticklabels(plot_labels, fontsize=10, fontweight='bold')
    ax.axvline(x=1, color='black', linestyle='-', lw=0.8, zorder=1)
    ax.set_xlabel("Odds Ratio (95% CI)")
    ax.set_xscale('symlog', linthresh=1, linscale=0.5); ax.set_xlim(0.03, 75)
    ax.set_title(f"F. {TARGET_GENE}-KO Perturbed ∩ scDEGs\n(Fisher exact test, one-sided, |log2FC|>0.5, P<0.01)",
                 fontweight='bold', fontsize=10)
else:
    ax.text(0.5, 0.5, "No overlap data", ha='center', va='center', transform=ax.transAxes)
    ax.set_title("F. Overlap Forest Plot", fontweight='bold', fontsize=10)

# ---- Final ----
plt.suptitle(f"scTenifoldKnk Virtual Knockout of {TARGET_GENE} Reveals\n"
             f"Perturbed Gene Regulatory Networks in MG-Relevant Myeloid and T Cells",
             fontweight='bold', fontsize=14, y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.96])

fig_path = os.path.join(fig_dir, f"Fig_AGER_virtual_KO")
for fmt, dpi_val, kw in [("png", 600, {}), ("pdf", None, {}),
                          ("tiff", 600, {"pil_kwargs": {"compression": "tiff_lzw"}})]:
    plt.savefig(f"{fig_path}.{fmt}", dpi=dpi_val, bbox_inches='tight',
                facecolor='white', edgecolor='none', **kw)
plt.close()
print(f"[OK] Figure saved: {fig_path}.{{png, pdf, tiff}}")

# ============================================================
# 8. FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print(f"FINAL RESULTS — AGER Virtual Knockout")
print("=" * 70)

summary_lines = [f"\n{TARGET_GENE} Expression:"]
for ct in ["Mono_C", "Mono_NC", "CD4_NC", "CD8_NC"]:
    if ct in cell_type_stats:
        s = cell_type_stats[ct]
        summary_lines.append(f"  {ct:12s}: mean={s['mean']:.4f}, detect={s['detect']:.3f}, "
                            f"HC={s['n_hc']}, donors={s['n_donors']}")

summary_lines.append(f"\nVIRTUAL KO (20 iterations):")
for ct, ko in ko_results.items():
    n35 = (ko["frequency"] >= 0.35).sum()
    n50 = (ko["frequency"] >= 0.50).sum()
    top5 = ko[ko["gene"] != TARGET_GENE].head(5)["gene"].tolist()
    summary_lines.append(f"  {ct}: {n35} recurrent (freq≥0.35), {n50} stable (freq≥0.50)")
    summary_lines.append(f"       Top: {', '.join(top5)}")

summary_lines.append(f"\nscDEGs (Wilcoxon, |log2FC|>0.5):")
for ct, deg in sc_degs.items():
    summary_lines.append(f"  {ct}: strict={deg['pass_sc'].sum()}, relaxed={deg['pass_relaxed'].sum()}")

summary_lines.append(f"\nOVERLAP (KO ∩ scDEG relaxed):")
for key, ov in all_overlaps.items():
    if "relaxed" in key:
        summary_lines.append(f"  {key}: {ov['n_overlap']}/{ov['n_ko']}, OR={ov['OR']:.1f}, P={ov['p_fisher']:.2e}")

summary_lines.append(f"\nGSEA (FDR < 0.10):")
for ct, gsea in gsea_results.items():
    n_sig = (gsea["FDR"] < 0.10).sum()
    top_pathways = gsea.head(3)["pathway"].tolist() if n_sig > 0 else ["(none)"]
    summary_lines.append(f"  {ct}: {n_sig} pathways → {', '.join(top_pathways)}")

summary_lines.append(f"\nNEGATIVE CONTROLS (n={len(control_genes)}):")
for ct_label in [primary_ct, "CD4_NC"]:
    ak = f"{TARGET_GENE}_{ct_label}"
    if ak in control_metrics:
        ck = [k for k in control_metrics if k.endswith(f"_{ct_label}") and not k.startswith(TARGET_GENE)]
        if ck:
            av = control_metrics[ak]["top10_mean_abs_z"]
            cv = [control_metrics[k]["top10_mean_abs_z"] for k in ck]
            nh = sum(1 for v in cv if v >= av)
            ep = (1 + nh) / (1 + len(cv))
            summary_lines.append(f"  {ct_label}: {TARGET_GENE} top10|Z|={av:.3f}, "
                                f"control median={np.median(cv):.3f}, P={ep:.4f}")

summary_lines.append(f"\nOutput: {fig_path}.png/pdf/tiff")
summary_lines.append(f"Results: {out_dir}/")

summary_text = "\n".join(summary_lines)
print(summary_text)

# Save summary
with open(os.path.join(out_dir, "AGER_KO_summary.txt"), 'w') as f:
    f.write(summary_text)

print("=" * 70)
print("AGER virtual KO analysis complete.")
print("=" * 70)
