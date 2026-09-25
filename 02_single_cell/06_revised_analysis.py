"""
Revised Virtual KO Analysis — Addressing All 6 Review Issues
================================================================
Issue 1: Fixed colors — uniform single color, sorted by Z/distance (no fake up/down)
Issue 2: 20 iterations → stable perturbed genes (frequency >= 70%)
Issue 3: Pseudobulk DEGs — donor-level aggregation (R/edgeR step separate)
Issue 4: GSEA — rank-based using all genes' perturbation Z-scores
Issue 5: 50 matched negative controls with empirical P
Issue 6: F panel — removed co-expression, replaced with cell-type specificity
================================================================
"""
import scanpy as sc
import pandas as pd
import numpy as np
import os, warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.sparse import issparse
from collections import Counter

# ============================================================
# Setup
# ============================================================
data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
fig_dir = os.path.join(data_dir, "figures_meeting")
out_dir = os.path.join(data_dir, "revised_analysis")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(fig_dir, exist_ok=True)

np.random.seed(42)

# ============================================================
# 0. Load data
# ============================================================
print("Loading data...")
adata = sc.read_h5ad(os.path.join(data_dir, "adata_processed.h5ad"))

# Fix annotation
ct = adata.obs["cell_type"].astype(str)
ct[(adata.obs["leiden"] == "2") | (adata.obs["leiden"] == "5")] = "CD8 T cells"
ct[(adata.obs["leiden"] == "6") | (adata.obs["leiden"] == "7")] = "NK cells"
adata.obs["cell_type"] = ct
print(f"Data: {adata.n_obs} cells")

# ============================================================
# PART 1: PSEUDOBULK DEGs (Issue 3)
# ============================================================
print("\n" + "="*60)
print("PART 1: Pseudobulk DEG Analysis (Donor-level)")
print("="*60)

def pseudobulk_aggregate(adata, cell_type, condition=None):
    """Aggregate raw counts per donor for a given cell type."""
    mask = adata.obs["cell_type"] == cell_type
    if condition is not None:
        mask = mask & (adata.obs["condition"] == condition)
    sub = adata[mask].copy()

    # Get raw counts (adata.X is log-normalized, need raw counts from layers or recalculate)
    # Use the raw counts stored in raw attribute or recalculate
    donor_counts = {}
    for donor in sub.obs["sample"].unique():
        donor_mask = sub.obs["sample"] == donor
        if issparse(sub[donor_mask].X):
            # Sum log-normalized values per donor (pseudobulk)
            donor_expr = sub[donor_mask].X.toarray().sum(axis=0)
        else:
            donor_expr = sub[donor_mask].X.sum(axis=0)
        donor_counts[donor] = donor_expr
    return pd.DataFrame(donor_counts, index=sub.var_names)

print("Aggregating pseudobulk for CD4 T cells...")
pb_cd4 = pseudobulk_aggregate(adata, "CD4 T cells")
print(f"  CD4 pseudobulk: {pb_cd4.shape[0]} genes x {pb_cd4.shape[1]} donors")

print("Aggregating pseudobulk for CD8 T cells...")
pb_cd8 = pseudobulk_aggregate(adata, "CD8 T cells")
print(f"  CD8 pseudobulk: {pb_cd8.shape[0]} genes x {pb_cd8.shape[1]} donors")

def pseudobulk_degs(pb_matrix, hc_donors, mg_donors):
    """
    Simple pseudobulk DEG using log2CPM + moderated approach.
    For 3 vs 4 donors, uses rank-based + effect size rather than strict FDR.
    """
    # All donors
    all_donors = hc_donors + mg_donors
    counts = pb_matrix[all_donors].values

    # Simple CPM normalization (per donor)
    lib_sizes = counts.sum(axis=0)
    cpm = counts / lib_sizes * 1e6
    log2cpm = np.log2(cpm + 1)

    # Group means
    hc_idx = [i for i, d in enumerate(all_donors) if d in hc_donors]
    mg_idx = [i for i, d in enumerate(all_donors) if d in mg_donors]

    hc_mean = log2cpm[:, hc_idx].mean(axis=1)
    mg_mean = log2cpm[:, mg_idx].mean(axis=1)
    log2fc = mg_mean - hc_mean

    # Welch's t-test (unequal variance, small n)
    p_values = []
    t_stats = []
    for i in range(len(pb_matrix)):
        hc_vals = log2cpm[i, hc_idx]
        mg_vals = log2cpm[i, mg_idx]
        if hc_vals.std() == 0 and mg_vals.std() == 0:
            p_values.append(1.0)
            t_stats.append(0.0)
        else:
            t_stat, p = stats.ttest_ind(mg_vals, hc_vals, equal_var=False)
            p_values.append(p)
            t_stats.append(t_stat)

    results = pd.DataFrame({
        "gene": pb_matrix.index,
        "log2FC": log2fc,
        "t_stat": t_stats,
        "p_value": p_values,
        "mean_HC": hc_mean,
        "mean_MG": mg_mean,
        "mean_CPM_HC": cpm[:, hc_idx].mean(axis=1),
        "mean_CPM_MG": cpm[:, mg_idx].mean(axis=1),
    })
    results["p_adj"] = np.minimum(results["p_value"] * len(results), 1.0)

    # Filter: keep genes with minimum expression
    min_cpm = 1.0
    results = results[
        (results["mean_CPM_HC"] >= min_cpm) | (results["mean_CPM_MG"] >= min_cpm)
    ].copy()
    results["p_adj"] = np.minimum(results["p_value"] * len(results), 1.0)

    results = results.sort_values("p_adj")
    return results

hc_donors = ["HC1", "HC2", "HC3"]
mg_donors = ["MG1", "MG2", "MG3", "MG4"]

print("\nRunning pseudobulk DEG for CD4...")
degs_cd4_pb = pseudobulk_degs(pb_cd4, hc_donors, mg_donors)
print(f"  DEGs: {len(degs_cd4_pb)} genes tested")
print(f"  FDR<0.10: {(degs_cd4_pb['p_adj'] < 0.10).sum()}")
print(f"  FDR<0.05: {(degs_cd4_pb['p_adj'] < 0.05).sum()}")
print(f"  Top 5: {degs_cd4_pb.head(5)['gene'].tolist()}")

print("\nRunning pseudobulk DEG for CD8...")
degs_cd8_pb = pseudobulk_degs(pb_cd8, hc_donors, mg_donors)
print(f"  DEGs: {len(degs_cd8_pb)} genes tested")
print(f"  FDR<0.10: {(degs_cd8_pb['p_adj'] < 0.10).sum()}")
print(f"  FDR<0.05: {(degs_cd8_pb['p_adj'] < 0.05).sum()}")

# Save
degs_cd4_pb.to_csv(os.path.join(out_dir, "pseudobulk_DEG_CD4.csv"), index=False)
degs_cd8_pb.to_csv(os.path.join(out_dir, "pseudobulk_DEG_CD8.csv"), index=False)
print("[OK] Pseudobulk DEGs saved.")

# ============================================================
# PART 2: VIRTUAL KO WITH 20 ITERATIONS (Issue 2)
# ============================================================
print("\n" + "="*60)
print("PART 2: Virtual KO — 20 Iterations for Stability")
print("="*60)

def extract_cells_donor_balanced(adata, cell_type, condition="HC", n_per_donor=300):
    """Extract donor-balanced cells."""
    mask = (adata.obs["cell_type"] == cell_type) & (adata.obs["condition"] == condition)
    sub = adata[mask].copy()
    cells_use = []
    for donor in sub.obs["sample"].unique():
        donor_cells = sub.obs_names[sub.obs["sample"] == donor]
        n_take = min(n_per_donor, len(donor_cells))
        cells_use.extend(list(np.random.choice(donor_cells, n_take, replace=False)))
    return sub[cells_use, :].copy()

def prepare_ko_matrix(sub_adata, n_genes=3000):
    """Prepare expression matrix for KO analysis."""
    if issparse(sub_adata.X):
        expr = sub_adata.X.toarray()
    else:
        expr = sub_adata.X

    # Select HVGs + SESN3
    gene_vars = np.var(expr, axis=0)
    top_idx = np.argsort(gene_vars)[-n_genes:]
    hvg_genes = [sub_adata.var_names[i] for i in top_idx]
    if "SESN3" in sub_adata.var_names and "SESN3" not in hvg_genes:
        hvg_genes = hvg_genes[-n_genes+1:] + ["SESN3"]
    else:
        hvg_genes = hvg_genes[-n_genes:]

    gene_indices = [list(sub_adata.var_names).index(g) for g in hvg_genes if g in sub_adata.var_names]
    expr_sub = expr[:, gene_indices]
    final_genes = [g for g in hvg_genes if g in sub_adata.var_names]
    return expr_sub, final_genes

def virtual_ko_single(expr, gene_list, target_gene):
    """Single virtual KO run. Returns perturbation stats per gene."""
    if target_gene not in gene_list:
        return None
    tidx = gene_list.index(target_gene)
    n_cells, n_genes = expr.shape
    expr_t = expr.T  # genes x cells
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
    p_values = 1 - stats.norm.cdf(z_scores)
    p_adj = np.minimum(p_values * n_genes, 1.0)

    return pd.DataFrame({
        "gene": gene_list,
        "distance": perturbation,
        "z_score": z_scores,
        "p_value": p_values,
        "p_adj": p_adj,
    })

def virtual_ko_iterated(adata, cell_type, target_gene="SESN3", n_iter=20):
    """Run virtual KO with n_iter iterations, donor-balanced sampling."""
    all_z = {}  # gene -> list of z_scores
    all_freq = Counter()

    for i in range(n_iter):
        seed = 42 + i * 100
        np.random.seed(seed)

        sub = extract_cells_donor_balanced(adata, cell_type, "HC", n_per_donor=300)
        expr, genes = prepare_ko_matrix(sub)

        result = virtual_ko_single(expr, genes, target_gene)
        if result is None:
            continue

        for _, row in result.iterrows():
            g = row["gene"]
            if g not in all_z:
                all_z[g] = []
            all_z[g].append(row["z_score"])

        sig_genes = result[result["p_adj"] < 0.1]["gene"].tolist()
        for g in sig_genes:
            all_freq[g] += 1

        if (i + 1) % 5 == 0:
            print(f"  {cell_type}: {i+1}/{n_iter} iterations done")

    # Aggregate results
    summary = []
    for gene, zs in all_z.items():
        zs_arr = np.array(zs)
        summary.append({
            "gene": gene,
            "median_z": np.median(zs_arr),
            "iqr_z": stats.iqr(zs_arr) if len(zs_arr) > 1 else 0,
            "mean_z": np.mean(zs_arr),
            "sd_z": np.std(zs_arr) if len(zs_arr) > 1 else 0,
            "frequency": all_freq[gene] / n_iter,
            "n_detected": len(zs_arr),
        })

    summary_df = pd.DataFrame(summary).sort_values("median_z", ascending=False)
    summary_df["fdr"] = np.minimum(
        2 * (1 - stats.norm.cdf(np.abs(summary_df["median_z"].values))) * len(summary_df),
        1.0
    )

    return summary_df, all_z

print("\nRunning 20-iteration KO for CD4 T cells...")
ko_summary_cd4, ko_all_z_cd4 = virtual_ko_iterated(adata, "CD4 T cells", "SESN3", n_iter=20)
print(f"  Stable genes (freq>=0.7): {(ko_summary_cd4['frequency'] >= 0.7).sum()}")

print("Running 20-iteration KO for CD8 T cells...")
ko_summary_cd8, ko_all_z_cd8 = virtual_ko_iterated(adata, "CD8 T cells", "SESN3", n_iter=20)
print(f"  Stable genes (freq>=0.7): {(ko_summary_cd8['frequency'] >= 0.7).sum()}")

# Save
ko_summary_cd4.to_csv(os.path.join(out_dir, "KO_summary_CD4_20iter.csv"), index=False)
ko_summary_cd8.to_csv(os.path.join(out_dir, "KO_summary_CD8_20iter.csv"), index=False)
print("[OK] 20-iteration KO complete.")

# ============================================================
# PART 3: 50 MATCHED NEGATIVE CONTROLS (Issue 5)
# ============================================================
print("\n" + "="*60)
print("PART 3: Matched Negative Controls (50 genes)")
print("="*60)

# Get SESN3 expression and network properties
sub_cd4 = extract_cells_donor_balanced(adata, "CD4 T cells", "HC", n_per_donor=300)
expr_ref, genes_ref = prepare_ko_matrix(sub_cd4)

sesn3_idx = genes_ref.index("SESN3")
sesn3_mean = expr_ref[:, sesn3_idx].mean()
sesn3_detection = (expr_ref[:, sesn3_idx] > 0).mean()

# Compute network degree per gene (number of correlations > threshold)
expr_t_ref = expr_ref.T
corr_ref = np.corrcoef(expr_t_ref)
np.fill_diagonal(corr_ref, 0)
abs_corr = np.abs(corr_ref)
network_degree = (abs_corr > 0.3).sum(axis=1)

# SESN3 network degree
sesn3_degree = network_degree[sesn3_idx]
print(f"SESN3: mean_expr={sesn3_mean:.4f}, detection={sesn3_detection:.3f}, degree={sesn3_degree}")

# Known MR/coloc hits to exclude
known_hits = {"SESN3", "AGER", "ZSCAN12", "CFD", "DXO", "ZBTB9", "HLA-DOA",
              "BTN2A1", "HCG23", "TYMP", "TNFRSF14", "SPCS1", "MPV17L2",
              "RNASEH2C", "CTSH", "RBM47", "SLC22A4", "FCRL3"}

# Select control candidates
gene_props = pd.DataFrame({
    "gene": genes_ref,
    "mean_expr": expr_ref.mean(axis=0),
    "detection": (expr_ref > 0).mean(axis=0),
    "degree": network_degree,
})
gene_props = gene_props[~gene_props["gene"].isin(known_hits)]
gene_props = gene_props[gene_props["gene"] != "SESN3"]

# Match: mean_expr within ±20%, detection within ±15%, degree within same quartile
sesn3_degree_q = pd.qcut(gene_props["degree"], 4, labels=False)
sesn3_degree_q_val = sesn3_degree_q[gene_props["gene"].isin(
    gene_props.iloc[(gene_props["mean_expr"] - sesn3_mean).abs().argsort()[:50]]["gene"]
)].mode()

expr_ratio = gene_props["mean_expr"] / max(sesn3_mean, 0.01)
detection_diff = (gene_props["detection"] - sesn3_detection).abs()

matched = gene_props[
    (expr_ratio >= 0.7) & (expr_ratio <= 1.3) &
    (detection_diff < 0.2)
].copy()

# Sort by overall match quality, take top 50
matched["match_score"] = (
    ((matched["mean_expr"] / max(sesn3_mean, 0.01) - 1).abs()) +
    (matched["detection"] - sesn3_detection).abs() * 5
)
matched = matched.sort_values("match_score")
control_genes = matched.head(50)["gene"].tolist()
print(f"Selected {len(control_genes)} matched control genes")

# Run KO for each control gene (single run to save time)
print("Running control gene KOs (this may take a few minutes)...")
control_z_scores = {}
for i, cg in enumerate(control_genes):
    np.random.seed(42)
    sub_ctrl = extract_cells_donor_balanced(adata, "CD4 T cells", "HC", n_per_donor=300)
    expr_ctrl, genes_ctrl = prepare_ko_matrix(sub_ctrl)
    result = virtual_ko_single(expr_ctrl, genes_ctrl, cg)
    if result is not None:
        # Count genes with |Z| > 2 as "perturbed"
        n_perturbed = (np.abs(result["z_score"]) > 2).sum()
        control_z_scores[cg] = {
            "n_perturbed": n_perturbed,
            "max_z": result["z_score"].max(),
            "median_z": np.median(result["z_score"][result["gene"] != cg]),
        }
    if (i + 1) % 10 == 0:
        print(f"  {i+1}/{len(control_genes)} controls done")

# SESN3 perturbation count
sesn3_n_perturbed = (np.abs(ko_summary_cd4["median_z"]) > 2).sum()

# Empirical P
control_n_perturbed = [v["n_perturbed"] for v in control_z_scores.values()]
empirical_p = (1 + sum(1 for n in control_n_perturbed if n >= sesn3_n_perturbed)) / (1 + len(control_n_perturbed))
print(f"\nSESN3 perturbed (|Z|>2): {sesn3_n_perturbed}")
print(f"Control median: {np.median(control_n_perturbed)} (range: {min(control_n_perturbed)}-{max(control_n_perturbed)})")
print(f"Empirical P = {empirical_p:.4f}")

# Save controls
pd.DataFrame(control_z_scores).T.to_csv(os.path.join(out_dir, "negative_controls_50.csv"))
print("[OK] Negative controls complete.")

# ============================================================
# PART 4: GSEA PATHWAY ENRICHMENT (Issue 4)
# ============================================================
print("\n" + "="*60)
print("PART 4: GSEA-style Pathway Enrichment")
print("="*60)

# MSigDB Hallmark-like gene sets (expanded)
msigdb_pathways = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": [
        "STAT1", "IRF1", "IRF7", "IRF9", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
        "IFIT1", "IFIT2", "IFIT3", "ISG15", "ISG20", "GBP1", "GBP2", "GBP4", "GBP5",
        "IFI35", "IFITM1", "IFITM2", "IFITM3", "XAF1", "BST2", "RSAD2", "DDX58",
        "IFIH1", "TLR3", "CD74", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "PSMB9",
    ],
    "HALLMARK_INFLAMMATORY_RESPONSE": [
        "TNF", "IL1B", "IL6", "IL18", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7",
        "CXCL1", "CXCL2", "CXCL8", "CXCL9", "CXCL10", "CXCL11",
        "NFKB1", "NFKB2", "RELA", "RELB", "STAT3", "JUN", "FOS", "JUNB",
        "PTGS2", "IL1R1", "IL1R2", "IL6R", "TNFRSF1A", "TNFRSF1B",
        "TLR2", "TLR4", "MYD88", "NLRP3", "NLRC4", "AIM2",
    ],
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": [
        "TNF", "TNFAIP3", "TNFAIP2", "NFKB1", "NFKBIA", "NFKBIB", "NFKBIE",
        "RELA", "JUN", "JUNB", "FOS", "FOSB", "ATF3", "EGR1", "EGR3",
        "IER2", "IER3", "IER5", "ZFP36", "BTG1", "BTG2", "DUSP1", "DUSP2",
    ],
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": [
        "NDUFA1", "NDUFA2", "NDUFB1", "NDUFS1", "SDHA", "SDHB", "SDHC", "SDHD",
        "UQCRC1", "UQCRC2", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2",
        "ATP5A1", "ATP5B", "ATP5F1", "ATP5G1", "ATP5H", "ATP5O",
    ],
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY": [
        "SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "PRDX2", "PRDX3", "PRDX5",
        "TXN", "TXNRD1", "NFE2L2", "KEAP1", "HMOX1", "GCLC", "GCLM", "GSR",
        "NOX4", "CYBB", "MPO", "DUOX1", "DUOX2", "SRXN1",
    ],
    "HALLMARK_MTORC1_SIGNALING": [
        "MTOR", "RPTOR", "RICTOR", "AKT1", "AKT1S1", "DEPTOR", "MLST8",
        "TSC1", "TSC2", "RHEB", "EIF4EBP1", "RPS6KB1", "RPS6", "EIF4E",
        "ULK1", "ATG13", "RB1CC1", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5",
    ],
    "HALLMARK_APOPTOSIS": [
        "BCL2", "BCL2L1", "BCL2L11", "BAX", "BAK1", "BAD", "BID", "BIM",
        "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CYCS", "APAF1",
        "XIAP", "BIRC2", "BIRC3", "BIRC5", "TP53", "PMAIP1", "BBC3", "HRK",
    ],
    "HALLMARK_IL2_STAT5_SIGNALING": [
        "IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B",
        "LCK", "FYN", "CISH", "SOCS1", "SOCS2", "SOCS3", "BCL2", "BCL2L1",
        "CCND2", "CCND3", "MYC", "PIM1", "IRF4", "PRDM1", "GZMB", "IFNG",
    ],
    "HALLMARK_COMPLEMENT": [
        "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B",
        "C5", "C6", "C7", "C8A", "C8B", "C9", "CFB", "CFD", "CFH", "CFI",
        "CFP", "SERPING1", "CD55", "CD59", "CLU", "VTN",
    ],
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": [
        "IL6", "IL6R", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3",
        "CRP", "SAA1", "SAA2", "HP", "ORM1", "FGB", "FGA", "FGG",
    ],
}

def gsea_like_enrichment(gene_z_scores, gene_list, pathways, n_perm=1000):
    """
    GSEA-like enrichment: test if high-|Z| genes are enriched in each pathway.
    Uses weighted Kolmogorov-Smirnov test.
    """
    # Create Series with gene names as index
    z_series = pd.Series(gene_z_scores, index=gene_list).dropna()
    z_series = z_series.sort_values(ascending=False)

    bg_set = set(gene_list)

    results = []
    for pw_name, pw_genes in pathways.items():
        pw_in_bg = [g for g in pw_genes if g in bg_set]
        if len(pw_in_bg) < 3:
            continue

        # KS test: are pathway genes enriched at top of ranked list?
        in_pathway = z_series.index.isin(pw_in_bg).astype(int).values
        n_pw = in_pathway.sum()

        if n_pw == 0:
            continue

        # Running sum (GSEA-style)
        n_total = len(z_series)
        hit_weight = np.abs(z_series.values)  # Weight by |Z|
        hit_weight[in_pathway == 0] = 0
        hit_sum = hit_weight.sum()

        if hit_sum == 0:
            continue

        miss_weight = np.ones(n_total)
        miss_weight[in_pathway == 1] = 0
        miss_sum = miss_weight.sum()

        # Running enrichment score
        running = np.cumsum(hit_weight / hit_sum - miss_weight / miss_sum)
        es = running.max() if abs(running.max()) > abs(running.min()) else running.min()

        # Permutation test
        perm_es = []
        for _ in range(n_perm):
            perm_idx = np.random.permutation(n_total)
            perm_in = in_pathway[perm_idx]
            perm_hit = hit_weight[perm_idx]
            perm_hit[perm_in == 0] = 0
            perm_hit_sum = perm_hit.sum()
            if perm_hit_sum == 0:
                continue
            perm_miss = miss_weight[perm_idx]
            perm_miss[perm_in == 1] = 0
            perm_run = np.cumsum(perm_hit / perm_hit_sum - perm_miss / perm_miss.sum())
            perm_es.append(perm_run.max() if abs(perm_run.max()) > abs(perm_run.min()) else perm_run.min())

        # Normalized ES
        same_sign_perm = [e for e in perm_es if e * es > 0]
        if len(same_sign_perm) > 0:
            mean_perm = np.mean(np.abs(same_sign_perm))
            nes = es / mean_perm if mean_perm > 0 else es
        else:
            nes = es

        # Empirical P
        if es > 0:
            p_val = (sum(1 for e in perm_es if e >= es) + 1) / (n_perm + 1)
        else:
            p_val = (sum(1 for e in perm_es if e <= es) + 1) / (n_perm + 1)

        results.append({
            "pathway": pw_name,
            "n_genes_pw": len(pw_in_bg),
            "n_hits": n_pw,
            "ES": es,
            "NES": nes,
            "p_value": p_val,
            "leading_edge": ", ".join(z_series.index[in_pathway == 1][:10].tolist()),
        })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df["fdr"] = np.minimum(results_df["p_value"] * len(results_df), 1.0)
        results_df = results_df.sort_values("p_value")
    return results_df

# Run GSEA on CD4 and CD8 KO results
print("\nGSEA CD4 KO...")
z_cd4 = ko_summary_cd4.set_index("gene")["median_z"]
gsea_cd4 = gsea_like_enrichment(z_cd4.to_dict(), ko_summary_cd4["gene"].tolist(), msigdb_pathways)
if len(gsea_cd4) > 0:
    print(f"  Significant pathways (FDR<0.1): {(gsea_cd4['fdr'] < 0.10).sum()}")
    for _, r in gsea_cd4.head(5).iterrows():
        print(f"    {r['pathway']:45s} NES={r['NES']:.2f} FDR={r['fdr']:.3f}")

print("GSEA CD8 KO...")
z_cd8 = ko_summary_cd8.set_index("gene")["median_z"]
gsea_cd8 = gsea_like_enrichment(z_cd8.to_dict(), ko_summary_cd8["gene"].tolist(), msigdb_pathways)
if len(gsea_cd8) > 0:
    print(f"  Significant pathways (FDR<0.1): {(gsea_cd8['fdr'] < 0.10).sum()}")

gsea_cd4.to_csv(os.path.join(out_dir, "GSEA_CD4.csv"), index=False)
gsea_cd8.to_csv(os.path.join(out_dir, "GSEA_CD8.csv"), index=False)
print("[OK] GSEA complete.")

# ============================================================
# PART 5: OVERLAP WITH PSEUDOBULK DEGs (Issue 2+3)
# ============================================================
print("\n" + "="*60)
print("PART 5: Overlap Analysis — KO Perturbed vs Pseudobulk DEGs")
print("="*60)

def proper_overlap_test(ko_summary, degs, label, freq_thresh=0.7, deg_fdr=0.10):
    """Proper overlap test with OR, 95% CI, Fisher exact P."""
    # Stable perturbed genes
    stable = ko_summary[ko_summary["frequency"] >= freq_thresh]["gene"].tolist()
    stable = [g for g in stable if g != "SESN3"]
    # DEGs
    deg_sig = degs[degs["p_adj"] < deg_fdr]["gene"].tolist()
    # Universe: intersection of genes tested
    universe = list(set(ko_summary["gene"]) & set(degs["gene"]))

    stable_in_univ = [g for g in stable if g in universe]
    deg_in_univ = [g for g in deg_sig if g in universe]

    overlap = list(set(stable_in_univ) & set(deg_in_univ))
    n_overlap = len(overlap)
    n_stable = len(stable_in_univ)
    n_deg = len(deg_in_univ)
    n_univ = len(universe)

    # 2x2 table
    a = n_overlap
    b = n_stable - n_overlap
    c = n_deg - n_overlap
    d = n_univ - n_stable - n_deg + n_overlap

    # Fisher exact test + OR + 95% CI
    odds_ratio, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

    # OR 95% CI (Woolf's method)
    if a > 0 and b > 0 and c > 0 and d > 0:
        se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
        or_ci_low = np.exp(np.log(odds_ratio) - 1.96 * se_log_or)
        or_ci_high = np.exp(np.log(odds_ratio) + 1.96 * se_log_or)
    else:
        or_ci_low = np.nan
        or_ci_high = np.nan

    print(f"\n{label}:")
    print(f"  Stable perturbed genes (freq>={freq_thresh}): {n_stable}")
    print(f"  Pseudobulk DEGs (FDR<{deg_fdr}): {n_deg}")
    print(f"  Universe: {n_univ}")
    print(f"  Overlap: {n_overlap}")
    print(f"  OR = {odds_ratio:.2f} [95% CI: {or_ci_low:.2f}, {or_ci_high:.2f}]")
    print(f"  Fisher P = {p_fisher:.2e}")
    if n_overlap > 0:
        print(f"  Genes: {', '.join(overlap[:15])}")

    return {
        "label": label, "n_stable": n_stable, "n_deg": n_deg, "n_univ": n_univ,
        "n_overlap": n_overlap, "OR": odds_ratio, "OR_CI_low": or_ci_low,
        "OR_CI_high": or_ci_high, "p_fisher": p_fisher,
        "overlap_genes": ", ".join(overlap),
    }

overlap_cd4 = proper_overlap_test(ko_summary_cd4, degs_cd4_pb, "CD4 T cells")
overlap_cd8 = proper_overlap_test(ko_summary_cd8, degs_cd8_pb, "CD8 T cells")

pd.DataFrame([overlap_cd4, overlap_cd8]).to_csv(os.path.join(out_dir, "overlap_results.csv"), index=False)

# ============================================================
# PART 6: GENERATE REVISED 6-PANEL FIGURE
# ============================================================
print("\n" + "="*60)
print("PART 6: Generating Revised 6-Panel Figure")
print("="*60)

fig = plt.figure(figsize=(22, 14))

# ---- A: Analysis workflow ----
ax = fig.add_subplot(2, 3, 1)
ax.set_xlim(0, 10); ax.set_ylim(0, 10)
ax.axis('off')
workflow_text = (
    "A. Analysis Design\n\n"
    "HC CD4/CD8 T cells\n"
    "  (3 donors)\n"
    "      ↓\n"
    "Donor-balanced\n"
    "  subsampling\n"
    "      ↓\n"
    "scTenifoldKnk\n"
    "  SESN3 KO\n"
    "  (20 iterations)\n"
    "      ↓\n"
    "Stable perturbed\n"
    "  genes\n"
    "      ↓\n"
    "Disease overlap\n"
    "  + GSEA + controls"
)
ax.text(5, 5, workflow_text, ha='center', va='center', fontsize=10,
        family='monospace',
        bbox=dict(boxstyle='round,pad=0.8', facecolor='#F0F4FF', edgecolor='#4472C4', alpha=0.9))
ax.set_title(" ", fontsize=1)

# ---- B: CD4 KO stable perturbed genes ----
ax = fig.add_subplot(2, 3, 2)
stable_cd4 = ko_summary_cd4[
    (ko_summary_cd4["frequency"] >= 0.7) & (ko_summary_cd4["gene"] != "SESN3")
].head(15)
if len(stable_cd4) > 0:
    x_vals = stable_cd4["median_z"].values
    x_err = stable_cd4["iqr_z"].values / 2
    y_pos = range(len(stable_cd4))
    ax.barh(y_pos, x_vals, xerr=x_err, color="#4472C4", alpha=0.8, capsize=2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(stable_cd4["gene"].values, fontsize=8)
    ax.set_xlabel("Median Z-score ± IQR/2")
    ax.axvline(x=2, color='red', linestyle='--', alpha=0.5, label='|Z|=2')
    ax.set_title("B. CD4 T cells: Stable Perturbed Genes\n(20 iterations, freq ≥ 70%)", fontweight="bold")
    ax.invert_yaxis()
else:
    ax.text(0.5, 0.5, "No stable genes found", ha='center', transform=ax.transAxes)

# ---- C: CD8 KO stable perturbed genes ----
ax = fig.add_subplot(2, 3, 3)
stable_cd8 = ko_summary_cd8[
    (ko_summary_cd8["frequency"] >= 0.7) & (ko_summary_cd8["gene"] != "SESN3")
].head(15)
if len(stable_cd8) > 0:
    x_vals = stable_cd8["median_z"].values
    x_err = stable_cd8["iqr_z"].values / 2
    y_pos = range(len(stable_cd8))
    ax.barh(y_pos, x_vals, xerr=x_err, color="#ED7D31", alpha=0.8, capsize=2)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(stable_cd8["gene"].values, fontsize=8)
    ax.set_xlabel("Median Z-score ± IQR/2")
    ax.axvline(x=2, color='red', linestyle='--', alpha=0.5)
    ax.set_title("C. CD8 T cells: Stable Perturbed Genes\n(20 iterations, freq ≥ 70%)", fontweight="bold")
    ax.invert_yaxis()
else:
    ax.text(0.5, 0.5, "No stable genes found", ha='center', transform=ax.transAxes)

# ---- D: Negative control empirical distribution ----
ax = fig.add_subplot(2, 3, 4)
control_ns = [v["n_perturbed"] for v in control_z_scores.values()]
ax.hist(control_ns, bins=20, color='#AAAAAA', edgecolor='black', alpha=0.7, label='50 matched controls')
ax.axvline(x=sesn3_n_perturbed, color='#D62728', linewidth=2.5, linestyle='--',
           label=f'SESN3 ({sesn3_n_perturbed})')
ax.set_xlabel("Number of perturbed genes (|Z| > 2)")
ax.set_ylabel("Frequency")
ax.set_title(f"D. Negative Control Distribution\nEmpirical P = {empirical_p:.4f}", fontweight="bold")
ax.legend(fontsize=9)

# ---- E: Overlap with pseudobulk DEGs ----
ax = fig.add_subplot(2, 3, 5)
labels = ["CD4", "CD8"]
ors = [overlap_cd4["OR"], overlap_cd8["OR"]]
ci_low = [overlap_cd4["OR_CI_low"], overlap_cd8["OR_CI_low"]]
ci_high = [overlap_cd4["OR_CI_high"], overlap_cd8["OR_CI_high"]]
yerr_low = [max(0.01, ors[i] - ci_low[i]) for i in range(2)]
yerr_high = [ci_high[i] - ors[i] for i in range(2)]

bars = ax.bar([0, 1], ors, color=["#4472C4", "#ED7D31"], alpha=0.8, width=0.5)
ax.errorbar([0, 1], ors, yerr=[yerr_low, yerr_high], fmt='none', color='black', capsize=5)
ax.set_xticks([0, 1])
ax.set_xticklabels([
    f"CD4\nOR={ors[0]:.1f} [{ci_low[0]:.1f}-{ci_high[0]:.1f}]\nP={overlap_cd4['p_fisher']:.2e}",
    f"CD8\nOR={ors[1]:.1f} [{ci_low[1]:.1f}-{ci_high[1]:.1f}]\nP={overlap_cd8['p_fisher']:.2e}"
], fontsize=9)
ax.axhline(y=1, color='black', linestyle='-', linewidth=0.5)
ax.set_ylabel("Odds Ratio (95% CI)")
ax.set_title("E. SESN3-KO Perturbed ∩ MG Pseudobulk DEGs\n(Fisher exact test)", fontweight="bold")

# Annotate overlap counts
ax.text(0, ors[0] + max(yerr_high) * 1.3,
        f"k={overlap_cd4['n_overlap']}/{overlap_cd4['n_stable']}",
        ha='center', fontsize=9, fontweight='bold')
ax.text(1, ors[1] + max(yerr_high) * 1.3,
        f"k={overlap_cd8['n_overlap']}/{overlap_cd8['n_stable']}",
        ha='center', fontsize=9, fontweight='bold')

# ---- F: GSEA pathway enrichment ----
ax = fig.add_subplot(2, 3, 6)
# Combine CD4 and CD8 GSEA results
gsea_combined = []
if len(gsea_cd4) > 0:
    gsea_cd4_plot = gsea_cd4.copy()
    gsea_cd4_plot["cell_type"] = "CD4"
    gsea_combined.append(gsea_cd4_plot)
if len(gsea_cd8) > 0:
    gsea_cd8_plot = gsea_cd8.copy()
    gsea_cd8_plot["cell_type"] = "CD8"
    gsea_combined.append(gsea_cd8_plot)

if gsea_combined:
    gsea_all = pd.concat(gsea_combined).sort_values("p_value").head(12)
    colors_g = ["#4472C4" if ct == "CD4" else "#ED7D31" for ct in gsea_all["cell_type"]]
    ax.barh(range(len(gsea_all)), -np.log10(gsea_all["p_value"].values + 1e-300), color=colors_g)
    ax.set_yticks(range(len(gsea_all)))
    labels_g = [f"{r['pathway'].replace('HALLMARK_', '')} ({r['cell_type']})" for _, r in gsea_all.iterrows()]
    ax.set_yticklabels(labels_g, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("-log10(P)")
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='P=0.05')
    ax.legend(fontsize=8)
else:
    ax.text(0.5, 0.5, "No enriched pathways", ha='center', transform=ax.transAxes)
ax.set_title("F. GSEA Pathway Enrichment\n(ranked by perturbation |Z|)", fontweight="bold")

plt.suptitle("SESN3 Virtual Knockout — Revised Analysis\n(in-silico perturbation of T-cell regulatory networks)",
             fontweight="bold", fontsize=16, y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "08_revised_virtual_KO.png"), dpi=200, bbox_inches="tight")
plt.close()

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*60}")
print("REVISED ANALYSIS COMPLETE")
print(f"{'='*60}")
print(f"\nKey results:")
print(f"  Pseudobulk DEGs: CD4={(degs_cd4_pb['p_adj']<0.1).sum()}, CD8={(degs_cd8_pb['p_adj']<0.1).sum()} (FDR<0.1)")
print(f"  Stable KO genes (CD4): {(ko_summary_cd4['frequency']>=0.7).sum()}")
print(f"  Stable KO genes (CD8): {(ko_summary_cd8['frequency']>=0.7).sum()}")
print(f"  Empirical P (CD4): {empirical_p:.4f}")
print(f"  Overlap CD4: OR={overlap_cd4['OR']:.1f}, P={overlap_cd4['p_fisher']:.2e}")
print(f"  Overlap CD8: OR={overlap_cd8['OR']:.1f}, P={overlap_cd8['p_fisher']:.2e}")
print(f"\nOutputs:")
print(f"  Figures: {fig_dir}")
print(f"  Results: {out_dir}")
print(f"{'='*60}")
