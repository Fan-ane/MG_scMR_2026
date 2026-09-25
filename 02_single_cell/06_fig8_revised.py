"""
Fig 8 Revised — SESN3 Virtual KO final figure
Changes per advisor:
  1. No main title (suptitle removed)
  2. Panel labels: A-F only, no descriptive text
  3. Larger fonts for readability
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

data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
fig_dir = os.path.join(data_dir, "figures_meeting")
out_dir = os.path.join(data_dir, "virtual_ko_output")
os.makedirs(out_dir, exist_ok=True)

adata = sc.read_h5ad(os.path.join(data_dir, "adata_processed.h5ad"))
print(f"Data: {adata.n_obs} cells, {adata.n_vars} genes")

# Fix cell type annotation
ct = adata.obs["cell_type"].astype(str)
ct[(adata.obs["leiden"] == "2") | (adata.obs["leiden"] == "5")] = "CD8 T cells"
ct[(adata.obs["leiden"] == "6") | (adata.obs["leiden"] == "7")] = "NK cells"
adata.obs["cell_type"] = ct

# ── Virtual KO function ──────────────────────────────────────────────────────
def virtual_knockout(expr_cells_x_genes, gene_list, target_gene):
    if target_gene not in gene_list:
        return None
    target_idx = gene_list.index(target_gene)
    n_cells, n_genes = expr_cells_x_genes.shape
    expr_t = expr_cells_x_genes.T
    print(f"  Computing {n_genes}x{n_genes} correlation matrix...")
    corr_orig = np.corrcoef(expr_t)
    corr_ko = corr_orig.copy()
    corr_ko[target_idx, :] = 0
    corr_ko[:, target_idx] = 0
    perturbation = np.array([
        np.sqrt(np.mean((corr_orig[i, :] - corr_ko[i, :]) ** 2))
        for i in range(n_genes)
    ])
    mask = np.ones(n_genes, dtype=bool)
    mask[target_idx] = False
    mean_p = perturbation[mask].mean()
    std_p = perturbation[mask].std()
    z_scores = np.zeros(n_genes)
    if std_p > 0:
        z_scores[mask] = (perturbation[mask] - mean_p) / std_p
    p_values = 1 - stats.norm.cdf(z_scores)
    p_adj = np.minimum(p_values * n_genes, 1.0)
    results = pd.DataFrame({
        "gene": gene_list, "distance": perturbation,
        "z_score": z_scores, "p_value": p_values, "p_adj": p_adj,
        "corr_with_target": corr_orig[target_idx, :]
    }).sort_values("p_adj")
    return results, corr_orig, corr_ko

# ── Extract data ──────────────────────────────────────────────────────────────
def extract_for_ko(adata, cell_type, condition="HC", n_cells_max=2000, n_genes=3000):
    mask = (adata.obs["cell_type"] == cell_type) & (adata.obs["condition"] == condition)
    sub = adata[mask].copy()
    np.random.seed(123)
    if sub.n_obs > n_cells_max:
        idx = np.random.choice(sub.n_obs, n_cells_max, replace=False)
        sub = sub[idx, :].copy()
    if issparse(sub.X):
        expr = sub.X.toarray()
    else:
        expr = sub.X
    gene_vars = np.var(expr, axis=0)
    top_idx = np.argsort(gene_vars)[-n_genes:]
    hvg_genes = []
    for i in top_idx:
        g = sub.var_names[i]
        if g not in hvg_genes:
            hvg_genes.append(g)
    if "SESN3" in sub.var_names and "SESN3" not in hvg_genes:
        hvg_genes = hvg_genes[:n_genes-1] + ["SESN3"]
    else:
        hvg_genes = hvg_genes[:n_genes]
    gene_indices = [list(sub.var_names).index(g) for g in hvg_genes if g in sub.var_names]
    expr_sub = expr[:, gene_indices]
    final_genes = [g for g in hvg_genes if g in sub.var_names]
    sesn3_mean = expr_sub[:, final_genes.index("SESN3")].mean() if "SESN3" in final_genes else 0
    print(f"  {cell_type}: {sub.n_obs} cells, {len(final_genes)} genes, SESN3 mean={sesn3_mean:.3f}")
    return expr_sub, final_genes

print("\nExtracting HC CD4 T cells...")
expr_cd4, genes_cd4 = extract_for_ko(adata, "CD4 T cells")
print("Extracting HC CD8 T cells...")
expr_cd8, genes_cd8 = extract_for_ko(adata, "CD8 T cells")

# ── Run virtual KO ────────────────────────────────────────────────────────────
print("\n=== SESN3 Virtual KO: CD4 T cells ===")
res_cd4, corr_cd4, ko_cd4 = virtual_knockout(expr_cd4, genes_cd4, "SESN3")
print("=== SESN3 Virtual KO: CD8 T cells ===")
res_cd8, corr_cd8, ko_cd8 = virtual_knockout(expr_cd8, genes_cd8, "SESN3")

print(f"\nCD4: {(res_cd4['p_adj'] < 0.05).sum()} genes significant")
print(f"CD8: {(res_cd8['p_adj'] < 0.05).sum()} genes significant")

# ── Negative controls ─────────────────────────────────────────────────────────
sesn3_idx = genes_cd4.index("SESN3")
sesn3_mean = expr_cd4[:, sesn3_idx].mean()
gene_means = expr_cd4.mean(axis=0)
nearby = np.abs(gene_means - sesn3_mean) / max(sesn3_mean, 0.01) < 0.2
known = ["SESN3", "AGER", "ZSCAN12", "CFD", "DXO", "ZBTB9", "HLA-DOA", "BTN2A1"]
nearby_genes = [g for g, m in zip(genes_cd4, nearby) if m and g not in known]
np.random.seed(42)
ctrl_genes = list(np.random.choice(nearby_genes, min(2, len(nearby_genes)), replace=False))
ctrl_results = {}
for cg in ctrl_genes:
    cres, _, _ = virtual_knockout(expr_cd4, genes_cd4, cg)
    if cres is not None:
        ctrl_results[cg] = cres

# ── MG vs HC DEGs ─────────────────────────────────────────────────────────────
def find_degs(adata, ct):
    sub = adata[adata.obs["cell_type"] == ct].copy()
    hc = sub.obs["condition"] == "HC"
    mg = sub.obs["condition"] == "MG"
    if issparse(sub.X):
        expr = sub.X.toarray()
    else:
        expr = sub.X
    results = []
    for i, gene in enumerate(sub.var_names):
        hc_expr = expr[hc.values, i]
        mg_expr = expr[mg.values, i]
        stat, p = stats.mannwhitneyu(mg_expr, hc_expr, alternative="two-sided")
        results.append({"gene": gene, "logFC": mg_expr.mean() - hc_expr.mean(), "p_value": p})
    deg_df = pd.DataFrame(results)
    deg_df["p_adj"] = np.minimum(deg_df["p_value"] * len(deg_df), 1.0)
    return deg_df.sort_values("p_adj")

degs_cd4 = find_degs(adata, "CD4 T cells")
degs_cd8 = find_degs(adata, "CD8 T cells")

# ── Overlap ────────────────────────────────────────────────────────────────────
def overlap_test(ko_res, degs, label, n_top=100):
    top = ko_res[ko_res["gene"] != "SESN3"].head(n_top)["gene"].tolist()
    sig = degs[degs["p_adj"] < 0.05]["gene"].tolist()
    ov = list(set(top) & set(sig))
    n_total = len(set(ko_res["gene"]) | set(degs["gene"]))
    p = stats.hypergeom.sf(len(ov) - 1, max(n_total, 1), len(top), max(len(sig), 1))
    print(f"  {label}: {len(ov)} overlapping genes (P={p:.4f})")
    return ov, p

ov_cd4, p_ov_cd4 = overlap_test(res_cd4, degs_cd4, "CD4")
ov_cd8, p_ov_cd8 = overlap_test(res_cd8, degs_cd8, "CD8")

# ── Pathway enrichment ────────────────────────────────────────────────────────
pathways = {
    "Oxidative Stress": ["SOD1", "SOD2", "CAT", "GPX1", "GPX4", "PRDX1", "TXN", "TXNRD1", "NFE2L2", "KEAP1", "HMOX1", "GCLC", "GCLM", "GSR", "CYBB", "NOX4", "MPO"],
    "T-cell Activation": ["CD3D", "CD3E", "LCK", "ZAP70", "LAT", "NFATC1", "NFKB1", "JUN", "FOS", "CD69", "IL2RA", "PRKCQ"],
    "mTOR Signaling": ["MTOR", "RPTOR", "RICTOR", "AKT1", "TSC1", "TSC2", "RHEB", "EIF4EBP1", "RPS6KB1", "RPS6", "ULK1", "DEPTOR"],
    "Inflammatory Response": ["TNF", "IL1B", "IL6", "CCL2", "CCL5", "CXCL8", "CXCL10", "NFKB1", "RELA", "STAT1", "STAT3", "TLR4", "MYD88"],
    "Apoptosis": ["BCL2", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "BID", "CYCS", "XIAP", "TP53", "PMAIP1", "BBC3"],
    "Autophagy": ["ULK1", "ATG5", "ATG7", "ATG12", "BECN1", "SQSTM1", "MAP1LC3B", "GABARAPL1"],
    "IFN Response": ["IFNG", "STAT1", "IRF1", "MX1", "OAS1", "IFIT1", "IFIT3", "ISG15", "GBP1", "GBP2"],
    "Complement": ["C1QA", "C1QB", "C2", "C3", "CFB", "CFD", "CFH", "CFI", "SERPING1", "CD55", "CD59"],
}
perturbed = res_cd4[res_cd4["p_adj"] < 0.1]["gene"].tolist()
bg_set = set(genes_cd4)
enrich = []
for pw, pw_genes in pathways.items():
    ov = set(perturbed) & set(pw_genes) & bg_set
    if len(ov) >= 1:
        a = len(ov); b = len(perturbed) - a
        c = len([g for g in pw_genes if g in bg_set]) - a
        d_val = len(bg_set) - len(perturbed) - c
        _, pf = stats.fisher_exact([[a, b], [c, d_val]], alternative="greater")
        enrich.append({"pathway": pw, "n_overlap": len(ov), "p_value": pf, "genes": ", ".join(sorted(ov))})
enrich_df = pd.DataFrame(enrich).sort_values("p_value")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE — Revised per advisor feedback
# ═══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(22, 14))

PANEL_FONT = 14      # panel letter
LABEL_FONT = 12      # axis labels
TICK_FONT  = 11      # tick labels
GENE_FONT  = 10      # gene/pathway names

# ── A: CD4 top perturbed genes ────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 1)
top = res_cd4[res_cd4["gene"] != "SESN3"].head(15)
colors_a = ["#D62728" if v > 0 else "#1F77B4" for v in top["corr_with_target"]]
ax.barh(range(len(top)), -np.log10(top["p_adj"].values + 1e-300), color=colors_a)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top["gene"].values, fontsize=GENE_FONT)
ax.set_xlabel("-log10(P.adj)", fontsize=LABEL_FONT)
ax.tick_params(axis='x', labelsize=TICK_FONT)
ax.set_title("A", loc="left", fontweight="bold", fontsize=PANEL_FONT)
ax.invert_yaxis()
ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.5)

# ── B: CD8 top perturbed genes ────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 2)
top8 = res_cd8[res_cd8["gene"] != "SESN3"].head(15)
colors_b = ["#D62728" if v > 0 else "#1F77B4" for v in top8["corr_with_target"]]
ax.barh(range(len(top8)), -np.log10(top8["p_adj"].values + 1e-300), color=colors_b)
ax.set_yticks(range(len(top8)))
ax.set_yticklabels(top8["gene"].values, fontsize=GENE_FONT)
ax.set_xlabel("-log10(P.adj)", fontsize=LABEL_FONT)
ax.tick_params(axis='x', labelsize=TICK_FONT)
ax.set_title("B", loc="left", fontweight="bold", fontsize=PANEL_FONT)
ax.invert_yaxis()
ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.5)

# ── C: Pathway enrichment ─────────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 3)
if len(enrich_df) > 0:
    ax.barh(range(len(enrich_df)), -np.log10(enrich_df["p_value"].values + 1e-300), color="#4472C4")
    ax.set_yticks(range(len(enrich_df)))
    ax.set_yticklabels([f"{r['pathway']} ({r['n_overlap']})" for _, r in enrich_df.iterrows()], fontsize=GENE_FONT)
    ax.set_xlabel("-log10(P)", fontsize=LABEL_FONT)
    ax.tick_params(axis='x', labelsize=TICK_FONT)
    ax.set_title("C", loc="left", fontweight="bold", fontsize=PANEL_FONT)
    ax.invert_yaxis()
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.5)
else:
    ax.text(0.5, 0.5, "No enrichment found", ha="center", transform=ax.transAxes)
    ax.set_title("C", loc="left", fontweight="bold", fontsize=PANEL_FONT)

# ── D: Negative control ───────────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 4)
n_sesn3 = (res_cd4["p_adj"] < 0.05).sum()
comp_data = {"SESN3": n_sesn3}
for cg, cres in ctrl_results.items():
    comp_data[cg] = (cres["p_adj"] < 0.05).sum()
colors_d = ["#D62728"] + ["#AAAAAA"] * len(ctrl_results)
ax.bar(range(len(comp_data)), list(comp_data.values()), color=colors_d)
ax.set_xticks(range(len(comp_data)))
ax.set_xticklabels(list(comp_data.keys()), fontsize=TICK_FONT)
ax.set_ylabel("# Perturbed (P<0.05)", fontsize=LABEL_FONT)
ax.tick_params(axis='y', labelsize=TICK_FONT)
ax.set_title("D", loc="left", fontweight="bold", fontsize=PANEL_FONT)

# ── E: MG DEG overlap ─────────────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 5)
ax.bar([0, 1], [len(ov_cd4), len(ov_cd8)], color=["#4472C4", "#ED7D31"])
ax.set_xticks([0, 1])
ax.set_xticklabels([f"CD4 (P={p_ov_cd4:.3f})", f"CD8 (P={p_ov_cd8:.3f})"], fontsize=LABEL_FONT)
ax.set_ylabel("# Overlapping Genes", fontsize=LABEL_FONT)
ax.tick_params(axis='y', labelsize=TICK_FONT)
ax.set_title("E", loc="left", fontweight="bold", fontsize=PANEL_FONT)

# ── F: SESN3 co-expression ────────────────────────────────────────────────────
ax = fig.add_subplot(2, 3, 6)
sesn3_corr = pd.DataFrame({"gene": genes_cd4, "corr": corr_cd4[genes_cd4.index("SESN3"), :]})
non_sesn3 = sesn3_corr[sesn3_corr["gene"] != "SESN3"]
top_pos = non_sesn3.nlargest(10, "corr")
top_neg = non_sesn3.nsmallest(10, "corr")
top_corr_all = pd.concat([top_pos, top_neg])
colors_f = ["#D62728" if v > 0 else "#1F77B4" for v in top_corr_all["corr"]]
ax.barh(range(len(top_corr_all)), top_corr_all["corr"].values, color=colors_f)
ax.set_yticks(range(len(top_corr_all)))
ax.set_yticklabels(top_corr_all["gene"].values, fontsize=GENE_FONT)
ax.set_xlabel("Correlation with SESN3", fontsize=LABEL_FONT)
ax.tick_params(axis='x', labelsize=TICK_FONT)
ax.set_title("F", loc="left", fontweight="bold", fontsize=PANEL_FONT)
ax.axvline(x=0, color="black", linewidth=0.5)

# ── No suptitle, tight layout ─────────────────────────────────────────────────
plt.tight_layout(pad=2.0)
plt.savefig(os.path.join(fig_dir, "Fig8_SESN3_virtual_KO_final.pdf"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(fig_dir, "Fig8_SESN3_virtual_KO_final.png"), dpi=300, bbox_inches="tight")
plt.close()
print("\nFig 8 saved (revised).")
