"""
Simplified Virtual Knockout Analysis for SESN3 in Python
Based on scTenifoldKnk methodology (Osorio et al., Patterns 2022)

Method:
1. Build gene-gene correlation network from HC CD4/CD8 T cells
2. Simulate SESN3 KO by removing all edges connected to SESN3
3. Compare original vs KO network to identify perturbed genes
4. GO enrichment of perturbed genes
5. Overlap with MG vs HC T-cell DEGs
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.sparse import issparse
import os
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=150, frameon=False, facecolor='white')

data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
fig_dir = os.path.join(data_dir, "figures_meeting")
out_dir = os.path.join(data_dir, "virtual_ko_output")
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# Load processed data
print("Loading processed data...")
adata = sc.read_h5ad(os.path.join(data_dir, "adata_processed.h5ad"))
print(f"Data: {adata.n_obs} cells, {adata.n_vars} genes")

# ============================================================
# 1. Extract HC CD4/CD8 T cells
# ============================================================

def extract_cell_population(adata, cell_type, condition='HC'):
    """Extract cells of a specific type and condition."""
    mask = (adata.obs['cell_type'] == cell_type) & (adata.obs['condition'] == condition)
    return adata[mask].copy()

cd4_hc = extract_cell_population(adata, 'CD4 T cells', 'HC')
cd8_hc = extract_cell_population(adata, 'CD8 T cells', 'HC')

print(f"HC CD4 T cells: {cd4_hc.n_obs}")
print(f"HC CD8 T cells: {cd8_hc.n_obs}")

# Subsample to 2000 cells max for computational efficiency
np.random.seed(123)
for sub in [cd4_hc, cd8_hc]:
    if sub.n_obs > 2000:
        idx = np.random.choice(sub.n_obs, 2000, replace=False)
        # Use .copy() and proper indexing
        sub._inplace_subset_obs(idx)

# ============================================================
# 2. Select HVGs + SESN3
# ============================================================

def select_genes_for_ko(adata_sub, n_hvg=3000, target_gene='SESN3'):
    """Select HVGs and ensure target gene is included."""
    # Get raw counts or log-normalized
    if issparse(adata_sub.X):
        expr = adata_sub.X.toarray()
    else:
        expr = adata_sub.X

    # Use highly variable genes based on variance
    gene_vars = np.var(expr, axis=0)
    top_hvg_idx = np.argsort(gene_vars)[-n_hvg:]
    hvg_genes = adata_sub.var_names[top_hvg_idx].tolist()

    # Ensure target gene is included
    if target_gene in adata_sub.var_names and target_gene not in hvg_genes:
        hvg_genes.append(target_gene)

    # Filter to HVGs
    expr_hvg = expr[:, [adata_sub.var_names.get_loc(g) for g in hvg_genes if g in adata_sub.var_names]]

    print(f"Selected {len(hvg_genes)} genes for KO analysis")
    return expr_hvg, hvg_genes

print("\nPreparing CD4 T-cell matrix...")
expr_cd4, genes_cd4 = select_genes_for_ko(cd4_hc)
print("Preparing CD8 T-cell matrix...")
expr_cd8, genes_cd8 = select_genes_for_ko(cd8_hc)

# ============================================================
# 3. Virtual KO: correlation-based approach
# ============================================================

def virtual_knockout(expr_matrix, gene_list, target_gene='SESN3', threshold=0.3):
    """
    Simulate gene knockout by:
    1. Computing gene-gene correlation network
    2. Removing all edges to/from target gene
    3. Computing perturbation score for each gene
    """
    if target_gene not in gene_list:
        print(f"WARNING: {target_gene} not in gene list!")
        return None

    target_idx = gene_list.index(target_gene)
    n_genes, n_cells = expr_matrix.shape

    # Compute gene-gene Spearman correlation matrix (original network)
    print(f"  Computing {n_genes}x{n_genes} correlation matrix...")
    corr_original = np.corrcoef(expr_matrix)  # Pearson on log-normalized data

    # Create KO network: zero out target gene edges
    corr_ko = corr_original.copy()
    corr_ko[target_idx, :] = 0
    corr_ko[:, target_idx] = 0

    # For each gene, compute perturbation = change in its correlation profile
    # Distance measure: Euclidean distance between original and KO edge weights
    perturbation_scores = np.zeros(n_genes)

    for i in range(n_genes):
        if i == target_idx:
            perturbation_scores[i] = 0  # Target gene itself
            continue
        # Change in correlation profile for gene i
        delta = corr_original[i, :] - corr_ko[i, :]
        perturbation_scores[i] = np.sqrt(np.mean(delta ** 2))

    # Z-score normalize perturbation scores (excluding target)
    mask = np.ones(n_genes, dtype=bool)
    mask[target_idx] = False
    mean_p = perturbation_scores[mask].mean()
    std_p = perturbation_scores[mask].std()
    if std_p > 0:
        z_scores = np.zeros(n_genes)
        z_scores[mask] = (perturbation_scores[mask] - mean_p) / std_p
    else:
        z_scores = perturbation_scores

    # P-values from Z-scores (one-tailed)
    from scipy.stats import norm
    p_values = 1 - norm.cdf(z_scores)

    # Multiple testing correction (Bonferroni)
    p_adj = np.minimum(p_values * n_genes, 1.0)

    # Results
    results = pd.DataFrame({
        'gene': gene_list,
        'distance': perturbation_scores,
        'z_score': z_scores,
        'p_value': p_values,
        'p_adj': p_adj
    })

    # Original correlation with target gene
    results['corr_with_target'] = corr_original[target_idx, :]

    results = results.sort_values('p_adj')
    return results, corr_original, corr_ko

print("\nRunning SESN3 virtual KO in CD4 T cells...")
res_cd4, corr_orig_cd4, corr_ko_cd4 = virtual_knockout(expr_cd4, genes_cd4, 'SESN3')

print("Running SESN3 virtual KO in CD8 T cells...")
res_cd8, corr_orig_cd8, corr_ko_cd8 = virtual_knockout(expr_cd8, genes_cd8, 'SESN3')

# ============================================================
# 4. Negative controls
# ============================================================

# Pick 2 control genes with similar mean expression to SESN3
sesn3_mean_cd4 = expr_cd4[:, genes_cd4.index('SESN3')].mean() if 'SESN3' in genes_cd4 else 0
gene_means_cd4 = expr_cd4.mean(axis=0)
# Find genes within 20% of SESN3 expression
nearby_mask = np.abs(gene_means_cd4 - sesn3_mean_cd4) / max(sesn3_mean_cd4, 0.01) < 0.2
known_hits = ['SESN3', 'AGER', 'ZSCAN12', 'CFD', 'DXO', 'ZBTB9',
              'HLA-DOA', 'BTN2A1', 'HCG23', 'TYMP', 'TNFRSF14']
nearby_genes = [g for g, m in zip(genes_cd4, nearby_mask) if m and g not in known_hits]
np.random.seed(42)
control_genes = list(np.random.choice(nearby_genes, min(2, len(nearby_genes)), replace=False))
print(f"\nNegative control genes: {control_genes}")

control_results = {}
for cg in control_genes:
    print(f"Running control KO: {cg} in CD4 T cells...")
    ctrl_res, _, _ = virtual_knockout(expr_cd4, genes_cd4, cg)
    if ctrl_res is not None:
        control_results[cg] = ctrl_res

# ============================================================
# 5. MG vs HC DEG analysis in T cells
# ============================================================

def find_degs(adata, cell_type):
    """Find DEGs between MG and HC for a given cell type."""
    sub = adata[adata.obs['cell_type'] == cell_type].copy()

    hc_mask = sub.obs['condition'] == 'HC'
    mg_mask = sub.obs['condition'] == 'MG'

    if hc_mask.sum() < 3 or mg_mask.sum() < 3:
        return None

    if issparse(sub.X):
        expr = sub.X.toarray()
    else:
        expr = sub.X

    hc_expr = expr[hc_mask.values, :]
    mg_expr = expr[mg_mask.values, :]

    deg_results = []
    for i, gene in enumerate(sub.var_names):
        stat, p = stats.mannwhitneyu(mg_expr[:, i], hc_expr[:, i], alternative='two-sided')
        fc = np.mean(mg_expr[:, i]) - np.mean(hc_expr[:, i])  # logFC

        deg_results.append({
            'gene': gene,
            'logFC': fc,
            'p_value': p,
            'mean_MG': np.mean(mg_expr[:, i]),
            'mean_HC': np.mean(hc_expr[:, i])
        })

    deg_df = pd.DataFrame(deg_results)
    deg_df['p_adj'] = np.minimum(deg_df['p_value'] * len(deg_df), 1.0)
    deg_df = deg_df.sort_values('p_adj')
    return deg_df

print("\nComputing MG vs HC DEGs in T cells...")
degs_cd4 = find_degs(adata, 'CD4 T cells')
degs_cd8 = find_degs(adata, 'CD8 T cells')

# ============================================================
# 6. Overlap analysis: SESN3-KO perturbed genes vs MG DEGs
# ============================================================

def overlap_analysis(ko_results, degs, label, n_top=100):
    """Test overlap between top KO-perturbed genes and MG DEGs."""
    top_ko = ko_results[ko_results['gene'] != 'SESN3'].head(n_top)['gene'].tolist()
    sig_degs = degs[degs['p_adj'] < 0.05]['gene'].tolist() if degs is not None else []

    overlap = set(top_ko) & set(sig_degs)
    print(f"\n{label}:")
    print(f"  Top {n_top} KO-perturbed: {len(top_ko)} genes")
    print(f"  Significant DEGs: {len(sig_degs)}")
    print(f"  Overlap: {len(overlap)} genes")

    # Hypergeometric test
    n_total = len(set(ko_results['gene']) | set(degs['gene']))
    n_overlap = len(overlap)
    n_ko = len(top_ko)
    n_deg = len(sig_degs)
    p_hyper = stats.hypergeom.sf(n_overlap - 1, n_total, n_ko, n_deg) if n_total > 0 else 1

    print(f"  Hypergeometric P = {p_hyper:.4f}")

    return list(overlap), p_hyper

overlap_cd4, p_cd4 = overlap_analysis(res_cd4, degs_cd4, "CD4 T cells")
overlap_cd8, p_cd8 = overlap_analysis(res_cd8, degs_cd8, "CD8 T cells")

# ============================================================
# 7. GO enrichment of perturbed genes
# ============================================================

def go_enrichment_simple(gene_list, background_genes, label):
    """
    Simple GO enrichment using MSigDB hallmark-like keyword matching.
    (Full GO requires enrichr API or local database; here we use keyword enrichment)
    """
    # Define pathway gene sets (MSigDB Hallmark-inspired)
    pathways = {
        'Oxidative Stress / ROS': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4', 'PRDX1', 'PRDX2', 'PRDX3',
                                    'PRDX5', 'PRDX6', 'TXN', 'TXNRD1', 'NFE2L2', 'KEAP1', 'HMOX1',
                                    'GCLC', 'GCLM', 'GSR', 'NOS2', 'CYBB', 'NOX4', 'MPO'],
        'T-cell Activation': ['CD3D', 'CD3E', 'CD3G', 'CD28', 'LCK', 'ZAP70', 'LAT', 'ITK',
                              'NFATC1', 'NFATC2', 'NFKB1', 'NFKBIA', 'JUN', 'FOS', 'CD69',
                              'IL2RA', 'IL2RB', 'IL2RG', 'PRKCQ', 'CARD11', 'MALT1'],
        'mTOR Signaling': ['MTOR', 'RPTOR', 'RICTOR', 'AKT1', 'AKT2', 'TSC1', 'TSC2', 'RHEB',
                           'EIF4EBP1', 'RPS6KB1', 'RPS6', 'ULK1', 'ATG13', 'DEPTOR', 'LAMTOR1'],
        'Inflammatory Response': ['TNF', 'IL1B', 'IL6', 'IL18', 'CCL2', 'CCL3', 'CCL4', 'CCL5',
                                   'CXCL8', 'CXCL10', 'NFKB1', 'RELA', 'STAT1', 'STAT3', 'IRF1',
                                   'IRF7', 'TLR2', 'TLR4', 'MYD88'],
        'Apoptosis': ['BCL2', 'BCL2L1', 'BAX', 'BAK1', 'CASP3', 'CASP8', 'CASP9', 'BID',
                      'CYCS', 'APAF1', 'XIAP', 'BIRC5', 'TP53', 'PMAIP1', 'BBC3'],
        'Interferon Response': ['IFNG', 'IFNGR1', 'STAT1', 'IRF1', 'IRF9', 'MX1', 'MX2', 'OAS1',
                                'OAS2', 'IFIT1', 'IFIT3', 'ISG15', 'GBP1', 'GBP2'],
        'Autophagy': ['ULK1', 'ULK2', 'ATG3', 'ATG5', 'ATG7', 'ATG12', 'ATG13', 'ATG16L1',
                      'BECN1', 'SQSTM1', 'MAP1LC3A', 'MAP1LC3B', 'GABARAP', 'GABARAPL1'],
        'Complement': ['C1QA', 'C1QB', 'C1QC', 'C2', 'C3', 'C4A', 'C4B', 'C5', 'CFB', 'CFD',
                       'CFH', 'CFI', 'SERPING1', 'CD55', 'CD59'],
    }

    enrichment = []
    for pathway, pw_genes in pathways.items():
        overlap = set(gene_list) & set(pw_genes)
        n_overlap = len(overlap)
        n_pw = len(pw_genes)
        n_list = len(gene_list)
        n_bg = len(background_genes)

        if n_overlap >= 1:
            # Fisher's exact test
            a = n_overlap
            b = n_list - n_overlap
            c = n_pw - n_overlap
            d = n_bg - n_list - n_pw + n_overlap
            _, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative='greater')

            enrichment.append({
                'pathway': pathway,
                'n_overlap': n_overlap,
                'n_pathway': n_pw,
                'overlap_genes': ', '.join(sorted(overlap)),
                'p_value': p_fisher
            })

    enrich_df = pd.DataFrame(enrichment)
    if len(enrich_df) > 0:
        enrich_df = enrich_df.sort_values('p_value')
    return enrich_df

# Get perturbed genes (p_adj < 0.1)
perturbed_cd4 = res_cd4[res_cd4['p_adj'] < 0.1]['gene'].tolist()
perturbed_cd8 = res_cd8[res_cd8['p_adj'] < 0.1]['gene'].tolist()

print(f"\nCD4: {len(perturbed_cd4)} genes perturbed by SESN3 KO (p_adj<0.1)")
print(f"CD8: {len(perturbed_cd8)} genes perturbed by SESN3 KO (p_adj<0.1)")

enrich_cd4 = go_enrichment_simple(perturbed_cd4, genes_cd4, "CD4")
enrich_cd8 = go_enrichment_simple(perturbed_cd8, genes_cd8, "CD8")

# ============================================================
# 8. Generate final figures
# ============================================================

fig = plt.figure(figsize=(20, 14))

# ---- A: Top perturbed genes in CD4 ----
ax_a = fig.add_subplot(2, 3, 1)
top_cd4 = res_cd4[res_cd4['gene'] != 'SESN3'].head(15)
colors_a = ['#4472C4' if v > 0 else '#ED7D31' for v in top_cd4['corr_with_target']]
ax_a.barh(range(len(top_cd4)), -np.log10(top_cd4['p_adj'].values), color=colors_a)
ax_a.set_yticks(range(len(top_cd4)))
ax_a.set_yticklabels(top_cd4['gene'].values, fontsize=8)
ax_a.set_xlabel('-log10(P.adj)')
ax_a.set_title('A. CD4 T cells: Top Perturbed Genes\n(SESN3 KO)', fontweight='bold')
ax_a.invert_yaxis()
ax_a.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='P=0.05')

# ---- B: Top perturbed genes in CD8 ----
ax_b = fig.add_subplot(2, 3, 2)
top_cd8 = res_cd8[res_cd8['gene'] != 'SESN3'].head(15)
colors_b = ['#4472C4' if v > 0 else '#ED7D31' for v in top_cd8['corr_with_target']]
ax_b.barh(range(len(top_cd8)), -np.log10(top_cd8['p_adj'].values), color=colors_b)
ax_b.set_yticks(range(len(top_cd8)))
ax_b.set_yticklabels(top_cd8['gene'].values, fontsize=8)
ax_b.set_xlabel('-log10(P.adj)')
ax_b.set_title('B. CD8 T cells: Top Perturbed Genes\n(SESN3 KO)', fontweight='bold')
ax_b.invert_yaxis()
ax_b.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)

# ---- C: GO enrichment ----
ax_c = fig.add_subplot(2, 3, 3)
all_enrich = pd.concat([enrich_cd4.assign(cell_type='CD4'), enrich_cd8.assign(cell_type='CD8')])
if len(all_enrich) > 0:
    all_enrich = all_enrich.sort_values('p_value').head(10)
    colors_c = ['#4472C4' if ct == 'CD4' else '#ED7D31' for ct in all_enrich['cell_type']]
    ax_c.barh(range(len(all_enrich)), -np.log10(all_enrich['p_value'].values), color=colors_c)
    ax_c.set_yticks(range(len(all_enrich)))
    labels_c = [f"{r['pathway']} ({r['cell_type']})" for _, r in all_enrich.iterrows()]
    ax_c.set_yticklabels(labels_c, fontsize=8)
    ax_c.set_xlabel('-log10(P)')
    ax_c.set_title('C. Pathway Enrichment of\nPerturbed Genes', fontweight='bold')
    ax_c.invert_yaxis()
    ax_c.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)

# ---- D: Negative control comparison ----
ax_d = fig.add_subplot(2, 3, 4)
n_sesn3_cd4 = (res_cd4['p_adj'] < 0.05).sum()
comparison_data = {'SESN3\n(CD4)': n_sesn3_cd4}
for cg, cres in control_results.items():
    comparison_data[f'{cg}\n(CD4)'] = (cres['p_adj'] < 0.05).sum()

bars = ax_d.bar(range(len(comparison_data)), list(comparison_data.values()),
                color=['#D62728'] + ['#AAAAAA'] * len(control_results))
ax_d.set_xticks(range(len(comparison_data)))
ax_d.set_xticklabels(list(comparison_data.keys()), fontsize=9)
ax_d.set_ylabel('# Perturbed Genes (P<0.05)')
ax_d.set_title('D. Negative Control Comparison', fontweight='bold')

# ---- E: MG DEG overlap ----
ax_e = fig.add_subplot(2, 3, 5)
overlap_data = {
    'CD4: KO ∩ DEG': len(overlap_cd4),
    'CD8: KO ∩ DEG': len(overlap_cd8),
}
ax_e.bar(range(len(overlap_data)), list(overlap_data.values()),
         color=['#4472C4', '#ED7D31'])
ax_e.set_xticks(range(len(overlap_data)))
ax_e.set_xticklabels(list(overlap_data.keys()), fontsize=10)
ax_e.set_ylabel('# Genes')
ax_e.set_title(f'E. Overlap with MG DEGs\nCD4 P={p_cd4:.3f} | CD8 P={p_cd8:.3f}', fontweight='bold')

# ---- F: SESN3 correlation network summary ----
ax_f = fig.add_subplot(2, 3, 6)
# Show top correlated genes with SESN3
corr_with_sesn3_cd4 = pd.DataFrame({
    'gene': genes_cd4,
    'corr': corr_orig_cd4[genes_cd4.index('SESN3'), :]
}).sort_values('corr', ascending=False)

top_pos = corr_with_sesn3_cd4[corr_with_sesn3_cd4['gene'] != 'SESN3'].head(10)
top_neg = corr_with_sesn3_cd4[corr_with_sesn3_cd4['gene'] != 'SESN3'].tail(10)
top_corr = pd.concat([top_pos, top_neg])

colors_f = ['#D62728' if v > 0 else '#1F77B4' for v in top_corr['corr']]
ax_f.barh(range(len(top_corr)), top_corr['corr'].values, color=colors_f)
ax_f.set_yticks(range(len(top_corr)))
ax_f.set_yticklabels(top_corr['gene'].values, fontsize=8)
ax_f.set_xlabel('Correlation with SESN3')
ax_f.set_title('F. SESN3 Co-expression Network\n(CD4 T cells, HC)', fontweight='bold')
ax_f.axvline(x=0, color='black', linewidth=0.5)

plt.suptitle('SESN3 Virtual Knockout Analysis — scTenifoldKnk-style in-silico Perturbation',
             fontweight='bold', fontsize=16, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, '08_virtual_KO_results.png'), dpi=200, bbox_inches='tight')
plt.close()
print("\n[OK] Virtual KO results figure saved.")

# ============================================================
# 9. Summary statistics
# ============================================================

print(f"\n{'='*60}")
print("SESN3 Virtual KO — Summary")
print(f"{'='*60}")
print(f"CD4 T cells — {len(perturbed_cd4)} perturbed genes (p_adj<0.1)")
print(f"  Top hit: {top_cd4.iloc[0]['gene']} (P={top_cd4.iloc[0]['p_adj']:.2e})")
print(f"CD8 T cells — {len(perturbed_cd8)} perturbed genes (p_adj<0.1)")
print(f"  Top hit: {top_cd8.iloc[0]['gene']} (P={top_cd8.iloc[0]['p_adj']:.2e})")
print(f"MG DEG overlap CD4: {len(overlap_cd4)} genes (P={p_cd4:.4f})")
print(f"MG DEG overlap CD8: {len(overlap_cd8)} genes (P={p_cd8:.4f})")
print(f"\nNegative controls: SESN3 KO induces {n_sesn3_cd4} perturbations vs")
for cg, cres in control_results.items():
    print(f"  {cg} KO: {(cres['p_adj'] < 0.05).sum()} perturbations")

# Save all results
res_cd4.to_csv(os.path.join(out_dir, "SESN3_KO_CD4_results.csv"), index=False)
res_cd8.to_csv(os.path.join(out_dir, "SESN3_KO_CD8_results.csv"), index=False)
if enrich_cd4 is not None:
    enrich_cd4.to_csv(os.path.join(out_dir, "SESN3_KO_CD4_enrichment.csv"), index=False)
if enrich_cd8 is not None:
    enrich_cd8.to_csv(os.path.join(out_dir, "SESN3_KO_CD8_enrichment.csv"), index=False)
if degs_cd4 is not None:
    degs_cd4.to_csv(os.path.join(out_dir, "MG_vs_HC_CD4_DEGs.csv"), index=False)
if degs_cd8 is not None:
    degs_cd8.to_csv(os.path.join(out_dir, "MG_vs_HC_CD8_DEGs.csv"), index=False)

print(f"\nAll results saved to: {out_dir}")
print("Analysis complete!")
