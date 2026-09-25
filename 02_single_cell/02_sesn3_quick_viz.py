"""
SESN3 Quick Visualization + Virtual KO for Advisor Meeting (Fixed)
"""
import scanpy as sc
import scipy.io.mmio
import anndata
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats, sparse
import os, warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=150, frameon=False, facecolor='white')

data_dir = r"c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis"
fig_dir = os.path.join(data_dir, "figures_meeting")
out_dir = os.path.join(data_dir, "virtual_ko_output")
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# 1. Custom 10X loader
# ============================================================
def read_10x_mtx_custom(path):
    mtx_file = os.path.join(path, 'matrix.mtx')
    barcode_file = os.path.join(path, 'barcodes.tsv')
    gene_file = os.path.join(path, 'genes.tsv')

    mtx = scipy.io.mmio.mmread(mtx_file).tocsr().T
    barcodes = pd.read_csv(barcode_file, sep='\t', header=None)[0].values
    genes_df = pd.read_csv(gene_file, sep='\t', header=None)
    genes_df.columns = ['gene_id', 'gene_name']
    adata = anndata.AnnData(X=mtx)
    adata.obs_names = barcodes
    adata.obs_names_make_unique()
    adata.var_names = genes_df['gene_name'].values
    adata.var['gene_id'] = genes_df['gene_id'].values
    adata.var_names_make_unique()
    return adata

# ============================================================
# 2. Load and find common genes
# ============================================================
samples_info = {
    "HC1": "H1_matrix", "HC2": "H2_matrix", "HC3": "H3_matrix",
    "MG1": "M0707-01_matrix_10X", "MG2": "M0709-02_matrix_10X",
    "MG3": "M0709-03_matrix_10X", "MG4": "M210715-01_matrix_10X",
}

adatas_raw = {}
all_gene_sets = []
for name, dirname in samples_info.items():
    path = os.path.join(data_dir, dirname)
    print(f"Loading {name}...")
    adata = read_10x_mtx_custom(path)
    adata.obs['sample'] = name
    adata.obs['condition'] = 'HC' if name.startswith('HC') else 'MG'
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adatas_raw[name] = adata
    all_gene_sets.append(set(adata.var_names))
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

# Find common genes
common_genes = sorted(set.intersection(*all_gene_sets))
print(f"\nCommon genes across ALL samples: {len(common_genes)}")
print(f"SESN3 in common genes: {'SESN3' in common_genes}")

# Filter each sample to common genes
adatas = {}
for name, adata in adatas_raw.items():
    gene_mask = [g in common_genes for g in adata.var_names]
    adata_filtered = adata[:, gene_mask].copy()
    adatas[name] = adata_filtered
    print(f"  {name}: {adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes (filtered)")

# ============================================================
# 3. QC summary
# ============================================================
qc_df = pd.DataFrame([
    {"Sample": name, "nCells": a.n_obs,
     "nGenes_median": a.obs['n_genes_by_counts'].median(),
     "pct_mt_median": a.obs['pct_counts_mt'].median(),
     "Condition": a.obs['condition'].iloc[0]}
    for name, a in adatas.items()
])

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
colors = ['#4472C4' if c == 'HC' else '#ED7D31' for c in qc_df['Condition']]
axes[0].bar(qc_df['Sample'], qc_df['nCells'], color=colors)
axes[0].set_title('Cell Count per Sample', fontweight='bold')
for i, (x, y) in enumerate(zip(range(len(qc_df)), qc_df['nCells'])):
    axes[0].text(x, y + max(qc_df['nCells'])*0.02, str(y), ha='center', fontsize=8)
axes[1].bar(qc_df['Sample'], qc_df['nGenes_median'], color=colors)
axes[1].set_title('Median Genes per Cell', fontweight='bold')
axes[2].bar(qc_df['Sample'], qc_df['pct_mt_median'], color=colors)
axes[2].set_title('Median MT% per Cell', fontweight='bold')
axes[2].axhline(y=10, color='red', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, '01_QC_summary.png'), dpi=200, bbox_inches='tight')
plt.close()
print("\n[OK] QC figure saved.")

# ============================================================
# 4. Concatenate
# ============================================================
print("\nConcatenating...")
adata_all = anndata.concat(
    [adatas['HC1'], adatas['HC2'], adatas['HC3'],
     adatas['MG1'], adatas['MG2'], adatas['MG3'], adatas['MG4']],
    join='inner',
    label='sample',
    keys=['HC1', 'HC2', 'HC3', 'MG1', 'MG2', 'MG3', 'MG4']
)
print(f"Combined: {adata_all.n_obs} cells, {adata_all.n_vars} genes")

# QC filter
sc.pp.filter_cells(adata_all, min_genes=200)
sc.pp.filter_cells(adata_all, max_genes=5000)
adata_all = adata_all[adata_all.obs['pct_counts_mt'] < 10, :].copy()
print(f"After QC: {adata_all.n_obs} cells, {adata_all.n_vars} genes")

# ============================================================
# 5. Normalize, HVG, PCA, Harmony, UMAP
# ============================================================
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, n_top_genes=3000, batch_key='sample')
sc.tl.pca(adata_all, n_comps=50, svd_solver='arpack')
print("Using BBKNN for batch correction...")
sc.pp.neighbors(adata_all, n_pcs=30)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, resolution=0.8)
print(f"Clusters: {adata_all.obs['leiden'].nunique()}")

# ============================================================
# 6. Cell type annotation
# ============================================================
marker_genes = {
    'CD4 T cells': ['CD3D', 'CD4', 'IL7R'],
    'CD8 T cells': ['CD3D', 'CD8A', 'CD8B'],
    'NK cells': ['NKG7', 'GNLY', 'KLRD1'],
    'B cells': ['MS4A1', 'CD79A', 'CD79B'],
    'CD14 Monocytes': ['CD14', 'LYZ', 'S100A9'],
    'FCGR3A Monocytes': ['FCGR3A', 'MS4A7'],
    'Dendritic cells': ['FCER1A', 'CST3', 'CLEC10A'],
    'Plasma cells': ['MZB1', 'SDC1'],
}

for ct, genes in marker_genes.items():
    valid = [g for g in genes if g in adata_all.var_names]
    if valid:
        sc.tl.score_genes(adata_all, gene_list=valid, score_name=f'score_{ct.replace(" ", "_")}')

ct_names = list(marker_genes.keys())
cluster_ct = {}
for cl in adata_all.obs['leiden'].unique():
    cl_mask = adata_all.obs['leiden'] == cl
    mean_scores = {}
    for ct in ct_names:
        col = f'score_{ct.replace(" ", "_")}'
        if col in adata_all.obs.columns:
            mean_scores[ct] = adata_all.obs.loc[cl_mask, col].mean()
    if mean_scores:
        cluster_ct[cl] = max(mean_scores, key=mean_scores.get)

adata_all.obs['cell_type'] = adata_all.obs['leiden'].map(cluster_ct).fillna('Unknown')
print("Cell types:")
for ct in sorted(adata_all.obs['cell_type'].unique()):
    n = (adata_all.obs['cell_type'] == ct).sum()
    print(f"  {ct}: {n} ({n/adata_all.n_obs*100:.1f}%)")

# ============================================================
# 7. UMAP figures
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(20, 5))
sc.pl.umap(adata_all, color='condition', ax=axes[0], title='Condition', show=False)
sc.pl.umap(adata_all, color='cell_type', ax=axes[1], title='Cell Types', show=False, legend_loc='right margin')
sc.pl.umap(adata_all, color='sample', ax=axes[2], title='Sample', show=False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, '02_UMAP_overview.png'), dpi=200, bbox_inches='tight')
plt.close()
print("[OK] UMAP saved.")

# ============================================================
# 8. SESN3 expression
# ============================================================
if 'SESN3' in adata_all.var_names:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(adata_all, color='SESN3', ax=axes[0], title='SESN3 Expression (UMAP)',
               cmap='Reds', vmax='p99', show=False)
    sc.pl.violin(adata_all, 'SESN3', groupby='cell_type', rotation=45, ax=axes[1], show=False)
    axes[1].set_title('SESN3 by Cell Type', fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, '03_SESN3_expression.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] SESN3 expression figure saved.")

    # HC vs MG per cell type
    ct_list = sorted([ct for ct in adata_all.obs['cell_type'].unique() if ct != 'Unknown'])
    ncols = min(4, len(ct_list))
    nrows = (len(ct_list) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3.5, nrows*3))
    axes = axes.flatten() if len(ct_list) > 1 else [axes]

    for i, ct in enumerate(ct_list):
        ct_data = adata_all[adata_all.obs['cell_type'] == ct]
        hc_mask = ct_data.obs['condition'] == 'HC'
        mg_mask = ct_data.obs['condition'] == 'MG'
        hc_vals = ct_data[hc_mask, 'SESN3'].X.toarray().flatten() if hc_mask.sum() > 0 else np.array([])
        mg_vals = ct_data[mg_mask, 'SESN3'].X.toarray().flatten() if mg_mask.sum() > 0 else np.array([])
        hc_pct = (hc_vals > 0).mean() * 100 if len(hc_vals) > 0 else 0
        mg_pct = (mg_vals > 0).mean() * 100 if len(mg_vals) > 0 else 0

        data_plot = [hc_vals[hc_vals > 0], mg_vals[mg_vals > 0]]
        bp = axes[i].boxplot(data_plot, positions=[1, 2], widths=0.5,
                              patch_artist=True, showfliers=False)
        bp['boxes'][0].set_facecolor('#4472C4'); bp['boxes'][0].set_alpha(0.7)
        bp['boxes'][1].set_facecolor('#ED7D31'); bp['boxes'][1].set_alpha(0.7)
        axes[i].set_xticks([1, 2]); axes[i].set_xticklabels(['HC', 'MG'])
        axes[i].set_title(ct, fontsize=10, fontweight='bold')
        axes[i].set_ylabel('log(Expr)')
        axes[i].text(0.5, 0.95, f'Det: HC={hc_pct:.1f}% MG={mg_pct:.1f}%',
                    ha='center', fontsize=7, transform=axes[i].transAxes, va='top')

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    plt.suptitle('SESN3 Expression: HC vs MG by Cell Type', fontweight='bold', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, '04_SESN3_HC_vs_MG_by_celltype.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] SESN3 HC vs MG figure saved.")

    # Dotplot
    genes_plot = ['SESN3', 'CD3D', 'CD4', 'CD8A', 'CD8B', 'CD14', 'NKG7', 'MS4A1']
    valid = [g for g in genes_plot if g in adata_all.var_names]
    sc.pl.dotplot(adata_all, valid, groupby='cell_type', show=False)
    plt.savefig(os.path.join(fig_dir, '05_SESN3_dotplot.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] Dotplot saved.")

# ============================================================
# 9. T-cell UMAP
# ============================================================
t_cts = [ct for ct in adata_all.obs['cell_type'].unique() if 'T cell' in ct]
if t_cts:
    t_cells = adata_all[adata_all.obs['cell_type'].isin(t_cts)].copy()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sc.pl.umap(t_cells, color='cell_type', ax=axes[0], title='T-cell Subtypes', show=False)
    if 'SESN3' in adata_all.var_names:
        sc.pl.umap(t_cells, color='SESN3', ax=axes[1], title='SESN3 in T cells',
                   cmap='Reds', vmax='p99', show=False)
    sc.pl.umap(t_cells, color='condition', ax=axes[2], title='Condition', show=False)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, '06_Tcell_UMAP_detail.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print("[OK] T-cell UMAP saved.")

# ============================================================
# 10. Save processed data
# ============================================================
adata_all.write(os.path.join(data_dir, "adata_processed.h5ad"), compression='gzip')
print(f"\n[OK] Data saved. {adata_all.n_obs} cells, {adata_all.n_vars} genes")
print("Quick viz script complete!")
