"""Restore original Fig8 (with titles, small fonts, suptitle)."""
import scanpy as sc; import pandas as pd; import numpy as np; import os, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy.sparse import issparse

data_dir = os.path.dirname(os.path.abspath(__file__))
fig_dir = os.path.join(data_dir, "figures_meeting")

adata = sc.read_h5ad(os.path.join(data_dir, "adata_processed.h5ad"))
ct = adata.obs["cell_type"].astype(str)
ct[(adata.obs["leiden"]=="2")|(adata.obs["leiden"]=="5")] = "CD8 T cells"
ct[(adata.obs["leiden"]=="6")|(adata.obs["leiden"]=="7")] = "NK cells"
adata.obs["cell_type"] = ct

def virtual_knockout(expr_cells_x_genes, gene_list, target_gene):
    if target_gene not in gene_list: return None
    target_idx = gene_list.index(target_gene)
    expr_t = expr_cells_x_genes.T
    corr_orig = np.corrcoef(expr_t)
    corr_ko = corr_orig.copy()
    corr_ko[target_idx,:]=0; corr_ko[:,target_idx]=0
    n_genes = expr_cells_x_genes.shape[1]
    perturbation = np.array([np.sqrt(np.mean((corr_orig[i,:]-corr_ko[i,:])**2)) for i in range(n_genes)])
    mask = np.ones(n_genes,dtype=bool); mask[target_idx]=False
    mean_p,std_p = perturbation[mask].mean(),perturbation[mask].std()
    z_scores = np.zeros(n_genes)
    if std_p>0: z_scores[mask] = (perturbation[mask]-mean_p)/std_p
    p_values = 1-stats.norm.cdf(z_scores)
    p_adj = np.minimum(p_values*n_genes,1.0)
    return pd.DataFrame({"gene":gene_list,"distance":perturbation,"z_score":z_scores,"p_value":p_values,"p_adj":p_adj,"corr_with_target":corr_orig[target_idx,:]}).sort_values("p_adj"), corr_orig, corr_ko

def extract_for_ko(adata, cell_type, condition="HC", n_cells_max=2000, n_genes=3000):
    mask = (adata.obs["cell_type"]==cell_type)&(adata.obs["condition"]==condition)
    sub = adata[mask].copy()
    np.random.seed(123)
    if sub.n_obs>n_cells_max: sub = sub[np.random.choice(sub.n_obs,n_cells_max,replace=False),:].copy()
    expr = sub.X.toarray() if issparse(sub.X) else sub.X
    gene_vars = np.var(expr,axis=0)
    top_idx = np.argsort(gene_vars)[-n_genes:]
    hvg_genes=[]
    for i in top_idx:
        g=sub.var_names[i]
        if g not in hvg_genes: hvg_genes.append(g)
    if "SESN3" in sub.var_names and "SESN3" not in hvg_genes: hvg_genes=hvg_genes[:n_genes-1]+["SESN3"]
    else: hvg_genes=hvg_genes[:n_genes]
    gene_indices = [list(sub.var_names).index(g) for g in hvg_genes if g in sub.var_names]
    return expr[:,gene_indices], [g for g in hvg_genes if g in sub.var_names]

print("Extracting...")
expr_cd4,genes_cd4 = extract_for_ko(adata,"CD4 T cells")
expr_cd8,genes_cd8 = extract_for_ko(adata,"CD8 T cells")
print("KO...")
res_cd4,corr_cd4,ko_cd4 = virtual_knockout(expr_cd4,genes_cd4,"SESN3")
res_cd8,corr_cd8,ko_cd8 = virtual_knockout(expr_cd8,genes_cd8,"SESN3")

# Controls
sesn3_idx = genes_cd4.index("SESN3")
sesn3_mean = expr_cd4[:,sesn3_idx].mean()
nearby = np.abs(expr_cd4.mean(axis=0)-sesn3_mean)/max(sesn3_mean,0.01)<0.2
known = ["SESN3","AGER","ZSCAN12","CFD","DXO","ZBTB9","HLA-DOA","BTN2A1"]
nearby_genes = [g for g,m in zip(genes_cd4,nearby) if m and g not in known]
np.random.seed(42)
ctrl_genes=list(np.random.choice(nearby_genes,min(2,len(nearby_genes)),replace=False))
ctrl_results={}
for cg in ctrl_genes:
    cres,_,_=virtual_knockout(expr_cd4,genes_cd4,cg)
    if cres is not None: ctrl_results[cg]=cres

def find_degs(adata,ct):
    sub=adata[adata.obs["cell_type"]==ct].copy()
    expr=sub.X.toarray() if issparse(sub.X) else sub.X
    hc_expr=expr[sub.obs["condition"]=="HC"]; mg_expr=expr[sub.obs["condition"]=="MG"]
    results=[]
    for i,g in enumerate(sub.var_names):
        s,p=stats.mannwhitneyu(mg_expr[:,i],hc_expr[:,i],alternative="two-sided")
        results.append({"gene":g,"logFC":mg_expr[:,i].mean()-hc_expr[:,i].mean(),"p_value":p})
    df=pd.DataFrame(results); df["p_adj"]=np.minimum(df["p_value"]*len(df),1.0)
    return df.sort_values("p_adj")

degs_cd4=find_degs(adata,"CD4 T cells"); degs_cd8=find_degs(adata,"CD8 T cells")

def overlap_test(ko_res,degs,n_top=100):
    top=ko_res[ko_res["gene"]!="SESN3"].head(n_top)["gene"].tolist()
    sig=degs[degs["p_adj"]<0.05]["gene"].tolist()
    ov=list(set(top)&set(sig))
    n_total=len(set(ko_res["gene"])|set(degs["gene"]))
    p=stats.hypergeom.sf(len(ov)-1,max(n_total,1),len(top),max(len(sig),1))
    return ov,p

ov_cd4,p_ov_cd4=overlap_test(res_cd4,degs_cd4)
ov_cd8,p_ov_cd8=overlap_test(res_cd8,degs_cd8)

pathways={
    "Oxidative Stress":["SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","TXN","TXNRD1","NFE2L2"],
    "T-cell Activation":["CD3D","CD3E","LCK","ZAP70","LAT","NFATC1","NFKB1","JUN","FOS","CD69"],
    "mTOR Signaling":["MTOR","RPTOR","RICTOR","AKT1","TSC1","TSC2","RHEB","EIF4EBP1","RPS6KB1"],
    "Inflammatory Response":["TNF","IL1B","IL6","CCL2","CCL5","NFKB1","RELA","STAT1","STAT3"],
    "Apoptosis":["BCL2","BAX","BAK1","CASP3","CASP8","CASP9","BID","CYCS","XIAP","TP53"],
    "Autophagy":["ULK1","ATG5","ATG7","ATG12","BECN1","SQSTM1","MAP1LC3B"],
    "IFN Response":["IFNG","STAT1","IRF1","MX1","OAS1","IFIT1","IFIT3","ISG15","GBP1","GBP2"],
    "Complement":["C1QA","C1QB","C2","C3","CFB","CFD","CFH","CFI","SERPING1","CD55","CD59"],
}
perturbed=res_cd4[res_cd4["p_adj"]<0.1]["gene"].tolist()
bg_set=set(genes_cd4)
enrich=[]
for pw,pw_genes in pathways.items():
    ov=set(perturbed)&set(pw_genes)&bg_set
    if len(ov)>=1:
        a=len(ov); b=len(perturbed)-a; c=len([g for g in pw_genes if g in bg_set])-a; d=len(bg_set)-len(perturbed)-c
        _,pf=stats.fisher_exact([[a,b],[c,d]],alternative="greater")
        enrich.append({"pathway":pw,"n_overlap":len(ov),"p_value":pf,"genes":", ".join(sorted(ov))})
enrich_df=pd.DataFrame(enrich).sort_values("p_value")

# ORIGINAL STYLE FIGURE
fig=plt.figure(figsize=(20,14))

ax=fig.add_subplot(2,3,1)
top=res_cd4[res_cd4["gene"]!="SESN3"].head(15)
colors_a=["#D62728" if v>0 else "#1F77B4" for v in top["corr_with_target"]]
ax.barh(range(len(top)),-np.log10(top["p_adj"].values+1e-300),color=colors_a)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top["gene"].values,fontsize=8)
ax.set_xlabel("-log10(P.adj)"); ax.set_title("A. CD4 T cells: Perturbed by SESN3 KO",fontweight="bold")
ax.invert_yaxis(); ax.axvline(x=-np.log10(0.05),color="red",linestyle="--",alpha=0.5)

ax=fig.add_subplot(2,3,2)
top8=res_cd8[res_cd8["gene"]!="SESN3"].head(15)
colors_b=["#D62728" if v>0 else "#1F77B4" for v in top8["corr_with_target"]]
ax.barh(range(len(top8)),-np.log10(top8["p_adj"].values+1e-300),color=colors_b)
ax.set_yticks(range(len(top8))); ax.set_yticklabels(top8["gene"].values,fontsize=8)
ax.set_xlabel("-log10(P.adj)"); ax.set_title("B. CD8 T cells: Perturbed by SESN3 KO",fontweight="bold")
ax.invert_yaxis(); ax.axvline(x=-np.log10(0.05),color="red",linestyle="--",alpha=0.5)

ax=fig.add_subplot(2,3,3)
if len(enrich_df)>0:
    ax.barh(range(len(enrich_df)),-np.log10(enrich_df["p_value"].values+1e-300),color="#4472C4")
    ax.set_yticks(range(len(enrich_df)))
    ax.set_yticklabels([f"{r['pathway']} ({r['n_overlap']})" for _,r in enrich_df.iterrows()],fontsize=8)
    ax.set_xlabel("-log10(P)"); ax.set_title("C. Pathway Enrichment (CD4 perturbed)",fontweight="bold")
    ax.invert_yaxis(); ax.axvline(x=-np.log10(0.05),color="red",linestyle="--",alpha=0.5)

n_sesn3=(res_cd4["p_adj"]<0.05).sum()
comp_data={"SESN3":n_sesn3}
for cg,cres in ctrl_results.items(): comp_data[cg]=(cres["p_adj"]<0.05).sum()

ax=fig.add_subplot(2,3,4)
ax.bar(range(len(comp_data)),list(comp_data.values()),color=["#D62728"]+["#AAAAAA"]*len(ctrl_results))
ax.set_xticks(range(len(comp_data))); ax.set_xticklabels(list(comp_data.keys()),fontsize=9)
ax.set_ylabel("# Perturbed (P<0.05)"); ax.set_title("D. Negative Control Comparison",fontweight="bold")

ax=fig.add_subplot(2,3,5)
ax.bar([0,1],[len(ov_cd4),len(ov_cd8)],color=["#4472C4","#ED7D31"])
ax.set_xticks([0,1]); ax.set_xticklabels([f"CD4 (P={p_ov_cd4:.3f})",f"CD8 (P={p_ov_cd8:.3f})"],fontsize=10)
ax.set_ylabel("# Overlapping Genes"); ax.set_title("E. SESN3-KO Perturbed ∩ MG DEGs",fontweight="bold")

ax=fig.add_subplot(2,3,6)
sesn3_corr=pd.DataFrame({"gene":genes_cd4,"corr":corr_cd4[genes_cd4.index("SESN3"),:]})
non_sesn3=sesn3_corr[sesn3_corr["gene"]!="SESN3"]
top_corr_all=pd.concat([non_sesn3.nlargest(10,"corr"),non_sesn3.nsmallest(10,"corr")])
colors_f=["#D62728" if v>0 else "#1F77B4" for v in top_corr_all["corr"]]
ax.barh(range(len(top_corr_all)),top_corr_all["corr"].values,color=colors_f)
ax.set_yticks(range(len(top_corr_all))); ax.set_yticklabels(top_corr_all["gene"].values,fontsize=8)
ax.set_xlabel("Correlation with SESN3"); ax.set_title("F. SESN3 Co-expression (CD4 T cells, HC)",fontweight="bold")
ax.axvline(x=0,color="black",linewidth=0.5)

plt.suptitle("SESN3 Virtual Knockout — In-silico Gene Perturbation in T-cell Regulatory Networks",fontweight="bold",fontsize=16)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,"Fig8_SESN3_virtual_KO_final.pdf"),dpi=300,bbox_inches="tight")
plt.savefig(os.path.join(fig_dir,"Fig8_SESN3_virtual_KO_final.png"),dpi=300,bbox_inches="tight")
plt.close()
print("Original Fig8 restored!")
