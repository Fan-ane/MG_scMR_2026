"""Prepare KO results for R visualization"""
import pandas as pd, numpy as np, os
from scipy import stats

out_r = r'c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis\for_R'
os.makedirs(out_r, exist_ok=True)

for ct in ['CD4_NC', 'CD8_NC']:
    ko = pd.read_csv(r'c:\Users\zyf\Desktop\课题一生信文章\验证部分工作\scRNA_analysis\final_results\KO_' + ct + '.csv')

    # Compute two-sided P from Z-score
    ko['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(ko['median_z'])))
    ko['neg_log10_p'] = -np.log10(ko['p_value'].clip(lower=1e-300))

    # Significance: freq >= 0.35
    ko['significant'] = (ko['frequency'] >= 0.35) & (ko['gene'] != 'SESN3')

    # Label: top genes with freq >= 0.50 or top 10 by Z
    top_z = ko[ko['gene'] != 'SESN3'].nlargest(15, 'median_z')
    top_labels = set(top_z['gene'])
    ko['label'] = ko.apply(lambda r: r['gene'] if r['gene'] in top_labels else '', axis=1)

    # Save
    cols = ['gene', 'median_z', 'iqr_z', 'frequency', 'p_value', 'neg_log10_p', 'significant', 'label']
    ko[cols].to_csv(os.path.join(out_r, 'ko_volcano_' + ct + '.csv'), index=False)

    # Perturbed gene list for GO
    perturbed = ko[ko['significant']]['gene'].tolist()
    with open(os.path.join(out_r, 'perturbed_genes_' + ct + '.txt'), 'w') as f:
        f.write('\n'.join(perturbed))

    n_sig = ko['significant'].sum()
    n_hi = (ko['frequency'] >= 0.50).sum()
    print(ct + ': ' + str(n_sig) + ' significant (freq>=0.35), ' + str(n_hi) + ' highly recurrent (freq>=0.50)')

print('Data ready in: ' + out_r)
