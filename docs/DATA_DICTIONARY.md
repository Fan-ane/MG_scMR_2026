# Data dictionary (selected files)

## 04_derived_tables/analysis_tables

| File | Rows | Content |
| --- | --- | --- |
| `Table01_postMHC_discovery_counts.csv` | | pairs / distinct lead variants / position-based 1-Mb groups for each discovery set |
| `Table02_MHC_before_after_UNRECONCILED.csv` | | pre- vs post-MHC counts; the pre-exclusion baseline is **unreconciled** with the submitted version |
| `Table03_MR_coloc_evidence.csv` | 7,328 | one row per post-exclusion gene-cell type pair: lead variant, MR P, global FDR, PP.H4, evidence class |
| `Table04_SESN3_donor_expression.csv` | 14 | SESN3 log2 CPM per donor and cell type |
| `Table05_SESN3_LODO.csv` | 16 | full model and seven leave-one-donor-out models |
| `Table06_SESN3_exact_inference.csv` | 2 | exact permutation P values and exact BH q values |
| `Table07_pseudobulk_DE.csv` | | donor-level pseudobulk differential expression across models |
| `Table08_Hallmark_50x2.csv` | 100 | 50 Hallmark sets x 2 populations, NES and within-population FDR |
| `Table09_sample_QC.csv` | 7 | per-sample QC, doublet rates and singlet counts |
| `Table10_per_donor_composition.csv` | | cell-type proportions per donor |
| `Table11_hierarchical_sensitivity.csv` | 3 | pair / lead-variant / physical-interval level correction |
| `Table12_PheWAS_descriptive_only.csv` | 6 | descriptive SNP-trait category counts (no enrichment test) |
| `Table13_DICE_historical_only.csv` | 2 | DICE overlap audit |
| `Table14_instrument_alignment.csv` | 7,328 | instrument strength, allele alignment, MR estimates |
| `Table15_cross_modality_evidence.csv` | 7,328 | evidence matrix across data layers |

## 03_revision_figures_tables/LD_locus_sensitivity

`01_locked_39_gene_cell_pairs.csv`, `02_locked_33_lead_variants.csv`,
`03_pairwise_LD_55_with_sensitivity_annotation.csv`, `04`/`05`/`06` locus membership, locus summary and
gene-cell mapping (observed support), `07`/`08`/`09` the same under the conservative forced-merge
scenario, `10_MHC_boundary_LD_summary.csv`, `11_unresolved_pair_EUR_subpopulation_audit.csv`,
`12_LD_locus_count_sensitivity_summary.csv`.

## Not redistributed

OneK1K cis-eQTL records, MG GWAS summary statistics, DICE records and Open Targets exports are public
resources; they are retrieved by `00_download_data.R` and linked here rather than re-hosted.
