# Analysis overview

## Stage 1 - genetic screen (`01_genetic_pipeline/`)

| Script | Purpose |
| --- | --- |
| `00_download_data.R` | retrieves OneK1K cis-eQTL records and the MG GWAS summary statistics |
| `01_data_preprocessing.R` | allele harmonisation, effect-direction alignment, exclusion of the extended MHC, conversion of OneK1K correlations to per-allele effects |
| `02_mr_analysis.R` | single-instrument Wald ratios, instrument-strength audit, global Benjamini-Hochberg correction |
| `03_colocalization.R` | `coloc` v5.2.3 with the converted per-allele effects; conditional shared-signal probability |
| `04_visualization.R` (+ `_v2`, `_v3`) | first-submission figures |
| `05_dice_validation.R` | DICE expression-QTL audit (reported as an audit only) |
| `06_phewas.R`, `06b_phewas_snp.R` | phenome-wide records for the post-exclusion lead variants |
| `07_drug_repurposing.R` | directionality-aware drug annotation |

## Stage 2 - single-cell disease context (`02_single_cell/`)

`01_load_and_qc.R` loads the seven 10x matrices, applies donor-aware adaptive QC and scDblFinder, and
writes the annotated object. The Python scripts (`06`-`11`, `ager_*`) run the donor-level comparisons
and the virtual knockout; `prepare_R_data.py` exports count matrices for the R analyses;
`R_figure_4panel.R` and `make_reference_style_validation_figures.R` draw the single-cell panels.

## Stage 3 - revision analyses (`03_revision_figures_tables/`)

`MG_V4_complete_figures_and_tables.R` rebuilds all main and supplementary figures plus the 15 analysis
tables (post-MHC discovery counts, MR/colocalisation evidence, SESN3 donor expression, leave-one-donor-out,
exact permutation, pseudobulk DE, Hallmark 50x2, sample QC, composition, hierarchical sensitivity,
descriptive PheWAS, DICE audit, instrument alignment, cross-modality matrix).

`LD_locus_sensitivity/` consolidates the lead variants of the discovery set into LD-based loci using
pairwise r2 from 1000 Genomes Phase 3 EUR (Ensembl GRCh37 REST API; r2 >= 0.1, 1-Mb window) and reports
the MHC-boundary LD of the three variants immediately outside the coordinate-defined extended MHC.
These are LD-based sensitivity analyses and are **not** PLINK clumping results.
