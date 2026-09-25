# Cell-type-resolved cis-eQTL prioritisation in myasthenia gravis

Code and derived data accompanying the revised manuscript
*Cell-type-resolved cis-eQTL prioritisation identifies SESN3 in naive/central memory T cells in
myasthenia gravis* (Frontiers in Immunology, revision).

This repository contains everything needed to **inspect and reproduce the reported numbers**:
the genetic pipeline, the single-cell pipeline, the revision-stage figure/table code, every derived
table, and the LD locus-consolidation analysis added during revision.

## What is here

| Folder | Contents |
| --- | --- |
| `01_genetic_pipeline/` | OneK1K cis-eQTL download, preprocessing, instrument selection, single-instrument Wald ratios, Bayesian colocalisation, DICE audit, PheWAS, drug annotation (`00`-`07`, R) |
| `02_single_cell/` | In-house PBMC scRNA-seq: loading and QC, annotation-based sensitivity, donor-level pseudobulk, virtual knockout, figure code (R + Python) |
| `03_revision_figures_tables/` | Revision-stage code: all main/supplementary figures and 15 analysis tables (`MG_V4_complete_figures_and_tables.R`) and the LD-based locus consolidation analysis |
| `04_derived_tables/` | The 15 analysis tables, the per-panel source tables, and the delivery status files |
| `05_reference_tables/` | Tables as first submitted (MR, colocalisation, DICE, PheWAS, drug annotation) |
| `06_environment/` | `sessionInfo.txt` of the analysis session used for the revision figures |
| `07_qPCR_TO_BE_ADDED/` | Reserved for the qPCR data (added separately) |
| `docs/` | Analysis overview, data dictionary, availability statement drafts |

## Analysis in one paragraph

Immune-cell-resolved cis-eQTLs from the OneK1K cohort (14 peripheral blood cell types) were combined
with myasthenia gravis GWAS summary statistics (5,708 cases / 432,028 controls). The extended MHC
region (chr6:25,726,063-33,400,644, GRCh37) was excluded before instrument selection and
colocalisation, and OneK1K correlations were converted to per-allele effects
(beta = rho / sqrt(2p(1-p))). Single-instrument Wald ratios and `coloc` (v5.2.3) were computed for
7,328 post-exclusion gene-cell type pairs; 39 pairs passed global Benjamini-Hochberg correction
(13 LD-based loci) and 31 reached PP.H4 > 0.80. SESN3 was prioritised in CD4+ and CD8+
naive/central memory T cells through the same lead variant (rs4409785). Disease-context evidence came
from donor-level pseudobulk analysis of an in-house PBMC scRNA-seq dataset (3 healthy controls,
4 patients) with exact donor-label permutation over all 35 allocations.

## How to reproduce

1. `01_genetic_pipeline/00_download_data.R` - fetches the public input resources (see
   `docs/DATA_DICTIONARY.md`; the GWAS summary statistics are not redistributed here).
2. `01_genetic_pipeline/00_run_pipeline.R` - runs `01` to `07` in order and writes the result tables.
3. `02_single_cell/` - the scRNA-seq scripts are numbered in execution order; the raw 10x matrices
   are deposited in GEO (accession in the manuscript).
4. `03_revision_figures_tables/MG_V4_complete_figures_and_tables.R` - regenerates all figures and the
   15 analysis tables of the revised manuscript from the intermediate objects produced in step 2.
5. `03_revision_figures_tables/LD_locus_sensitivity/` - pairwise LD and locus consolidation
   (1000 Genomes Phase 3 EUR through the Ensembl GRCh37 REST API).

R version and package versions are recorded in `06_environment/sessionInfo.txt`.

## Data availability

Public resources used (not redistributed here): OneK1K cis-eQTLs (https://onek1k.org/), the MG GWAS
meta-analysis (GWAS Catalog GCST90432156), DICE (https://dice-database.org/), Open Targets,
NHGRI-EBI GWAS Catalog, DSigDB and Enrichr. The in-house PBMC scRNA-seq dataset (3 healthy controls,
4 patients) is deposited in GEO (accession to be inserted).

## Licence

Code: MIT (`LICENSE`). Derived tables and documentation: CC BY 4.0 (`LICENSE_DATA.md`).
