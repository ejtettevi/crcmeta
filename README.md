# crcmeta

Reproducible R workflow for **study-adjusted meta-analysis of multicohort shotgun metagenomes**, demonstrated on public colorectal cancer (CRC) cohorts from [curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData).

Naïve pooling of samples across studies treats batch, geography and protocol differences as biology. `crcmeta` instead:

1. Exports species/genus MetaPhlAn profiles + metadata from *curatedMetagenomicData*
2. CLR-transforms abundances **within each study**
3. Runs pooled differential abundance (Wilcoxon + FDR + Cliff’s δ)
4. Combines within-study Hedges’ *g* by **DerSimonian–Laird random-effects meta-analysis**
5. Reports **leave-one-study-out (LOSO)** stability and **cross-study effect concordance**
6. Regenerates all manuscript figures at 750 dpi from a single script

## Validation case study

On 11 public CRC shotgun cohorts (1,650 faecal metagenomes: 701 CRC, 694 controls, 209 adenomas):

- **Disease signal is robust:** 40 of 161 species significant after study adjustment; all top-20 species stable in 100% of LOSO refits; median cross-study effect correlation ρ = 0.61
- **Young- vs old-onset CRC is correctly null after study adjustment:** 0 of 161 species; continuous age correlations weak (max |ρ| = 0.19)

See `results/RESULTS_SUMMARY.md` and `results/tables/`.

## Requirements

- R ≥ 4.0
- Bioconductor ≥ 3.14 with **curatedMetagenomicData**
- CRAN: `ggplot2`, `patchwork`, `vegan`, `readr`, `dplyr`, `tidyr`, `forcats`, `scales` (and optionally `metafor`)
- ~4 GB RAM; ~2 GB disk after profile export

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedMetagenomicData")
install.packages(c("ggplot2", "patchwork", "vegan", "readr",
                   "dplyr", "tidyr", "forcats", "scales"))
```

## Pipeline

Run from the repository root:

```bash
bash  analysis/01_fetch_metadata.sh          # cache ENA/SRA run reports (optional)
Rscript analysis/02_fetch_profiles.R         # export MetaPhlAn profiles + metadata
Rscript analysis/03_differential_abundance.R # pooled Wilcoxon + Cliff's delta
Rscript analysis/05_meta_analysis.R          # study-adjusted RE meta + LOSO + cross-study ρ
Rscript analysis/06_genus_level.R            # genus-level pooled contrasts
Rscript analysis/07_extended_analyses.R      # continuous age, overlap, clinical, permanova inputs
Rscript analysis/08_figures_publication.R    # all nine figures at 750 dpi
```

Outputs land in `results/tables/` and `results/figures/`.

### Cohort selection

Edit `config/cohorts.tsv` to change which *curatedMetagenomicData* study IDs are included. Defaults match the 11-cohort CRC analysis (FengQ_2015 … ZellerG_2014).

## Layout

```
config/     cohort list
analysis/   pipeline scripts 00–08
results/    tables/ and figures/ (regenerable)
logs/       run logs (not committed)
data/       profile cache (not committed; regenerate with 02)
manuscript/ submission materials (not part of this release)
```

## Software

- **License:** MIT
- **Language:** R
- **Versioning:** tag releases as `vX.Y.Z`; archive snapshots to Zenodo for a citable DOI

## Citation

When using this workflow, please cite the software release (DOI assigned on archival) and the original cohort studies listed in `results/tables/study_characteristics.tsv`.

## Acknowledgements

Profiles and curated metadata are distributed by the *curatedMetagenomicData* project (Waldron Lab). Underlying primary data belong to the original cohort investigators.
