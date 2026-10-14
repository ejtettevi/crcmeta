# CRC Shotgun Project — Results Summary

Analysis route: pre-profiled MetaPhlAn3 taxonomic profiles from `curatedMetagenomicData`
(fast, < 1 GB) instead of re-processing raw FASTQ (200 metagenomes ≈ 1–2 TB, exceeds local disk).

## Data pooled
11 CRC/adenoma studies (species-level profiles), pooled on shared species keys:

| Study | Taxa | Samples |
|-------|-----:|--------:|
| YachidaS_2019 | 718 | 616 |
| ZellerG_2014 | 652 | 156 |
| FengQ_2015 | 606 | 154 |
| YuJ_2015 | 575 | 128 |
| WirbelJ_2018 | 537 | 125 |
| VogtmannE_2016 | 540 | 110 |
| HanniganGD_2017 | 292 | 81 |
| ThomasAM_2018a | 477 | 80 |
| ThomasAM_2019_c | 519 | 80 |
| GuptaA_2019 | 308 | 60 |
| ThomasAM_2018b | 503 | 60 |

**Shared species-level taxa across all 11 studies: 161**

## Differential abundance (CLR, Wilcoxon, BH-FDR, Cliff's delta)

| Comparison | nA | nB | Sig. taxa (q<0.05) |
|------------|----|----|-------------------:|
| CRC vs Control        | 701 | 694 | **88 / 161** |
| Adenoma vs Control    | —   | —   | see `diffabund_Adenoma_vs_control.tsv` |
| Young (<50) vs Old (≥50) CRC | 82 | 619 | **3 / 161** |
| Taxon vs age (CRC, Spearman) | 701 | — | **21 taxa** |

## Key interpretation (relevant to revision)
1. **Power from pooling (Reviewer 4, point 1):** pooling raises CRC/control to ~700/group;
   the disease contrast recovers **88/161** significant species — a strong, well-powered signal,
   versus the manuscript's 8/group design.
2. **Age effects are small (Reviewer 4, point 3):** young-vs-old CRC yields only **3/161**
   significant species, and species-vs-age associations (n=701) have modest effect sizes
   (max |Spearman rho| ≈ 0.18). This is consistent with age-related microbiota differences
   being real but weak relative to the disease signal — supporting a tempered "novelty" claim.
3. Drivers: *Escherichia coli*↑, *Bifidobacterium*/*Haemophilus* patterns with age;
   classic CRC taxa (e.g. *F. prausnitzii*, *Roseburia* spp.) shift with disease.

## Figures
- `results/figures/fig_diffabund_CRC_vs_control.png`
- `results/figures/fig_diffabund_Young_vs_Old.png`

## Tables
- `results/tables/diffabund_CRC_vs_control.tsv`
- `results/tables/diffabund_Adenoma_vs_control.tsv`
- `results/tables/diffabund_Young_vs_Old_CRC.tsv`
- `results/tables/assoc_taxon_vs_age_CRC.tsv`

## Study-adjusted analysis (per-cohort Hedges' g + random-effects meta-analysis)

`analysis/05_meta_analysis.R` combines studies as *effect sizes* (DerSimonian–Laird
random effects, tau^2, I^2) instead of pooling samples, so batch effects become
between-study heterogeneity rather than signal contamination.

| Comparison | k studies | Sig. taxa (q<0.05) | median I^2 |
|------------|----------:|-------------------:|-----------:|
| CRC vs Control     | 11 | **40 / 161** | 31% |
| Adenoma vs Control | 5  | 0 / 161      | 6%  |
| Young vs Old CRC   | 10 | **0 / 161**  | 0%  |

**Top CRC vs Control (direction-consistent across studies):**
- Enriched: *Peptostreptococcus stomatis* (g=+0.65, 11/11 studies), *Gemella morbillorum*
  (+0.65, 10/11), *Solobacterium moorei* (+0.54, 11/11), *Clostridium symbiosum* (+0.43, 11/11),
  *Porphyromonas asaccharolytica* (+0.41, 11/11).
- Depleted: *Streptococcus salivarius* (g=−0.38, 11/11), *Roseburia intestinalis*
  (−0.38, 10/11), *Faecalibacterium prausnitzii* (−0.39, 10/11), *Eubacterium eligens* (−0.29, 10/11).

This reproduces the **known CRC microbiome signature** with strong per-study consistency,
and after study adjustment the **young-vs-old contrast is null** (0/161; best Fisher p≈0.003,
not FDR-significant).

## Bottom line for Reviewer 4
- Point 1 (power): study-adjusted CRC/control contrast is robust at k=11 studies / ~700 per arm.
- Point 3 ("novel" age differences): the age effect **does not survive** study adjustment —
  it is weak/heterogeneous, not a robust novel finding.
- Data/scripts/tables/figures are reproducible in `analysis/` + `results/`.

> Notes: single-chunk CLR + Wilcoxon = `03`; study-adjusted meta-analysis = `05`.
> Figures: `fig_forest_*` (forest), `fig_diffabund_*` (effect bars).
