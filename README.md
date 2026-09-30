# Code for: Genetic architecture of Bell's palsy in 5,062 Finnish cases

Analysis code accompanying the manuscript

> **Genetic architecture of Bell's palsy in 5,062 Finnish cases: two non-MHC risk loci
> and no classical HLA association**

This repository contains the Python and R code used for harmonisation, meta-analysis,
SuSiE credible-set handling, colocalisation, credible-set interrogation, cross-trait
Mendelian randomisation and every figure in the paper. See `MANIFEST.md` for the mapping
between each analysis and the script that performs it.

## Layout

```
.
├── scripts/          analysis and figure code
├── results/          small result tables produced by the scripts (for checking values)
├── figures.R         drawing code for all figures (base R)
├── MANIFEST.md       analysis -> script -> output mapping
├── requirements.txt
└── LICENSE
```

The scripts expect a workspace with three subdirectories next to `scripts/`:

```
<workspace>/
├── data/     public summary statistics you download (see below)
├── work/     intermediate caches (created automatically)
└── out/      result tables
```

By default `<workspace>` is the parent of `scripts/`. Set the environment variable
`BELLSPALSY_ROOT` to override it.

## Environment

- Python 3.14, standard library plus **numpy** only. No pandas, no scipy: the special
  functions and the Mendelian randomisation estimators are implemented from scratch in
  `mrlib.py`, and the colocalisation and fine-mapping routines in `coloclib.py` and
  `susie.py`.
- R 4.6.1 with **base graphics only** (no ggplot2 or other packages) for `figures.R`.

## Data

All input data are public. Nothing here redistributes them. To re-run:

| Data | Source | Used for |
|---|---|---|
| FinnGen release 13 summary statistics (`G6_BELLPA`, `G6_NERPLEX`, `G6_CARPTU`, `G6_EPLEPSY`, `AUTOIMMUNE_NONTHYROID`, `E4_HYTHY_AI_STRICT`, `T1D_WIDE`, `M13_RHEUMA`, `L12_PSORIASIS`, `G6_MS`, `T2D`) | <https://r13.finngen.fi/> | main analyses |
| FinnGen R13 SuSiE credible sets | FinnGen R13 `finemap/summary/` | fine-mapping |
| FinnGen R13 HLA imputation results | FinnGen R13 | HLA allele analysis |
| FinnGen R13 pheweb variant API | `https://r13.finngen.fi/api/variant/<chrom>-<pos>-<ref>-<alt>` | phenome-wide scan |
| Bell's palsy GWAS, GWAS Catalog **GCST90018843** (Sakaue et al.) | <https://www.ebi.ac.uk/gwas/> | meta-analysis, cross-trait comparison |
| eQTLGen phase 1 cis-eQTL (full) | <https://www.eqtlgen.org/> | cis-eQTL annotation (Fig. 2c) |
| Ensembl release 112 GTF, GRCh38 | <https://www.ensembl.org/> | gene coordinates |
| 1000 Genomes phase 3 genotypes (for LD) | <https://www.internationalgenome.org/> | LD clumping, checks |

By default the scripts look for these under `data/`.

## Running

```bash
export BELLSPALSY_ROOT=/path/to/workspace     # optional
python scripts/phenome_scan.py                # writes out/phenome_scan_rs9357446.tsv
python scripts/mhc_snp_scan.py
python scripts/hla_alleles.py
python scripts/colocalisation_and_credible_set.py
python scripts/cross_trait_mr.py
python scripts/extract_figure_data.py
python scripts/build_figure_data.py
Rscript figures.R                             # writes figures/Fig*.pdf and .png
```

Scripts are ordered so that later ones consume earlier outputs. `colocalisation_and_credible_set.py`
caches per-trait regional data in `work/trait_cache.pkl` on first run.

**Figure pipeline.** `figures.R` reads `figures/data/*.tsv`, which are shipped here so that
the figures can be redrawn without re-running the genetics. To rebuild them from scratch:

1. `python scripts/extract_figure_data.py` - scans the FinnGen summary statistics and writes
   `work/fig_manhattan.tsv` (~24 MB) and `work/fig_region.tsv`. These two intermediates are
   **not** shipped because they are large and fully regenerable from the public data.
2. `python scripts/build_figure_data.py` - writes `figures/data/*.tsv`.
3. `Rscript figures.R` - writes `figures/Fig*.pdf` and `.png`, plus the graphical abstract.

## Notes on two implementation choices

- **Colocalisation marginal likelihoods use means, not sums.** Summing the per-variant
  approximate Bayes factors makes the H3 likelihood grow with the square of the number of
  variants, so every region is called H3. The implementation was validated by colocalising
  Bell's palsy against itself, which returns PP.H4 = 1.000.
- **Distance pruning for cross-trait Mendelian randomisation walks the genome in order**
  and retains the most significant variant per window. An earlier version sorted variants
  by P value first, which left thousands of correlated instruments.

## Licence

MIT (see `LICENSE`). Please change this if your institution requires otherwise.

## Citation

If you use this code, please cite the manuscript above.
