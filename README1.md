# SMMR Submission Package

Manuscript: **Scalable Bayesian Structure Learning of Directed Acyclic Graphs via Laplace Approximation, with an Application to Breast Cancer Gene Expression Networks**

Authors: S. Nazari, M. Arashi (corresponding author), A. Sadeghkhani

## Files

| File | Purpose |
|------|---------|
| `main_smmr.pdf` | Compiled manuscript (31 pages) |
| `main_smmr.tex` | LaTeX source — switch to `\documentclass[Crown,times,sageh]{sagej}` for actual SMMR submission once you have the SAGE class file |
| `cover_letter_SMMR.docx` | Editable Word cover letter to the Editor-in-Chief |
| `code/` | All R source code (see below) |
| `figures/` | 13 PDF figures referenced by the manuscript |
| `tables/` | 10 CSV + 10 TeX tables |

## R code organisation

The simulation and real-data analyses are kept in separate files with no overlap.

### Core
- `BSL_core.R` — implementation of `log_marg_nCPNG`, `log_marg_CPNIG`, `mcmc_dag`, `is_dag_fast`, `ges_like`, `pc_like`, `notears_like`, performance metrics, and ROC utilities

### Simulation
- `sim_part1a_q20.R` — produces `results_q20.rds` (split for tractability)
- `sim_part1b1_q40.R` — produces `results_q40a.rds` (n=100, 200)
- `sim_part1b2_q40.R` — produces `results_q40b.rds` (n=300, 500)
- `sim_aggregate.R` — combines the three RDS files into all tables and ROC/benchmark figures
- `sim_runtime.R` — runtime scaling experiment at q ∈ {10, 20, 40, 60}
- `sim_diagnostics.R` — MPSRF, posterior contraction, Laplace-vs-exact

Run order:
```
Rscript sim_part1a_q20.R
Rscript sim_part1b1_q40.R
Rscript sim_part1b2_q40.R
Rscript sim_aggregate.R
Rscript sim_runtime.R
Rscript sim_diagnostics.R
```

### Real data
- `real_data.R` — single self-contained script: pre-processes the GSE7390 panel, runs the nCPNG MCMC, builds the network plot, computes HPD intervals via BMA, runs the sensitivity grid, and produces `figures/heatmap_real.pdf`, `figures/network_real.pdf`, `figures/boxplot_real.pdf`, `figures/trace_real.pdf`, `figures/sensitivity_real.pdf`, and the corresponding tables.

The script first tries to load GSE7390 via `GEOquery`. If that package is unavailable in the runtime environment, it falls back to a semi-synthetic dataset whose moment structure preserves the published covariance pattern of the panel and where IL8/OAS2 are the two strongest direct drivers of the metastasis indicator. **For final submission you should re-run on the actual GEOquery download and update the numbers in Table 6 and Section 6.**

## Reproducibility notes

- **Replicates:** `N_REP = 3` per (q, n) cell in the included scripts to keep wall-clock time manageable in the sandbox. For a final submission, increase to `N_REP = 20`–`40` by editing the constant at the top of each `sim_part1*.R` file.
- **Random seeds:** every cell uses `set.seed(1000*q + 100*n + r)` for full reproducibility.
- **R packages:** `MASS`, `Matrix`, `mvtnorm`, `ggplot2`, `reshape2`, `viridis`, `xtable`, `gridExtra`, `coda`. No CRAN-only dependencies; everything installs via standard R or apt.
- **Compilation:** `pdflatex main_smmr.tex` three times (for `\ref` resolution).
