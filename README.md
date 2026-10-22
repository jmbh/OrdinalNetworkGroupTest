# Bayesian ordinal-network group comparison

Reproducibility archive for *Testing Group Differences in Network Models using
Bayesian Model Selection: A Tutorial* by Jonas Haslbeck, Mona Klau, Tessa
Blanken, and Maarten Marsman.

The application compares the conditional associations among six problematic
pornography use (PPU) symptoms in women and men using `bgms::bgmCompare()`.
This repository contains only the raw PPU file, the scripts needed for this
application, an archived model fit, and the resulting paper outputs.

## Reproduce

Use R 4.6.1 and `bgms` 0.2.0.0, and run from the repository root:

```sh
Rscript run_all.R
```

This rebuilds the analysis data, numerical summaries, tables, and all figures
from the archived fit. To rerun the Bayesian sampler from the raw data as well:

```sh
Rscript run_all.R --refit
```

The refit uses four chains with 10,000 warmup and 20,000 retained draws per
chain, so it takes longer. The random seed is fixed at 1. Posterior summaries can still
differ slightly across software versions or platforms; the archived fit is
included to reproduce the exact reported output.

The standard reproduction also constructs 500 posterior predictive replicated
datasets per group. This takes a few minutes and uses fixed simulation seeds.

Install the required packages if needed:

```r
install.packages(c(
  "bgms", "ggraph", "ggplot2", "patchwork", "readxl", "tidygraph", "xtable"
))
```

The archived fit was created with `bgms` 0.2.0.0 using its Normal(0, 1)
defaults for baseline interactions and nonzero group differences, a Bernoulli(.5)
prior on difference inclusion, and a beta-prime(.5, .5) threshold prior. Package
versions used by the last successful run are recorded in
`results/diagnostics/session_info.txt`.

## Workflow

| Script | Purpose |
|---|---|
| `R/01_prepare_data.R` | Select the six PPU items, retain women and men, validate the data, and save group-specific analysis files. |
| `R/02_fit_model.R` | Fit the Bayesian ordinal MRF comparison, or load the archived fit unless `--refit` is supplied. |
| `R/03_summarize_results.R` | Export edge estimates, Bayes factors, threshold tables, diagnostics, and the Withdrawal--Conflict interpretation. |
| `R/04_plot_descriptives.R` | Plot the six observed marginal distributions. |
| `R/05_plot_model.R` | Create the main four-panel network figure and its heatmap counterpart. |
| `R/06_plot_priors.R` | Visualize the default priors and alternatives discussed in the paper. |
| `R/07_posterior_predictive_checks.R` | Generate joint-distribution posterior predictive checks and chain-specific inclusion-indicator diagnostics. |
| `R/08_prior_sensitivity.R` | Check edge-difference Bayes factors across Normal slab scales 0.5--2 and create the appendix summary and figure. |

Shared constants and validation live in `R/00_config.R`; reusable extraction
functions live in `R/functions.R`. `run_all.R` executes the scripts in order.

## Inputs and outputs

- `data/raw/PPU Data.xlsx`: supplied source data.
- `data/derived/`: analysis-ready female and male data, rebuilt on every run.
- `results/model/bgms_fit.rds`: archived `bgmCompare` fit; replaced by `--refit`.
- `results/model/prior_sensitivity.rds`: archived scale-sensitivity check; recomputed by `--refit`.
- `results/tables/`: complete edge results, key interpreted quantities, and
  appendix threshold tables in CSV and LaTeX formats.
- `results/diagnostics/`: MCMC diagnostics, reproducibility checks, and session
  information, including complete posterior predictive summaries and
  chain-specific diagnostics for all inclusion indicators.
- `figures/`: the nine PDFs used in the manuscript: two main-text figures,
  four prior figures, the appendix network heatmaps, the posterior
  predictive-check figure, and the prior-sensitivity figure.

All generated files use paths relative to the repository root. No working
directory or user-specific path is embedded in the code.

## Reproducibility checks

The workflow stops if data preparation does not yield 520 women and 387 men, if
the six selected items contain missing or invalid values, or if the archived fit
has an unexpected variable set, sample size, or `bgms` version. It also checks
the reported convergence and sampler-diagnostic thresholds. A summary is
exported to `results/diagnostics/reproducibility_checks.csv`.

The pipeline contains automated consistency checks and stops when core data or
fit assumptions are violated, but most failures currently produce basic
assertion messages rather than fully tailored diagnostic guidance. These checks
protect the analysis archive against accidental changes; the corresponding
values in the manuscript must still be kept aligned manually.

## Data and citation

The data originate from Klau et al. (2026), *Problematic pornography use:
uncovering symptom importance and differences between women and men through
Bayesian analysis of networks*, *The Journal of Sexual Medicine*, 23(6),
qdag155. <https://doi.org/10.1093/jsxmed/qdag155>

The raw workbook contains study data beyond the seven columns used here. Treat
it as research-participant data and follow the original study's consent,
ethical, and reuse conditions. See `LICENSE.md`.
