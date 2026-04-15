# causaldosemc

`causaldosemc` provides matrix-completion-based tools for causal estimation in
panel data with a continuous treatment that can switch on and off over time.
The package combines zero-dose masking, low-rank baseline estimation,
distributed-lag treatment modeling, weighting, diagnostics, bootstrap
uncertainty, and downstream dose-response summaries.

## Status

This is an experimental research package intended for empirical and
methodological workflows. The public API is usable, but the package is still
best treated as an actively developing GitHub release rather than a frozen
long-term interface.

Current scope:

- control state fixed at `0`
- balanced panels and unbalanced panels padded internally onto the full
  unit-time grid
- staged and joint `cdmc_fit()` estimation paths
- blocked cross-validation or heuristic lambda selection
- cross-fitted DR estimation with internal GPS or balancing-based weighting
- post-fit dose-response and dynamic/pathwise estimands
- placebo, carryover, joint placebo, and SCIA diagnostics
- bootstrap inference and bounded-bias or replay-based sensitivity analysis

## Installation

From a local clone:

```r
install.packages("remotes")
remotes::install_local(".")
```

From GitHub after publication:

```r
remotes::install_github("<owner>/<repo>")
```

## Quick start

```r
library(causaldosemc)

panel <- simulate_cdmc_data(
  n_units = 20,
  n_times = 12,
  beta = 1,
  lag_beta = 0.3,
  seed = 42
)

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  washout = 1,
  lag_order = 1,
  rank_max = 3,
  lambda = NULL,
  lambda_selection = "cv",
  cv_rounds = 3L,
  cv_block_size = 2L,
  seed = 42
)

summary(fit)

dose_response <- cdmc_dose_response(
  fit,
  model = "spline",
  lag_order = 1,
  df = 4
)

predict(
  dose_response,
  dose = seq(0, 1, length.out = 5),
  type = "response"
)
```

For empirical work, prefer blocked cross-validation whenever the data contain
enough eligible zero-dose support to hold out contiguous control blocks safely.
The heuristic lambda rule remains a fallback when support is too thin for a
stable blocked tuning exercise.

## Main entry points

- `cdmc_fit()`: staged or joint matrix-completion estimator for the baseline
  counterfactual surface and treatment-effect layer.
- `cdmc_dr_fit()`: cross-fitted doubly robust effect estimator with internal
  Gaussian or kernel GPS weighting, fold-local balancing weights, or external
  weights.
- `cdmc_dose_response()`: post-fit dose-response modeling using linear,
  spline, GAM, tree, or forest specifications.
- `cdmc_dynamic_estimand()`: response, contrast, and slope summaries for
  user-supplied treatment paths.
- `cdmc_bootstrap()`: unit-level bootstrap inference for fits, downstream
  estimands, and supported diagnostics.
- `cdmc_sensitivity_scan()` and `cdmc_sensitivity_bounds()`: replay-based and
  bounded-bias sensitivity workflows.

## Vignettes

- [Getting started](vignettes/getting-started.Rmd)
- [Empirical workflow](vignettes/empirical-workflow.Rmd)
- [DR weighting workflows](vignettes/dr-weighting.Rmd)
- [Dose-response and dynamic estimands](vignettes/dose-response-dynamics.Rmd)
- [Manuscript-style case study](vignettes/policy-case-study.Rmd)

After installation, run `browseVignettes("causaldosemc")` to open the rendered
articles locally.

## Development

- Run package tests with `testthat::test_local()`.
- Build the source package from the repository root with `R CMD build .`.

This repository tracks package source only. Local tarballs, ad hoc debugging
scripts, test-output dumps, and editor settings are intentionally excluded from
version control.