# causaldosemc

`causaldosemc` is a research-use R package for matrix-completion-based causal
estimation in panel data with a continuous treatment and `0` as the control
state. Active doses may be positive or negative.

Version `0.1` is intentionally narrow:

- balanced or unbalanced panels, with unbalanced inputs padded internally onto
  the full unit-time grid
- zero-dose masking with optional washout removal
- low-rank baseline imputation via `softImpute` or a weighted proximal SVT path
- optional blocked cross-validation for lambda selection
- optional observation weights in the main estimator
- built-in linear or spline contemporaneous or distributed-lag effect fitting
- optional staged or joint objective estimation in `cdmc_fit()`
- cross-fitted doubly robust-style linear effect estimation with either internal
  Gaussian or kernel GPS weights with linear or spline mean models, or
  external weights
- residual-based carryover screening on post-exit zero-dose periods
- SCIA-oriented dynamic assignment screening via lagged outcome history tests
- unit-level bootstrap inference for main fits, post-fit dose-response objects,
  DR fits, dynamic/pathwise estimands, and diagnostic summaries, including
  max-t simultaneous bands for multi-statistic dynamic estimands
- replay-based sensitivity scans over washout, zero-dose tolerance, and
  weighting choices for main fits, DR fits, dose-response objects, and dynamic
  targets
- formal bounded-bias sensitivity bounds for staged or joint main fits, DR
  fits, dose-response objects, and compatible dynamic/pathwise summaries, with
  exact stage-2 and baseline-counterfactual perturbation layers plus a local
  optimization-refit layer for cdmc_fit- and cdmc_dr_fit-backed targets
- refit-based placebo and carryover diagnostics for held-out control windows
- linear or spline-based post-fit dose-response estimation with optional weights

## Quick start

```r
devtools::load_all("causaldosemc")

panel <- simulate_cdmc_data(
  n_units = 20,
  n_times = 12,
  beta = 1.2,
  lag_beta = 0.4,
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
  verbose = TRUE
)

print(fit)
summary(fit)
head(predict(fit, type = "y0"))
```

The main estimator can also use a spline basis directly in stage 2:

```r
fit_spline <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  washout = 1,
  lag_order = 1,
  rank_max = 3,
  lambda = 0.2,
  effect_model = "spline",
  effect_df = 3,
  seed = 42
)

head(predict(fit_spline, type = "tau_model"))
```

If you want the manuscript-style global objective instead of the default staged
baseline-then-effect pipeline, `cdmc_fit()` now supports a joint mode:

```r
fit_joint <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lag_order = 1,
  objective = "joint",
  effect_model = "linear",
  rank_max = 3,
  lambda = 0.2,
  seed = 42
)

head(predict(fit_joint, type = "tau_model"))
```

## Main function

`cdmc_fit()` returns a `cdmc_fit` object with:

- the internally padded full unit-time grid plus `.cdmc_observed` markers for
  originally observed rows
- the imputed baseline matrix `y0_hat`
- residual treatment effects `tau_hat = y - y0_hat`
- eligible control masks after washout handling
- the optimization mask used by the chosen staged or joint objective
- normalized observation weights when a weighted fit is requested
- fitted stage-2 coefficients for either linear lags or spline bases

By default, `cdmc_fit()` uses the staged reversible-treatment estimator: it
fits the baseline on eligible zero-dose cells and then projects the residual
treatment component on the requested linear or spline dose basis. With
`objective = "joint"`, the same low-rank/additive structure is fit jointly with
the dose basis over the observed estimation sample. The staged mode remains the
default because it is still the more conservative choice when zero-dose support
and reversibility are the primary identification anchors.

`cdmc_bootstrap()` provides unit-level bootstrap inference for `cdmc_fit()`,
`cdmc_dynamic_estimand()`, `cdmc_dose_response()`, `cdmc_dr_fit()`, and
supported diagnostic objects. The main-fit path supports coefficient and
average residual-treatment summaries, the post-fit dose-response path supports
coefficient and prediction summaries, the dynamic-estimand path supports
user-supplied history responses, contrasts, and slopes with automatic max-t
simultaneous bands when multiple pathwise summaries are requested, the DR path
additionally supports lag-path averages and dose-path contrasts, and the
diagnostic path can resample placebo, carryover, equivalence, SCIA, and joint
placebo summaries directly.

`cdmc_carryover_test()` provides a simple residual-based carryover diagnostic on
post-exit zero-dose periods. It is intended as an initial screening tool rather
than a full refit-based inferential procedure.

`cdmc_placebo_test()` and `cdmc_carryover_refit_test()` remove targeted control
windows from the main baseline objective, refit the baseline surface, and summarize the
resulting pseudo effects out of sample.

`cdmc_equivalence_test()` applies a TOST-style equivalence check to placebo or
carryover diagnostics so the package can test whether those pseudo effects are
small enough to be treated as practically negligible.

`cdmc_joint_placebo_test()` runs multiple placebo windows and aggregates them
into a single pretrend diagnostic, either through Fisher's method or through a
Bonferroni-adjusted equivalence check.

`cdmc_scia_test()` provides a research-use screen for dynamic assignment bias by
testing whether lagged observed outcomes or lagged residual treatment effects
still improve prediction of the current dose after conditioning on lagged doses,
optional covariates, and optional unit or time effects.

`cdmc_dr_fit()` adds a unit-fold cross-fitted doubly robust-style linear effect
layer on top of the weighted estimator. When `weights` is omitted it estimates
fold-specific stabilized generalized propensity score weights internally from
the supplied `weight_covariates` and optional time fixed effects. The internal
path supports either a Gaussian density approximation, a kernel-smoothed
assignment density, or a fold-local CBPS balancing fit. The Gaussian and kernel
paths can use linear terms, natural splines, GAM smooths, regression trees, or
regression forests for the assignment mean model, while the CBPS path currently
supports linear or spline covariate formulas. When richer learner-based
assignment models are needed, the same function can still consume externally
supplied weights such as the output of `cdmc_cbps_weights()`.

`cdmc_dynamic_estimand()` turns any fitted main-stage, DR, or post-fit
dose-response model into a bootstrapable dynamic or pathwise estimand. It can
evaluate a user-supplied dose history as a fitted response, a contrast against
another path, or a local slope with respect to the contemporaneous dose.
Bootstrapping those objects now supports joint uncertainty over the requested
summary vector through max-t simultaneous bands.

## Cross-fitted DR layer

Internal Gaussian GPS weighting is the simplest parametric way to run the DR
layer:

```r
dr_fit <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  weight_covariates = "x1",
  gps_time_effects = TRUE,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If you already have external balancing weights, pass them directly instead:

```r
dr_weights <- cdmc_cbps_weights(panel, dose = "dose", covariates = "x1")

dr_fit_external <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weights = dr_weights,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If you want a built-in balancing-oriented DR nuisance path without precomputing
weights up front, the same function can estimate fold-local CBPS weights
internally:

```r
dr_fit_cbps <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "cbps",
  cbps_iterations = 250L,
  gps_time_effects = TRUE,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

When treatment assignment is visibly nonlinear in the observed covariates, the
same internal path can use spline-expanded Gaussian GPS nuisance models:

```r
dr_fit_spline <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  gps_model = "spline",
  gps_df = 3,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If the assignment surface is smooth but not well captured by a low-order linear
or spline basis, the Gaussian GPS path can also use a GAM mean model:

```r
dr_fit_gam <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  gps_model = "gam",
  gps_df = 4,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If assignment appears to change by thresholds or regime splits rather than by a
smooth surface, the same Gaussian GPS path can use a regression tree mean
model:

```r
dr_fit_tree <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  gps_model = "tree",
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If the assignment rule looks interaction-heavy and unstable across individual
tree splits, the Gaussian GPS path can instead use a regression forest mean
model:

```r
dr_fit_forest <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = c("x1", "x2"),
  weight_method = "gaussian_gps",
  gps_model = "forest",
  gps_forest_trees = 200L,
  gps_forest_mtry = 2L,
  gps_forest_min_node_size = 5L,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

If the conditional dose distribution looks visibly non-Gaussian, the same mean
model can be paired with a kernel-smoothed density estimate instead:

```r
dr_fit_kernel <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "kernel_gps",
  gps_model = "spline",
  gps_df = 3,
  gps_bandwidth = 0.3,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 42
)
```

The internal nuisance path is still a pragmatic research-use family. Gaussian
and kernel variants now cover linear, spline, GAM, tree, and forest mean-model
choices, while the new `cbps` mode supplies fold-local balancing weights. Use
external weights when your design needs richer boosting, stacking, or other
assignment models beyond that built-in menu.

Bootstrap inference works on DR fits as well:

```r
dr_boot <- cdmc_bootstrap(
  dr_fit,
  n_boot = 50,
  statistics = c("coefficients", "average_tau_dr", "average_tau_linear"),
  seed = 42
)

print(dr_boot)
```

You can also bootstrap richer DR functionals, such as lag-path contributions or
a custom dose-path contrast:

```r
dr_boot_paths <- cdmc_bootstrap(
  dr_fit,
  n_boot = 50,
  statistics = c("lag_average_tau_dr", "dose_contrast_dr"),
  contrast_history = c(dose_lag0 = 1, dose_lag1 = 0.5),
  reference_history = c(dose_lag0 = 0, dose_lag1 = 0),
  seed = 42
)
```

Diagnostic objects can be bootstrapped through the same interface:

```r
placebo <- cdmc_placebo_test(fit, periods = -2:0)
placebo_boot <- cdmc_bootstrap(placebo, n_boot = 50, seed = 42)

joint_placebo <- cdmc_joint_placebo_test(fit, periods = -2:-1)
joint_boot <- cdmc_bootstrap(joint_placebo, n_boot = 50, seed = 42)
```

The SCIA screen can be run directly on a fitted estimator:

```r
scia <- cdmc_scia_test(
  fit,
  lags = 1,
  outcome_proxy = "tau"
)

print(scia)
```

`cdmc_dose_response()` fits a post-baseline dose-response model to the residual
treatment surface. It currently supports linear and natural-spline models, with
optional observation weights so external continuous-treatment balancing methods
can be plugged into the second stage. When the parent `cdmc_fit()` call used
weights, `cdmc_dose_response()` reuses those weights by default unless they are
overridden explicitly.

Those post-fit dose-response objects can also be bootstrapped directly when you
want uncertainty for nonlinear fitted values or local slopes:

```r
dose_response <- cdmc_dose_response(fit, model = "spline", lag_order = 1, df = 3)

dose_response_boot <- cdmc_bootstrap(
  dose_response,
  n_boot = 50,
  statistics = c("coefficients", "prediction"),
  prediction_dose = c(-0.5, 0.5),
  prediction_type = "response",
  seed = 42
)
```

The same bootstrap workflow now supports dynamic or pathwise estimands:

```r
path_effect <- cdmc_dynamic_estimand(
  fit,
  history = data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1)
  ),
  type = "contrast",
  aggregate = "both"
)

path_boot <- cdmc_bootstrap(
  path_effect,
  n_boot = 50,
  statistics = "estimate",
  seed = 42
)

print(path_boot)
```

`cdmc_sensitivity_scan()` provides a replay-based sensitivity layer for
`cdmc_fit()`, `cdmc_dr_fit()`, `cdmc_dose_response()`, and
`cdmc_dynamic_estimand()`. It refits the stored source model over a grid of
washout choices, zero-dose tolerances, and optional weighting scenarios,
records support failures instead of aborting the whole run, and returns a
reference scenario together with explicit support, washout-dependence, and
weight-dependence summaries:

```r
sensitivity <- cdmc_sensitivity_scan(
  fit,
  washout_grid = 0:2,
  include_unweighted = TRUE
)

print(sensitivity)
summary(sensitivity)
sensitivity$washout_summary
```

For staged `cdmc_fit()` sources, those support summaries track the eligible
zero-dose optimization sample directly. For joint `cdmc_fit()` sources, the
scan now distinguishes between the zero-dose anchor summary and the actual
observed-panel optimization sample, so `fit_success` follows the optimization
support columns rather than the zero-dose-anchor columns.

The same helper can also replay downstream objects through their stored parent
fits:

```r
path_sensitivity <- cdmc_sensitivity_scan(path_effect, washout_grid = 0:2)

summary(path_sensitivity)
```

`cdmc_sensitivity_bounds()` complements those replay scans with a formal
bias-budget model for the stored post-fit layers. The default
`perturbation_layer = "stage2"` perturbs the stage-2 response or pseudo-outcome
sample directly, while `perturbation_layer = "baseline"` perturbs the stored
baseline-counterfactual layer instead. Those two layers use the stored linear
post-fit map and therefore return exact worst-case intervals for any supported
linear summary. By default the perturbation can hit every sensitivity-sample
cell, but `contamination_fraction_grid` can also restrict the bias budget to a
sparse fraction of the largest-loading cells. The default
`perturbation_constraint = "cellwise"` uses the original coordinatewise bias
budget, while `perturbation_constraint = "energy"` switches to an exact
Euclidean energy budget over the contaminated cells:

```r
bounds <- cdmc_sensitivity_bounds(
  path_effect,
  gamma_grid = c(0, 0.25, 0.5, 1),
  scale = "manual",
  scale_value = 1
)

sparse_bounds <- cdmc_sensitivity_bounds(
  path_effect,
  gamma_grid = c(0, 0.25, 0.5, 1),
  contamination_fraction_grid = c(0.1, 0.25, 1),
  scale = "manual",
  scale_value = 1
)

energy_bounds <- cdmc_sensitivity_bounds(
  path_effect,
  gamma_grid = c(0, 0.25, 0.5, 1),
  perturbation_constraint = "energy",
  contamination_fraction_grid = c(0.1, 1),
  scale = "manual",
  scale_value = 1
)

summary(bounds)

dr_baseline_bounds <- cdmc_sensitivity_bounds(
  dr_fit,
  perturbation_layer = "baseline",
  gamma_grid = c(0, 0.25, 0.5, 1),
  scale = "manual",
  scale_value = 1
)

joint_refit_bounds <- cdmc_sensitivity_bounds(
  fit_joint,
  perturbation_layer = "optimization",
  gamma_grid = c(0, 0.25, 0.5),
  scale = "manual",
  scale_value = 1,
  refit_step = 1e-3
)
```

That formal helper is exact for staged and joint `cdmc_fit()` objects,
`cdmc_dr_fit()`, `cdmc_dose_response()`, and compatible
`cdmc_dynamic_estimand()` objects built from those stored linear post-fit maps.
For targets backed by `cdmc_fit()` or `cdmc_dr_fit()`,
`perturbation_layer = "optimization"` now adds a local refit mode that
perturbs the source optimization sample itself, rebuilds the staged baseline,
joint objective, or cross-fitted DR pipeline with fixed tuning, and estimates
the resulting sensitivity multiplier from a one-step numerical Jacobian. For
DR fits, the helper now reuses the stored fold assignments so these same-panel
refits are deterministic. That mode is local rather than exact, but it closes
the most immediate gap for perturbations that move through the fitting problem
instead of only through the stored post-fit maps. The remaining gap is broader
identification-robust sensitivity modeling beyond these cellwise or energy
bias-budget bounds.

`cdmc_cbps_weights()` wraps `CBPS::CBPS()` for the continuous-treatment setting
and returns observation weights that can be supplied directly to
`cdmc_fit()`, `cdmc_dr_fit()`, or `cdmc_dose_response()`.

`cdmc_bootstrap()` performs unit-level bootstrap resampling, preserving whole
time paths for each sampled unit and optionally rerunning lambda selection in
each bootstrap replicate.

## Project roadmap

The current implementation status, design choices, limitations, and next steps
are tracked in [ROADMAP.md](ROADMAP.md).
