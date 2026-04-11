# Project Roadmap

## Original goal

The library goal remains the same as the manuscript goal: build an R library for
causal estimation with matrix completion under a continuous, reversible
treatment, where treatment can switch on and off multiple times and may have
lagging effects. The implementation should stay aligned with that target even
when early versions intentionally narrow scope for correctness and research-use
stability.

## Current status

The package currently exists as a research-use v0.1 core engine in
[README.md](README.md) with the main estimator entry point in [R/fit.R](R/fit.R).

Implemented now:

- Balanced and unbalanced panel support through internal padding onto the full
  unit-time grid.
- Continuous dose support with 0 as the only control state, including signed
  active doses.
- Exact unit-time uniqueness checks.
- Zero-dose control masking with optional washout removal after treatment.
- Optional blocked cross-validation for lambda selection on eligible zero-dose
  observations.
- Stage-1 baseline estimation using additive unit and time effects, optional
  covariates, and low-rank matrix completion through softImpute in the
  unweighted case or a weighted proximal SVT routine when observation weights
  are supplied.
- Stage-2 residual treatment-effect estimation using either a linear or spline
  contemporaneous or distributed-lag dose model, with optional observation
  weights in the main fit.
- An optional joint objective for `cdmc_fit()` that estimates the low-rank
  surface, additive nuisance terms, and a global linear or spline treatment
  basis together over the observed estimation sample.
- A research-use cross-fitted doubly robust-style linear effect estimator using
  either internally estimated Gaussian or kernel GPS weights with linear,
  spline, GAM, tree, or forest assignment mean models, fold-local CBPS balancing
  weights with linear or spline covariate formulas, or externally supplied
  observation weights.
- A residual-based carryover diagnostic for post-exit zero-dose periods.
- Refit-based placebo and carryover diagnostics for held-out control windows.
- A research-use SCIA-oriented diagnostic screen that tests whether lagged
  outcomes still predict current dose assignment after conditioning on observed
  history.
- Unit-level bootstrap inference for linear effect coefficients and mean residual
  treatment effects in `cdmc_fit()`, plus coefficient, mean augmented effect,
  mean fitted DR-effect, lag-path contribution, and dose-path contrast
  summaries for `cdmc_dr_fit()`, and bootstrap support for placebo,
  carryover, equivalence, SCIA, and joint placebo diagnostic objects, with
  optional rerunning of lambda selection where the stored fit supports it.
- Bootstrap support for post-fit `cdmc_dose_response()` coefficients and
  user-specified response or slope predictions.
- A bootstrapable `cdmc_dynamic_estimand()` helper for user-supplied pathwise
  responses, contrasts, and slopes built on top of `cdmc_fit()`,
  `cdmc_dr_fit()`, or `cdmc_dose_response()`.
- Max-t simultaneous bootstrap bands for multi-statistic
  `cdmc_dynamic_estimand()` summaries.
- A replay-based sensitivity helper, `cdmc_sensitivity_scan()`, for washout,
  zero-dose tolerance, and weighting scenarios in `cdmc_fit()`,
  `cdmc_dr_fit()`, `cdmc_dose_response()`, and
  `cdmc_dynamic_estimand()`, including explicit zero-dose anchor summaries,
  objective-aware optimization-support summaries, and aggregated washout- and
  weight-dependence tables relative to the stored source fit.
- A formal sensitivity helper, `cdmc_sensitivity_bounds()`, for staged or joint
  `cdmc_fit()`, `cdmc_dr_fit()`, `cdmc_dose_response()`, and compatible
  `cdmc_dynamic_estimand()` objects, reporting exact worst-case bounds for
  supported linear summaries under bounded additive perturbations to the
  stored stage-2 response layer or the baseline-counterfactual layer, plus a
  local optimization-refit layer for `cdmc_fit`- and `cdmc_dr_fit`-backed
  targets that numerically linearizes perturbations requiring the staged
  baseline, joint objective, or cross-fitted DR pipeline to be re-solved.
- Optional CBPS-based continuous-treatment balancing weights that can be passed
  directly to the main fit, used internally inside the DR layer, or passed to
  the post-fit dose-response layer.
- Linear and spline-based post-fit dose-response estimation with optional
  observation weights in the second stage.
- Basic simulation helper for smoke tests and development.
- Basic tests for panel validation, washout masking, and approximate recovery of
  a simple lagged effect.
- Manual pages for the current exported API.

## Design choices already made

### Scope choices

These were explicitly fixed for v0.1:

- 0 is the only control state, while active doses may be positive or negative.
- The package targets research use rather than CRAN-ready completeness.

Why:

- This removes ambiguity about control eligibility and keeps masking logic
  tractable.
- Research-use scope allows the core estimator to be stress-tested before adding
  broader API and inferential commitments.

### Estimation choices

- The main estimator defaults to a two-stage reversible-treatment fit, but now
  also supports an optional joint objective when a global parametric treatment
  basis is acceptable.
- The staged path uses eligible zero-dose cells to estimate the untreated
  baseline.
- The unweighted low-rank step uses softImpute with type = "svd".
- Weighted fits use a proximal singular-value-thresholding update on the masked
  weighted loss rather than the unweighted softImpute call.
- Lambda can now be chosen either heuristically or by blocked cross-validation.
- The heuristic default, when lambda = NULL and lambda_selection = "heuristic",
  is based on a fraction of softImpute::lambda0() after partialling out additive
  nuisance terms.
- The blocked cross-validation path masks contiguous within-unit control blocks
  while preserving at least one eligible control per unit and time period.
- The built-in second-stage effect model does not include an intercept.
- Lagging effects are represented by dose_lag0, dose_lag1, and so on up to
  lag_order.

Why:

- A staged estimator remains easier to validate against the reversible-treatment
  identification logic in the paper and is therefore still the default, while
  the optional joint mode covers the manuscript's global regularized objective
  when that specification is acceptable.
- softImpute is an established and stable implementation of the SVT or
  Soft-Impute strategy described in the paper.
- The built-in stage-2 layer stays intercept-free so its fitted component
  remains anchored at the zero-dose control state.

### Safety and validation choices

- The package fails fast when a unit or time period loses all eligible zero-dose
  observations after washout masking.
- The package requires complete data in the modeled variables.
- Tiny numerical doses are coerced to 0 using zero_tolerance.
- Prediction is currently restricted to the training panel.
- The current carryover diagnostic is residual-based and does not refit the
  baseline model under masked exit periods.
- The bootstrap resamples whole unit trajectories and reassigns synthetic unit
  identifiers inside each bootstrap sample to preserve each sampled unit's
  observed time path before internal panel padding is rebuilt.

Why:

- These checks prevent silent misuse of the matrix completion step.
- Restricting prediction avoids pretending that out-of-sample panel completion is
  implemented when it is not.

## What is not implemented yet

The package is not yet a full implementation of the manuscript target.

Missing core features:

- Richer learner-based or balancing-oriented fold-specific
  treatment-assignment estimation inside the DR layer beyond the current
  Gaussian, kernel, forest, and CBPS options, especially boosting-, stacking-,
  or balancing-oriented nuisance learners.
- Formal sensitivity analysis beyond the current bounded-bias and replay-based
  helpers, which now include dense and sparse contamination-bias bounds under
  both cellwise and energy budgets plus a local optimization-refit layer,
  especially broader identification-robust sensitivity models.

## Known limitations of the current implementation

- The default `cdmc_fit()` path is still sequential, while the optional joint
  path is limited to the built-in global linear or spline treatment basis and
  therefore remains less flexible than the staged-plus-postfit workflow for
  complex reversible-treatment response shapes.
- The carryover diagnostic is a screening tool, not a full inferential routine
  that re-estimates the baseline under masked exit periods.
- The refit-based placebo and carryover diagnostics currently hold lambda fixed
  at the value chosen in the original fit instead of rerunning tuning inside the
  diagnostic refit.
- The new `cdmc_dr_fit()` layer is limited to the built-in linear stage-2 target.
  Its internal treatment-assignment nuisance model now supports Gaussian and
  kernel generalized propensity score variants with linear, spline, GAM, tree,
  or forest mean models plus a fold-local CBPS balancing option, all with
  optional time fixed effects, but more flexible boosting or stacking
  assignment models still require externally supplied weights.
- The new `cdmc_scia_test()` helper is a screening diagnostic, not a formal DAG
  verification routine or a proof of sequential conditional independence.
- The current bootstrap reports percentile intervals for a small set of summary
  statistics only; it now adds max-t simultaneous bands for finite collections
  of dynamic pathwise summaries, but it still does not provide continuous
  functional confidence bands or asymptotic uncertainty theory.
- The current effect model uses only treated-history cells with fully observed
  lag history and nonzero cumulative lag exposure.
- The current sensitivity tooling now combines replay-based scan summaries with
  formal bounded-bias helpers for the stored stage-2 response layer, the
  baseline-counterfactual layer, and a local optimization-refit layer for
  `cdmc_fit`- and `cdmc_dr_fit`-backed targets, including direct
  joint-objective `cdmc_fit()` objects through either their stored post-fit
  treatment maps or a numerical refit Jacobian and cross-fitted DR objects
  through stored fold reuse, but it still does not cover broader
  identification-robust sensitivity models.
- The default lambda rule is heuristic, not validated by blocked
  cross-validation.
- Convergence is based on relative change in the low-rank matrix only.
- The package currently exposes bootstrap percentile uncertainty summaries, but
  not asymptotic inference or joint uncertainty for diagnostic estimands.
- Manual pages exist now, but there are still no vignettes or worked empirical
  examples.

## Important alignment note

The original goal is broader than the current code. To stay aligned with the
goal, future work should keep the following separation clear:

- The stage-1 baseline engine should continue to target untreated potential
  outcomes under 0 dose.
- The stage-2 layer should grow toward a flexible continuous dose-response model
  with lag structure, not drift toward a binary-treatment interface.
- Diagnostics should explicitly test the assumptions created by reversibility and
  lagging effects, especially carryover and insufficient zero-dose support.

## Recommended next implementation order

1. Extend the formal sensitivity layer beyond the current exact stored-map and
  local optimization-refit helpers toward richer identification-robust
  sensitivity analyses.
2. Add richer learner-based treatment-assignment estimation inside the DR layer
  beyond the current Gaussian, kernel, forest, and CBPS options, especially
  boosting- or stacking-style nuisance learners.

## Files to know first

- [R/fit.R](R/fit.R): public estimator API and S3 methods.
- [R/panel.R](R/panel.R): panel validation, eligible-control masking, lag-array
  construction.
- [R/nuisance.R](R/nuisance.R): additive nuisance-term estimation.
- [R/baseline.R](R/baseline.R): low-rank baseline estimation and default penalty.
- [R/effect.R](R/effect.R): linear lagged second-stage effect model.
- [R/simulate.R](R/simulate.R): synthetic data generator.
- [tests/testthat](tests/testthat): current automated checks.

## Handoff note for another agent

If another agent continues this work, the first question should always be:

"Does the next change make the library more faithful to matrix completion for
continuous, reversible treatment with possible lagging effects, or is it just API
surface growth?"

Prefer changes that strengthen the estimator contract, diagnostics, and
identification story before expanding convenience features.
