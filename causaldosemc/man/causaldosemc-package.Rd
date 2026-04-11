\name{causaldosemc-package}
\alias{causaldosemc}
\alias{causaldosemc-package}
\title{causaldosemc: Matrix Completion for Continuous Treatment Panels}
\description{
Research-use tools for causal estimation in panel data with a
continuous treatment and 0 as the control state.
}
\details{
The package is organized around a staged estimator with an optional joint
objective for the main \code{cdmc_fit()} path.

In stage 1, the untreated baseline surface is estimated from eligible zero-dose
observations using additive unit and time effects, optional covariates, and a
low-rank matrix completion step. Unweighted fits use \code{softImpute}, while
weighted fits switch to a weighted proximal singular-value-thresholding routine.

In stage 2, residual treatment effects are estimated from the difference between
observed outcomes and the imputed baseline. Version 0.1 provides either a
linear or spline contemporaneous or distributed-lag dose-response model, can
optionally reuse observation weights from the main fit, and provides a
residual-based carryover diagnostic for post-exit zero-dose periods.

When the main fit is run with \code{objective = "joint"}, the package instead
fits the low-rank surface, additive nuisance terms, and the chosen linear or
spline treatment basis jointly over the observed estimation sample while still
returning the implied baseline counterfactual and fitted treatment component.

The package also provides unit-level bootstrap inference that preserves each
sampled unit's full time path and can optionally rerun lambda selection inside
each bootstrap replicate. The same bootstrap entry point can now summarize both
the main estimator, post-fit dose-response objects, the cross-fitted DR layer,
and the refit-based diagnostic objects.

For diagnostics, the package now includes both a residual carryover screen and
refit-based placebo and carryover tests that hold out targeted control windows
from stage-1 fitting, a joint placebo pretrend diagnostic across multiple
windows, a SCIA-oriented lagged-outcome screen for dynamic assignment, and a
TOST-style equivalence helper for assessing whether those diagnostic effects
are practically negligible.

For robustness checks, the package provides both a replay-based sensitivity
helper, \code{cdmc_sensitivity_scan()}, which varies washout, zero-dose
tolerance, and weighting choices for a fitted main estimator, DR fit,
post-fit dose-response object, or dynamic estimand, records support failures
without aborting the full scan, and summarizes zero-dose support, washout
dependence, and weight dependence relative to the stored source fit, and a
formal bounded-bias helper, \code{cdmc_sensitivity_bounds()}, which reports
exact worst-case intervals for supported linear post-fit summaries under either
a stage-2-response or baseline-counterfactual bounded perturbation model,
optionally restricted to sparse contamination fractions and evaluated under
either cellwise or energy perturbation budgets, plus a local
optimization-refit layer for \code{cdmc_fit}- and \code{cdmc_dr_fit}-backed
targets that numerically linearizes perturbations requiring the staged
baseline, joint objective, or cross-fitted DR pipeline to be re-solved.

For dose-response estimation beyond the built-in main-fit stage-2 model, the package
also provides a post-fit linear or spline dose-response layer with optional
observation weights. Continuous-treatment balancing weights can be estimated
through \code{cdmc_cbps_weights()} when \code{CBPS} is available.

For dynamic or pathwise summaries, the package also provides
\code{cdmc_dynamic_estimand()}, which evaluates user-supplied dose histories on
top of a fitted \code{cdmc_fit()}, \code{cdmc_dr_fit()}, or
\code{cdmc_dose_response()} object and can be passed directly to
\code{cdmc_bootstrap()} for pointwise percentile inference and finite-sample
max-t simultaneous bands over multi-statistic summary vectors.

The package also includes a cross-fitted doubly robust-style linear estimator,
\code{cdmc_dr_fit()}, which combines out-of-fold baseline prediction with
either internally estimated Gaussian or kernel GPS weights, fold-local CBPS
balancing weights, or externally supplied weights. The internal Gaussian and
kernel GPS paths can use raw linear terms, natural-spline expansions, GAM
smooths, regression trees, or regression forests in the assignment mean model,
and the kernel path replaces the Gaussian density assumption with
kernel-smoothed residual and marginal dose densities. The CBPS path currently
supports linear or spline covariate formulas.

The current implementation is intentionally narrow and aligned to the first
research milestone of the broader library goal: continuous, reversible treatment
panels with possible lagging effects. It supports both balanced and unbalanced
panels by padding missing unit-time combinations internally, still requires
complete data for the modeled variables on the observed rows, and uses zero to
define the control state.
}
\seealso{
\code{\link{cdmc_fit}} for estimation, \code{\link{cdmc_bootstrap}} for the
current bootstrap workflow, \code{\link{cdmc_carryover_test}} for the current
carryover screen, \code{\link{cdmc_placebo_test}} and
\code{\link{cdmc_carryover_refit_test}}, \code{\link{cdmc_joint_placebo_test}},
and \code{\link{cdmc_equivalence_test}} for refit-based diagnostics,
aggregated placebo checks, and equivalence checks, and
\code{\link{cdmc_dynamic_estimand}}, \code{\link{cdmc_dose_response}}, \code{\link{cdmc_cbps_weights}},
\code{\link{cdmc_dr_fit}}, \code{\link{cdmc_sensitivity_scan}},
\code{\link{cdmc_sensitivity_bounds}}, and
\code{\link{simulate_cdmc_data}} for post-fit dose-response work, the
cross-fitted DR layer, replay-based and formal sensitivity checks, and small
simulated examples.
}
\keyword{package}

