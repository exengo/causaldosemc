\name{cdmc_bootstrap}
\alias{cdmc_bootstrap}
\alias{print.cdmc_bootstrap}
\title{Unit-level bootstrap inference for causaldosemc fits}
\description{
Resample whole unit trajectories with replacement, refit the estimator on each
bootstrap sample, and compute percentile intervals for selected summary
statistics, with optional max-t simultaneous bands for dynamic estimands.
}
\usage{
cdmc_bootstrap(
  object,
  n_boot = 100L,
  statistics = c("coefficients", "average_tau"),
  rerun_tuning = NULL,
  conf_level = 0.95,
  simultaneous = NULL,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  workers = 1L,
  seed = NULL,
  verbose = FALSE
)

\method{print}{cdmc_bootstrap}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit}, \code{cdmc_dynamic_estimand},
  \code{cdmc_dose_response}, or \code{cdmc_dr_fit} object, or a
  supported diagnostic object such as \code{cdmc_placebo_test},
  \code{cdmc_carryover_test}, \code{cdmc_carryover_refit_test},
  \code{cdmc_joint_placebo_test}, \code{cdmc_equivalence_test}, or
  \code{cdmc_scia_test}.}

  \item{n_boot}{Number of bootstrap replications.}

  \item{statistics}{Which summaries to collect from each bootstrap fit.
  For \code{cdmc_fit} objects, \code{"coefficients"} collects the linear
  second-stage effect coefficients and \code{"average_tau"} collects the
  average residual treatment effect over active nonzero-dose observations. For
  \code{cdmc_dynamic_estimand} objects, \code{"estimate"} collects the named
  dynamic or pathwise summaries stored on the object, such as individual path
  effects, path contrasts, slopes, or their weighted mean. For
  \code{cdmc_dose_response} objects, \code{"coefficients"} collects the
  post-fit dose-response coefficients and \code{"prediction"} collects fitted
  response or slope values at user-supplied dose histories. For
  \code{cdmc_dr_fit} objects, \code{"coefficients"} collects the final DR
  linear coefficients, \code{"average_tau_dr"} collects the mean augmented DR
  pseudo outcome over the DR sample cells, and \code{"average_tau_linear"}
  collects the mean fitted DR linear effect over those same cells.
  \code{"lag_average_tau_dr"} collects lag-specific average fitted
  contributions, and \code{"dose_contrast_dr"} collects a linear contrast
  between two dose histories. For placebo and carryover diagnostics,
  \code{"mean_tau"} collects the mean pseudo effect. For joint placebo
  diagnostics, \code{"period_mean_tau"} collects the placebo means when the
  diagnostic was built from one window per period, \code{"window_mean_tau"}
  collects the placebo means for named or multi-period windows, and
  \code{"joint_p_value"} collects the aggregated joint p value. For
  \code{cdmc_equivalence_test} objects, \code{"mean_tau"} collects the mean
  diagnostic effect and \code{"equivalence_p_value"} collects the TOST-style
  equivalence p value. For \code{cdmc_scia_test} objects, \code{"f_statistic"}
  and \code{"p_value"} collect the global incremental lagged-outcome
  statistics, while \code{"restriction_p_value"} and
  \code{"restriction_adj_p_value"} collect the block-level raw and adjusted
  SCIA restriction p values. When
  \code{statistics} is omitted, a class-specific default set is used.}

  \item{rerun_tuning}{Logical flag controlling whether lambda selection should be
  rerun in each bootstrap replicate when the original fit used heuristic or
  cross-validated lambda selection. If \code{NULL}, the default is to rerun
  tuning unless the original fit used a fixed lambda. For empirical workflows
  built on automatic tuning, keeping this at its default behavior or setting it
  explicitly to \code{TRUE} is recommended so the bootstrap reflects retuning
  rather than treating the stored lambda as fixed.}

  \item{conf_level}{Confidence level used for the reported percentile interval.}

  \item{simultaneous}{Logical flag controlling whether a max-t simultaneous
  confidence band should also be reported. This is currently supported for
  \code{cdmc_dynamic_estimand} bootstrap objects only. If \code{NULL}, the
  default is to report simultaneous bands automatically when bootstrapping a
  dynamic estimand with more than one requested summary.}

  \item{contrast_history}{Optional single-row dose history used when
  \code{"dose_contrast_dr"} is requested for a \code{cdmc_dr_fit} object.
  Supply either a numeric vector, a matrix, or a data frame with lag columns.
  When omitted, the contrast defaults to a unit dose at every included lag.}

  \item{reference_history}{Optional reference dose history used when
  \code{"dose_contrast_dr"} is requested for a \code{cdmc_dr_fit} object. When
  omitted, the reference defaults to the all-zero history.}

  \item{prediction_dose}{Optional vector of contemporaneous doses used when
  \code{"prediction"} is requested for a \code{cdmc_dose_response} object.
  When \code{lag_order > 0}, unspecified lagged doses are set to 0.}

  \item{prediction_history}{Optional matrix or data frame of explicit lagged
  dose histories used when \code{"prediction"} is requested for a
  \code{cdmc_dose_response} object.}

  \item{prediction_type}{Prediction target used when \code{"prediction"} is
  requested for a \code{cdmc_dose_response} object. \code{"response"} returns
  the fitted nonlinear response and \code{"slope"} returns the finite-difference
  derivative with respect to \code{dose_lag0}.}

  \item{workers}{Number of workers used for bootstrap refits. Set
  \code{workers = 1} for sequential execution (default). Values greater than 1
  enable parallel bootstrap refits via multicore execution on supported
  platforms.}

  \item{seed}{Optional random seed used before bootstrapping.}

  \item{verbose}{Logical flag indicating whether replicate progress messages
  should be printed.}

  \item{x}{A \code{cdmc_bootstrap} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The bootstrap operates at the unit level. Entire unit histories are sampled with
replacement, and each sampled history is assigned a synthetic unit identifier so
that the balanced-panel contract remains valid in the bootstrap dataset.

When the source estimator used automatic lambda selection, empirical workflows
should generally keep \code{rerun_tuning = TRUE} so bootstrap uncertainty
propagates the retuning step instead of conditioning on one stored penalty.

For \code{cdmc_dr_fit()} objects, the bootstrap refits the full cross-fitted DR
pipeline inside each replicate. When the DR object used the built-in internal
weighting nuisance model, whether Gaussian, kernel-smoothed, or fold-local
CBPS, that weighting model is refit within each bootstrap sample. When the DR object used external
weights, the stored normalized weight column is resampled with the panel and
passed back into the DR refit.

For supported diagnostic objects, the bootstrap resamples the parent fitted
panel, refits the original \code{cdmc_fit()} model, and reruns the requested
diagnostic on that bootstrap replicate. This keeps placebo, carryover,
equivalence, SCIA, and joint placebo summaries aligned with the estimator
configuration that generated the original diagnostic, including whether the
diagnostic refits reuse the selected lambda or rerun the stored automatic
tuning rule, and for joint placebo objects it preserves the original named
window protocol.

For \code{cdmc_dose_response()} objects, the bootstrap resamples the parent
\code{cdmc_fit()} panel, refits that parent model, recreates the stored
post-fit dose-response object, and then evaluates either coefficients or the
requested response or slope functionals on each bootstrap replicate.

For \code{cdmc_dynamic_estimand()} objects, the bootstrap resamples the parent
fitted model, recreates the stored dynamic estimand on the bootstrap fit, and
collects the stored pathwise summaries. This provides pointwise percentile
intervals for user-supplied dynamic responses, contrasts, or slopes. When
\code{simultaneous = TRUE}, or when the automatic default is triggered for a
multi-statistic dynamic estimand, the bootstrap also reports max-t simultaneous
bands across the full requested summary vector.

When \code{workers > 1}, bootstrap refits are distributed across workers and
then reassembled in the original replicate order. If \code{seed} is supplied,
the bootstrap uses deterministic per-replicate seeds so repeated calls with the
same inputs return the same draws under a fixed worker configuration.

The current implementation returns percentile intervals for a small set of
summary statistics. It is intended as a first inference layer for the research
workflow and now supports simultaneous dynamic uncertainty control over finite
sets of user-supplied path summaries, but it does not yet provide continuous
functional confidence bands or asymptotic inference.
}
\value{
An object of class \code{cdmc_bootstrap} containing:

\item{estimates}{The statistic values from the original fitted object.}

\item{draws}{A matrix of bootstrap draws, with one row per bootstrap replicate.}

\item{summary}{A data frame containing the original estimate, bootstrap mean,
bootstrap standard deviation, percentile interval, and number of successful draws
for each statistic. When simultaneous bands are requested, this table also
includes \code{simultaneous_conf_low} and \code{simultaneous_conf_high}.}

\item{simultaneous}{Either \code{NULL} when no simultaneous bands were
requested, or a list describing the max-t band calculation, including the
critical value and the number of jointly successful bootstrap draws.}

\item{n_boot, conf_level, workers, parallel, rerun_tuning, lambda_method, n_failures}{Bootstrap
metadata describing the resampling run.}
}
\examples{
\donttest{
panel <- simulate_cdmc_data(
  n_units = 10,
  n_times = 8,
  beta = 1,
  lag_beta = 0.3,
  seed = 1
)

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lambda = NULL,
  lambda_selection = "cv",
  cv_rounds = 3L,
  cv_block_size = 2L,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

cdmc_bootstrap(fit, n_boot = 10, rerun_tuning = TRUE, seed = 1)

dr_fit <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  lambda = 0.2,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

cdmc_bootstrap(
  dr_fit,
  n_boot = 10,
  statistics = c("coefficients", "average_tau_dr", "average_tau_linear"),
  seed = 1
)

cdmc_bootstrap(
  dr_fit,
  n_boot = 10,
  statistics = c("lag_average_tau_dr", "dose_contrast_dr"),
  contrast_history = c(dose_lag0 = 1, dose_lag1 = 0.5),
  reference_history = c(dose_lag0 = 0, dose_lag1 = 0),
  seed = 1
)

placebo <- cdmc_placebo_test(fit, periods = -2:0)
cdmc_bootstrap(placebo, n_boot = 10, seed = 1)

dose_response <- cdmc_dose_response(fit, model = "spline", lag_order = 1, df = 3)
cdmc_bootstrap(
  dose_response,
  n_boot = 10,
  statistics = c("coefficients", "prediction"),
  prediction_dose = c(-0.5, 0.5),
  prediction_type = "response",
  seed = 1
)

dynamic_path <- cdmc_dynamic_estimand(
  fit,
  history = data.frame(dose_lag0 = c(1, 0.5), dose_lag1 = c(0.5, 1)),
  type = "contrast",
  aggregate = "both"
)
cdmc_bootstrap(
  dynamic_path,
  n_boot = 10,
  statistics = "estimate",
  simultaneous = TRUE,
  seed = 1
)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dynamic_estimand}},
\code{\link{cdmc_dose_response}},
\code{\link{cdmc_dr_fit}},
\code{\link{cdmc_placebo_test}}, \code{\link{cdmc_carryover_test}},
\code{\link{cdmc_carryover_refit_test}}, \code{\link{cdmc_equivalence_test}},
\code{\link{cdmc_scia_test}}, and
\code{\link{cdmc_joint_placebo_test}} for model fitting and diagnostics.
}