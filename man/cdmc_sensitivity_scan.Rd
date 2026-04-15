\name{cdmc_sensitivity_scan}
\alias{cdmc_sensitivity_scan}
\alias{print.cdmc_sensitivity_scan}
\alias{summary.cdmc_sensitivity_scan}
\alias{print.summary.cdmc_sensitivity_scan}
\title{Replay a fitted estimator over washout, support, and weight scenarios}
\description{
Run a replay-based sensitivity-analysis sweep for a fitted
\code{cdmc_fit()}, \code{cdmc_dr_fit()}, \code{cdmc_dose_response()}, or
\code{cdmc_dynamic_estimand()} object by replaying the stored source model over
alternative washout values, zero-dose tolerances, and weighting scenarios,
while summarizing zero-dose support, objective-specific optimization support,
washout dependence, and weight dependence relative to the stored fit.
}
\usage{
cdmc_sensitivity_scan(
  object,
  washout_grid = cdmc_default_washout_grid(object),
  zero_tolerance_grid = cdmc_default_zero_tolerance_grid(object),
  weight_specs = NULL,
  include_unweighted = NULL,
  statistics = NULL,
  rerun_tuning = NULL,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  verbose = FALSE
)

\method{print}{cdmc_sensitivity_scan}(x, ...)
\method{summary}{cdmc_sensitivity_scan}(object, ...)

\method{print}{summary.cdmc_sensitivity_scan}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit}, \code{cdmc_dr_fit},
  \code{cdmc_dose_response}, or \code{cdmc_dynamic_estimand} object. For
  downstream objects, the helper replays the stored parent fit and then
  reconstructs the downstream target on each replayed scenario.}

  \item{washout_grid}{One or more nonnegative integers describing the washout
  values to replay. By default, the function evaluates the stored source-fit
  washout value together with its immediate neighbors when feasible.}

  \item{zero_tolerance_grid}{One or more nonnegative numeric values describing
  the zero-dose tolerance values to replay. By default, the function evaluates
  the stored source-fit tolerance together with one smaller and one larger
  order-of-magnitude value.}

  \item{weight_specs}{Optional weighting scenarios to compare against the
  stored source model. For \code{cdmc_fit} sources, supply either a named list
  or a single weight specification that would also be accepted by
  \code{cdmc_fit(weights = ...)}. For \code{cdmc_dr_fit} sources, each entry may
  either be a raw external weight specification or a named list overriding
  fields such as \code{weights}, \code{weight_method}, \code{gps_model},
  \\code{gps_bandwidth}, \\code{gps_forest_trees}, \\code{gps_forest_mtry},
  \\code{gps_forest_min_node_size}, \\code{cbps_method},
  \\code{adaptive_balance_methods}, \\code{entropy_balance_degree},
  \\code{kernel_balance_centers}, or related internal weighting replay
  settings.}

  \item{include_unweighted}{Logical indicating whether an explicit unweighted
  replay should be included alongside the stored weighting setup. By default,
  this is enabled only when a \code{cdmc_fit} source used explicit weights. For
  \code{cdmc_dr_fit} sources with \code{weight_method = "external"}, setting
  \code{include_unweighted = TRUE} adds a replay under the built-in
  \code{"gaussian_gps"} weighting path.}

  \item{statistics}{Optional replay summaries to collect from each scenario.
  When omitted, the helper uses the same class-specific defaults as
  \code{cdmc_bootstrap()}. This means coefficient and average treatment-effect
  summaries for \code{cdmc_fit()} or \code{cdmc_dr_fit()}, coefficient summaries
  for \code{cdmc_dose_response()}, and stored pathwise summaries for
  \code{cdmc_dynamic_estimand()}.}

  \item{rerun_tuning}{Logical flag controlling whether lambda selection should
  be rerun in each replay when the original fit used heuristic or
  cross-validated tuning. If \code{NULL}, the helper follows the same replay
  logic used by \code{cdmc_bootstrap()}.}

  \item{contrast_history}{Optional single-row dose history used when replaying
  \code{"dose_contrast_dr"} summaries for a \code{cdmc_dr_fit} object.}

  \item{reference_history}{Optional reference dose history used when replaying
  \code{"dose_contrast_dr"} summaries for a \code{cdmc_dr_fit} object.}

  \item{prediction_dose}{Optional vector of contemporaneous doses used when
  replaying \code{"prediction"} summaries for a \code{cdmc_dose_response}
  object.}

  \item{prediction_history}{Optional matrix or data frame of explicit lagged
  dose histories used when replaying \code{"prediction"} summaries for a
  \code{cdmc_dose_response} object.}

  \item{prediction_type}{Prediction target used when replaying
  \code{"prediction"} summaries for a \code{cdmc_dose_response} object.}

  \item{verbose}{Logical flag indicating whether replay progress messages should
  be printed.}

  \item{x}{A \code{cdmc_sensitivity_scan} object.}

  \item{...}{Unused additional arguments for the print or summary methods.}
}
\details{
The helper reuses the stored public replay specification from the fitted source
model and refits that source over the requested grid. For each scenario it
records observed and eligible zero-dose support, the optimization sample used by
the replayed objective, whether the optimization-support checks passed, and
class-appropriate summary statistics from the replayed target.

For \code{cdmc_fit()} objects, the target is the replayed fit itself.

When a \code{cdmc_fit()} source uses the default staged objective, the reported
support metrics match the eligible zero-dose optimization sample. When the
source uses \code{objective = "joint"}, the scan distinguishes between the
eligible zero-dose anchor summary and the actual observed-panel optimization
sample through columns such as \code{objective}, \code{optimization_sample},
\code{optimization_cells}, and \code{optimization_fraction}. In that joint
case, replay viability follows the optimization sample rather than the zero-dose
anchor summary.

For \code{cdmc_dr_fit()} objects, the helper replays the full cross-fitted DR
pipeline and then summarizes the replayed DR coefficients or DR treatment
effect summaries.

For \code{cdmc_dose_response()} objects, the helper replays the stored parent
fit and reconstructs the stored post-fit dose-response target before collecting
coefficient or prediction summaries.

For \code{cdmc_dynamic_estimand()} objects, the helper replays the stored source
fit, reconstructs any nested post-fit dose-response object when necessary, and
then reevaluates the stored pathwise dynamic summaries.

If a scenario loses support, the helper records the support failure and keeps
going instead of stopping the whole scan. This makes it useful for checking how
fragile the estimator is to more aggressive washout choices, alternative
zero-dose thresholds, or alternative weighting configurations.

In addition to the row-level replay table, the returned object stores the
stored-fit reference scenario, a compact support table with both zero-dose and
optimization-sample summaries plus reference deltas, and aggregated summaries
that report how much support and fit statistics vary across washout and
weighting dimensions.

This is still a replay-based first sensitivity layer rather than a formal
identification-robust sensitivity framework.
}
\value{
An object of class \code{cdmc_sensitivity_scan} containing:

\item{results}{A data frame with one row per replay scenario. The table
includes scenario labels, objective and optimization-sample labels, zero-dose
anchor summaries, optimization-support summaries, fit-success indicators, any
captured fit error, replayed effect summaries, and reference-relative deltas.}

\item{reference}{A named list describing the stored fitted scenario that serves
as the default reference point for the scan summaries.}

\item{support_summary}{A compact scenario-level table focused on zero-dose
support and optimization-sample support, including support fractions, support
losses, and deltas relative to the stored fit.}

\item{washout_summary}{An aggregated table that summarizes how zero-dose
support and fitted effect summaries vary across washout choices within each
fixed weighting and zero-tolerance configuration.}

\item{weight_summary}{An aggregated table that summarizes how zero-dose
support and fitted effect summaries vary across weighting scenarios within each
fixed washout and zero-tolerance configuration.}

\item{base_object, source_object}{The original target object passed to the
helper and the replayed source fit used to generate sensitivity scenarios.}

\item{summary(object)}{Returns a list containing the same reference,
support, washout, and weight summaries in a compact summary object.}
}
\examples{
\donttest{
panel <- simulate_cdmc_data(
  n_units = 12,
  n_times = 10,
  beta = 1,
  lag_beta = 0.2,
  seed = 1
)
panel$w <- 0.5 + abs(panel$x1)

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weights = "w",
  lambda = 0.2,
  rank_max = 3,
  washout = 1,
  lag_order = 1,
  seed = 1
)

scan <- cdmc_sensitivity_scan(
  fit,
  washout_grid = 0:2,
  include_unweighted = TRUE
)

print(scan)
summary(scan)

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

dr_scan <- cdmc_sensitivity_scan(dr_fit, washout_grid = 0:1)
summary(dr_scan)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}},
\code{\link{cdmc_dose_response}}, \code{\link{cdmc_dynamic_estimand}},
\code{\link{cdmc_sensitivity_bounds}}, and \code{\link{cdmc_bootstrap}}.
}