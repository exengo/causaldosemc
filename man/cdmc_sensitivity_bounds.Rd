\name{cdmc_sensitivity_bounds}
\alias{cdmc_sensitivity_bounds}
\alias{print.cdmc_sensitivity_bounds}
\alias{summary.cdmc_sensitivity_bounds}
\alias{print.summary.cdmc_sensitivity_bounds}
\title{Compute formal bounded-bias sensitivity intervals for fitted summaries}
\description{
Report formal bounded-bias intervals for supported scalar summaries under a
bounded additive perturbation model on the stored stage-2 response layer, the
baseline-counterfactual layer, or, for \code{cdmc_fit}- or
\code{cdmc_dr_fit}-backed targets, the source optimization sample itself
through a local refit Jacobian.
}
\usage{
cdmc_sensitivity_bounds(
  object,
  statistics = NULL,
  gamma_grid = c(0, 0.25, 0.5, 1),
  contamination_fraction_grid = 1,
  perturbation_constraint = c("cellwise", "energy"),
  perturbation_scope = c("cell", "unit", "time"),
  perturbation_layer = c("stage2", "baseline", "optimization"),
  scale = c("residual_sd", "response_sd", "manual"),
  scale_value = NULL,
  null = 0,
  contrast_history = NULL,
  reference_history = NULL,
  prediction_dose = NULL,
  prediction_history = NULL,
  prediction_type = c("response", "slope"),
  refit_step = NULL
)

\method{print}{cdmc_sensitivity_bounds}(x, ...)
\method{summary}{cdmc_sensitivity_bounds}(object, ...)
\method{print}{summary.cdmc_sensitivity_bounds}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit}, \code{cdmc_dr_fit},
  \code{cdmc_dose_response}, or \code{cdmc_dynamic_estimand} object.}

  \item{statistics}{Optional scalar summaries to bound. When omitted, the
  helper uses the same class-specific defaults as \code{cdmc_bootstrap()}.}

  \item{gamma_grid}{One or more nonnegative sensitivity-parameter values.
  Each \code{gamma} value is converted to a raw perturbation bound of
  \code{gamma * scale_value}.}

  \item{contamination_fraction_grid}{One or more contamination fractions in
  \eqn{[0, 1]}. Each value limits the adversarial perturbation to that share of
  the underlying sensitivity sample, so \code{1} recovers the dense worst-case
  \eqn{\ell_\infty} model and smaller values report sparse contamination
  bounds based on the largest loading entries only.}

  \item{perturbation_constraint}{How the adversarial perturbation is bounded on
  the selected contamination set. \code{"cellwise"} uses the original
  coordinatewise bound \eqn{|u_i| \le \Gamma s} and therefore produces
  \eqn{\ell_1}-loading radii. \code{"energy"} instead bounds the Euclidean
  energy of the perturbation over the contaminated cells and produces
  \eqn{\ell_2}-loading radii.}

  \item{perturbation_scope}{Which grouped sample elements the adversary may
  perturb. \code{"cell"} treats each unit-time cell separately. \code{"unit"}
  aggregates the loading vector within unit before applying the chosen dense or
  sparse perturbation budget, so the perturbation is interpreted as a
  unit-level shock. \code{"time"} analogously aggregates within time period and
  reports period-level shocks.}

  \item{perturbation_layer}{Which stored post-fit layer is perturbed.
  \code{"stage2"} perturbs the stage-2 response or pseudo-outcome sample used
  by the stored linear map. \code{"baseline"} perturbs the underlying
  baseline-counterfactual layer instead. For staged and joint main fits and for
  post-fit dose-response objects, the resulting bound geometry is often the
  same because the fitted treatment summaries are linear in
  \eqn{Y - \hat{Y}(0)}. For \code{cdmc_dr_fit()}, the baseline layer is
  distinct because baseline errors are reweighted before the pseudo-outcome
  regression. \code{"optimization"} perturbs the source optimization sample
  itself and recomputes the downstream summary through a local numerical refit
  Jacobian. This layer is currently available only for targets backed by
  \code{cdmc_fit()} or \code{cdmc_dr_fit()}, including post-fit
  \code{cdmc_dose_response()} objects and compatible
  \code{cdmc_dynamic_estimand()} targets that can be rebuilt from those source
  fits. For nonlinear post-fit dose-response models such as \code{"gam"},
  \code{"tree"}, and \code{"forest"}, this is the supported formal
  sensitivity layer.}

  \item{scale}{How to calibrate \code{scale_value}. \code{"residual_sd"}
  uses the standard deviation of the relevant residual sample for the chosen
  perturbation layer, \code{"response_sd"} uses the standard deviation of the
  corresponding perturbed sample itself, and \code{"manual"} uses the
  user-supplied \code{scale_value}.}

  \item{scale_value}{Optional positive scalar used when
  \code{scale = "manual"}. Ignored otherwise.}

  \item{null}{Null value used when reporting breakdown points. Supply either a
  single scalar, an unnamed vector matching the number of reported statistics,
  or a named numeric vector keyed by statistic name.}

  \item{contrast_history}{Optional single-row dose history used when bounding
  \code{"dose_contrast_dr"} summaries for a \code{cdmc_dr_fit} object.}

  \item{reference_history}{Optional reference dose history used when bounding
  \code{"dose_contrast_dr"} summaries for a \code{cdmc_dr_fit} object.}

  \item{prediction_dose}{Optional vector of contemporaneous doses used when
  bounding \code{"prediction"} summaries for a \code{cdmc_dose_response}
  object.}

  \item{prediction_history}{Optional matrix or data frame of explicit lagged
  dose histories used when bounding \code{"prediction"} summaries for a
  \code{cdmc_dose_response} object.}

  \item{prediction_type}{Prediction target used when bounding
  \code{"prediction"} summaries for a \code{cdmc_dose_response} object.}

  \item{refit_step}{Optional positive step size used when
  \code{perturbation_layer = "optimization"}. When omitted, the helper uses a
  default proportional to the resolved perturbation scale. Ignored for the
  exact stored-map layers.}

  \item{x}{A \code{cdmc_sensitivity_bounds} object or a
  \code{summary.cdmc_sensitivity_bounds} object.}

  \item{...}{Unused additional arguments for the print or summary methods.}
}
\details{
The helper imposes a bounded additive perturbation model on the relevant
sample. Under \code{perturbation_constraint = "cellwise"}, the constraint is
coordinatewise,
\deqn{|u_i| \le \Gamma s,}
where \eqn{\Gamma} is taken from \code{gamma_grid} and \eqn{s} is the resolved
\code{scale_value}. For any supported scalar summary that can be written as a
linear functional of that perturbed sample,
\deqn{\hat{\theta} = h^\top y,}
the helper reports worst-case bounds of the form
\deqn{\hat{\theta} \pm \Gamma s \|h\|_1.}

Under \code{perturbation_constraint = "energy"}, the perturbation is instead
bounded by an \eqn{\ell_2} energy budget over the contaminated cells, so the
same linear summary receives worst-case bounds of the form
\deqn{\hat{\theta} \pm \Gamma s \|h\|_2.}

When \code{contamination_fraction_grid} contains values below \code{1}, the
helper also reports sparse contamination bounds where the adversary may perturb
only the largest-loading fraction of cells. In that case the interval radius is
computed from the largest-loading subset rather than the full loading vector,
using either the corresponding \eqn{\ell_1} or \eqn{\ell_2} geometry implied by
\code{perturbation_constraint}. \code{contamination_fraction_grid = 1}
recovers the original dense worst-case model exactly.

When \code{perturbation_scope = "unit"} or \code{"time"}, the helper first
aggregates the loading vector within each unit or period and then applies the
chosen dense or sparse cellwise or energy budget to those grouped loadings.
This lets the same formal bound machinery represent coherent unit-level or
period-level perturbations instead of purely cell-level shocks.

With \code{perturbation_layer = "stage2"}, the perturbed sample is the stored
stage-2 response or pseudo-outcome layer. With
\code{perturbation_layer = "baseline"}, the perturbed sample is the stored
baseline-counterfactual layer, which then propagates through the same stored
linear map. In the DR path this changes the operator because baseline errors
are first multiplied by the holdout weight used in the pseudo outcome.

With \code{perturbation_layer = "optimization"}, the perturbed sample is the
source optimization outcome sample used by the parent \code{cdmc_fit()} or
\code{cdmc_dr_fit()} object. The helper then perturbs each optimization cell by
\code{refit_step}, rebuilds the source fit with fixed tuning, re-evaluates the
requested target object, and uses the resulting one-step numerical Jacobian as
the loading vector \eqn{h}. In the DR path, the refit reuses the stored fold
assignments so the same-panel local Jacobian is deterministic. This mode is
therefore local rather than exact, but it propagates perturbations through
re-solving the staged low-rank baseline, the joint penalized objective, or the
cross-fitted DR pipeline.

This is exact for the stored linear post-fit maps used by staged and joint
\code{cdmc_fit()} objects, \code{cdmc_dr_fit()}, \code{cdmc_dose_response()},
and compatible \code{cdmc_dynamic_estimand()} objects built from those linear
layers. The optimization layer extends that contract for \code{cdmc_fit()}-
and \code{cdmc_dr_fit()}-backed targets through a local refit approximation.
The remaining limitations are broader identification-robust sensitivity models.

The output includes both the reference estimate table and a grid-expanded bound
table. The reference table reports the effective sensitivity multiplier under
the chosen perturbation constraint and scope for each statistic together with
the breakdown value at which the bounded-bias interval first covers the
requested null.
}
\value{
An object of class \code{cdmc_sensitivity_bounds} containing:

\item{reference}{A data frame with one row per reported statistic, including
the point estimate, null value, perturbation layer, perturbation constraint,
perturbation scope, basis label, effective sensitivity multiplier, the
resolved scale information, and the breakdown bias bound or breakdown
\code{gamma}.}

\item{bounds}{A data frame with one row per statistic, contamination fraction,
and \code{gamma} value, including the implied number of contaminated cells, the
effective multiplier under that sparse or dense perturbation budget, the raw
perturbation bound, interval radius, lower and upper bounds, and an indicator
for whether the bounded interval excludes the requested null.}

\item{summary(object)}{Returns a compact summary containing the reference table
and the bound row associated with the largest requested \code{gamma} and
contamination fraction for each statistic.}
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

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lambda = 0.2,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

bounds <- cdmc_sensitivity_bounds(
  fit,
  statistics = c("coefficients", "average_tau"),
  gamma_grid = c(0, 0.25, 0.5, 1),
  scale = "residual_sd"
)

summary(bounds)

summary(cdmc_sensitivity_bounds(
  fit,
  statistics = "coefficients",
  perturbation_scope = "unit",
  scale = "manual",
  scale_value = 1
))

joint_fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  objective = "joint",
  lambda = 0.2,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

summary(cdmc_sensitivity_bounds(joint_fit, scale = "manual", scale_value = 1))

summary(cdmc_sensitivity_bounds(
  joint_fit,
  perturbation_layer = "optimization",
  scale = "manual",
  scale_value = 1,
  refit_step = 1e-3
))

dose_response <- cdmc_dose_response(
  fit,
  model = "linear",
  lag_order = 1,
  include_zero_dose = TRUE
)

prediction_bounds <- cdmc_sensitivity_bounds(
  dose_response,
  statistics = "prediction",
  prediction_history = data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1)
  ),
  gamma_grid = c(0, 0.5, 1),
  scale = "manual",
  scale_value = 1
)

summary(prediction_bounds)

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

summary(cdmc_sensitivity_bounds(
  dr_fit,
  perturbation_layer = "baseline",
  scale = "manual",
  scale_value = 1
))

summary(cdmc_sensitivity_bounds(
  dr_fit,
  perturbation_layer = "optimization",
  scale = "manual",
  scale_value = 1,
  refit_step = 1e-3
))
}
}
\seealso{
\code{\link{cdmc_sensitivity_scan}}, \code{\link{cdmc_fit}},
\code{\link{cdmc_dr_fit}}, \code{\link{cdmc_dose_response}},
\code{\link{cdmc_dynamic_estimand}}, and \code{\link{cdmc_bootstrap}}.
}