\name{cdmc_entropy_balance_weights}
\alias{cdmc_entropy_balance_weights}
\alias{print.cdmc_entropy_balance_weights}
\title{Estimate continuous-treatment entropy-balance weights}
\description{
Estimate exponential-tilting balancing weights for a continuous treatment and
return them in a form that can be supplied directly to \code{cdmc_fit()},
\code{cdmc_dr_fit()}, or \code{cdmc_dose_response()}.
}
\usage{
cdmc_entropy_balance_weights(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  degree = 1L,
  standardize = TRUE,
  iterations = 1000L,
  reltol = 1e-8,
  sample_weights = NULL,
  max_weight = "adaptive"
)

\method{print}{cdmc_entropy_balance_weights}(x, ...)
}
\arguments{
  \item{data}{A data frame containing the continuous treatment and balancing covariates.}

  \item{dose}{Name of the treatment column.}

  \item{covariates}{Character vector of covariate column names to balance, or
  \code{NULL} when balancing should rely on \code{time_effects} only.}

  \item{time}{Optional time-index column used when \code{time_effects = TRUE}.}

  \item{time_effects}{Logical indicating whether time fixed effects should be
  included in the marginal balance basis.}

  \item{model}{Formula expansion used for the balancing covariates.
  \code{"linear"} uses raw covariate terms, while \code{"spline"} expands
  selected numeric covariates with natural splines.}

  \item{df}{Positive integer controlling the spline basis size when
  \code{model = "spline"}.}

  \item{spline_covariates}{Character vector of covariate column names to
  expand with natural splines when \code{model = "spline"}. Defaults to
  \code{covariates}.}

  \item{degree}{Positive integer controlling the highest power of the centered
  dose included in the balance moments. The default \code{1L} balances
  first-order dose-feature moments.}

  \item{standardize}{Logical indicating whether the balance constraints should
  be scaled to unit weighted RMS before the dual optimization.}

  \item{iterations}{Maximum number of damped Newton iterations used by the
  entropy-balance solver.}

  \item{reltol}{Positive convergence target used by the entropy-balance solver.}

  \item{sample_weights}{Optional sampling weights, supplied either as a column
  name or a numeric vector of length \code{nrow(data)}.}

  \item{max_weight}{Optional clipping rule applied after the entropy-balance
  optimization. The default, \code{"adaptive"}, caps the returned weights at
  the smaller of the sample 99th percentile and \eqn{\sqrt{n}}. Supply
  \code{NULL} to disable clipping or a single positive numeric value to impose
  a fixed cap.}

  \item{x}{A \code{cdmc_entropy_balance_weights} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
This helper fits a continuous-treatment exponential-tilting problem on the
observed rows of \code{data}. The optimization preserves the weighted means of
its linear or spline balance basis and tilts the weights to reduce centered
moment conditions between the treatment and the selected covariate basis.

When \code{time_effects = TRUE}, time fixed effects enter the marginal balance
basis. The dose-balance moments focus on the covariate basis itself, or on the
full main basis when no covariates are supplied.

When the input includes padded panel rows flagged by \code{.cdmc_observed}, the
helper estimates the balancing weights on observed rows only and returns zero
weights on the padded rows so the result can be reused directly inside the
panel-fitting helpers.

These weights can be supplied directly to \code{cdmc_fit(weights = ...)},
\code{cdmc_dr_fit(weights = ...)} for the cross-fitted augmented linear effect
estimator, used internally through
\code{cdmc_dr_fit(weight_method = "entropy_balance")}, or passed to
\code{cdmc_dose_response(weights = ...)}.
}
\value{
An object of class \code{cdmc_entropy_balance_weights} with the estimated
weight vector, the balancing formula, optimization diagnostics, and a table of
pre- and post-weighted balance moments.
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
w <- cdmc_entropy_balance_weights(
  panel,
  dose = "dose",
  covariates = "x1",
  time = "time",
  time_effects = TRUE
)

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weights = w,
  lambda = 0.2,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

print(fit)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}},
\code{\link{cdmc_dose_response}}, \code{\link{cdmc_cbps_weights}}, and
\code{\link{cdmc_kernel_balance_weights}}.
}
