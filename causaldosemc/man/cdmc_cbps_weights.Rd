\name{cdmc_cbps_weights}
\alias{cdmc_cbps_weights}
\alias{print.cdmc_cbps_weights}
\title{Estimate continuous-treatment CBPS weights}
\description{
Estimate covariate-balancing propensity weights for a continuous treatment using
\code{CBPS::CBPS()} and return them in a form that can be supplied directly to
\code{cdmc_fit()}, \code{cdmc_dr_fit()}, or \code{cdmc_dose_response()}.
}
\usage{
cdmc_cbps_weights(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  standardize = TRUE,
  method = c("over", "exact"),
  iterations = 1000L,
  twostep = TRUE,
  sample_weights = NULL,
  ...
)

\method{print}{cdmc_cbps_weights}(x, ...)
}
\arguments{
  \item{data}{A data frame containing the continuous treatment and balancing covariates.}

  \item{dose}{Name of the treatment column.}

  \item{covariates}{Character vector of covariate column names to balance, or
  \code{NULL} when balancing should rely on \code{time_effects} only.}

  \item{time}{Optional time-index column used when \code{time_effects = TRUE}.}

  \item{time_effects}{Logical indicating whether time fixed effects should be
  added to the CBPS balancing formula.}

  \item{model}{Formula expansion used for the balancing covariates.
  \code{"linear"} uses raw covariate terms, while \code{"spline"} expands
  selected numeric covariates with natural-spline terms.}

  \item{df}{Positive integer controlling the spline basis size when
  \code{model = "spline"}.}

  \item{spline_covariates}{Character vector of covariate column names to
  expand with natural splines when \code{model = "spline"}. Defaults to
  \code{covariates}.}

  \item{standardize}{Passed to \code{CBPS::CBPS()}. When \code{TRUE}, weights are
  normalized to sum to 1 over the sample for continuous treatments.}

  \item{method}{Optimization mode passed to \code{CBPS::CBPS()}.}

  \item{iterations}{Maximum optimizer iterations.}

  \item{twostep}{Logical flag forwarded to \code{CBPS::CBPS()}.}

  \item{sample_weights}{Optional sampling weights, supplied either as a column
  name or a numeric vector of length \code{nrow(data)}.}

  \item{...}{Additional arguments passed through to \code{CBPS::CBPS()}.}

  \item{x}{A \code{cdmc_cbps_weights} object.}
}
\details{
This helper is intentionally narrow: it targets the continuous-treatment case
that matches the current \code{causaldosemc} scope. The returned object contains
both the weight vector and the underlying \code{CBPS} fit object.

When the input includes padded panel rows flagged by \code{.cdmc_observed}, the
helper estimates the CBPS fit on observed rows only and returns zero weights on
the padded rows so the result can be reused directly inside the panel-fitting
helpers.

When the default over-identified continuous-treatment fit fails because the GMM
weighting matrix is singular, the helper automatically retries the just-identified
\code{method = "exact"} variant.

These weights can now be supplied directly to \code{cdmc_fit(weights = ...)} so
that both the main baseline fit and the built-in treatment-effect fit
use the same balancing weights. They can also be supplied to
\code{cdmc_dr_fit(weights = ...)} for the cross-fitted augmented linear effect
estimator, used internally through \code{cdmc_dr_fit(weight_method = "cbps")},
or passed to \code{cdmc_dose_response(weights = ...)}.
}
\value{
An object of class \code{cdmc_cbps_weights} with the estimated weight vector,
the matching formula, and the underlying \code{CBPS} fit object.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 10, beta = 1, lag_beta = 0.2, seed = 1)
w <- cdmc_cbps_weights(panel, dose = "dose", covariates = "x1")

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

dr <- cdmc_dose_response(fit, model = "linear", lag_order = 1)
print(dr)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}}, and
\code{\link{cdmc_dose_response}}.
}