\name{cdmc_adaptive_balance_weights}
\alias{cdmc_adaptive_balance_weights}
\alias{print.cdmc_adaptive_balance_weights}
\title{Estimate continuous-treatment adaptive balance weights}
\description{
Estimate continuous-treatment balancing weights by fitting several balancing
helpers on the same data and selecting the candidate with the best common
fold-local balance score.
}
\usage{
cdmc_adaptive_balance_weights(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  methods = NULL,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  entropy_balance_degree = 1L,
  entropy_balance_standardize = TRUE,
  entropy_balance_iterations = 1000L,
  entropy_balance_reltol = 1e-8,
  kernel_balance_degree = 1L,
  kernel_balance_centers = 25L,
  kernel_balance_bandwidth = NULL,
  kernel_balance_standardize = TRUE,
  kernel_balance_iterations = 1000L,
  kernel_balance_reltol = 1e-8,
  sample_weights = NULL,
  max_weight = "adaptive"
)

\method{print}{cdmc_adaptive_balance_weights}(x, ...)
}
\arguments{
  \item{data}{A data frame containing the continuous treatment and balancing covariates.}

  \item{dose}{Name of the treatment column.}

  \item{covariates}{Character vector of covariate column names to balance, or
  \code{NULL} when balancing should rely on \code{time_effects} only.}

  \item{time}{Optional time-index column used when \code{time_effects = TRUE}.}

  \item{time_effects}{Logical indicating whether time fixed effects should be
  added to the balance basis.}

  \item{model}{Formula expansion used by the candidate balancing helpers.
  \code{"linear"} uses raw covariate terms, while \code{"spline"} expands
  selected numeric covariates with natural-spline terms.}

  \item{df}{Positive integer controlling the spline basis size when
  \code{model = "spline"}.}

  \item{spline_covariates}{Character vector of covariate column names to
  expand with natural splines when \code{model = "spline"}. Defaults to
  \code{covariates}.}

  \item{methods}{Optional character vector choosing which balancing helpers to
  compare. Allowed values are \code{"cbps"}, \code{"entropy_balance"}, and
  \code{"kernel_balance"}. When left as \code{NULL}, the helper compares
  entropy-balance and kernel-balance by default and also includes \code{"cbps"}
  when \code{CBPS} is installed.}

  \item{cbps_standardize, cbps_method, cbps_iterations, cbps_twostep}{Tuning
  arguments forwarded to the \code{CBPS} candidate whenever
  \code{"cbps"} is included in \code{methods}.}

  \item{entropy_balance_degree, entropy_balance_standardize,
  entropy_balance_iterations, entropy_balance_reltol}{Tuning arguments
  forwarded to the entropy-balance candidate whenever
  \code{"entropy_balance"} is included in \code{methods}.}

  \item{kernel_balance_degree, kernel_balance_centers,
  kernel_balance_bandwidth, kernel_balance_standardize,
  kernel_balance_iterations, kernel_balance_reltol}{Tuning arguments forwarded
  to the reduced-kernel balancing candidate whenever
  \code{"kernel_balance"} is included in \code{methods}.}

  \item{sample_weights}{Optional sampling weights, supplied either as a column
  name or a numeric vector of length \code{nrow(data)}.}

  \item{max_weight}{Optional clipping rule applied to each candidate after
  fitting. The default, \code{"adaptive"}, caps the realized weights at the
  smaller of the 99th percentile and \eqn{\sqrt{n}}. Supply \code{NULL} to
  disable clipping or a single positive numeric value to impose a fixed cap.}

  \item{x}{A \code{cdmc_adaptive_balance_weights} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
This helper compares several balancing-oriented weighting paths on the same
sample and then selects one automatically. Each candidate is fit on the same
observed rows and then scored on a common balance basis. When the candidate set
includes \code{"kernel_balance"}, the common score basis uses the same reduced
radial-basis representation used by the kernel-balance helper; otherwise it
falls back to the centered linear or spline moment basis used by the
entropy-balance helper.

The selected candidate is the one with the smallest maximum absolute
standardized imbalance on that common score basis, with effective sample size as
a tie-breaker. The returned object stores the full candidate score table, the
selected method, and the selected weight vector.

When the input includes padded panel rows flagged by \code{.cdmc_observed}, the
helper estimates all candidate fits on observed rows only and returns zero
weights on the padded rows so the result can be reused directly inside the
panel-fitting helpers.

These weights can be supplied directly to \code{cdmc_fit(weights = ...)} so
that both the main baseline fit and the built-in treatment-effect fit use the
same selected balancing weights. They can also be supplied to
\code{cdmc_dr_fit(weights = ...)} for the cross-fitted augmented linear effect
estimator, used internally through
\code{cdmc_dr_fit(weight_method = "adaptive_balance")}, or passed to
\code{cdmc_dose_response(weights = ...)}.
}
\value{
An object of class \code{cdmc_adaptive_balance_weights} with the selected weight
vector, the selected candidate method, the candidate score table, and the
selected candidate fit.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 20, n_times = 10, beta = 1, lag_beta = 0.2, seed = 1)
w <- cdmc_adaptive_balance_weights(
  panel,
  dose = "dose",
  covariates = "x1",
  time = "time",
  time_effects = TRUE,
  methods = c("entropy_balance", "kernel_balance"),
  kernel_balance_centers = 8L
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
\code{\link{cdmc_dose_response}}, \code{\link{cdmc_entropy_balance_weights}},
\code{\link{cdmc_kernel_balance_weights}}, and
\code{\link{cdmc_cbps_weights}}.
}
