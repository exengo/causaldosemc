\name{cdmc_kernel_balance_weights}
\alias{cdmc_kernel_balance_weights}
\alias{print.cdmc_kernel_balance_weights}
\title{Estimate continuous-treatment kernel-balance weights}
\description{
Estimate nonlinear continuous-treatment balancing weights through a reduced
kernel basis and return them in a form that can be supplied directly to
\code{cdmc_fit()}, \code{cdmc_dr_fit()}, or \code{cdmc_dose_response()}.
}
\usage{
cdmc_kernel_balance_weights(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  degree = 1L,
  n_centers = 25L,
  bandwidth = NULL,
  standardize = TRUE,
  iterations = 1000L,
  reltol = 1e-8,
  sample_weights = NULL,
  max_weight = "adaptive"
)

\method{print}{cdmc_kernel_balance_weights}(x, ...)
}
\arguments{
  \item{data}{A data frame containing the continuous treatment and balancing covariates.}

  \item{dose}{Name of the treatment column.}

  \item{covariates}{Character vector of covariate column names to balance, or
  \code{NULL} when balancing should rely on \code{time_effects} only.}

  \item{time}{Optional time-index column used when \code{time_effects = TRUE}.}

  \item{time_effects}{Logical indicating whether time fixed effects should be
  added to the linear or spline balance basis.}

  \item{model}{Formula expansion used for the balance basis.
  \code{"linear"} uses raw covariate terms, while \code{"spline"} expands
  selected numeric covariates with natural-spline terms.}

  \item{df}{Positive integer controlling the spline basis size when
  \code{model = "spline"}.}

  \item{spline_covariates}{Character vector of covariate column names to
  expand with natural splines when \code{model = "spline"}. Defaults to
  \code{covariates}.}

  \item{degree}{Positive integer controlling the highest centered dose power
  used in the kernel-balance moment constraints. The default \code{1L}
  balances first-order dose-feature moments.}

  \item{n_centers}{Positive integer upper bound on the number of kernel centers
  used in the reduced radial-basis balance map. When the requested basis is too
  aggressive for the sample, the helper may fall back to fewer centers and
  records the realized count in the returned object.}

  \item{bandwidth}{Optional positive kernel bandwidth. When left as
  \code{NULL}, the helper chooses a data-dependent bandwidth and may
  automatically widen it when narrower kernels fail to achieve acceptable
  approximate balance.}

  \item{standardize}{Logical indicating whether the kernel-balance constraint
  columns should be scaled to unit weighted RMS before the dual optimization.}

  \item{iterations}{Positive integer giving the maximum number of damped Newton
  iterations used by the kernel-balance solver.}

  \item{reltol}{Positive convergence target used by the kernel-balance solver.}

  \item{sample_weights}{Optional sampling weights, supplied either as a column
  name or a numeric vector of length \code{nrow(data)}.}

  \item{max_weight}{Optional clipping rule applied after the kernel-balance
  fit. The default, \code{"adaptive"}, caps the realized weights at the
  smaller of the 99th percentile and \eqn{\sqrt{n}}. Supply \code{NULL} to
  disable clipping or a single positive numeric value to impose a fixed cap.}

  \item{x}{A \code{cdmc_kernel_balance_weights} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
This helper fits a positive reweighting scheme by exponential tilting, but it
targets a nonlinear balance basis rather than only raw linear moments. The
balance basis contains the same centered linear or spline covariate and time
terms used by the entropy-balance helper together with a reduced radial-basis
feature map built from the observed covariates. The resulting weights therefore
act like an approximate kernel-balancing layer for the continuous-treatment
setting.

Because the kernel basis is intentionally reduced through a finite center set,
this helper treats balance as approximate rather than exact. When the initial
auto-selected kernel bandwidth or center count is too aggressive for the sample,
the helper automatically widens the bandwidth and, if necessary, simplifies the
center set before declaring failure.

When the input includes padded panel rows flagged by \code{.cdmc_observed}, the
helper estimates the kernel-balance fit on observed rows only and returns zero
weights on the padded rows so the result can be reused directly inside the
panel-fitting helpers.

These weights can be supplied directly to \code{cdmc_fit(weights = ...)} so
that both the main baseline fit and the built-in treatment-effect fit use the
same balancing weights. They can also be supplied to
\code{cdmc_dr_fit(weights = ...)} for the cross-fitted augmented linear effect
estimator, used internally through
\code{cdmc_dr_fit(weight_method = "kernel_balance")}, or passed to
\code{cdmc_dose_response(weights = ...)}.
}
\value{
An object of class \code{cdmc_kernel_balance_weights} with the estimated weight
vector, the balance formula, the realized kernel bandwidth, the realized center
count, and the optimization summary.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 10, beta = 1, lag_beta = 0.2, seed = 1)
w <- cdmc_kernel_balance_weights(
  panel,
  dose = "dose",
  covariates = "x1",
  time = "time",
  time_effects = TRUE,
  n_centers = 8L
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

dr <- cdmc_dose_response(fit, model = "linear", lag_order = 1)
print(dr)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}}, and
\code{\link{cdmc_dose_response}}.
}