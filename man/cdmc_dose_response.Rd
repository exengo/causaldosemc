\name{cdmc_dose_response}
\alias{cdmc_dose_response}
\alias{print.cdmc_dose_response}
\alias{predict.cdmc_dose_response}
\title{Fit a post-baseline dose-response model}
\description{
Estimate a linear or natural-spline dose-response model from the residual
treatment surface obtained after baseline matrix completion.
}
\usage{
cdmc_dose_response(
  object,
  model = c("linear", "spline"),
  lag_order = object$lag_order,
  df = 4L,
  include_zero_dose = TRUE,
  weights = NULL
)

\method{print}{cdmc_dose_response}(x, ...)

\method{predict}{cdmc_dose_response}(object, dose = NULL, history = NULL, type = c("response", "slope"), eps = 1e-4, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object for \code{cdmc_dose_response()},
  or a fitted \code{cdmc_dose_response} object for the prediction method.}

  \item{model}{Dose-response model type. \code{"linear"} fits an additive linear
  model in the contemporaneous and lagged doses. \code{"spline"} fits additive
  natural-spline terms centered so that the predicted effect at all-zero dose
  history is 0.}

  \item{lag_order}{Lag order to include in the post-fit dose-response model.}

  \item{df}{Degrees of freedom for the natural-spline basis when
  \code{model = "spline"}.}

  \item{include_zero_dose}{Logical flag indicating whether all valid histories,
  including all-zero histories, should be included in the dose-response fit.
  Including zero-dose histories helps anchor the response near 0.}

  \item{weights}{Optional observation weights for the second-stage dose-response
  fit. This can be a column name in the fitted panel data, a numeric vector of
  length equal to the panel row count, a numeric matrix matching the panel
  dimensions, or an object with a numeric \code{weights} component such as the
  output of \code{cdmc_cbps_weights()}. When omitted after a weighted
  \code{cdmc_fit()} call, the stored fit weights are reused automatically.}

  \item{x}{A \code{cdmc_dose_response} object.}

  \item{dose}{For prediction, a vector of contemporaneous doses. When
  \code{lag_order > 0}, unspecified lagged doses are set to 0.}

  \item{history}{For prediction, an explicit data frame or matrix containing the
  required lag columns \code{dose_lag0}, \code{dose_lag1}, and so on.}

  \item{type}{Prediction target. \code{"response"} returns the fitted dose-response
  value. \code{"slope"} returns a finite-difference derivative with respect to
  \code{dose_lag0}.}

  \item{eps}{Finite-difference step size used when \code{type = "slope"}.}

  \item{...}{Unused additional arguments for the print and prediction methods.}
}
\details{
The function works on the residual treatment surface
\eqn{\hat{\tau}_{it} = Y_{it} - \hat{Y}_{it}(0)}
from a fitted \code{cdmc_fit} object.

The current implementation is post-fit rather than joint: it does not change the
stage-1 matrix-completion baseline. This makes it useful for flexible dose-response
exploration while keeping the core estimator stable.

If the parent \code{cdmc_fit()} object was estimated with observation weights and
\code{weights = NULL}, the stored fit weights are used by default.

The fitted \code{cdmc_dose_response} object stores its parent \code{cdmc_fit()}
object so it can be supplied directly to \code{cdmc_bootstrap()} for
coefficient and prediction inference.
}
\value{
An object of class \code{cdmc_dose_response} containing the fitted history sample,
design information, coefficients, fitted values, residuals, and the stored
parent fit object.

The prediction method returns a data frame with the requested dose history and the
predicted response or slope.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 10, beta = 1, lag_beta = 0.2, seed = 1)
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

dr <- cdmc_dose_response(fit, model = "spline", lag_order = 1, df = 3)
predict(dr, dose = c(0, 0.5, 1, 1.5))
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_cbps_weights}}.
}