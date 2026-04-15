\name{cdmc_dynamic_estimand}
\alias{cdmc_dynamic_estimand}
\alias{print.cdmc_dynamic_estimand}
\title{Evaluate dynamic or pathwise causal estimands}
\description{
Construct a bootstrapable dynamic estimand object from a fitted
\code{cdmc_fit}, \code{cdmc_dr_fit}, or \code{cdmc_dose_response} object by
evaluating user-supplied dose histories.
}
\usage{
cdmc_dynamic_estimand(
  object,
  history,
  reference_history = NULL,
  type = c("response", "contrast", "slope"),
  aggregate = c("individual", "mean", "both"),
  path_weights = NULL,
  eps = 1e-4,
  labels = NULL
)

\method{print}{cdmc_dynamic_estimand}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit}, \code{cdmc_dr_fit}, or
  \code{cdmc_dose_response} object. For \code{cdmc_fit} objects, the stored
  built-in second-stage effect model must be \code{"linear"} or
  \code{"spline"}; \code{effect_model = "none"} is not supported.}

  \item{history}{The requested dose history or histories. Supply either a
  numeric scalar, a named numeric vector, or a matrix or data frame with the
  lag columns \code{dose_lag0}, \code{dose_lag1}, and so on up to the stored
  lag order. A scalar is interpreted as a contemporaneous dose with all lagged
  doses set to 0.}

  \item{reference_history}{Optional reference dose history used when
  \code{type = "contrast"}. When omitted, the reference defaults to the all-zero
  history. If either \code{history} or \code{reference_history} has one row and
  the other has multiple rows, the single-row input is recycled across the
  requested paths.}

  \item{type}{Dynamic estimand type. \code{"response"} returns the fitted effect
  of each supplied history, \code{"contrast"} returns the fitted effect relative
  to \code{reference_history}, and \code{"slope"} returns a finite-difference or
  analytic derivative with respect to \code{dose_lag0}.}

  \item{aggregate}{How the object should summarize the requested paths.
  \code{"individual"} stores one estimate per supplied history,
  \code{"mean"} stores only the weighted mean across those histories, and
  \code{"both"} stores both the individual path estimates and the weighted
  mean.}

  \item{path_weights}{Optional positive weights used when
  \code{aggregate = "mean"} or \code{aggregate = "both"}. Supply either a scalar
  or a vector matching the number of requested paths.}

  \item{eps}{Positive finite-difference step size used when
  \code{type = "slope"} and the source model does not admit an analytic slope
  shortcut.}

  \item{labels}{Optional character vector of labels for the requested paths. If
  omitted, labels are taken from the history row names when those row names are
  non-default; otherwise the object uses \code{path1}, \code{path2}, and so on.}

  \item{x}{A \code{cdmc_dynamic_estimand} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The helper formalizes dynamic or pathwise effect summaries without introducing a
new estimation layer.

For \code{cdmc_fit} objects, the function reuses the stored built-in second-stage
linear or spline effect model.

For \code{cdmc_dr_fit} objects, the function evaluates the final linear DR lag
coefficients on the requested histories.

For \code{cdmc_dose_response} objects, the function reuses the stored post-fit
dose-response model and prediction logic.

The returned object stores the requested histories, the evaluated path-level
estimates, and a named numeric vector of summary statistics. That vector is the
bootstrap target used by \code{cdmc_bootstrap()}, so the same dynamic estimand
definition can be resampled directly without manual reconstruction. When the
stored summary vector contains multiple pathwise statistics,
\code{cdmc_bootstrap()} now reports max-t simultaneous bands by default.
}
\value{
An object of class \code{cdmc_dynamic_estimand} containing:

\item{estimate_table}{A data frame of path-level estimates, including the
requested labels, the fitted path estimate, and path weights. Contrast objects
also include the fitted treated and reference effects.}

\item{estimates}{A named numeric vector containing the requested bootstrapable
summary statistics, such as \code{dynamic_contrast_path1} or
\code{mean_dynamic_slope}.}

\item{history, reference_history}{The resolved history inputs used in the
evaluation.}

\item{source_object, source_class}{The fitted source model used to generate the
dynamic estimand.}
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

path_effect <- cdmc_dynamic_estimand(
  fit,
  history = data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1)
  ),
  type = "contrast",
  aggregate = "both"
)

print(path_effect)
cdmc_bootstrap(path_effect, n_boot = 10, simultaneous = TRUE, seed = 1)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_dr_fit}},
\code{\link{cdmc_dose_response}}, and \code{\link{cdmc_bootstrap}}.
}