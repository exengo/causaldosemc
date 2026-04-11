\name{cdmc_carryover_test}
\alias{cdmc_carryover_test}
\alias{print.cdmc_carryover_test}
\title{Residual-based carryover diagnostic for post-exit zero-dose periods}
\description{
Screen for potential carryover by summarizing residual treatment effects on
zero-dose observations that occur shortly after a positive treatment period.
}
\usage{
cdmc_carryover_test(object, periods = 1L)

\method{print}{cdmc_carryover_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{One or more positive integers identifying the post-exit zero-dose
  periods to include in the screen. For example, \code{periods = 1:2} examines the
  first and second zero-dose periods after exposure ends.}

  \item{x}{A \code{cdmc_carryover_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The function computes the exit distance for each unit-time observation, defined
as the number of periods since the most recent nonzero dose, conditional on the
current dose being 0. It then extracts the residual treatment effects from the
fitted object,
\eqn{\hat{\tau}_{it} = Y_{it} - \hat{Y}_{it}(0)},
for the requested post-exit periods.

The output is best interpreted as a residual-based screening diagnostic. It does
not refit the baseline model with masked exit periods and should therefore not be
treated as a full inferential carryover procedure.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the mean residual carryover effect.
}
\value{
An object of class \code{cdmc_carryover_test} containing:

\item{periods}{The requested post-exit periods.}

\item{n}{The number of matched zero-dose exit cells.}

\item{mean_tau, sd_tau, standard_error, t_statistic, p_value}{Summary statistics
for the residual treatment effects on the matched cells.}

\item{cells}{A data frame with the matched unit-time observations, exit distance,
residual effect, and control-eligibility flag.}
}
\examples{
\donttest{
panel <- simulate_cdmc_data(
  n_units = 12,
  n_times = 10,
  beta = 1,
  lag_beta = 0.5,
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

cdmc_carryover_test(fit, periods = 1)
}
}
\seealso{
\code{\link{cdmc_fit}} and \code{\link{cdmc_bootstrap}} for model fitting and
bootstrap inference.
}