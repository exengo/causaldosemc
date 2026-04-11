\name{cdmc_carryover_refit_test}
\alias{cdmc_carryover_refit_test}
\alias{print.cdmc_carryover_refit_test}
\title{Refit-based carryover diagnostic for post-exit zero-dose periods}
\description{
Hold out selected post-exit zero-dose periods, refit the main baseline objective,
and summarize the resulting pseudo effects on the held-out exit cells.
}
\usage{
cdmc_carryover_refit_test(object, periods = 1L, verbose = FALSE)

\method{print}{cdmc_carryover_refit_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{One or more positive integers identifying zero-dose periods after
  treatment exits. For example, \code{periods = 1:2} tests the first and second
  zero-dose periods after the most recent nonzero dose.}

  \item{verbose}{Logical flag controlling progress messages during the diagnostic
  refit.}

  \item{x}{A \code{cdmc_carryover_refit_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The function identifies eligible post-exit zero-dose observations by computing
the number of periods since the most recent nonzero dose. The requested cells are
removed from the baseline objective, and the refitted baseline surface is used
to compute pseudo effects on the held-out exit periods.

This is more aligned with the diagnostic logic discussed in the project roadmap
than the residual-only carryover screen, though it still reuses the original
lambda value instead of rerunning tuning inside the diagnostic refit.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the mean refit carryover pseudo effect.
}
\value{
An object containing the requested carryover periods, pseudo-effect summary
statistics, the matched exit cells, and the refitted baseline object used to
generate the out-of-sample predictions.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 10, beta = 1, lag_beta = 0.5, seed = 1)
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

cdmc_carryover_refit_test(fit, periods = 1)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}},
\code{\link{cdmc_carryover_test}}, and \code{\link{cdmc_placebo_test}}.
}