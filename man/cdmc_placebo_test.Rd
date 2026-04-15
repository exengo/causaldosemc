\name{cdmc_placebo_test}
\alias{cdmc_placebo_test}
\alias{print.cdmc_placebo_test}
\title{Refit-based placebo diagnostic for pre-exposure control periods}
\description{
Hold out selected zero-dose periods before the first positive exposure, refit the
main baseline objective, and summarize the resulting pseudo effects on the held-out
placebo cells.
}
\usage{
cdmc_placebo_test(object, periods = -2:0, rerun_tuning = FALSE, verbose = FALSE)

\method{print}{cdmc_placebo_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{One or more nonpositive integers indexing periods relative to
  the first positive exposure for each unit. For example, \code{-2:0} requests the
  two periods before first exposure and the immediately preceding period.}

  \item{rerun_tuning}{Logical indicating whether the diagnostic refit should
  rerun the source fit's automatic lambda-selection rule on the trimmed support
  mask. This requires the source \code{cdmc_fit} object to have been created
  with \code{lambda = NULL}; fixed-lambda source fits must use the default
  \code{FALSE}.}

  \item{verbose}{Logical flag controlling progress messages during the diagnostic
  refit.}

  \item{x}{A \code{cdmc_placebo_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The function identifies eligible zero-dose periods prior to the first positive
treatment for each unit, removes the requested placebo cells from the baseline
objective, and then computes pseudo effects on the held-out cells using the
refitted baseline surface.

By default, the refit reuses the original selected lambda. When
\code{rerun_tuning = TRUE}, the refit instead reruns the stored automatic
heuristic or blocked cross-validation lambda-selection rule on the trimmed
support mask before rebuilding the baseline.

The diagnostic is designed to check whether the baseline engine can reproduce
untreated outcomes in windows that should not contain treatment effects.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the mean placebo pseudo effect.
}
\value{
An object containing the requested placebo periods, pseudo-effect summary
statistics, the matched cells, and the refitted baseline object used to generate
the out-of-sample placebo predictions.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 10, beta = 1, lag_beta = 0.3, seed = 1)
fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lambda = NULL,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

cdmc_placebo_test(fit, periods = -2:0, rerun_tuning = TRUE)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_carryover_refit_test}}.
}
