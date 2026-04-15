\name{cdmc_joint_placebo_test}
\alias{cdmc_joint_placebo_test}
\alias{print.cdmc_joint_placebo_test}
\title{Aggregate placebo diagnostics across multiple pre-exposure windows}
\description{
Run placebo diagnostics separately for multiple negative periods relative to the
first treated observation and aggregate them into a joint pretrend diagnostic.
}
\usage{
cdmc_joint_placebo_test(
  object,
  periods = -3:-1,
  alpha = 0.05,
  equivalence_margin = NULL,
  verbose = FALSE
)

\method{print}{cdmc_joint_placebo_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{Strictly negative placebo periods relative to the first treated
  time for each unit.}

  \item{alpha}{Target significance level.}

  \item{equivalence_margin}{Optional positive equivalence margin. When supplied,
  the joint diagnostic aggregates period-specific equivalence tests rather than
  period-specific null tests.}

  \item{verbose}{Logical flag passed through to the underlying refit diagnostics.}

  \item{x}{A \code{cdmc_joint_placebo_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
Without an equivalence margin, the function combines period-specific placebo
diagnostic p values using Fisher's method. This targets the question of whether
any pre-exposure placebo window shows evidence against a zero pseudo effect.

When \code{equivalence_margin} is supplied, the function instead applies
\code{cdmc_equivalence_test()} period by period and aggregates the resulting
equivalence p values using a Bonferroni correction.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the period-specific placebo means and the joint
aggregated p value.
}
\value{
An object of class \code{cdmc_joint_placebo_test} containing the per-period test
table, the aggregation method, the joint p value, and a logical \code{passed}
indicator.
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 12, n_times = 12, beta = 1, lag_beta = 0.1, seed = 927)
fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lambda = 0.2,
  rank_max = 3,
  washout = 1,
  lag_order = 1,
  seed = 927
)

cdmc_joint_placebo_test(fit, periods = -2:-1, equivalence_margin = 0.2)
}
}
\seealso{
\code{\link{cdmc_placebo_test}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_equivalence_test}}.
}