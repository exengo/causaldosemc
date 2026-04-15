\name{cdmc_equivalence_test}
\alias{cdmc_equivalence_test}
\alias{print.cdmc_equivalence_test}
\title{Test placebo or carryover equivalence}
\description{
Apply a two one-sided tests equivalence procedure to a diagnostic effect sample
from \code{cdmc_placebo_test()}, \code{cdmc_carryover_refit_test()}, or
\code{cdmc_carryover_test()}.
}
\usage{
cdmc_equivalence_test(object, margin, alpha = 0.05)

\method{print}{cdmc_equivalence_test}(x, ...)
}
\arguments{
  \item{object}{A fitted diagnostic object from \code{cdmc_placebo_test()},
  \code{cdmc_carryover_refit_test()}, or \code{cdmc_carryover_test()}.}

  \item{margin}{Positive equivalence margin. The null hypothesis is that the
  mean diagnostic effect lies outside \eqn{[-margin, margin]}.}

  \item{alpha}{Target significance level for the two one-sided tests procedure.}

  \item{x}{A \code{cdmc_equivalence_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
This helper formalizes the manuscript's recommendation to treat placebo and
carryover diagnostics as equivalence questions rather than relying only on a
failure to reject zero. It computes one lower-tail and one upper-tail t test
for the mean diagnostic effect and reports the larger of the two one-sided
p values.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile summaries for the mean diagnostic effect and the equivalence
p value.
}
\value{
An object of class \code{cdmc_equivalence_test} containing the diagnostic source,
sample size, mean diagnostic effect, equivalence margin, one-sided t statistics,
one-sided p values, and a logical \code{equivalent} flag.
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

placebo <- cdmc_placebo_test(fit, periods = -2:0)
cdmc_equivalence_test(placebo, margin = 0.2)
}
}
\seealso{
\code{\link{cdmc_placebo_test}}, \code{\link{cdmc_carryover_refit_test}}, and
\code{\link{cdmc_carryover_test}}, and \code{\link{cdmc_bootstrap}}.
}