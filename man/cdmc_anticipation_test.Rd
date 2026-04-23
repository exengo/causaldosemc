\name{cdmc_anticipation_test}
\alias{cdmc_anticipation_test}
\alias{print.cdmc_anticipation_test}
\title{Residual-based anticipation diagnostic for pre-entry zero-dose periods}
\description{
Screen for anticipation by summarizing residual treatment effects on zero-dose
observations that occur shortly before a positive treatment period. This is the
forward-looking analogue of \code{cdmc_carryover_test()}.
}
\usage{
cdmc_anticipation_test(object, periods = 1L)

\method{print}{cdmc_anticipation_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{One or more positive integers identifying the pre-entry
  zero-dose periods (lead lengths) to include in the screen. For example,
  \code{periods = 1:2} examines the first and second zero-dose periods before
  exposure begins.}

  \item{x}{A \code{cdmc_anticipation_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The function computes the entry distance for each unit-time observation,
defined as the number of periods until the next nonzero dose, conditional on
the current dose being 0. It then extracts the residual treatment effects
from the fitted object for the requested pre-entry periods.

A non-trivial mean residual on these "lead" cells is evidence of anticipation
or non-parallel trends that the matrix-completion baseline did not absorb. The
test reports both a naive standard error and a unit-clustered standard error
(CR1), defaulting to the cluster-robust value for inference.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the mean residual anticipation effect.
}
\value{
A list with class \code{cdmc_anticipation_test} containing the matched cells,
mean residual effect, naive and cluster-robust standard errors, t statistic,
p-value, and the underlying \code{cdmc_fit}.
}
\seealso{
\code{\link{cdmc_carryover_test}}, \code{\link{cdmc_placebo_test}}
}
