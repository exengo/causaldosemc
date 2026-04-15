\name{cdmc_scia_test}
\alias{cdmc_scia_test}
\alias{print.cdmc_scia_test}
\title{Run a SCIA-style restriction diagnostic for dynamic assignment}
\description{
Run a research-use restriction diagnostic for sequential conditional
independence by testing whether named blocks of lagged observed outcomes or
lagged residual treatment effects still predict the current dose after
conditioning on observed treatment history and optional covariate history.
}
\usage{
cdmc_scia_test(
  object,
  lags = 1L,
  outcome_proxy = c("tau", "y"),
  covariates = object$covariates,
  include_current_covariates = TRUE,
  include_covariate_lags = TRUE,
  include_unit_effects = TRUE,
  include_time_effects = TRUE,
  restriction_blocks = NULL,
  p_adjust_method = c("holm", "none"),
  alpha = 0.05
)

\method{print}{cdmc_scia_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{lags}{Positive integer giving the number of lag periods to screen.}

  \item{outcome_proxy}{Whether to use lagged residual treatment effects
  \code{"tau"} or lagged observed outcomes \code{"y"} as the dynamic outcome
  history proxy in the augmented assignment model.}

  \item{covariates}{Character vector of covariate column names to include in the
  assignment screen. Defaults to the baseline-model covariates in
  \code{object}.}

  \item{include_current_covariates}{Logical indicating whether current-period
  covariate values should be included in the restricted assignment model.}

  \item{include_covariate_lags}{Logical indicating whether lagged covariate
  values should be included in the restricted assignment model.}

  \item{include_unit_effects}{Logical indicating whether unit fixed effects
  should be included in the assignment screen.}

  \item{include_time_effects}{Logical indicating whether time fixed effects
  should be included in the assignment screen.}

  \item{restriction_blocks}{Optional named list of lag-index blocks defining
  which lagged outcome histories should be checked as separate sequential
  restrictions. Each element should contain one or more integers between
  \code{1} and \code{lags}. When left as \code{NULL}, the helper creates one
  restriction per lagged outcome term.}

  \item{p_adjust_method}{Multiplicity adjustment used when summarizing the
  block-level restriction p values. Allowed values are \code{"holm"} and
  \code{"none"}.}

  \item{alpha}{Significance level used when flagging the screen as passed or
  failed.}

  \item{x}{A \code{cdmc_scia_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The diagnostic fits two linear models for the current continuous dose. The
restricted model conditions on lagged dose history, optional covariates, and
optional unit and time effects. The augmented model adds the requested lagged
outcome proxy terms. The reported F statistic and p value summarize whether the
lagged outcome terms add incremental predictive power.

In addition to that global comparison, the helper also evaluates a set of
named restriction blocks. Each block compares the restricted assignment model
against a version augmented only with the lagged outcome terms assigned to that
block. The returned object stores a block-level table with raw and adjusted
p values, making it easier to report which parts of the lagged outcome history
still appear to inform treatment assignment.

This is a screening diagnostic, not a formal proof of sequential conditional
independence and not a replacement for substantive DAG reasoning. A small p
value should be interpreted as evidence that dynamic outcome history may still
be informing treatment assignment beyond the observed conditioning set.

The resulting object can be supplied directly to \code{cdmc_bootstrap()} to
obtain percentile intervals for the global incremental F statistic and p value,
and optionally for the block-level restriction p values.
}
\value{
An object of class \code{cdmc_scia_test} containing the screened lag order, the
assignment-model formulas, fitted restricted and augmented regressions, the
global incremental F statistic and p value, a compact coefficient table for the
lagged outcome terms, an explicit conditioning-set summary, a block-level
restriction table with adjusted p values, and the sample used in the screen.
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

scia <- cdmc_scia_test(fit, lags = 1, outcome_proxy = "tau")
print(scia)

scia_blocks <- cdmc_scia_test(
  fit,
  lags = 2,
  outcome_proxy = "tau",
  restriction_blocks = list(recent_history = 1, deeper_history = 2)
)

print(scia_blocks)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_carryover_test}}.
}