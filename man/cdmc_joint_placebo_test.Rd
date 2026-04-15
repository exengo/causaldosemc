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
  placebo_windows = NULL,
  alpha = 0.05,
  equivalence_margin = NULL,
  rerun_tuning = FALSE,
  verbose = FALSE
)

\method{print}{cdmc_joint_placebo_test}(x, ...)
}
\arguments{
  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{periods}{Strictly negative placebo periods relative to the first treated
  time for each unit. This convenience argument is used when
  \code{placebo_windows = NULL}, in which case the helper creates one placebo
  window per requested period.}

  \item{placebo_windows}{Optional named list of placebo windows. Each element
  should contain one or more strictly negative integers relative to the first
  treated time for each unit. When supplied, the helper runs one placebo refit
  per named window and aggregates those window-specific diagnostics.}

  \item{alpha}{Target significance level.}

  \item{equivalence_margin}{Optional positive equivalence margin. When supplied,
  the joint diagnostic aggregates period-specific equivalence tests rather than
  period-specific null tests.}

  \item{rerun_tuning}{Logical indicating whether each period-specific placebo
  refit should rerun the source fit's automatic lambda-selection rule on its
  trimmed support mask. This requires the source \code{cdmc_fit} object to have
  been created with \code{lambda = NULL}; fixed-lambda source fits must use the
  default \code{FALSE}.}

  \item{verbose}{Logical flag passed through to the underlying refit diagnostics.}

  \item{x}{A \code{cdmc_joint_placebo_test} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
Without an equivalence margin, the function combines window-specific placebo
diagnostic p values using Fisher's method. This targets the question of whether
any pre-exposure placebo window shows evidence against a zero pseudo effect.

When \code{equivalence_margin} is supplied, the function instead applies
\code{cdmc_equivalence_test()} window by window and aggregates the resulting
equivalence p values using a Bonferroni correction.

By default, the component placebo refits reuse the original selected lambda.
When \code{rerun_tuning = TRUE}, each component placebo refit reruns the stored
automatic heuristic or blocked cross-validation lambda-selection rule on its
trimmed support mask before rebuilding the baseline.

The resulting object stores a clearer repeated-placebo protocol table with one
row per placebo window, including the window name, the pooled placebo periods,
the window-specific mean pseudo effect, and the window-specific test result.
It can be supplied directly to \code{cdmc_bootstrap()} to obtain percentile
intervals for the window-specific placebo means and the joint aggregated p
value.
}
\value{
An object of class \code{cdmc_joint_placebo_test} containing the per-window
test table, the aggregation method, the joint p value, and a logical
\code{passed} indicator.
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
  lambda = NULL,
  rank_max = 3,
  washout = 1,
  lag_order = 1,
  seed = 927
)

cdmc_joint_placebo_test(
  fit,
  placebo_windows = list(early_window = -3:-2, immediate_window = -1),
  equivalence_margin = 0.2,
  rerun_tuning = TRUE
)
}
}
\seealso{
\code{\link{cdmc_placebo_test}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_equivalence_test}}.
}