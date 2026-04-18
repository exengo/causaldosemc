\name{cdmc_fit}
\alias{cdmc_fit}
\alias{print.cdmc_fit}
\alias{summary.cdmc_fit}
\alias{print.summary.cdmc_fit}
\alias{predict.cdmc_fit}
\title{Fit a matrix-completion estimator for continuous treatment panels}
\description{
Estimate baseline untreated outcomes and residual treatment effects in a
panel with a continuous dose and 0 as the control state.
}
\usage{
cdmc_fit(
  data,
  outcome,
  dose,
  unit,
  time,
  covariates = NULL,
  weights = NULL,
  lambda = NULL,
  rank_max = NULL,
  lambda_fraction = 0.25,
  lambda_selection = c("heuristic", "cv"),
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  cv_workers = 1L,
  washout = 0L,
  lag_order = 0L,
  effect_model = c("linear", "spline", "none"),
  effect_df = 4L,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  zero_tolerance = 1e-8,
  seed = NULL,
  verbose = FALSE,
  objective = c("staged", "joint")
)

\method{print}{cdmc_fit}(x, ...)

\method{summary}{cdmc_fit}(object, ...)

\method{print}{summary.cdmc_fit}(x, ...)

\method{predict}{cdmc_fit}(object, newdata = NULL, type = c("y0", "tau", "tau_model", "tau_linear"), ...)
}
\arguments{
  \item{data}{A data frame containing one row per observed unit-time
  combination. Missing unit-time combinations are allowed and are padded
  internally onto the full unit-time grid.}

  \item{outcome}{Name of the numeric outcome column.}

  \item{dose}{Name of the numeric dose column. Active doses may be positive or
  negative; 0 defines the control state.}

  \item{unit}{Name of the unit identifier column.}

  \item{time}{Name of the time identifier column. Each unit must share the same
  ordered time support.}

  \item{covariates}{Optional character vector of numeric covariate columns to
  include in the stage-1 baseline model.}

  \item{weights}{Optional observation weights used in both the stage-1 baseline
  fit and the built-in stage-2 effect fit. This can be a column name in
  \code{data}, a numeric vector of length \code{nrow(data)}, a numeric matrix
  matching the panel dimensions, or an object with a numeric \code{weights}
  component such as the output of \code{cdmc_cbps_weights()} or
  \code{cdmc_entropy_balance_weights()} or
  \code{cdmc_kernel_balance_weights()} or
  \code{cdmc_adaptive_balance_weights()}. The internal fit normalizes supplied
  weights to have mean 1 over the active optimization sample: eligible zero-
  dose control cells under \code{objective = "staged"} and observed estimation
  cells under \code{objective = "joint"}.}

  \item{lambda}{Nonnegative nuclear-norm penalty passed to
  the low-rank solver. If \code{NULL}, a default value is computed from the
  residualized optimization sample after partialling out additive nuisance
  terms.}

  \item{rank_max}{Maximum rank supplied to the low-rank solver. If
  \code{NULL}, a conservative default is chosen from the panel dimensions.}

  \item{lambda_fraction}{Scaling factor used when \code{lambda = NULL}. Must lie
  in the interval \code{(0, 1]}.}

  \item{lambda_selection}{How \code{lambda} should be chosen when
  \code{lambda = NULL}. \code{"heuristic"} uses a fraction of
  a reference singular-value threshold after partialling out additive nuisance
  terms. \code{"cv"} uses blocked cross-validation on the active optimization
  sample and is the recommended path for empirical workflows when the panel has
  enough eligible zero-dose support to hold out contiguous control blocks
  safely.}

  \item{lambda_grid}{Optional numeric grid of candidate lambda values used when
  \code{lambda_selection = "cv"}. If \code{NULL}, a decreasing logarithmic grid is
  constructed automatically from the same residual reference threshold used by
  the heuristic rule.}

  \item{nlambda}{Number of automatically generated candidate lambda values used
  when \code{lambda_selection = "cv"} and \code{lambda_grid = NULL}.}

  \item{lambda_min_ratio}{Smallest automatic lambda as a fraction of the largest
  candidate lambda when \code{lambda_selection = "cv"}.}

  \item{cv_rounds}{Number of blocked cross-validation rounds used when
  \code{lambda_selection = "cv"}.}

  \item{cv_block_size}{Requested within-unit contiguous holdout block size used
  when \code{lambda_selection = "cv"}. The implementation may fall back to smaller
  blocks if the requested size cannot produce a valid holdout set.}

  \item{cv_workers}{Number of workers used when evaluating candidate lambda
  values inside each blocked cross-validation round. Set
  \code{cv_workers = 1} for sequential evaluation (default). Values greater
  than 1 enable parallel lambda scoring via multicore execution on supported
  platforms.}

  \item{washout}{Nonnegative integer giving the number of post-treatment periods
  to remove from the eligible zero-dose control set after each nonzero dose.
  This is a simple guard against carryover contamination in stage 1.}

  \item{lag_order}{Nonnegative integer giving the maximum lag order in the
  built-in second-stage dose-response model.}

  \item{effect_model}{Second-stage effect specification. Version 0.1 supports
  \code{"linear"}, \code{"spline"}, and \code{"none"}.}

  \item{effect_df}{Positive integer controlling the spline basis size when
  \code{effect_model = "spline"}. Ignored otherwise.}

  \item{objective}{Estimation objective. \code{"staged"} first fits the
  untreated baseline on eligible zero-dose cells and then estimates the effect
  model on the residual treatment component. \code{"joint"} fits the low-rank
  surface, additive nuisance terms, and the chosen linear or spline dose basis
  jointly over the observed estimation sample.}

  \item{outer_maxit}{Maximum number of outer alternations between nuisance-term
  estimation and low-rank matrix completion.}

  \item{fe_maxit}{Maximum number of iterations for the additive nuisance-term
  solver.}

  \item{soft_maxit}{Maximum number of iterations for
  the low-rank inner solver. This is passed to \code{softImpute} in the
  unweighted case and to the weighted proximal SVT routine in the weighted case.}

  \item{tol}{Convergence tolerance for the outer low-rank loop and the
  low-rank inner solver.}

  \item{fe_tol}{Convergence tolerance for the additive nuisance-term solver.}

  \item{zero_tolerance}{Absolute tolerance used to coerce very small doses to 0
  before building the eligible-control mask.}

  \item{seed}{Optional random seed used before fitting.}

  \item{verbose}{Logical flag indicating whether the outer low-rank iterations
  should print progress messages.}

  \item{x}{A fitted \code{cdmc_fit} object or a \code{summary.cdmc_fit} object.}

  \item{object}{A fitted \code{cdmc_fit} object.}

  \item{newdata}{Not used in version 0.1. Prediction is supported only on the
  training panel.}

  \item{type}{Prediction target. \code{"y0"} returns the imputed baseline,
  \code{"tau"} returns the residual treatment effect matrix, and
  \code{"tau_model"} returns the fitted second-stage component. The legacy
  \code{"tau_linear"} alias remains available for backward compatibility.}

  \item{...}{Unused additional arguments for the S3 methods.}
}
\details{
The package supports two related objectives.

Under the default staged objective, the package validates the observed panel
rows, orders the data by unit and time, pads any missing unit-time combinations
onto the full unit-time grid, and builds an eligible-control mask from observed
zero-dose cells. When \code{washout > 0}, zero-dose periods immediately
following nonzero treatment are removed from the control set. The low-rank
baseline is then estimated from that eligible zero-dose sample by alternating
between additive nuisance-term updates and low-rank matrix completion.
Residual treatment effects,
\eqn{\hat{\tau}_{it} = Y_{it} - \hat{Y}_{it}(0)}, are then projected on the
contemporaneous and lagged doses using either an unpenalized linear model or a
spline basis without an intercept.

Under \code{objective = "joint"}, the package keeps the same low-rank and
additive nuisance structure but augments the optimization problem with the
chosen dose basis and fits all components jointly over the observed estimation
sample. This mirrors the manuscript's global regularized objective more closely
than the staged estimator. When lagged doses are requested in the joint path,
pre-sample lag values are filled with 0 when constructing the global treatment
basis.

If \code{lambda = NULL}, the penalty can be chosen either heuristically or by
blocked cross-validation on the objective-specific optimization sample.
For empirical workflows, prefer blocked cross-validation when support allows
and treat the heuristic rule as a fallback for thin-support or exploratory
settings.
Unweighted fits use \code{softImpute} with \code{type = "svd"}, while weighted
fits use a weighted proximal singular-value-thresholding update on the same
objective-specific loss.

This version is aligned with the original library goal only at the core-engine
level. It now supports both balanced and unbalanced panels, but still treats
missing unit-time combinations as structurally unobserved rows with complete
data required on the rows that are present.
}
\value{
An object of class \code{cdmc_fit}. Important components include:

\item{data}{The full unit-time grid reordered by unit and time, augmented with
\code{.cdmc_observed}, \code{.cdmc_eligible_control}, \code{.cdmc_y0_hat},
\code{.cdmc_tau}, \code{.cdmc_tau_model}, the legacy alias
\code{.cdmc_tau_linear}, and \code{.cdmc_weight} when weighted fitting is
used. Rows introduced by internal padding retain \code{NA} values in the
original modeled variables.}

\item{eligible_mask}{Logical matrix identifying the zero-dose observations used
for stage-1 estimation after washout removal.}

\item{optimization_mask}{Logical matrix identifying the cells used in the
active optimization objective. This equals \code{eligible_mask} under the
staged estimator and the observed weighted sample under the joint estimator.}

\item{weight_matrix}{Numeric matrix of normalized observation weights when
weighted fitting is used, or \code{NULL} otherwise.}

\item{baseline}{A list containing the low-rank component, additive nuisance-term
fit, final baseline estimate, penalty value, the solver name, effective rank,
and convergence diagnostics.}

\item{effect}{A list containing the second-stage coefficients, fitted values,
residual treatment effects, and the sample mask used in the chosen effect model.}

\item{lambda_tuning}{A list describing how \code{lambda} was chosen. For
cross-validation fits, this includes the candidate grid, round-by-round scores,
mean scores, and holdout-set sizes.}

\item{lambda, rank_max, washout, lag_order, effect_model, effect_df, objective, weights}{The main
tuning and modeling choices used in the fit.}
}
\examples{
\donttest{
panel <- simulate_cdmc_data(
  n_units = 10,
  n_times = 8,
  beta = 1,
  lag_beta = 0.3,
  seed = 1
)

fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  washout = 1,
  lag_order = 1,
  lambda = NULL,
  lambda_selection = "cv",
  cv_rounds = 3L,
  cv_block_size = 2L,
  rank_max = 3,
  seed = 1
)

print(fit)
summary(fit)
head(predict(fit, type = "y0"))

joint_fit <- cdmc_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  lag_order = 1,
  effect_model = "linear",
  objective = "joint",
  lambda = 0.2,
  rank_max = 3,
  seed = 1
)

print(joint_fit)
}
}
\seealso{
\code{\link{cdmc_carryover_test}} for the current residual carryover screen and
\code{\link{simulate_cdmc_data}} to generate a small example panel.
}

