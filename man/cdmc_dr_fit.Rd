\name{cdmc_dr_fit}
\alias{cdmc_dr_fit}
\alias{print.cdmc_dr_fit}
\title{Fit a cross-fitted doubly robust linear effect estimator}
\description{
Estimate a research-use cross-fitted doubly robust-style linear dose-response
model by combining out-of-fold baseline predictions with either externally
supplied observation weights or internally estimated generalized propensity
score or balancing weights.
}
\usage{
cdmc_dr_fit(
  data,
  outcome,
  dose,
  unit,
  time,
  covariates = NULL,
  weights = NULL,
  weight_method = NULL,
  weight_covariates = covariates,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = weight_covariates,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  stabilize_weights = TRUE,
  max_weight = NULL,
  n_folds = 2L,
  fold_assignments = NULL,
  lambda = NULL,
  rank_max = NULL,
  lambda_fraction = 0.25,
  lambda_selection = c("heuristic", "cv"),
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  washout = 0L,
  lag_order = 0L,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  zero_tolerance = 1e-8,
  seed = NULL,
  verbose = FALSE
)

\method{print}{cdmc_dr_fit}(x, ...)
}
\arguments{
  \item{data}{A balanced panel data frame containing one row per unit-time combination.}

  \item{outcome}{Name of the numeric outcome column.}

  \item{dose}{Name of the numeric continuous-dose column. Active doses may be
  positive or negative; 0 defines the control state.}

  \item{unit}{Name of the unit identifier column.}

  \item{time}{Name of the time identifier column.}

  \item{covariates}{Optional character vector of covariate column names used in
  the stage-1 baseline model.}

  \item{weights}{Optional externally supplied observation weights used in the
  augmentation step and in the weighted stage-1 fits inside each training fold.
  This can be a column name, a numeric vector, a matrix matching the panel
  dimensions, or an object with a numeric \code{weights} component such as the
  output of \code{cdmc_cbps_weights()}. Leave this as \code{NULL} to use the
  built-in internal GPS nuisance model.}

  \item{weight_method}{Optional weighting mode. When omitted, the function uses
  \code{"external"} if \code{weights} is supplied and \code{"gaussian_gps"}
  otherwise. Supported internal choices are \code{"gaussian_gps"},
  \code{"kernel_gps"}, and \code{"cbps"}.}

  \item{weight_covariates}{Character vector of column names, or \code{NULL},
  used in the built-in internal GPS treatment-assignment model. Defaults to
  \code{covariates}. When \code{NULL}, the internal weighting model falls back
  to time fixed effects only if \code{gps_time_effects = TRUE}.}

  \item{gps_time_effects}{Logical indicating whether the internal GPS
  model should include time fixed effects.}

  \item{gps_model}{Built-in internal GPS mean-model specification.
  \code{"linear"} uses raw covariate terms, \code{"spline"} expands selected
  numeric covariates with natural-spline terms, \code{"gam"} fits smooth
  terms through \code{mgcv::gam()}, \code{"tree"} fits a regression tree
  through \code{rpart::rpart()}, and \code{"forest"} fits a regression
  forest through \code{ranger::ranger()}. The \code{"gam"}, \code{"tree"},
  and \code{"forest"} options are available for the Gaussian and kernel GPS
  paths, but not for
  \code{weight_method = "cbps"}.}

  \item{gps_df}{Positive integer controlling the nonlinear basis size when
  \code{gps_model = "spline"} or \code{"gam"}. For GAM fits this is passed as
  the smooth basis size \code{k}, so it must be at least 3. It is ignored when
  \code{gps_model = "linear"}, \code{"tree"}, or \code{"forest"}.}

  \item{gps_spline_covariates}{Character vector of covariate column names to
  expand with natural splines or GAM smooth terms when
  \code{gps_model = "spline"} or \code{"gam"}. Defaults to
  \code{weight_covariates}. Columns that are not numeric or that do not have
  enough unique values fall back to linear terms automatically. It is ignored
  when \code{gps_model = "tree"} or \code{"forest"}.}

  \item{gps_forest_trees}{Positive integer giving the number of trees used when
  \code{gps_model = "forest"}. It is ignored otherwise.}

  \item{gps_forest_mtry}{Optional positive integer giving the number of
  candidate predictors evaluated at each split when
  \code{gps_model = "forest"}. When left as \code{NULL}, the underlying
  \code{ranger::ranger()} default is used.}

  \item{gps_forest_min_node_size}{Optional positive integer giving the minimum
  terminal node size when \code{gps_model = "forest"}. When left as
  \code{NULL}, the underlying \code{ranger::ranger()} default is used.}

  \item{gps_bandwidth}{Optional positive scalar bandwidth used when
  \code{weight_method = "kernel_gps"}. When left as \code{NULL}, the kernel
  density estimator chooses a rule-of-thumb bandwidth separately within each
  fold.}

  \item{cbps_standardize}{Logical flag passed to the internal CBPS weighting
  fit when \code{weight_method = "cbps"}. This uses the same standardization
  convention as \code{cdmc_cbps_weights()}.}

  \item{cbps_method}{Optimization mode passed to the internal CBPS weighting
  fit when \code{weight_method = "cbps"}.}

  \item{cbps_iterations}{Positive integer giving the maximum optimizer
  iterations for the internal CBPS weighting fit.}

  \item{cbps_twostep}{Logical flag forwarded to the internal CBPS weighting
  fit when \code{weight_method = "cbps"}.}

  \item{stabilize_weights}{Logical indicating whether the internal GPS
  weights should use a marginal numerator density. This applies to the
  \code{"gaussian_gps"} and \code{"kernel_gps"} paths; \code{"cbps"} uses
  \code{cbps_standardize} instead.}

  \item{max_weight}{Optional positive scalar used to cap internally estimated
  GPS or CBPS weights.}

  \item{n_folds}{Number of unit folds used for cross-fitting.}

  \item{fold_assignments}{Optional data frame with columns \code{unit} and
  \code{fold} that fixes the unit-level cross-fitting split. When supplied,
  the rows must match the current unit set exactly and use consecutive fold
  labels \code{1, 2, ..., K}. This is mainly useful when deterministic refits
  are needed for replay or formal sensitivity routines.}

  \item{lambda, rank_max, lambda_fraction, lambda_selection, lambda_grid,
  nlambda, lambda_min_ratio, cv_rounds, cv_block_size, washout, lag_order,
  outer_maxit, fe_maxit, soft_maxit, tol, fe_tol, zero_tolerance, seed,
  verbose}{Passed through to the fold-specific \code{cdmc_fit()} calls. The
  built-in second stage is always linear in the contemporaneous and lagged
  doses.}

  \item{x}{A \code{cdmc_dr_fit} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\details{
The function partitions units into folds, or reuses the supplied
\code{fold_assignments}, fits the weighted matrix-completion
baseline and linear second-stage model on the training folds, and then predicts
the untreated baseline for each holdout unit using that unit's own eligible
zero-dose observations together with the time-side low-rank basis learned on the
training units.

If \code{weights} is not supplied, each training fold first fits an internal
generalized propensity score nuisance model for the observed dose using
\code{weight_covariates} and, by default, time fixed effects. The built-in
GPS mean model can use either raw linear terms, natural-spline expansions, GAM
smooths of selected numeric covariates, a regression tree, or a regression
forest on the same raw covariate and time-effect inputs. When
\code{weight_method = "gaussian_gps"}, the
conditional and marginal assignment densities use Gaussian approximations.
When \code{weight_method = "kernel_gps"}, the same mean model is paired with
kernel-smoothed residual and marginal dose densities; \code{gps_bandwidth}
controls that smoother when it is supplied. When
\code{weight_method = "cbps"}, the function reuses the same linear or spline
covariate formula inside fold-local \code{CBPS::CBPS()} fits. Because that
path returns realized balancing weights rather than an out-of-sample density
predictor, the training and holdout samples are weighted by separate fold-local
CBPS fits.

For holdout treated histories, the function forms an augmented pseudo outcome of
the form
\deqn{\hat{\tau}^{DR} = \hat{\tau}^{model} + w (\hat{\tau}^{obs} - \hat{\tau}^{model})}
and regresses that pseudo outcome on the lagged dose history without an
intercept.

This is a research-use implementation of a cross-fitted doubly robust-style
linear estimator. The built-in internal weighting path is intentionally simple:
it currently offers Gaussian and kernel GPS nuisance variants with linear,
spline, GAM, tree, or forest mean models plus a fold-local CBPS balancing
variant built on linear or spline covariate formulas. External weights remain
preferable when the design calls for richer balancing or for assignment models
that need boosting, stacking, or other learner-based nuisance estimation
beyond that built-in menu.

When the fitted object is supplied to \code{cdmc_bootstrap()}, the bootstrap can
summarize not only coefficients and mean DR effects but also lag-specific
average fitted contributions and user-defined linear dose-path contrasts.
}
\value{
An object of class \code{cdmc_dr_fit}. Important components include:

\item{data}{The reordered panel augmented with out-of-fold baseline,
out-of-fold residual treatment effects, augmented pseudo outcomes, fold ids, and
the normalized weights used in each holdout fold.}

\item{baseline}{A list containing the out-of-fold baseline matrix and a compact
summary of the fold-specific nuisance fits.}

\item{effect}{A list containing the final linear coefficients, the out-of-fold
residual treatment matrix, the out-of-fold model-based treatment component, the
augmented pseudo outcome matrix, and the sample mask used in the DR stage.}

\item{fold_assignments}{A data frame mapping each unit to its cross-fitting fold.}
}
\examples{
\donttest{
panel <- simulate_cdmc_data(n_units = 20, n_times = 12, beta = 1, lag_beta = 0.3, seed = 1)

fit <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "gaussian_gps",
  n_folds = 2,
  lambda = 0.2,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

print(fit)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}}, and
\code{\link{cdmc_cbps_weights}}.
}