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
  gps_model = c("linear", "spline", "gam", "tree", "forest", "stack", "boost"),
  gps_df = 4L,
  gps_spline_covariates = weight_covariates,
  gps_stack_models = NULL,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  gps_boost_trees = 200L,
  gps_boost_depth = 2L,
  gps_boost_shrinkage = 0.05,
  gps_boost_min_obs_node = 10L,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  adaptive_balance_methods = NULL,
  entropy_balance_degree = 1L,
  entropy_balance_standardize = TRUE,
  entropy_balance_iterations = 1000L,
  entropy_balance_reltol = 1e-8,
  kernel_balance_degree = 1L,
  kernel_balance_centers = 25L,
  kernel_balance_bandwidth = NULL,
  kernel_balance_standardize = TRUE,
  kernel_balance_iterations = 1000L,
  kernel_balance_reltol = 1e-8,
  stabilize_weights = TRUE,
  max_weight = "adaptive",
  baseline_weighting = c("gps", "none"),
  dr_score = c("aipw", "plr"),
  n_folds = 2L,
  dr_workers = 1L,
  fold_assignments = NULL,
  lambda = NULL,
  rank_max = NULL,
  lambda_fraction = 0.25,
  lambda_selection = c("cv", "heuristic"),
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  cv_workers = 1L,
  cv_top_k = NULL,
  cv_coarse_to_fine = FALSE,
  cv_coarse_nlambda = NULL,
  cv_warm_starts = FALSE,
  lambda_control = NULL,
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
  output of \code{cdmc_cbps_weights()} or
  \code{cdmc_entropy_balance_weights()} or
  \code{cdmc_kernel_balance_weights()} or
  \code{cdmc_adaptive_balance_weights()}. Leave this as \code{NULL} to use the
  built-in internal weighting nuisance model.}

  \item{weight_method}{Optional weighting mode. When omitted, the function uses
  \code{"external"} if \code{weights} is supplied and \code{"gaussian_gps"}
  otherwise. Supported internal choices are \code{"gaussian_gps"},
  \code{"kernel_gps"}, \code{"cbps"}, \code{"entropy_balance"}, and
  \code{"kernel_balance"}, and \code{"adaptive_balance"}.}

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
  forest through \code{ranger::ranger()}. \code{"stack"} fits a convex
  combination of one or more of those base learners using in-fold prediction
  stacking. \code{"boost"} fits a gradient boosting regression model through
  \code{gbm::gbm()}. The \code{"gam"}, \code{"tree"}, \code{"forest"},
  \code{"stack"}, and \code{"boost"} options are available for the Gaussian
  and kernel GPS paths, but not for balancing-only paths such as
  \code{weight_method = "cbps"}, \code{"entropy_balance"}, or
  \code{"kernel_balance"}, or \code{"adaptive_balance"}.}

  \item{gps_df}{Positive integer controlling the nonlinear basis size when
  \code{gps_model = "spline"} or \code{"gam"}. For GAM fits this is passed as
  the smooth basis size \code{k}, so it must be at least 3. When
  \code{gps_model = "stack"}, the same value is used by any spline or GAM
  components in \code{gps_stack_models}. It is otherwise ignored when
  \code{gps_model = "linear"}, \code{"tree"}, or \code{"forest"}.}

  \item{gps_spline_covariates}{Character vector of covariate column names to
  expand with natural splines or GAM smooth terms when
  \code{gps_model = "spline"} or \code{"gam"}. Defaults to
  \code{weight_covariates}. Columns that are not numeric or that do not have
  enough unique values fall back to linear terms automatically. It is ignored
  when \code{gps_model = "tree"} or \code{"forest"}, and is used only by the
  spline or GAM components when \code{gps_model = "stack"}.}

  \item{gps_stack_models}{Optional character vector giving the base learner
  set used when \code{gps_model = "stack"}. Allowed entries are
  \code{"linear"}, \code{"spline"}, \code{"gam"}, \code{"tree"}, and
  \code{"forest"}, and \code{"boost"}. When left as \code{NULL}, the stack
  uses the available installed learners from that list. Ignored otherwise.}

  \item{gps_forest_trees}{Positive integer giving the number of trees used when
  \code{gps_model = "forest"}. When \code{gps_model = "stack"}, the same
  value is passed to any forest component included in \code{gps_stack_models}.
  It is ignored otherwise.}

  \item{gps_forest_mtry}{Optional positive integer giving the number of
  candidate predictors evaluated at each split when
  \code{gps_model = "forest"}. When \code{gps_model = "stack"}, the same
  value is passed to any forest component included in \code{gps_stack_models}.
  When left as \code{NULL}, the underlying \code{ranger::ranger()} default is
  used.}

  \item{gps_forest_min_node_size}{Optional positive integer giving the minimum
  terminal node size when \code{gps_model = "forest"}. When left as
  \code{NULL}, the underlying \code{ranger::ranger()} default is used. When
  \code{gps_model = "stack"}, the same value is passed to any forest
  component included in \code{gps_stack_models}.}

  \item{gps_boost_trees}{Positive integer giving the number of boosting
  iterations used when \code{gps_model = "boost"}. When
  \code{gps_model = "stack"}, the same value is passed to any boosting
  component included in \code{gps_stack_models}.}

  \item{gps_boost_depth}{Positive integer giving the interaction depth passed
  to \code{gbm::gbm()} when \code{gps_model = "boost"}. When
  \code{gps_model = "stack"}, the same value is passed to any boosting
  component included in \code{gps_stack_models}.}

  \item{gps_boost_shrinkage}{Positive scalar giving the learning-rate
  parameter passed to \code{gbm::gbm()} when \code{gps_model = "boost"}. When
  \code{gps_model = "stack"}, the same value is passed to any boosting
  component included in \code{gps_stack_models}.}

  \item{gps_boost_min_obs_node}{Positive integer giving the minimum terminal
  node size passed to \code{gbm::gbm()} when \code{gps_model = "boost"}. When
  \code{gps_model = "stack"}, the same value is passed to any boosting
  component included in \code{gps_stack_models}.}

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

  \item{adaptive_balance_methods}{Optional character vector choosing which
  built-in balancing helpers the adaptive selector should compare when
  \code{weight_method = "adaptive_balance"}. Allowed entries are
  \code{"cbps"}, \code{"entropy_balance"}, and \code{"kernel_balance"}.
  When left as \code{NULL}, the selector compares entropy-balance and
  kernel-balance by default and also includes \code{"cbps"} when
  \code{CBPS} is installed.}

  \item{entropy_balance_degree}{Positive integer controlling the highest power
  of the centered dose used in the entropy-balance moment constraints when
  \code{weight_method = "entropy_balance"}. The default \code{1L} balances
  first-order dose-feature moments.}

  \item{entropy_balance_standardize}{Logical indicating whether the entropy-
  balance constraint columns should be scaled to unit weighted RMS before the
  dual optimization when \code{weight_method = "entropy_balance"}.}

  \item{entropy_balance_iterations}{Positive integer giving the maximum number
  of damped Newton iterations used by the internal entropy-balance solver when
  \code{weight_method = "entropy_balance"}.}

  \item{entropy_balance_reltol}{Positive numeric convergence target used by
  the internal entropy-balance solver when
  \code{weight_method = "entropy_balance"}.}

  \item{kernel_balance_degree}{Positive integer controlling the highest power
  of the centered dose used in the kernel-balance moment constraints when
  \code{weight_method = "kernel_balance"}. The default \code{1L} balances
  first-order dose-feature moments.}

  \item{kernel_balance_centers}{Positive integer upper bound on the number of
  radial-basis kernel centers used when \code{weight_method = "kernel_balance"}.
  When the requested kernel basis is too aggressive for a fold, the internal
  helper may fall back to fewer centers and records the realized count in the
  fold summaries.}

  \item{kernel_balance_bandwidth}{Optional positive kernel bandwidth used when
  \code{weight_method = "kernel_balance"}. When left as \code{NULL}, the
  internal helper chooses a data-dependent bandwidth and may automatically
  widen it when narrower kernels fail to achieve acceptable approximate
  balance.}

  \item{kernel_balance_standardize}{Logical indicating whether the
  kernel-balance constraint columns should be scaled to unit weighted RMS
  before the dual optimization when
  \code{weight_method = "kernel_balance"}.}

  \item{kernel_balance_iterations}{Positive integer giving the maximum number
  of damped Newton iterations used by the internal kernel-balance solver when
  \code{weight_method = "kernel_balance"}.}

  \item{kernel_balance_reltol}{Positive numeric convergence target used by the
  internal kernel-balance solver when
  \code{weight_method = "kernel_balance"}.}

  \item{stabilize_weights}{Logical indicating whether the internal GPS
  weights should use a marginal numerator density. This applies to the
  \code{"gaussian_gps"} and \code{"kernel_gps"} paths; balancing-only paths
  such as \code{"cbps"}, \code{"entropy_balance"}, and
  \code{"kernel_balance"}, and \code{"adaptive_balance"} use their own
  direct weight-fitting criteria instead.}

  \item{max_weight}{Optional clipping rule used for internally estimated GPS,
  CBPS, entropy-balance, kernel-balance, or adaptive-balance weights. The default,
  \code{"adaptive"}, caps each fold's internal weights at the smaller of the
  99th percentile of the fold-local training weights and \eqn{\sqrt{n}} for
  that training sample. Supply \code{NULL} to disable clipping or a single
  positive numeric value to impose a fixed cap.}

  \item{baseline_weighting}{Controls whether the train-fold matrix-completion
  baseline is fit with GPS/balancing weights. \code{"gps"} (default, matches
  historical behaviour) weights the train baseline by the same per-cell
  weights used in the augmentation term, tilting the imputed counterfactual
  surface toward the treatment-relevant covariate distribution. \code{"none"}
  fits the baseline unweighted, which is the canonical AIPW configuration:
  asymptotic orthogonality only requires the augmentation term to carry the
  weights. The unweighted choice can reduce variance when the dose-response
  surface is well-modeled across the full covariate distribution; the
  weighted choice can reduce bias when the surface is mainly smooth in the
  treatment-relevant region.}

  \item{dr_score}{Selects the orthogonal score used to combine the train-fold
  baseline with the holdout dose history. \code{"aipw"} (default, matches
  historical behaviour) computes the weighted augmented IPW pseudo-outcome
  \eqn{\hat\tau_{it} + w_{it}(Y_{it} - \hat Y_{it}^{(0)} - \hat\tau_{it})}
  per cell and regresses it on the dose-lag design. \code{"plr"} implements
  the Robinson partially-linear (Frisch-Waugh) score: a separate
  conditional-mean model \eqn{\hat g(X_{it}) = E[D|X,t]} is fit on each
  train fold using the same \code{gps_model} family, the holdout dose lags
  are residualized as \eqn{D_{it} - \hat g(X_{it})}, and the outcome
  residuals are regressed on the residualized lags without GPS reweighting.
  PLR provides a useful second leg for robustness checks because its
  consistency does not require the GPS weights to be correctly normalised.}

  \item{n_folds}{Number of unit folds used for cross-fitting.}

  \item{dr_workers}{Number of workers used to evaluate DR folds. Set
  \code{dr_workers = 1} for sequential fold execution (default). Values
  greater than 1 enable parallel fold execution via multicore on supported
  platforms.}

  \item{fold_assignments}{Optional data frame with columns \code{unit} and
  \code{fold} that fixes the unit-level cross-fitting split. When supplied,
  the rows must match the current unit set exactly and use consecutive fold
  labels \code{1, 2, ..., K}. This is mainly useful when deterministic refits
  are needed for replay or formal sensitivity routines.}

  \item{lambda, rank_max, lambda_fraction, lambda_selection, lambda_grid,
  nlambda, lambda_min_ratio, cv_rounds, cv_block_size, cv_workers, cv_top_k,
  cv_coarse_to_fine, cv_coarse_nlambda, cv_warm_starts, washout, lag_order,
  outer_maxit, fe_maxit, soft_maxit, tol, fe_tol, zero_tolerance, seed,
  verbose}{Passed through to the fold-specific \code{cdmc_fit()} calls. The
  built-in second stage is always linear in the contemporaneous and lagged
  doses. For empirical workflows with enough zero-dose support inside each
  training fold, prefer \code{lambda = NULL} together with
  \code{lambda_selection = "cv"} so each fold reuses the blocked
  cross-validation tuning rule rather than the heuristic fallback. Use
  \code{cv_workers > 1} to parallelize lambda-candidate scoring within each
  fold's blocked cross-validation routine on supported platforms. Set
  \code{cv_top_k} to a positive integer to enable approximate candidate
  screening after each round (all candidates are still evaluated in the first
  round), set \code{cv_coarse_to_fine = TRUE} to run an additional two-stage
  coarse-to-fine search over lambda candidates, use
  \code{cv_coarse_nlambda} to control coarse-stage grid size, and set
  \code{cv_warm_starts = TRUE} to allow sequential
  \code{cv_workers = 1} tuning to reuse baseline fit state across lambda
  candidates within each fold.}

  \item{lambda_control}{Optional \code{cdmc_lambda_control} object bundling
  every lambda / CV scalar (see \code{\link{cdmc_lambda_control}}). When
  supplied, its fields override the corresponding scalar arguments listed
  above. The legacy scalars remain on the public API for back-compatibility
  and continue to work unchanged when \code{lambda_control = NULL} (the
  default).}

  \item{x}{A \code{cdmc_dr_fit} object.}

  \item{...}{Unused additional arguments for the print method.}
}
\section{Three orthogonal design flags}{

\code{cdmc_dr_fit()} exposes three independent design choices that should not
be conflated. They can be set in any combination and are useful as robustness
contrasts:

\describe{
  \item{\code{dr_score} (\code{"aipw"} default, or \code{"plr"})}{Selects
    the cross-fit score. \code{"aipw"} forms the augmented pseudo-outcome
    \eqn{\hat\tau_{it} + w_{it}(Y_{it} - \hat Y_{it}^{(0)} - \hat\tau_{it})}
    and regresses it on the dose-lag design (doubly robust: bias vanishes if
    either the baseline or the GPS weights are correct). \code{"plr"} uses
    the Robinson partialling-out score: it fits a separate conditional-mean
    model \eqn{\hat g(X_{it}) = E[D_{it} \mid X_{it}, t]} on each train fold
    via the same \code{gps_model} family, residualises the holdout dose lags
    as \eqn{D_{it} - \hat g(X_{it})}, and regresses outcome residuals on the
    residualised lags without GPS reweighting.}
  \item{\code{baseline_weighting} (\code{"gps"} default, or \code{"none"})}{
    Controls whether the train-fold matrix-completion baseline is fit weighted
    by the per-cell GPS / balancing weights (\code{"gps"}, preserves the AIPW
    orthogonality argument under fixed train-fold weights) or unweighted
    (\code{"none"}, fastest and closer to a classical AIPW argument in which
    the baseline nuisance is estimated independently of the propensity
    nuisance).}
  \item{\code{weight_method} (\code{"gaussian_gps"}, \code{"kernel_gps"},
    \code{"cbps"}, \code{"entropy_balance"}, \code{"kernel_balance"},
    \code{"adaptive_balance"}, or \code{"external"} when \code{weights} are
    supplied)}{Selects the per-cell weighting scheme used by both the train
    baseline (when \code{baseline_weighting = "gps"}) and the AIPW
    augmentation. Has no effect on \code{dr_score = "plr"} consistency, only
    on efficiency.}
}

Use \code{\link{cdmc_dr_score_contrast}} to compare two fits that differ in
\code{dr_score} (typically \code{"aipw"} vs.\ \code{"plr"}) on the same panel
as a robustness check; substantial disagreement, particularly at
\code{dose_lag0}, signals that one of the nuisances is poorly specified. The
\code{summary()} method reports which combination of flags was used.
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
smooths of selected numeric covariates, a regression tree, a regression
forest, a gradient boosting model, or a convex stack of those learners on the
same raw covariate and time-effect inputs. When
\code{weight_method = "gaussian_gps"}, the
conditional and marginal assignment densities use Gaussian approximations.
When \code{weight_method = "kernel_gps"}, the same mean model is paired with
kernel-smoothed residual and marginal dose densities; \code{gps_bandwidth}
controls that smoother when it is supplied. When
\code{weight_method = "cbps"}, the function reuses the same linear or spline
covariate formula inside fold-local \code{CBPS::CBPS()} fits. Because that
path returns realized balancing weights rather than an out-of-sample density
predictor, the training and holdout samples are weighted by separate fold-local
CBPS fits. When \code{weight_method = "entropy_balance"}, the function fits a
fold-local exponential-tilting problem that preserves the weighted means of the
chosen linear or spline balance basis while targeting first-order centered
dose-feature moments. As with CBPS, this produces realized balancing weights,
so the training and holdout samples are weighted by separate fold-local fits.
When \code{weight_method = "kernel_balance"}, the function fits a reduced-
kernel version of that same positive balancing problem, replacing the raw
interaction basis with a radial-basis feature map over the covariates. Because
the kernel basis is intentionally reduced, balance is treated as approximate and
the helper may simplify the kernel basis further when a fold is too small for
the requested nonlinear basis.
When \code{weight_method = "adaptive_balance"}, the function fits several of
the built-in balancing helpers on each fold, scores them on a common balance
basis, and then uses the best-scoring candidate separately on the training and
holdout samples for that fold.

For holdout treated histories, the function forms an augmented pseudo outcome of
the form
\deqn{\hat{\tau}^{DR} = \hat{\tau}^{model} + w (\hat{\tau}^{obs} - \hat{\tau}^{model})}
and regresses that pseudo outcome on the lagged dose history without an
intercept.

This is a research-use implementation of a cross-fitted doubly robust-style
linear estimator. The built-in internal weighting path is intentionally simple:
it currently offers Gaussian and kernel GPS nuisance variants with linear,
spline, GAM, tree, forest, boosted, or stacked mean models plus fold-local
CBPS, entropy-balance, kernel-balance, and adaptive-balance variants built on
linear or spline covariate formulas. External weights remain preferable when
the design calls for richer balancing or for assignment models beyond that
built-in menu.

When \code{max_weight = "adaptive"}, the function computes the clipping rule
inside each training fold and reuses the resulting cap for the realized
training and holdout weights from that fold's internal nuisance fit.

When the fitted object is supplied to \code{cdmc_bootstrap()}, the bootstrap can
summarize not only coefficients and mean DR effects but also lag-specific
average fitted contributions and user-defined linear dose-path contrasts.
For empirical workflows, the same tuning advice as \code{cdmc_fit()} applies:
prefer fold-local blocked cross-validation when the fold support is rich enough
to hold out contiguous control blocks safely.
}
\value{
An object of class \code{cdmc_dr_fit}. Important components include:

\item{data}{The reordered panel augmented with out-of-fold baseline,
out-of-fold residual treatment effects, augmented pseudo outcomes, fold ids, and
the normalized weights used in each holdout fold.}

\item{baseline}{A list containing the out-of-fold baseline matrix and a compact
summary of the fold-specific nuisance fits.}
\item{weight_diagnostics}{A list containing overall and fold-level summaries of
the realized out-of-fold weights, including effective sample size and, for
internal weighting methods, clipping counts and clipping shares.}

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
  lambda = NULL,
  lambda_selection = "cv",
  cv_rounds = 3L,
  cv_block_size = 2L,
  rank_max = 3,
  lag_order = 1,
  seed = 1
)

print(fit)

fit_balance <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "entropy_balance",
  gps_model = "linear",
  entropy_balance_iterations = 500L,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 1
)

print(fit_balance)

fit_kernel_balance <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "kernel_balance",
  gps_model = "linear",
  kernel_balance_centers = 8L,
  kernel_balance_iterations = 500L,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 1
)

print(fit_kernel_balance)

fit_adaptive_balance <- cdmc_dr_fit(
  data = panel,
  outcome = "y",
  dose = "dose",
  unit = "unit",
  time = "time",
  covariates = "x1",
  weight_method = "adaptive_balance",
  adaptive_balance_methods = c("entropy_balance", "kernel_balance"),
  kernel_balance_centers = 8L,
  kernel_balance_iterations = 500L,
  lag_order = 1,
  lambda = 0.2,
  rank_max = 3,
  seed = 1
)

print(fit_adaptive_balance)
}
}
\seealso{
\code{\link{cdmc_fit}}, \code{\link{cdmc_bootstrap}},
\code{\link{cdmc_dr_score_contrast}} for AIPW vs.\ PLR robustness contrasts,
and the weighting helpers \code{\link{cdmc_adaptive_balance_weights}},
\code{\link{cdmc_cbps_weights}}, \code{\link{cdmc_entropy_balance_weights}},
\code{\link{cdmc_kernel_balance_weights}}. The \code{dr-weighting} vignette
(\code{vignette("dr-weighting", package = "causaldosemc")}) walks through the
three orthogonal flags with worked examples.
}