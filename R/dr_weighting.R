cdmc_resolve_dr_weight_method <- function(weights, weight_method) {
  if (is.null(weight_method)) {
    return(if (is.null(weights)) "gaussian_gps" else "external")
  }

  resolved <- match.arg(weight_method, c("external", "gaussian_gps", "kernel_gps", "cbps", "entropy_balance", "kernel_balance", "adaptive_balance"))
  if (identical(resolved, "external") && is.null(weights)) {
    stop("weights must be supplied when weight_method = 'external'.", call. = FALSE)
  }
  if (!identical(resolved, "external") && !is.null(weights)) {
    stop("Do not supply weights when weight_method uses an internal weighting nuisance model.", call. = FALSE)
  }

  resolved
}

cdmc_resolve_gps_bandwidth <- function(gps_bandwidth) {
  cdmc_resolve_scalar(gps_bandwidth, "gps_bandwidth", type = "numeric",
                      min = 0, strict_lower = TRUE, allow_null = TRUE)
}

cdmc_resolve_cbps_iterations <- function(cbps_iterations) {
  cdmc_resolve_scalar(cbps_iterations, "cbps_iterations", type = "integer", min = 1)
}

cdmc_resolve_gps_forest_trees <- function(gps_forest_trees) {
  cdmc_resolve_scalar(gps_forest_trees, "gps_forest_trees", type = "integer", min = 1)
}

cdmc_resolve_gps_forest_mtry <- function(gps_forest_mtry) {
  cdmc_resolve_scalar(gps_forest_mtry, "gps_forest_mtry", type = "integer",
                      min = 1, allow_null = TRUE)
}

cdmc_resolve_gps_forest_min_node_size <- function(gps_forest_min_node_size) {
  cdmc_resolve_scalar(gps_forest_min_node_size, "gps_forest_min_node_size",
                      type = "integer", min = 1, allow_null = TRUE)
}

cdmc_resolve_gps_boost_trees <- function(gps_boost_trees) {
  cdmc_resolve_scalar(gps_boost_trees, "gps_boost_trees", type = "integer", min = 1)
}

cdmc_resolve_gps_boost_depth <- function(gps_boost_depth) {
  cdmc_resolve_scalar(gps_boost_depth, "gps_boost_depth", type = "integer", min = 1)
}

cdmc_resolve_gps_boost_shrinkage <- function(gps_boost_shrinkage) {
  cdmc_resolve_scalar(gps_boost_shrinkage, "gps_boost_shrinkage", type = "numeric",
                      min = 0, strict_lower = TRUE)
}

cdmc_resolve_gps_boost_min_obs_node <- function(gps_boost_min_obs_node) {
  cdmc_resolve_scalar(gps_boost_min_obs_node, "gps_boost_min_obs_node",
                      type = "integer", min = 1)
}

cdmc_resolve_max_weight_spec <- function(max_weight) {
  if (is.null(max_weight)) {
    return(NULL)
  }

  if (is.character(max_weight)) {
    if (length(max_weight) != 1L || !identical(max_weight, "adaptive")) {
      stop("max_weight must be NULL, 'adaptive', or a single positive numeric value.", call. = FALSE)
    }
    return(max_weight)
  }

  if (!is.numeric(max_weight) || length(max_weight) != 1L || !is.finite(max_weight) || max_weight <= 0) {
    stop("max_weight must be NULL, 'adaptive', or a single positive numeric value.", call. = FALSE)
  }

  as.numeric(max_weight)
}

cdmc_compute_adaptive_max_weight <- function(weights) {
  finite_weights <- as.numeric(weights)
  finite_weights <- finite_weights[is.finite(finite_weights) & finite_weights > 0]
  if (length(finite_weights) < 1L) {
    return(NULL)
  }

  quantile_cap <- as.numeric(stats::quantile(finite_weights, probs = 0.99, names = FALSE, type = 8))
  size_cap <- sqrt(length(finite_weights))
  adaptive_cap <- min(quantile_cap, size_cap)

  if (!is.finite(adaptive_cap) || adaptive_cap <= 0) {
    return(NULL)
  }

  as.numeric(adaptive_cap)
}

cdmc_resolve_internal_max_weight <- function(weights, max_weight = "adaptive") {
  max_weight <- cdmc_resolve_max_weight_spec(max_weight)
  if (is.null(max_weight)) {
    return(NULL)
  }

  if (identical(max_weight, "adaptive")) {
    return(cdmc_compute_adaptive_max_weight(weights))
  }

  max_weight
}

cdmc_cap_internal_weights <- function(weights, max_weight = NULL) {
  resolved_max_weight <- cdmc_resolve_internal_max_weight(weights, max_weight = max_weight)
  if (is.null(resolved_max_weight)) {
    return(as.numeric(weights))
  }

  as.numeric(pmin(weights, resolved_max_weight))
}

cdmc_effective_sample_size <- function(weights) {
  weights <- as.numeric(weights)
  weights <- weights[is.finite(weights) & weights > 0]
  if (length(weights) < 1L) {
    return(NA_real_)
  }

  total_weight <- sum(weights)
  squared_weight <- sum(weights ^ 2)
  if (!is.finite(total_weight) || total_weight <= 0 || !is.finite(squared_weight) || squared_weight <= 0) {
    return(NA_real_)
  }

  total_weight ^ 2 / squared_weight
}

cdmc_weight_diagnostic_summary <- function(weights, reference_weights = NULL) {
  weights <- as.numeric(weights)
  if (!is.null(reference_weights)) {
    reference_weights <- as.numeric(reference_weights)
    if (length(reference_weights) != length(weights)) {
      stop("reference_weights must have the same length as weights.", call. = FALSE)
    }
  }

  valid <- is.finite(weights) & weights > 0
  if (!is.null(reference_weights)) {
    valid <- valid & is.finite(reference_weights) & reference_weights >= 0
  }

  if (!any(valid)) {
    return(list(
      n_weights = 0L,
      total_weight = NA_real_,
      mean_weight = NA_real_,
      sd_weight = NA_real_,
      min_weight = NA_real_,
      max_weight = NA_real_,
      ess = NA_real_,
      ess_fraction = NA_real_,
      clipped_count = if (is.null(reference_weights)) NA_integer_ else 0L,
      clipped_share = if (is.null(reference_weights)) NA_real_ else NA_real_
    ))
  }

  weights <- weights[valid]
  ess <- cdmc_effective_sample_size(weights)
  clipped_count <- if (is.null(reference_weights)) {
    NA_integer_
  } else {
    reference_weights <- reference_weights[valid]
    as.integer(sum(reference_weights > weights + 1e-10))
  }

  list(
    n_weights = length(weights),
    total_weight = sum(weights),
    mean_weight = mean(weights),
    sd_weight = if (length(weights) > 1L) stats::sd(weights) else 0,
    min_weight = min(weights),
    max_weight = max(weights),
    ess = ess,
    ess_fraction = if (is.finite(ess)) ess / length(weights) else NA_real_,
    clipped_count = clipped_count,
    clipped_share = if (is.na(clipped_count)) NA_real_ else clipped_count / length(weights)
  )
}

cdmc_bind_rows_simple <- function(rows) {
  if (length(rows) < 1L) {
    return(data.frame())
  }

  do.call(
    rbind,
    lapply(rows, function(row) as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE))
  )
}

cdmc_fold_weight_diagnostic_row <- function(
  fold,
  weight_method,
  requested_max_weight,
  applied_max_weight,
  observed_weights,
  observed_reference = NULL,
  dr_sample_weights,
  dr_sample_reference = NULL
) {
  observed_summary <- cdmc_weight_diagnostic_summary(
    weights = observed_weights,
    reference_weights = observed_reference
  )
  dr_sample_summary <- cdmc_weight_diagnostic_summary(
    weights = dr_sample_weights,
    reference_weights = dr_sample_reference
  )

  c(
    list(
      fold = fold,
      weight_method = weight_method,
      requested_max_weight = if (is.null(requested_max_weight)) NA_character_ else as.character(requested_max_weight),
      applied_max_weight = applied_max_weight %||% NA_real_
    ),
    stats::setNames(observed_summary, paste0("observed_", names(observed_summary))),
    stats::setNames(dr_sample_summary, paste0("dr_sample_", names(dr_sample_summary)))
  )
}

cdmc_build_dr_weight_diagnostics <- function(fitted_panel, by_fold, weight_method) {
  observed_weights <- fitted_panel$.cdmc_weight[!is.na(fitted_panel$.cdmc_observed) & fitted_panel$.cdmc_observed]
  dr_sample_weights <- fitted_panel$.cdmc_weight[!is.na(fitted_panel$.cdmc_dr_sample) & fitted_panel$.cdmc_dr_sample]

  observed_summary <- cdmc_weight_diagnostic_summary(observed_weights)
  dr_sample_summary <- cdmc_weight_diagnostic_summary(dr_sample_weights)

  if (nrow(by_fold) > 0L) {
    if (all(is.finite(by_fold$observed_clipped_count))) {
      observed_summary$clipped_count <- sum(by_fold$observed_clipped_count)
      observed_summary$clipped_share <- if (observed_summary$n_weights > 0L) {
        observed_summary$clipped_count / observed_summary$n_weights
      } else {
        NA_real_
      }
    }
    if (all(is.finite(by_fold$dr_sample_clipped_count))) {
      dr_sample_summary$clipped_count <- sum(by_fold$dr_sample_clipped_count)
      dr_sample_summary$clipped_share <- if (dr_sample_summary$n_weights > 0L) {
        dr_sample_summary$clipped_count / dr_sample_summary$n_weights
      } else {
        NA_real_
      }
    }
  }

  list(
    weight_method = weight_method,
    by_fold = by_fold,
    overall = list(
      observed = observed_summary,
      dr_sample = dr_sample_summary
    )
  )
}

cdmc_fit_kernel_density <- function(x, bandwidth = NULL, n = 512L) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) == 0L) {
    stop("Kernel density estimation requires at least one finite observation.", call. = FALSE)
  }

  resolved_bandwidth <- bandwidth %||% suppressWarnings(stats::bw.nrd0(x))
  if (!is.finite(resolved_bandwidth) || resolved_bandwidth <= sqrt(.Machine$double.eps)) {
    fallback_scale <- stats::sd(x)
    if (!is.finite(fallback_scale) || fallback_scale <= sqrt(.Machine$double.eps)) {
      fallback_scale <- 1
    }
    resolved_bandwidth <- fallback_scale * length(x)^(-1 / 5)
  }

  density_fit <- stats::density(
    x = x,
    bw = resolved_bandwidth,
    n = max(256L, as.integer(n)),
    from = min(x) - 4 * resolved_bandwidth,
    to = max(x) + 4 * resolved_bandwidth
  )

  list(
    bandwidth = resolved_bandwidth,
    x = density_fit$x,
    y = density_fit$y,
    support_min = min(density_fit$x),
    support_max = max(density_fit$x)
  )
}

# Evaluate the KDE at `values`, falling back to a Gaussian-tail extrapolation
# anchored at the boundary density when `values` lies outside the support
# captured by the training KDE grid. Returns the predicted density and a
# logical vector flagging out-of-support evaluations so callers can surface
# diagnostics about positivity violations in cross-fit folds.
cdmc_predict_kernel_density <- function(object, values) {
  values <- as.numeric(values)
  approximated <- stats::approx(
    x = object$x,
    y = object$y,
    xout = values,
    yleft = NA_real_,
    yright = NA_real_,
    rule = 1L
  )$y

  out_of_support <- is.na(approximated) & is.finite(values)
  if (any(out_of_support)) {
    bandwidth <- object$bandwidth
    if (!is.finite(bandwidth) || bandwidth <= sqrt(.Machine$double.eps)) {
      bandwidth <- sqrt(.Machine$double.eps)
    }
    boundary_y_left <- object$y[1L]
    boundary_y_right <- object$y[length(object$y)]
    support_min <- object$support_min %||% min(object$x)
    support_max <- object$support_max %||% max(object$x)

    left_idx <- out_of_support & values < support_min
    right_idx <- out_of_support & values > support_max

    if (any(left_idx)) {
      delta <- (support_min - values[left_idx]) / bandwidth
      approximated[left_idx] <- boundary_y_left * exp(-0.5 * delta * delta)
    }
    if (any(right_idx)) {
      delta <- (values[right_idx] - support_max) / bandwidth
      approximated[right_idx] <- boundary_y_right * exp(-0.5 * delta * delta)
    }
  }

  density_values <- pmax(as.numeric(approximated), .Machine$double.xmin)
  attr(density_values, "out_of_support") <- out_of_support
  density_values
}

cdmc_fit_gaussian_gps <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest", "stack", "boost"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_stack_models = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  gps_boost_trees = 200L,
  gps_boost_depth = 2L,
  gps_boost_shrinkage = 0.05,
  gps_boost_min_obs_node = 10L,
  stabilize_weights = TRUE,
  max_weight = "adaptive"
) {
  mean_model <- cdmc_fit_gps_mean_model(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_stack_models = gps_stack_models,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_boost_trees = gps_boost_trees,
    gps_boost_depth = gps_boost_depth,
    gps_boost_shrinkage = gps_boost_shrinkage,
    gps_boost_min_obs_node = gps_boost_min_obs_node
  )
  mean_fit <- mean_model$fit
  fitted_mean <- cdmc_predict_gps_mean_model(mean_fit, newdata = mean_model$prepare_newdata(data))
  residuals <- data[[dose]] - fitted_mean
  sigma_cond <- sqrt(mean(residuals ^ 2))
  if (!is.finite(sigma_cond) || sigma_cond <= sqrt(.Machine$double.eps)) {
    sigma_cond <- sqrt(.Machine$double.eps)
  }

  observed_dose <- stats::na.omit(data[[dose]])
  marginal_mean <- mean(observed_dose)
  marginal_sd <- stats::sd(observed_dose)
  if (!is.finite(marginal_sd) || marginal_sd <= sqrt(.Machine$double.eps)) {
    marginal_sd <- sigma_cond
  }

  gps_fit <- list(
    weight_method = "gaussian_gps",
    formula = mean_model$formula,
    lm_fit = mean_fit,
    sigma_cond = sigma_cond,
    marginal_mean = marginal_mean,
    marginal_sd = marginal_sd,
    covariates = covariates,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_stack_models = mean_model$gps_stack_models %||% NULL,
    gps_bandwidth = NULL,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_boost_trees = gps_boost_trees,
    gps_boost_depth = gps_boost_depth,
    gps_boost_shrinkage = gps_boost_shrinkage,
    gps_boost_min_obs_node = gps_boost_min_obs_node,
    gps_prepare_newdata = mean_model$prepare_newdata,
    requested_max_weight = cdmc_resolve_max_weight_spec(max_weight),
    applied_max_weight = NULL
  )

  training_weights <- cdmc_predict_gaussian_gps_weights(
    object = gps_fit,
    newdata = data,
    dose = dose,
    stabilize = stabilize_weights,
    max_weight = NULL
  )
  gps_fit$applied_max_weight <- cdmc_resolve_internal_max_weight(training_weights, max_weight = gps_fit$requested_max_weight)
  gps_fit
}

cdmc_fit_kernel_gps <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest", "stack", "boost"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_stack_models = NULL,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  gps_boost_trees = 200L,
  gps_boost_depth = 2L,
  gps_boost_shrinkage = 0.05,
  gps_boost_min_obs_node = 10L,
  stabilize_weights = TRUE,
  max_weight = "adaptive"
) {
  mean_model <- cdmc_fit_gps_mean_model(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_stack_models = gps_stack_models,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_boost_trees = gps_boost_trees,
    gps_boost_depth = gps_boost_depth,
    gps_boost_shrinkage = gps_boost_shrinkage,
    gps_boost_min_obs_node = gps_boost_min_obs_node
  )
  mean_fit <- mean_model$fit
  fitted_mean <- cdmc_predict_gps_mean_model(mean_fit, newdata = mean_model$prepare_newdata(data))
  residuals <- data[[dose]] - fitted_mean
  residual_density <- cdmc_fit_kernel_density(residuals, bandwidth = gps_bandwidth)
  marginal_density <- cdmc_fit_kernel_density(data[[dose]], bandwidth = gps_bandwidth)

  gps_fit <- list(
    weight_method = "kernel_gps",
    formula = mean_model$formula,
    lm_fit = mean_fit,
    residual_density = residual_density,
    marginal_density = marginal_density,
    covariates = covariates,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_stack_models = mean_model$gps_stack_models %||% NULL,
    gps_bandwidth = residual_density$bandwidth,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_boost_trees = gps_boost_trees,
    gps_boost_depth = gps_boost_depth,
    gps_boost_shrinkage = gps_boost_shrinkage,
    gps_boost_min_obs_node = gps_boost_min_obs_node,
    gps_prepare_newdata = mean_model$prepare_newdata,
    requested_max_weight = cdmc_resolve_max_weight_spec(max_weight),
    applied_max_weight = NULL
  )

  training_weights <- cdmc_predict_kernel_gps_weights(
    object = gps_fit,
    newdata = data,
    dose = dose,
    stabilize = stabilize_weights,
    max_weight = NULL
  )
  gps_fit$applied_max_weight <- cdmc_resolve_internal_max_weight(training_weights, max_weight = gps_fit$requested_max_weight)
  gps_fit
}

cdmc_predict_gaussian_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  conditional_mean <- cdmc_predict_gps_mean_model(object$lm_fit, newdata = object$gps_prepare_newdata(newdata))
  conditional_density <- stats::dnorm(
    newdata[[dose]],
    mean = conditional_mean,
    sd = object$sigma_cond
  )
  numerator <- if (isTRUE(stabilize)) {
    stats::dnorm(
      newdata[[dose]],
      mean = object$marginal_mean,
      sd = object$marginal_sd
    )
  } else {
    rep(1, nrow(newdata))
  }

  weights <- numerator / pmax(conditional_density, .Machine$double.xmin)
  cdmc_cap_internal_weights(weights, max_weight = max_weight)
}

cdmc_predict_kernel_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  conditional_mean <- cdmc_predict_gps_mean_model(object$lm_fit, newdata = object$gps_prepare_newdata(newdata))
  residuals <- newdata[[dose]] - conditional_mean
  conditional_density <- cdmc_predict_kernel_density(object$residual_density, residuals)
  numerator <- if (isTRUE(stabilize)) {
    cdmc_predict_kernel_density(object$marginal_density, newdata[[dose]])
  } else {
    rep(1, nrow(newdata))
  }

  weights <- as.numeric(numerator) / pmax(as.numeric(conditional_density), .Machine$double.xmin)

  cond_oos <- attr(conditional_density, "out_of_support")
  marg_oos <- attr(numerator, "out_of_support")
  out_of_support_flags <- (if (is.null(cond_oos)) FALSE else cond_oos) |
    (if (is.null(marg_oos)) FALSE else marg_oos)

  capped <- cdmc_cap_internal_weights(weights, max_weight = max_weight)
  attr(capped, "out_of_support") <- out_of_support_flags
  capped
}

cdmc_fit_cbps_weights <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  max_weight = "adaptive"
) {
  prepare_cbps_internal_fit <- get("cdmc_prepare_cbps_internal_fit", mode = "function")
  prepare_cbps_internal_fit(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    cbps_standardize = cbps_standardize,
    cbps_method = cbps_method,
    cbps_iterations = cbps_iterations,
    cbps_twostep = cbps_twostep,
    sample_weights = NULL,
    max_weight = max_weight
  )
}

cdmc_fit_internal_gps <- function(
  weight_method,
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest", "stack", "boost"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_stack_models = NULL,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  gps_boost_trees = 200L,
  gps_boost_depth = 2L,
  gps_boost_shrinkage = 0.05,
  gps_boost_min_obs_node = 10L,
  stabilize_weights = TRUE,
  max_weight = "adaptive"
) {
  weight_method <- match.arg(weight_method, c("gaussian_gps", "kernel_gps"))

  if (identical(weight_method, "gaussian_gps")) {
    return(cdmc_fit_gaussian_gps(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = gps_time_effects,
      gps_model = gps_model,
      gps_df = gps_df,
      gps_spline_covariates = gps_spline_covariates,
      gps_stack_models = gps_stack_models,
      gps_forest_trees = gps_forest_trees,
      gps_forest_mtry = gps_forest_mtry,
      gps_forest_min_node_size = gps_forest_min_node_size,
      gps_boost_trees = gps_boost_trees,
      gps_boost_depth = gps_boost_depth,
      gps_boost_shrinkage = gps_boost_shrinkage,
      gps_boost_min_obs_node = gps_boost_min_obs_node,
      stabilize_weights = stabilize_weights,
      max_weight = max_weight
    ))
  }

  cdmc_fit_kernel_gps(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_stack_models = gps_stack_models,
    gps_bandwidth = gps_bandwidth,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_boost_trees = gps_boost_trees,
    gps_boost_depth = gps_boost_depth,
    gps_boost_shrinkage = gps_boost_shrinkage,
    gps_boost_min_obs_node = gps_boost_min_obs_node,
    stabilize_weights = stabilize_weights,
    max_weight = max_weight
  )
}

# Generic cross-fit transport for balancing weights.
#
# Balancing solvers (CBPS, entropy_balance, kernel_balance, adaptive_balance)
# do not expose a closed-form predict step that maps a learned tilt to new
# observations: their constraint design is centered on the training sample and
# their parameters are not directly identified out-of-sample. To preserve the
# Chernozhukov-style cross-fitting property we therefore fit a *weight
# regression* on the training fold -- regressing log(train weights) onto the
# same covariate basis used by the GPS nuisance -- and apply that predictor to
# the holdout fold. The resulting holdout weights depend on holdout features
# only through the train-fold model, so the second-stage estimator remains
# Neyman-orthogonal under standard nuisance-rate conditions.
#
# This replaces the previous in-sample re-fit on the holdout fold, which was
# silently anti-conservative (the holdout score reused information already
# absorbed into the holdout weights).
cdmc_fit_balance_weight_transport <- function(
  train_data,
  train_weights,
  dose,
  weight_covariates,
  time,
  gps_time_effects,
  gps_model,
  gps_df,
  gps_spline_covariates
) {
  if (length(train_weights) != nrow(train_data)) {
    stop("train_weights must have one entry per row of train_data.", call. = FALSE)
  }

  observed_rows <- if (".cdmc_observed" %in% names(train_data)) {
    !is.na(train_data$.cdmc_observed) & train_data$.cdmc_observed
  } else {
    rep(TRUE, nrow(train_data))
  }
  fit_data <- train_data[observed_rows, , drop = FALSE]
  fit_weights <- train_weights[observed_rows]

  positive <- is.finite(fit_weights) & fit_weights > 0
  if (sum(positive) < 2L) {
    stop("Cross-fit weight transport requires at least two positive training weights.", call. = FALSE)
  }
  fit_data <- fit_data[positive, , drop = FALSE]
  fit_weights <- fit_weights[positive]

  log_response_column <- ".cdmc_log_train_weight"
  fit_data[[log_response_column]] <- log(fit_weights)

  # The GPS model family for the transport mirrors the user's chosen GPS
  # specification, except that we treat dose as an additional covariate so the
  # transport can capture the weight's dose-dependence (which is exactly the
  # balancing model's structural mechanism). Tree/forest/boost/stack models
  # require dose as a feature; spline/linear/gam expand it as a covariate.
  transport_covariates <- unique(c(weight_covariates %||% character(0), dose))
  transport_spline_covariates <- intersect(
    transport_covariates,
    c(gps_spline_covariates %||% character(0), dose)
  )

  mean_model <- cdmc_fit_gps_mean_model(
    data = fit_data,
    dose = log_response_column,
    covariates = transport_covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = transport_spline_covariates,
    gps_stack_models = NULL,
    gps_forest_trees = 200L,
    gps_forest_mtry = NULL,
    gps_forest_min_node_size = NULL,
    gps_boost_trees = 200L,
    gps_boost_depth = 2L,
    gps_boost_shrinkage = 0.05,
    gps_boost_min_obs_node = 10L
  )

  # Calibrate so the transported predictor reproduces the training weight scale
  # exactly: the residual between mean(exp(log_pred)) and mean(train_weights)
  # is absorbed into a multiplicative constant. This avoids drift in the
  # AIPW augmentation when the transport is biased on the log scale.
  fitted_log <- cdmc_predict_gps_mean_model(
    mean_model$fit,
    newdata = mean_model$prepare_newdata(fit_data)
  )
  fitted_weight <- exp(fitted_log)
  scale_factor <- mean(fit_weights) / mean(fitted_weight)
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    scale_factor <- 1
  }

  list(
    mean_model = mean_model,
    scale_factor = scale_factor,
    log_response_column = log_response_column
  )
}

cdmc_predict_balance_weight_transport <- function(transport, newdata) {
  log_pred <- cdmc_predict_gps_mean_model(
    transport$mean_model$fit,
    newdata = transport$mean_model$prepare_newdata(newdata)
  )
  weights <- exp(as.numeric(log_pred)) * transport$scale_factor
  weights[!is.finite(weights) | weights < 0] <- 0
  weights
}

cdmc_predict_internal_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  if (missing(max_weight)) {
    max_weight <- object$applied_max_weight %||% NULL
  }

  if (identical(object$weight_method, "gaussian_gps")) {
    return(cdmc_predict_gaussian_gps_weights(
      object = object,
      newdata = newdata,
      dose = dose,
      stabilize = stabilize,
      max_weight = max_weight
    ))
  }

  if (identical(object$weight_method, "kernel_gps")) {
    return(cdmc_predict_kernel_gps_weights(
      object = object,
      newdata = newdata,
      dose = dose,
      stabilize = stabilize,
      max_weight = max_weight
    ))
  }

  stop("Internal GPS prediction requires a fitted gaussian_gps or kernel_gps nuisance model.", call. = FALSE)
}
