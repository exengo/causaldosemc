cdmc_assign_unit_folds <- function(unit_levels, n_folds) {
  n_folds <- as.integer(n_folds)
  if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
    stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
  }
  if (n_folds > length(unit_levels)) {
    stop("n_folds cannot exceed the number of units.", call. = FALSE)
  }

  shuffled_units <- sample(unit_levels)
  fold_sizes <- rep(length(unit_levels) %/% n_folds, n_folds)
  fold_sizes[seq_len(length(unit_levels) %% n_folds)] <- fold_sizes[seq_len(length(unit_levels) %% n_folds)] + 1L
  fold_id <- integer(length(unit_levels))
  names(fold_id) <- unit_levels
  start_index <- 1L

  for (fold_index in seq_len(n_folds)) {
    end_index <- start_index + fold_sizes[fold_index] - 1L
    current_units <- shuffled_units[seq.int(start_index, end_index)]
    fold_id[current_units] <- fold_index
    start_index <- end_index + 1L
  }

  data.frame(
    unit = unit_levels,
    fold = unname(fold_id[unit_levels]),
    stringsAsFactors = FALSE
  )
}

cdmc_resolve_unit_folds <- function(unit_levels, n_folds, fold_assignments = NULL) {
  if (is.null(fold_assignments)) {
    return(cdmc_assign_unit_folds(unit_levels, n_folds = n_folds))
  }

  if (!is.data.frame(fold_assignments)) {
    stop("fold_assignments must be NULL or a data frame with columns 'unit' and 'fold'.", call. = FALSE)
  }

  required_columns <- c("unit", "fold")
  missing_columns <- setdiff(required_columns, names(fold_assignments))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "fold_assignments is missing required columns: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  provided_units <- as.character(fold_assignments$unit)
  if (anyNA(provided_units) || anyDuplicated(provided_units)) {
    stop("fold_assignments must contain each unit exactly once with no missing values.", call. = FALSE)
  }

  unit_keys <- as.character(unit_levels)
  missing_units <- setdiff(unit_keys, provided_units)
  extra_units <- setdiff(provided_units, unit_keys)
  if (length(missing_units) > 0L || length(extra_units) > 0L) {
    message_parts <- character(0)
    if (length(missing_units) > 0L) {
      message_parts <- c(message_parts, sprintf("missing units: %s", paste(missing_units, collapse = ", ")))
    }
    if (length(extra_units) > 0L) {
      message_parts <- c(message_parts, sprintf("unknown units: %s", paste(extra_units, collapse = ", ")))
    }

    stop(
      sprintf("fold_assignments must match the current unit set exactly (%s).", paste(message_parts, collapse = "; ")),
      call. = FALSE
    )
  }

  unit_index <- match(unit_keys, provided_units)
  fold_values_raw <- fold_assignments$fold[unit_index]
  if (anyNA(fold_values_raw) || any(!is.finite(fold_values_raw))) {
    stop("fold_assignments$fold must be finite for every unit.", call. = FALSE)
  }
  if (any(abs(fold_values_raw - round(fold_values_raw)) > sqrt(.Machine$double.eps))) {
    stop("fold_assignments$fold must contain integer-valued fold labels.", call. = FALSE)
  }

  fold_values <- as.integer(round(fold_values_raw))
  if (any(fold_values < 1L)) {
    stop("fold_assignments$fold must contain positive fold labels.", call. = FALSE)
  }

  unique_folds <- sort(unique(fold_values))
  if (!identical(unique_folds, seq_along(unique_folds))) {
    stop("fold_assignments$fold must use consecutive fold labels 1, 2, ..., K.", call. = FALSE)
  }
  if (length(unique_folds) < 2L) {
    stop("fold_assignments must define at least two folds.", call. = FALSE)
  }

  if (!is.null(n_folds)) {
    n_folds <- as.integer(n_folds)
    if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
      stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
    }
    if (length(unique_folds) != n_folds) {
      stop(
        sprintf(
          "Provided fold_assignments define %d folds, but n_folds = %d.",
          length(unique_folds),
          n_folds
        ),
        call. = FALSE
      )
    }
  }

  data.frame(
    unit = unit_levels,
    fold = fold_values,
    stringsAsFactors = FALSE
  )
}

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
  if (is.null(gps_bandwidth)) {
    return(NULL)
  }

  if (!is.numeric(gps_bandwidth) || length(gps_bandwidth) != 1L || !is.finite(gps_bandwidth) || gps_bandwidth <= 0) {
    stop("gps_bandwidth must be NULL or a single positive numeric value.", call. = FALSE)
  }

  as.numeric(gps_bandwidth)
}

cdmc_resolve_cbps_iterations <- function(cbps_iterations) {
  if (!is.numeric(cbps_iterations) || length(cbps_iterations) != 1L || !is.finite(cbps_iterations) || cbps_iterations < 1) {
    stop("cbps_iterations must be a positive integer.", call. = FALSE)
  }

  as.integer(cbps_iterations)
}

cdmc_resolve_gps_forest_trees <- function(gps_forest_trees) {
  if (!is.numeric(gps_forest_trees) || length(gps_forest_trees) != 1L || !is.finite(gps_forest_trees) || gps_forest_trees < 1) {
    stop("gps_forest_trees must be a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_trees)
}

cdmc_resolve_gps_forest_mtry <- function(gps_forest_mtry) {
  if (is.null(gps_forest_mtry)) {
    return(NULL)
  }

  if (!is.numeric(gps_forest_mtry) || length(gps_forest_mtry) != 1L || !is.finite(gps_forest_mtry) || gps_forest_mtry < 1) {
    stop("gps_forest_mtry must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_mtry)
}

cdmc_resolve_gps_forest_min_node_size <- function(gps_forest_min_node_size) {
  if (is.null(gps_forest_min_node_size)) {
    return(NULL)
  }

  if (!is.numeric(gps_forest_min_node_size) || length(gps_forest_min_node_size) != 1L || !is.finite(gps_forest_min_node_size) || gps_forest_min_node_size < 1) {
    stop("gps_forest_min_node_size must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_min_node_size)
}

cdmc_resolve_gps_boost_trees <- function(gps_boost_trees) {
  if (!is.numeric(gps_boost_trees) || length(gps_boost_trees) != 1L || !is.finite(gps_boost_trees) || gps_boost_trees < 1) {
    stop("gps_boost_trees must be a positive integer.", call. = FALSE)
  }

  as.integer(gps_boost_trees)
}

cdmc_resolve_gps_boost_depth <- function(gps_boost_depth) {
  if (!is.numeric(gps_boost_depth) || length(gps_boost_depth) != 1L || !is.finite(gps_boost_depth) || gps_boost_depth < 1) {
    stop("gps_boost_depth must be a positive integer.", call. = FALSE)
  }

  as.integer(gps_boost_depth)
}

cdmc_resolve_gps_boost_shrinkage <- function(gps_boost_shrinkage) {
  if (!is.numeric(gps_boost_shrinkage) || length(gps_boost_shrinkage) != 1L || !is.finite(gps_boost_shrinkage) || gps_boost_shrinkage <= 0) {
    stop("gps_boost_shrinkage must be a single positive numeric value.", call. = FALSE)
  }

  as.numeric(gps_boost_shrinkage)
}

cdmc_resolve_gps_boost_min_obs_node <- function(gps_boost_min_obs_node) {
  if (!is.numeric(gps_boost_min_obs_node) || length(gps_boost_min_obs_node) != 1L || !is.finite(gps_boost_min_obs_node) || gps_boost_min_obs_node < 1) {
    stop("gps_boost_min_obs_node must be a positive integer.", call. = FALSE)
  }

  as.integer(gps_boost_min_obs_node)
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

cdmc_prepare_dr_fold_sample <- function(object, prepared, baseline_hat, weight_matrix = NULL, fold_id = NA_integer_,
                                        dr_score = c("aipw", "plr"), dose_residual_matrix = NULL) {
  dr_score <- match.arg(dr_score)
  tau_matrix <- prepared$y_matrix - baseline_hat
  # In PLR mode, the lag design is built from residualized doses
  # D_t - E[D_t | X_t]. The orthogonal score reduces to OLS of the outcome
  # residual on residualized lags (no augmentation; weights are not needed
  # for orthogonality, only for efficiency, so we drop them here).
  if (identical(dr_score, "plr")) {
    if (is.null(dose_residual_matrix)) {
      stop("dose_residual_matrix must be supplied when dr_score = 'plr'.", call. = FALSE)
    }
    lag_array <- cdmc_build_lagged_doses(dose_residual_matrix, lag_order = object$lag_order)
  } else {
    lag_array <- cdmc_build_lagged_doses(prepared$dose_matrix, lag_order = object$lag_order)
  }
  valid_history <- apply(!is.na(lag_array), c(1, 2), all)
  # PLR exposure magnitude is computed on the *original* doses so the
  # zero-dose support condition is unchanged across the two scores.
  exposure_lag_array <- if (identical(dr_score, "plr")) {
    cdmc_build_lagged_doses(prepared$dose_matrix, lag_order = object$lag_order)
  } else {
    lag_array
  }
  exposure_magnitude <- cdmc_lag_exposure_magnitude(exposure_lag_array)

  sample_mask <- valid_history & !is.na(tau_matrix) & exposure_magnitude > 0
  if (!is.null(weight_matrix)) {
    sample_mask <- sample_mask & weight_matrix > 0
  }

  sample_indices <- which(sample_mask, arr.ind = TRUE)
  if (nrow(sample_indices) == 0L) {
    return(list(
      sample_mask = sample_mask,
      design = matrix(0, nrow = 0L, ncol = object$lag_order + 1L),
      pseudo_tau = numeric(0),
      tau_observed = numeric(0),
      tau_model = numeric(0),
      sample_table = data.frame()
    ))
  }

  design <- do.call(
    cbind,
    lapply(seq_len(dim(lag_array)[3L]), function(index) lag_array[, , index][sample_indices])
  )
  colnames(design) <- dimnames(lag_array)[[3L]]

  coefficient_vector <- numeric(ncol(design))
  names(coefficient_vector) <- colnames(design)
  if (length(object$effect$coefficients) > 0L) {
    common_columns <- intersect(names(object$effect$coefficients), names(coefficient_vector))
    coefficient_vector[common_columns] <- object$effect$coefficients[common_columns]
  }

  tau_observed <- tau_matrix[sample_indices]
  tau_model <- as.vector(design %*% coefficient_vector)
  sample_weights <- if (is.null(weight_matrix)) rep(1, length(tau_observed)) else weight_matrix[sample_indices]
  if (identical(dr_score, "plr")) {
    # Robinson PLR pseudo-outcome: y_resid = Y - m(X) only. Orthogonality
    # is achieved by using residualized lags in `design`. Per-cell weights
    # default to 1 so the downstream OLS is canonical FWL; users can still
    # supply weights for efficiency, in which case we honour them.
    pseudo_tau <- tau_observed
  } else {
    pseudo_tau <- tau_model + sample_weights * (tau_observed - tau_model)
  }

  linear_indices <- (sample_indices[, 1L] - 1L) * prepared$n_times + sample_indices[, 2L]
  sample_table <- prepared$data[linear_indices, c(names(prepared$data)), drop = FALSE]
  sample_table[, colnames(design)] <- design
  sample_table$.cdmc_fold_id <- fold_id
  sample_table$.cdmc_tau_observed <- tau_observed
  sample_table$.cdmc_tau_model <- tau_model
  sample_table$.cdmc_tau_dr <- pseudo_tau
  sample_table$.cdmc_weight_used <- sample_weights

  list(
    sample_mask = sample_mask,
    design = design,
    pseudo_tau = pseudo_tau,
    tau_observed = tau_observed,
    tau_model = tau_model,
    sample_weights = sample_weights,
    sample_indices = sample_indices,
    sample_table = sample_table
  )
}

cdmc_dr_fit <- function(
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
) {
  n_folds_requested <- if (missing(n_folds)) NULL else n_folds
  lambda_selection <- match.arg(lambda_selection)
  gps_model <- match.arg(gps_model)
  cbps_method <- match.arg(cbps_method)
  baseline_weighting <- match.arg(baseline_weighting)
  dr_score <- match.arg(dr_score)
  weight_method <- cdmc_resolve_dr_weight_method(weights = weights, weight_method = weight_method)
  max_weight <- if (identical(weight_method, "external")) {
    NULL
  } else {
    cdmc_resolve_max_weight_spec(max_weight)
  }

  if (!is.null(weight_covariates) && !is.character(weight_covariates)) {
    stop("weight_covariates must be NULL or a character vector of column names.", call. = FALSE)
  }

  if (!is.numeric(gps_df) || length(gps_df) != 1L || !is.finite(gps_df) || gps_df < 1) {
    stop("gps_df must be a positive integer.", call. = FALSE)
  }
  gps_df <- as.integer(gps_df)
  gps_stack_models <- if (identical(gps_model, "stack")) {
    cdmc_resolve_gps_stack_models(gps_stack_models)
  } else {
    NULL
  }
  if ((identical(gps_model, "gam") || (identical(gps_model, "stack") && "gam" %in% gps_stack_models)) && gps_df < 3L) {
    stop("gps_df must be at least 3 when gps_model = 'gam'.", call. = FALSE)
  }
  gps_spline_covariates <- cdmc_resolve_gps_spline_covariates(weight_covariates, gps_spline_covariates)
  gps_bandwidth <- cdmc_resolve_gps_bandwidth(gps_bandwidth)
  gps_forest_trees <- cdmc_resolve_gps_forest_trees(gps_forest_trees)
  gps_forest_mtry <- cdmc_resolve_gps_forest_mtry(gps_forest_mtry)
  gps_forest_min_node_size <- cdmc_resolve_gps_forest_min_node_size(gps_forest_min_node_size)
  gps_boost_trees <- cdmc_resolve_gps_boost_trees(gps_boost_trees)
  gps_boost_depth <- cdmc_resolve_gps_boost_depth(gps_boost_depth)
  gps_boost_shrinkage <- cdmc_resolve_gps_boost_shrinkage(gps_boost_shrinkage)
  gps_boost_min_obs_node <- cdmc_resolve_gps_boost_min_obs_node(gps_boost_min_obs_node)
  cbps_iterations <- cdmc_resolve_cbps_iterations(cbps_iterations)
  resolve_adaptive_balance_methods <- get("cdmc_resolve_adaptive_balance_methods", mode = "function")
  resolve_entropy_balance_degree <- get("cdmc_resolve_entropy_balance_degree", mode = "function")
  resolve_entropy_balance_iterations <- get("cdmc_resolve_entropy_balance_iterations", mode = "function")
  resolve_entropy_balance_reltol <- get("cdmc_resolve_entropy_balance_reltol", mode = "function")
  resolve_kernel_balance_centers <- get("cdmc_resolve_kernel_balance_centers", mode = "function")
  resolve_kernel_balance_bandwidth <- get("cdmc_resolve_kernel_balance_bandwidth", mode = "function")
  adaptive_balance_methods <- if (identical(weight_method, "adaptive_balance")) {
    resolve_adaptive_balance_methods(adaptive_balance_methods)
  } else {
    NULL
  }
  entropy_balance_degree <- resolve_entropy_balance_degree(entropy_balance_degree)
  entropy_balance_iterations <- resolve_entropy_balance_iterations(entropy_balance_iterations)
  entropy_balance_reltol <- resolve_entropy_balance_reltol(entropy_balance_reltol)
  kernel_balance_degree <- resolve_entropy_balance_degree(kernel_balance_degree)
  kernel_balance_centers <- resolve_kernel_balance_centers(kernel_balance_centers)
  kernel_balance_bandwidth <- resolve_kernel_balance_bandwidth(kernel_balance_bandwidth)
  kernel_balance_iterations <- resolve_entropy_balance_iterations(kernel_balance_iterations)
  kernel_balance_reltol <- resolve_entropy_balance_reltol(kernel_balance_reltol)

  if (!is.logical(cbps_standardize) || length(cbps_standardize) != 1L || is.na(cbps_standardize)) {
    stop("cbps_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(cbps_twostep) || length(cbps_twostep) != 1L || is.na(cbps_twostep)) {
    stop("cbps_twostep must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(entropy_balance_standardize) || length(entropy_balance_standardize) != 1L || is.na(entropy_balance_standardize)) {
    stop("entropy_balance_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(kernel_balance_standardize) || length(kernel_balance_standardize) != 1L || is.na(kernel_balance_standardize)) {
    stop("kernel_balance_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(dr_workers) || length(dr_workers) != 1L || !is.finite(dr_workers) || dr_workers < 1 || dr_workers != floor(dr_workers)) {
    stop("dr_workers must be a positive integer.", call. = FALSE)
  }
  dr_workers <- as.integer(dr_workers)
  if (!is.numeric(cv_workers) || length(cv_workers) != 1L || !is.finite(cv_workers) || cv_workers < 1 || cv_workers != floor(cv_workers)) {
    stop("cv_workers must be a positive integer.", call. = FALSE)
  }
  cv_workers <- as.integer(cv_workers)
  if (!is.null(cv_top_k)) {
    if (!is.numeric(cv_top_k) || length(cv_top_k) != 1L || !is.finite(cv_top_k) || cv_top_k < 1 || cv_top_k != floor(cv_top_k)) {
      stop("cv_top_k must be NULL or a positive integer.", call. = FALSE)
    }
    cv_top_k <- as.integer(cv_top_k)
  }
  if (!is.logical(cv_coarse_to_fine) || length(cv_coarse_to_fine) != 1L || is.na(cv_coarse_to_fine)) {
    stop("cv_coarse_to_fine must be TRUE or FALSE.", call. = FALSE)
  }
  cv_coarse_to_fine <- isTRUE(cv_coarse_to_fine)
  if (!is.null(cv_coarse_nlambda)) {
    if (!is.numeric(cv_coarse_nlambda) || length(cv_coarse_nlambda) != 1L || !is.finite(cv_coarse_nlambda) || cv_coarse_nlambda < 3 || cv_coarse_nlambda != floor(cv_coarse_nlambda)) {
      stop("cv_coarse_nlambda must be NULL or an integer >= 3.", call. = FALSE)
    }
    cv_coarse_nlambda <- as.integer(cv_coarse_nlambda)
  }
  if (!is.logical(cv_warm_starts) || length(cv_warm_starts) != 1L || is.na(cv_warm_starts)) {
    stop("cv_warm_starts must be TRUE or FALSE.", call. = FALSE)
  }
  cv_warm_starts <- isTRUE(cv_warm_starts)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  full_prepared <- cdmc_prepare_panel(
    data = data,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    zero_tolerance = zero_tolerance
  )
  full_data <- full_prepared$data
  full_data$.cdmc_global_unit_index <- full_data$.cdmc_unit_index
  full_data$.cdmc_global_time_index <- full_data$.cdmc_time_index

  full_eligible_mask <- cdmc_build_eligible_mask(
    dose_matrix = full_prepared$dose_matrix,
    zero_tolerance = zero_tolerance,
    washout = washout
  )
  if (identical(weight_method, "external")) {
    weight_seed_info <- cdmc_prepare_panel_weights(
      weights = weights,
      data = full_data,
      n_units = full_prepared$n_units,
      n_times = full_prepared$n_times,
      eligible_mask = full_eligible_mask,
      column_name = ".cdmc_dr_weight_raw"
    )
    full_data <- weight_seed_info$data
    full_eligible_mask <- full_eligible_mask & weight_seed_info$matrix > 0
  } else {
    missing_weight_covariates <- setdiff(weight_covariates, names(full_data))
    if (length(missing_weight_covariates) > 0L) {
      stop(
        sprintf("Missing weight_covariates columns: %s.", paste(missing_weight_covariates, collapse = ", ")),
        call. = FALSE
      )
    }
    if (weight_method %in% c("cbps", "entropy_balance", "kernel_balance", "adaptive_balance") && !gps_model %in% c("linear", "spline")) {
      stop(
        "gps_model is not available when weight_method uses a balancing-only path. Use 'linear' or 'spline', or switch to gaussian_gps or kernel_gps.",
        call. = FALSE
      )
    }
    if (weight_method %in% c("cbps", "entropy_balance", "kernel_balance", "adaptive_balance") &&
        length(weight_covariates %||% character(0)) == 0L &&
        !isTRUE(gps_time_effects)) {
      stop(
        "Internal balancing weights require at least one weight_covariate or gps_time_effects = TRUE.",
        call. = FALSE
      )
    }
  }
  cdmc_validate_support(full_eligible_mask)

  fold_assignments <- cdmc_resolve_unit_folds(
    unit_levels = full_prepared$unit_levels,
    n_folds = if (is.null(fold_assignments)) n_folds else n_folds_requested,
    fold_assignments = fold_assignments
  )
  n_folds <- length(unique(fold_assignments$fold))
  fold_lookup <- stats::setNames(fold_assignments$fold, fold_assignments$unit)

  baseline_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_model_oof <- matrix(0, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_dr_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  fitted_dr_oof <- matrix(0, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  sample_mask_oof <- matrix(FALSE, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  weight_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)

  design_parts <- vector("list", n_folds)
  pseudo_parts <- vector("list", n_folds)
  fold_samples <- vector("list", n_folds)
  fold_summaries <- vector("list", n_folds)
  weight_diagnostic_rows <- vector("list", n_folds)
  sample_counter <- 0L

  use_fold_parallel <- dr_workers > 1L
  if (use_fold_parallel && identical(.Platform$OS.type, "windows")) {
    warning(
      "Parallel DR fold execution currently uses multicore execution and is not available on Windows. Falling back to sequential execution.",
      call. = FALSE
    )
    use_fold_parallel <- FALSE
    dr_workers <- 1L
  }
  fold_worker_count <- if (use_fold_parallel) min(dr_workers, n_folds) else 1L
  fold_cv_workers <- if (use_fold_parallel && cv_workers > 1L) {
    if (verbose) {
      message("dr_workers > 1 detected; using cv_workers = 1 inside each fold to avoid nested oversubscription")
    }
    1L
  } else {
    cv_workers
  }
  fold_seeds <- if (use_fold_parallel) sample.int(.Machine$integer.max, n_folds, replace = TRUE) else integer(0)

  run_dr_fold <- function(fold_index) {
    if (use_fold_parallel) {
      set.seed(fold_seeds[[fold_index]])
    }

    holdout_units <- fold_assignments$unit[fold_assignments$fold == fold_index]
    train_data <- full_data[!full_data[[unit]] %in% holdout_units, , drop = FALSE]
    holdout_data <- full_data[full_data[[unit]] %in% holdout_units, , drop = FALSE]

    if (verbose) {
      message(sprintf("cross-fit fold %d/%d: %d holdout units", fold_index, n_folds, length(holdout_units)))
    }

    gps_fit <- NULL
    holdout_gps_fit <- NULL
    weight_transport <- NULL
    adaptive_balance_weights <- get("cdmc_adaptive_balance_weights", mode = "function")
    entropy_balance_weights <- get("cdmc_entropy_balance_weights", mode = "function")
    kernel_balance_weights <- get("cdmc_kernel_balance_weights", mode = "function")
    if (identical(weight_method, "cbps")) {
      gps_fit <- cdmc_fit_cbps_weights(
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates,
        cbps_standardize = cbps_standardize,
        cbps_method = cbps_method,
        cbps_iterations = cbps_iterations,
        cbps_twostep = cbps_twostep,
        max_weight = max_weight
      )
      train_data$.cdmc_dr_weight_raw <- gps_fit$weights
      weight_transport <- cdmc_fit_balance_weight_transport(
        train_data = train_data,
        train_weights = gps_fit$weights,
        dose = dose,
        weight_covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_cap_internal_weights(
        cdmc_predict_balance_weight_transport(weight_transport, holdout_data),
        max_weight = gps_fit$applied_max_weight
      )
    } else if (identical(weight_method, "entropy_balance")) {
      gps_fit <- entropy_balance_weights(
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        time_effects = gps_time_effects,
        model = gps_model,
        df = gps_df,
        spline_covariates = gps_spline_covariates,
        degree = entropy_balance_degree,
        standardize = entropy_balance_standardize,
        iterations = entropy_balance_iterations,
        reltol = entropy_balance_reltol,
        max_weight = max_weight
      )
      train_data$.cdmc_dr_weight_raw <- gps_fit$weights
      weight_transport <- cdmc_fit_balance_weight_transport(
        train_data = train_data,
        train_weights = gps_fit$weights,
        dose = dose,
        weight_covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_cap_internal_weights(
        cdmc_predict_balance_weight_transport(weight_transport, holdout_data),
        max_weight = gps_fit$applied_max_weight
      )
    } else if (identical(weight_method, "kernel_balance")) {
      gps_fit <- kernel_balance_weights(
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        time_effects = gps_time_effects,
        model = gps_model,
        df = gps_df,
        spline_covariates = gps_spline_covariates,
        degree = kernel_balance_degree,
        n_centers = kernel_balance_centers,
        bandwidth = kernel_balance_bandwidth,
        standardize = kernel_balance_standardize,
        iterations = kernel_balance_iterations,
        reltol = kernel_balance_reltol,
        max_weight = max_weight
      )
      train_data$.cdmc_dr_weight_raw <- gps_fit$weights
      weight_transport <- cdmc_fit_balance_weight_transport(
        train_data = train_data,
        train_weights = gps_fit$weights,
        dose = dose,
        weight_covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_cap_internal_weights(
        cdmc_predict_balance_weight_transport(weight_transport, holdout_data),
        max_weight = gps_fit$applied_max_weight
      )
    } else if (identical(weight_method, "adaptive_balance")) {
      gps_fit <- adaptive_balance_weights(
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        time_effects = gps_time_effects,
        model = gps_model,
        df = gps_df,
        spline_covariates = gps_spline_covariates,
        methods = adaptive_balance_methods,
        cbps_standardize = cbps_standardize,
        cbps_method = cbps_method,
        cbps_iterations = cbps_iterations,
        cbps_twostep = cbps_twostep,
        entropy_balance_degree = entropy_balance_degree,
        entropy_balance_standardize = entropy_balance_standardize,
        entropy_balance_iterations = entropy_balance_iterations,
        entropy_balance_reltol = entropy_balance_reltol,
        kernel_balance_degree = kernel_balance_degree,
        kernel_balance_centers = kernel_balance_centers,
        kernel_balance_bandwidth = kernel_balance_bandwidth,
        kernel_balance_standardize = kernel_balance_standardize,
        kernel_balance_iterations = kernel_balance_iterations,
        kernel_balance_reltol = kernel_balance_reltol,
        max_weight = max_weight
      )
      train_data$.cdmc_dr_weight_raw <- gps_fit$weights
      weight_transport <- cdmc_fit_balance_weight_transport(
        train_data = train_data,
        train_weights = gps_fit$weights,
        dose = dose,
        weight_covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_cap_internal_weights(
        cdmc_predict_balance_weight_transport(weight_transport, holdout_data),
        max_weight = gps_fit$applied_max_weight
      )
    } else if (!identical(weight_method, "external")) {
      gps_fit <- cdmc_fit_internal_gps(
        weight_method = weight_method,
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
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
      train_data$.cdmc_dr_weight_raw <- cdmc_predict_internal_gps_weights(
        object = gps_fit,
        newdata = train_data,
        dose = dose,
        stabilize = stabilize_weights
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_predict_internal_gps_weights(
        object = gps_fit,
        newdata = holdout_data,
        dose = dose,
        stabilize = stabilize_weights
      )
    }

    # Train baseline weighting choice. The standard AIPW orthogonal score is
    # asymptotically valid with an unweighted outcome regression (the GPS
    # weights enter only through the augmentation term). Weighting the train
    # baseline by the same weights tilts the imputed counterfactual toward
    # the treatment-relevant covariate distribution, which can reduce bias
    # when the dose-response surface is smoother in that distribution but is
    # not required for orthogonality. Default "gps" preserves the historical
    # behaviour; pass `baseline_weighting = "none"` for canonical AIPW.
    train_weight_arg <- if (identical(baseline_weighting, "none")) NULL else ".cdmc_dr_weight_raw"

    fold_fit <- cdmc_fit(
      data = train_data,
      outcome = outcome,
      dose = dose,
      unit = unit,
      time = time,
      covariates = covariates,
      weights = train_weight_arg,
      lambda = lambda,
      rank_max = rank_max,
      lambda_fraction = lambda_fraction,
      lambda_selection = lambda_selection,
      lambda_grid = lambda_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_workers = fold_cv_workers,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = cv_coarse_to_fine,
      cv_coarse_nlambda = cv_coarse_nlambda,
      cv_warm_starts = cv_warm_starts,
      washout = washout,
      lag_order = lag_order,
      effect_model = "linear",
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      zero_tolerance = zero_tolerance,
      seed = NULL,
      verbose = FALSE
    )

    holdout_prepared <- cdmc_prepare_panel(
      data = holdout_data,
      outcome = outcome,
      dose = dose,
      unit = unit,
      time = time,
      covariates = covariates,
      zero_tolerance = zero_tolerance
    )
    holdout_eligible_mask <- cdmc_build_eligible_mask(
      dose_matrix = holdout_prepared$dose_matrix,
      zero_tolerance = zero_tolerance,
      washout = washout
    )
    holdout_weight_info <- cdmc_prepare_panel_weights(
      weights = ".cdmc_dr_weight_raw",
      data = holdout_prepared$data,
      n_units = holdout_prepared$n_units,
      n_times = holdout_prepared$n_times,
      eligible_mask = holdout_eligible_mask,
      column_name = ".cdmc_weight"
    )
    holdout_prepared$data <- holdout_weight_info$data
    holdout_weight_matrix <- holdout_weight_info$matrix
    holdout_eligible_mask <- holdout_eligible_mask & holdout_weight_matrix > 0

    holdout_reference_normalized <- NULL
    if (!identical(weight_method, "external")) {
      if (weight_method %in% c("cbps", "entropy_balance", "kernel_balance", "adaptive_balance")) {
        # Reference: the transported holdout weights (from the train-fold
        # weight regression). This is the same predictor used to populate
        # `.cdmc_dr_weight_raw` above; we reattach it on the prepared (padded)
        # data frame so the diagnostic summary aligns with downstream rows.
        reference_lookup <- stats::setNames(
          holdout_data$.cdmc_dr_weight_raw,
          paste(holdout_data$.cdmc_global_unit_index, holdout_data$.cdmc_global_time_index, sep = "\r")
        )
        holdout_reference_raw <- unname(reference_lookup[
          paste(holdout_prepared$data$.cdmc_global_unit_index, holdout_prepared$data$.cdmc_global_time_index, sep = "\r")
        ])
      } else {
        holdout_reference_raw <- cdmc_predict_internal_gps_weights(
          object = gps_fit,
          newdata = holdout_prepared$data,
          dose = dose,
          stabilize = stabilize_weights,
          max_weight = NULL
        )
      }
      holdout_reference_normalized <- holdout_reference_raw / holdout_weight_info$scale
    }

    holdout_baseline <- cdmc_predict_holdout_baseline(
      object = fold_fit,
      prepared = holdout_prepared,
      eligible_mask = holdout_eligible_mask,
      weight_matrix = holdout_weight_matrix
    )

    # Robinson PLR: fit a conditional dose mean model E[D | X, t] on the
    # train fold and predict on the holdout. Residualizing dose lags by
    # this conditional mean gives the orthogonal Frisch-Waugh score.
    dose_residual_matrix <- NULL
    plr_dose_model <- NULL
    if (identical(dr_score, "plr")) {
      plr_train <- train_data
      plr_train[[dose]] <- as.numeric(plr_train[[dose]])
      plr_dose_model <- cdmc_fit_gps_mean_model(
        data = plr_train,
        dose = dose,
        covariates = weight_covariates,
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
      holdout_dose_pred_data <- plr_dose_model$prepare_newdata(holdout_prepared$data)
      holdout_dose_hat <- as.numeric(cdmc_predict_gps_mean_model(plr_dose_model$fit, newdata = holdout_dose_pred_data))
      dose_residual_vector <- holdout_prepared$data[[dose]] - holdout_dose_hat
      dose_residual_matrix <- matrix(
        dose_residual_vector,
        nrow = holdout_prepared$n_units,
        ncol = holdout_prepared$n_times,
        byrow = TRUE
      )
      dose_residual_matrix[is.na(holdout_prepared$dose_matrix)] <- NA_real_
    }

    fold_sample <- cdmc_prepare_dr_fold_sample(
      object = fold_fit,
      prepared = holdout_prepared,
      baseline_hat = holdout_baseline$baseline_hat,
      weight_matrix = holdout_weight_matrix,
      fold_id = fold_index,
      dr_score = dr_score,
      dose_residual_matrix = dose_residual_matrix
    )

    global_indices <- cbind(
      holdout_prepared$data$.cdmc_global_unit_index,
      holdout_prepared$data$.cdmc_global_time_index
    )
    baseline_flat <- cdmc_flatten_matrix(holdout_baseline$baseline_hat)
    tau_flat <- cdmc_flatten_matrix(holdout_prepared$y_matrix - holdout_baseline$baseline_hat)
    weight_flat <- holdout_prepared$data$.cdmc_weight

    sample_global_indices <- matrix(integer(0), ncol = 2)
    sample_tau_model <- numeric(0)
    sample_pseudo_tau <- numeric(0)
    design_part <- NULL
    sample_table <- NULL

    if (nrow(fold_sample$design) > 0L) {
      sample_linear_indices <- (fold_sample$sample_indices[, 1L] - 1L) * holdout_prepared$n_times + fold_sample$sample_indices[, 2L]
      sample_global_indices <- cbind(
        holdout_prepared$data$.cdmc_global_unit_index[sample_linear_indices],
        holdout_prepared$data$.cdmc_global_time_index[sample_linear_indices]
      )
      sample_tau_model <- fold_sample$tau_model
      sample_pseudo_tau <- fold_sample$pseudo_tau
      design_part <- fold_sample$design
      sample_table <- fold_sample$sample_table
    }

    sample_linear_indices <- if (nrow(fold_sample$sample_indices) > 0L) {
      (fold_sample$sample_indices[, 1L] - 1L) * holdout_prepared$n_times + fold_sample$sample_indices[, 2L]
    } else {
      integer(0)
    }
    observed_rows <- !is.na(holdout_prepared$data$.cdmc_observed) & holdout_prepared$data$.cdmc_observed
    weight_diag_row <- cdmc_fold_weight_diagnostic_row(
      fold = fold_index,
      weight_method = weight_method,
      requested_max_weight = if (is.null(gps_fit)) NULL else gps_fit$requested_max_weight %||% NULL,
      applied_max_weight = if (is.null(gps_fit)) NULL else gps_fit$applied_max_weight %||% NULL,
      observed_weights = holdout_prepared$data$.cdmc_weight[observed_rows],
      observed_reference = if (is.null(holdout_reference_normalized)) NULL else holdout_reference_normalized[observed_rows],
      dr_sample_weights = fold_sample$sample_weights,
      dr_sample_reference = if (is.null(holdout_reference_normalized) || length(sample_linear_indices) < 1L) {
        NULL
      } else {
        holdout_reference_normalized[sample_linear_indices]
      }
    )

    fold_summary <- list(
      fold = fold_index,
      holdout_units = holdout_units,
      weight_method = weight_method,
      lambda = fold_fit$lambda,
      lambda_method = fold_fit$lambda_tuning$method,
      baseline_solver = fold_fit$baseline$solver,
      gps_formula = if (is.null(gps_fit)) NULL else gps_fit$formula,
      gps_weight_method = if (is.null(gps_fit)) NULL else gps_fit$weight_method,
      gps_model = if (is.null(gps_fit)) NULL else gps_fit$gps_model,
      gps_df = if (is.null(gps_fit)) NULL else gps_fit$gps_df,
      gps_stack_models = if (is.null(gps_fit)) NULL else gps_fit$gps_stack_models,
      gps_bandwidth = if (is.null(gps_fit)) NULL else gps_fit$gps_bandwidth,
      gps_forest_trees = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_trees,
      gps_forest_mtry = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_mtry,
      gps_forest_min_node_size = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_min_node_size,
      gps_boost_trees = if (is.null(gps_fit)) NULL else gps_fit$gps_boost_trees,
      gps_boost_depth = if (is.null(gps_fit)) NULL else gps_fit$gps_boost_depth,
      gps_boost_shrinkage = if (is.null(gps_fit)) NULL else gps_fit$gps_boost_shrinkage,
      gps_boost_min_obs_node = if (is.null(gps_fit)) NULL else gps_fit$gps_boost_min_obs_node,
      requested_max_weight = if (is.null(gps_fit)) NULL else gps_fit$requested_max_weight %||% NULL,
      applied_max_weight = if (is.null(gps_fit)) NULL else gps_fit$applied_max_weight %||% NULL,
      gps_sigma_cond = if (is.null(gps_fit)) NULL else gps_fit$sigma_cond %||% NULL,
      cbps_requested_method = if (is.null(gps_fit) || !identical(weight_method, "cbps")) NULL else gps_fit$requested_method,
      cbps_train_method = if (is.null(gps_fit) || !identical(weight_method, "cbps")) NULL else gps_fit$method,
      cbps_holdout_method = NULL,
      weight_transport_model = if (is.null(weight_transport)) NULL else gps_model,
      weight_transport_scale_factor = if (is.null(weight_transport)) NULL else weight_transport$scale_factor,
      baseline_weighting = baseline_weighting,
      dr_score = dr_score,
      cbps_standardize = if (identical(weight_method, "cbps")) cbps_standardize else NULL,
      cbps_iterations = if (identical(weight_method, "cbps")) cbps_iterations else NULL,
      cbps_twostep = if (identical(weight_method, "cbps")) cbps_twostep else NULL,
      entropy_balance_degree = if (identical(weight_method, "entropy_balance")) entropy_balance_degree else NULL,
      entropy_balance_standardize = if (identical(weight_method, "entropy_balance")) entropy_balance_standardize else NULL,
      entropy_balance_iterations = if (identical(weight_method, "entropy_balance")) entropy_balance_iterations else NULL,
      entropy_balance_reltol = if (identical(weight_method, "entropy_balance")) entropy_balance_reltol else NULL,
      entropy_balance_max_abs_balance = if (is.null(gps_fit) || !identical(weight_method, "entropy_balance")) NULL else gps_fit$max_abs_balance,
      kernel_balance_degree = if (identical(weight_method, "kernel_balance")) kernel_balance_degree else NULL,
      kernel_balance_centers = if (identical(weight_method, "kernel_balance")) kernel_balance_centers else NULL,
      kernel_balance_bandwidth = if (is.null(gps_fit) || !identical(weight_method, "kernel_balance")) NULL else gps_fit$applied_bandwidth,
      kernel_balance_standardize = if (identical(weight_method, "kernel_balance")) kernel_balance_standardize else NULL,
      kernel_balance_iterations = if (identical(weight_method, "kernel_balance")) kernel_balance_iterations else NULL,
      kernel_balance_reltol = if (identical(weight_method, "kernel_balance")) kernel_balance_reltol else NULL,
      kernel_balance_max_abs_balance = if (is.null(gps_fit) || !identical(weight_method, "kernel_balance")) NULL else gps_fit$max_abs_balance,
      adaptive_balance_methods = if (identical(weight_method, "adaptive_balance")) gps_fit$methods %||% adaptive_balance_methods else NULL,
      adaptive_balance_selected_method = if (is.null(gps_fit) || !identical(weight_method, "adaptive_balance")) NULL else gps_fit$selected_method,
      adaptive_balance_candidate_scores = if (is.null(gps_fit) || !identical(weight_method, "adaptive_balance")) NULL else gps_fit$candidate_scores,
      adaptive_balance_max_abs_standardized_balance = if (is.null(gps_fit) || !identical(weight_method, "adaptive_balance")) NULL else gps_fit$max_abs_standardized_balance,
      sample_size = nrow(fold_sample$design)
    )

    list(
      fold_index = fold_index,
      global_indices = global_indices,
      baseline_flat = baseline_flat,
      tau_flat = tau_flat,
      weight_flat = weight_flat,
      sample_global_indices = sample_global_indices,
      sample_tau_model = sample_tau_model,
      sample_pseudo_tau = sample_pseudo_tau,
      design_part = design_part,
      sample_table = sample_table,
      weight_diagnostic_row = weight_diag_row,
      fold_summary = fold_summary
    )
  }

  fold_results <- if (use_fold_parallel) {
    parallel::mclapply(
      seq_len(n_folds),
      run_dr_fold,
      mc.cores = fold_worker_count,
      mc.set.seed = FALSE
    )
  } else {
    lapply(seq_len(n_folds), run_dr_fold)
  }

  for (fold_result in fold_results) {
    baseline_oof[fold_result$global_indices] <- fold_result$baseline_flat
    tau_oof[fold_result$global_indices] <- fold_result$tau_flat
    weight_oof[fold_result$global_indices] <- fold_result$weight_flat

    if (nrow(fold_result$sample_global_indices) > 0L) {
      sample_mask_oof[fold_result$sample_global_indices] <- TRUE
      tau_model_oof[fold_result$sample_global_indices] <- fold_result$sample_tau_model
      tau_dr_oof[fold_result$sample_global_indices] <- fold_result$sample_pseudo_tau

      sample_counter <- sample_counter + 1L
      design_parts[[sample_counter]] <- fold_result$design_part
      pseudo_parts[[sample_counter]] <- fold_result$sample_pseudo_tau
      fold_samples[[sample_counter]] <- fold_result$sample_table
    }

    weight_diagnostic_rows[[fold_result$fold_index]] <- fold_result$weight_diagnostic_row
    fold_summaries[[fold_result$fold_index]] <- fold_result$fold_summary
  }

  design_parts <- design_parts[seq_len(sample_counter)]
  pseudo_parts <- pseudo_parts[seq_len(sample_counter)]
  fold_samples <- fold_samples[seq_len(sample_counter)]

  if (sample_counter == 0L) {
    stop("No nonzero-dose histories were available for the doubly robust stage.", call. = FALSE)
  }

  final_design <- do.call(rbind, design_parts)
  final_pseudo_tau <- unlist(pseudo_parts, use.names = FALSE)
  final_effect_fit <- cdmc_fit_weighted_regression(
    design = final_design,
    response = final_pseudo_tau,
    weights = NULL
  )

  for (sample_table in fold_samples) {
    design <- as.matrix(sample_table[, paste0("dose_lag", seq.int(0L, lag_order)), drop = FALSE])
    fitted_values <- as.vector(design %*% final_effect_fit$coefficients)
    sample_global_indices <- cbind(
      sample_table$.cdmc_global_unit_index,
      sample_table$.cdmc_global_time_index
    )
    fitted_dr_oof[sample_global_indices] <- fitted_values
  }

  fitted_panel <- full_data[, setdiff(names(full_data), c(".cdmc_global_unit_index", ".cdmc_global_time_index")), drop = FALSE]
  fitted_panel$.cdmc_fold_id <- unname(fold_lookup[fitted_panel[[unit]]])
  fitted_panel$.cdmc_y0_oof <- cdmc_flatten_matrix(baseline_oof)
  fitted_panel$.cdmc_tau_oof <- cdmc_flatten_matrix(tau_oof)
  fitted_panel$.cdmc_tau_model_oof <- cdmc_flatten_matrix(tau_model_oof)
  fitted_panel$.cdmc_tau_dr <- cdmc_flatten_matrix(tau_dr_oof)
  fitted_panel$.cdmc_tau_dr_linear <- cdmc_flatten_matrix(fitted_dr_oof)
  fitted_panel$.cdmc_dr_sample <- cdmc_flatten_matrix(sample_mask_oof)
  fitted_panel$.cdmc_weight <- cdmc_flatten_matrix(weight_oof)
  fitted_panel$.cdmc_tau_model_oof[!fitted_panel$.cdmc_observed] <- NA_real_
  fitted_panel$.cdmc_tau_dr_linear[!fitted_panel$.cdmc_observed] <- NA_real_

  weight_diagnostics <- cdmc_build_dr_weight_diagnostics(
    fitted_panel = fitted_panel,
    by_fold = cdmc_bind_rows_simple(weight_diagnostic_rows),
    weight_method = weight_method
  )

  result <- list(
    call = match.call(),
    data = fitted_panel,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    weights = ".cdmc_weight",
    n_units = full_prepared$n_units,
    n_times = full_prepared$n_times,
    unit_levels = full_prepared$unit_levels,
    time_levels = full_prepared$time_levels,
    observed_mask = full_prepared$observed_mask,
    y_matrix = full_prepared$y_matrix,
    dose_matrix = full_prepared$dose_matrix,
    weight_matrix = weight_oof,
    fold_assignments = fold_assignments,
    n_folds = n_folds,
    dr_workers = fold_worker_count,
    dr_parallel = use_fold_parallel,
    weight_method = weight_method,
    weight_covariates = weight_covariates,
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
    cbps_standardize = cbps_standardize,
    cbps_method = cbps_method,
    cbps_iterations = cbps_iterations,
    cbps_twostep = cbps_twostep,
    adaptive_balance_methods = adaptive_balance_methods,
    entropy_balance_degree = entropy_balance_degree,
    entropy_balance_standardize = entropy_balance_standardize,
    entropy_balance_iterations = entropy_balance_iterations,
    entropy_balance_reltol = entropy_balance_reltol,
    kernel_balance_degree = kernel_balance_degree,
    kernel_balance_centers = kernel_balance_centers,
    kernel_balance_bandwidth = kernel_balance_bandwidth,
    kernel_balance_standardize = kernel_balance_standardize,
    kernel_balance_iterations = kernel_balance_iterations,
    kernel_balance_reltol = kernel_balance_reltol,
    stabilize_weights = stabilize_weights,
    max_weight = max_weight,
    baseline_weighting = baseline_weighting,
    dr_score = dr_score,
    rank_max = rank_max,
    washout = as.integer(washout),
    lag_order = as.integer(lag_order),
    lambda = lambda,
    lambda_selection = lambda_selection,
    zero_tolerance = zero_tolerance,
    fit_control = list(
      weights = if (identical(weight_method, "external")) ".cdmc_weight" else NULL,
      weight_method = weight_method,
      weight_covariates = weight_covariates,
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
      cbps_standardize = cbps_standardize,
      cbps_method = cbps_method,
      cbps_iterations = cbps_iterations,
      cbps_twostep = cbps_twostep,
      adaptive_balance_methods = adaptive_balance_methods,
      entropy_balance_degree = entropy_balance_degree,
      entropy_balance_standardize = entropy_balance_standardize,
      entropy_balance_iterations = entropy_balance_iterations,
      entropy_balance_reltol = entropy_balance_reltol,
      kernel_balance_degree = kernel_balance_degree,
      kernel_balance_centers = kernel_balance_centers,
      kernel_balance_bandwidth = kernel_balance_bandwidth,
      kernel_balance_standardize = kernel_balance_standardize,
      kernel_balance_iterations = kernel_balance_iterations,
      kernel_balance_reltol = kernel_balance_reltol,
      stabilize_weights = stabilize_weights,
      max_weight = max_weight,
      baseline_weighting = baseline_weighting,
      dr_score = dr_score,
      n_folds = n_folds,
      dr_workers = fold_worker_count,
      lambda = lambda,
      rank_max = rank_max,
      lambda_fraction = lambda_fraction,
      lambda_selection = lambda_selection,
      lambda_grid = lambda_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_workers = cv_workers,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = cv_coarse_to_fine,
      cv_coarse_nlambda = cv_coarse_nlambda,
      cv_warm_starts = cv_warm_starts,
      washout = as.integer(washout),
      lag_order = as.integer(lag_order),
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      zero_tolerance = zero_tolerance
    ),
    weight_diagnostics = weight_diagnostics,
    baseline = list(
      baseline_hat = baseline_oof,
      fold_summaries = fold_summaries
    ),
    effect = list(
      coefficients = final_effect_fit$coefficients,
      fitted = fitted_dr_oof,
      tau = tau_oof,
      tau_model = tau_model_oof,
      tau_dr = tau_dr_oof,
      sample_mask = sample_mask_oof,
      design_columns = colnames(final_design),
      residuals = final_effect_fit$residuals
    )
  )

  class(result) <- "cdmc_dr_fit"
  result
}

print.cdmc_dr_fit <- function(x, ...) {
  cat("causaldosemc cross-fitted DR fit\n")
  cat(sprintf("  units x times: %d x %d\n", x$n_units, x$n_times))
  cat(sprintf("  folds: %d\n", x$n_folds))
  cat(sprintf("  weight method: %s\n", x$weight_method))
  fold_lambda_methods <- vapply(
    x$baseline$fold_summaries,
    function(summary) summary$lambda_method %||% NA_character_,
    character(1)
  )
  fold_lambda_methods <- unique(fold_lambda_methods[!is.na(fold_lambda_methods)])
  if (length(fold_lambda_methods) == 1L) {
    cat(sprintf("  lambda selection: %s\n", fold_lambda_methods[[1L]]))
    cdmc_print_empirical_lambda_note(fold_lambda_methods[[1L]])
  } else if (length(fold_lambda_methods) > 1L) {
    cat(sprintf(
      "  lambda selection: mixed (%s)\n",
      paste(fold_lambda_methods, collapse = ", ")
    ))
  }
  fold_lambdas <- vapply(
    x$baseline$fold_summaries,
    function(summary) summary$lambda %||% NA_real_,
    numeric(1)
  )
  fold_lambdas <- fold_lambdas[is.finite(fold_lambdas)]
  if (length(fold_lambdas) > 0L) {
    if (stats::sd(fold_lambdas) <= sqrt(.Machine$double.eps)) {
      cat(sprintf("  fold lambda: %.6g\n", fold_lambdas[[1L]]))
    } else {
      cat(sprintf(
        "  fold lambda range: [%.6g, %.6g]\n",
        min(fold_lambdas),
        max(fold_lambdas)
      ))
    }
  }
  if (identical(x$weight_method, "cbps")) {
    cat(sprintf("  cbps method: %s\n", x$cbps_method))
    cat(sprintf("  cbps standardized: %s\n", if (isTRUE(x$cbps_standardize)) "yes" else "no"))
  }
  if (identical(x$weight_method, "entropy_balance")) {
    cat(sprintf("  entropy-balance model: %s\n", x$gps_model))
    cat(sprintf("  entropy-balance degree: %d\n", x$entropy_balance_degree))
    cat(sprintf("  entropy-balance standardized: %s\n", if (isTRUE(x$entropy_balance_standardize)) "yes" else "no"))
    if (identical(x$gps_model, "spline")) {
      cat(sprintf("  entropy-balance spline df: %d\n", x$gps_df))
    }
  }
  if (identical(x$weight_method, "kernel_balance")) {
    cat(sprintf("  kernel-balance model: %s\n", x$gps_model))
    cat(sprintf("  kernel-balance degree: %d\n", x$kernel_balance_degree))
    cat(sprintf("  kernel-balance centers: %d\n", x$kernel_balance_centers))
    if (!is.null(x$kernel_balance_bandwidth)) {
      cat(sprintf("  kernel-balance bandwidth: %.6g\n", x$kernel_balance_bandwidth))
    }
    cat(sprintf("  kernel-balance standardized: %s\n", if (isTRUE(x$kernel_balance_standardize)) "yes" else "no"))
    if (identical(x$gps_model, "spline")) {
      cat(sprintf("  kernel-balance spline df: %d\n", x$gps_df))
    }
  }
  if (identical(x$weight_method, "adaptive_balance")) {
    selected_methods <- vapply(
      x$baseline$fold_summaries,
      function(summary) summary$adaptive_balance_selected_method %||% NA_character_,
      character(1)
    )
    selected_methods <- selected_methods[!is.na(selected_methods)]
    cat(sprintf("  adaptive-balance candidates: %s\n", paste(x$adaptive_balance_methods, collapse = ", ")))
    if (length(selected_methods) > 0L) {
      selected_counts <- table(selected_methods)
      cat(sprintf(
        "  adaptive-balance selected folds: %s\n",
        paste(sprintf("%s x%d", names(selected_counts), as.integer(selected_counts)), collapse = ", ")
      ))
    }
    cat(sprintf("  adaptive-balance model: %s\n", x$gps_model))
    if (identical(x$gps_model, "spline")) {
      cat(sprintf("  adaptive-balance spline df: %d\n", x$gps_df))
    }
  }
  if (identical(x$weight_method, "gaussian_gps") || identical(x$weight_method, "kernel_gps")) {
    cat(sprintf("  gps model: %s\n", x$gps_model))
    if (identical(x$max_weight, "adaptive")) {
      applied_caps <- vapply(
        x$baseline$fold_summaries,
        function(summary) summary$applied_max_weight %||% NA_real_,
        numeric(1)
      )
      applied_caps <- applied_caps[is.finite(applied_caps)]
      if (length(applied_caps) > 0L) {
        cat(sprintf("  gps weight cap: adaptive [%.6g, %.6g]\n", min(applied_caps), max(applied_caps)))
      } else {
        cat("  gps weight cap: adaptive\n")
      }
    } else if (is.null(x$max_weight)) {
      cat("  gps weight cap: none\n")
    } else {
      cat(sprintf("  gps weight cap: %.6g\n", x$max_weight))
    }
    if (identical(x$gps_model, "spline")) {
      cat(sprintf("  gps spline df: %d\n", x$gps_df))
    }
    if (identical(x$gps_model, "gam")) {
      cat(sprintf("  gps smooth basis k: %d\n", x$gps_df))
    }
    if (identical(x$gps_model, "stack")) {
      cat(sprintf("  gps stack models: %s\n", paste(x$gps_stack_models, collapse = ", ")))
    }
    if (identical(x$gps_model, "forest")) {
      cat(sprintf("  gps forest trees: %d\n", x$gps_forest_trees))
    }
    if (identical(x$gps_model, "boost")) {
      cat(sprintf("  gps boost trees: %d\n", x$gps_boost_trees))
      cat(sprintf("  gps boost depth: %d\n", x$gps_boost_depth))
      cat(sprintf("  gps boost shrinkage: %.6g\n", x$gps_boost_shrinkage))
    }
    if (identical(x$weight_method, "kernel_gps")) {
      cat(sprintf("  gps kernel bandwidth: %.6g\n", x$gps_bandwidth))
    }
  }
  if (identical(x$weight_method, "cbps")) {
    if (identical(x$max_weight, "adaptive")) {
      applied_caps <- vapply(
        x$baseline$fold_summaries,
        function(summary) summary$applied_max_weight %||% NA_real_,
        numeric(1)
      )
      applied_caps <- applied_caps[is.finite(applied_caps)]
      if (length(applied_caps) > 0L) {
        cat(sprintf("  cbps weight cap: adaptive [%.6g, %.6g]\n", min(applied_caps), max(applied_caps)))
      } else {
        cat("  cbps weight cap: adaptive\n")
      }
    } else if (is.null(x$max_weight)) {
      cat("  cbps weight cap: none\n")
    } else {
      cat(sprintf("  cbps weight cap: %.6g\n", x$max_weight))
    }
  }
  if (identical(x$weight_method, "entropy_balance")) {
    if (identical(x$max_weight, "adaptive")) {
      applied_caps <- vapply(
        x$baseline$fold_summaries,
        function(summary) summary$applied_max_weight %||% NA_real_,
        numeric(1)
      )
      applied_caps <- applied_caps[is.finite(applied_caps)]
      if (length(applied_caps) > 0L) {
        cat(sprintf("  entropy-balance weight cap: adaptive [%.6g, %.6g]\n", min(applied_caps), max(applied_caps)))
      } else {
        cat("  entropy-balance weight cap: adaptive\n")
      }
    } else if (is.null(x$max_weight)) {
      cat("  entropy-balance weight cap: none\n")
    } else {
      cat(sprintf("  entropy-balance weight cap: %.6g\n", x$max_weight))
    }
  }
  if (identical(x$weight_method, "kernel_balance")) {
    if (identical(x$max_weight, "adaptive")) {
      applied_caps <- vapply(
        x$baseline$fold_summaries,
        function(summary) summary$applied_max_weight %||% NA_real_,
        numeric(1)
      )
      applied_caps <- applied_caps[is.finite(applied_caps)]
      if (length(applied_caps) > 0L) {
        cat(sprintf("  kernel-balance weight cap: adaptive [%.6g, %.6g]\n", min(applied_caps), max(applied_caps)))
      } else {
        cat("  kernel-balance weight cap: adaptive\n")
      }
    } else if (is.null(x$max_weight)) {
      cat("  kernel-balance weight cap: none\n")
    } else {
      cat(sprintf("  kernel-balance weight cap: %.6g\n", x$max_weight))
    }
  }
  if (identical(x$weight_method, "adaptive_balance")) {
    if (identical(x$max_weight, "adaptive")) {
      applied_caps <- vapply(
        x$baseline$fold_summaries,
        function(summary) summary$applied_max_weight %||% NA_real_,
        numeric(1)
      )
      applied_caps <- applied_caps[is.finite(applied_caps)]
      if (length(applied_caps) > 0L) {
        cat(sprintf("  adaptive-balance weight cap: adaptive [%.6g, %.6g]\n", min(applied_caps), max(applied_caps)))
      } else {
        cat("  adaptive-balance weight cap: adaptive\n")
      }
    } else if (is.null(x$max_weight)) {
      cat("  adaptive-balance weight cap: none\n")
    } else {
      cat(sprintf("  adaptive-balance weight cap: %.6g\n", x$max_weight))
    }
  }

  diagnostics <- x$weight_diagnostics %||% NULL
  if (!is.null(diagnostics) && !is.null(diagnostics$overall)) {
    observed_diag <- diagnostics$overall$observed %||% NULL
    dr_sample_diag <- diagnostics$overall$dr_sample %||% NULL
    if (!is.null(observed_diag) && is.finite(observed_diag$ess) && observed_diag$n_weights > 0L) {
      cat(sprintf("  observed weight ESS: %.6g / %d\n", observed_diag$ess, observed_diag$n_weights))
    }
    if (!is.null(dr_sample_diag) && is.finite(dr_sample_diag$ess) && dr_sample_diag$n_weights > 0L) {
      cat(sprintf("  DR sample weight ESS: %.6g / %d\n", dr_sample_diag$ess, dr_sample_diag$n_weights))
    }
    if (identical(x$weight_method, "gaussian_gps") || identical(x$weight_method, "kernel_gps") || identical(x$weight_method, "cbps") || identical(x$weight_method, "adaptive_balance")) {
      by_fold <- diagnostics$by_fold %||% data.frame()
      if (nrow(by_fold) > 0L) {
        observed_clip <- by_fold$observed_clipped_share
        observed_clip <- observed_clip[is.finite(observed_clip)]
        if (length(observed_clip) > 0L) {
          cat(sprintf(
            "  fold observed clipped share: [%.2f%%, %.2f%%]\n",
            100 * min(observed_clip),
            100 * max(observed_clip)
          ))
        }
        dr_sample_ess_fraction <- by_fold$dr_sample_ess_fraction
        dr_sample_ess_fraction <- dr_sample_ess_fraction[is.finite(dr_sample_ess_fraction)]
        if (length(dr_sample_ess_fraction) > 0L) {
          cat(sprintf(
            "  fold DR-sample ESS fraction: [%.3f, %.3f]\n",
            min(dr_sample_ess_fraction),
            max(dr_sample_ess_fraction)
          ))
        }
      }
    }
  }
  cat(sprintf("  DR sample cells: %d\n", sum(x$effect$sample_mask, na.rm = TRUE)))

  if (length(x$effect$coefficients) > 0L) {
    cat("  DR effect coefficients:\n")
    for (coefficient_name in names(x$effect$coefficients)) {
      cat(sprintf("    %s = %.6g\n", coefficient_name, x$effect$coefficients[[coefficient_name]]))
    }
  }

  invisible(x)
}
#' Summarise a cross-fitted DR fit
#'
#' Returns a structured summary of [cdmc_dr_fit()] output: the dose-lag
#' coefficients, fold-level lambda choices, baseline rank usage, weight
#' diagnostics, and the score (`aipw` / `plr`) plus train baseline weighting
#' setting. For inferential intervals, pair this with [cdmc_bootstrap()].
#'
#' @param object A `cdmc_dr_fit` object.
#' @param ... Unused.
#'
#' @return A `summary.cdmc_dr_fit` object (a list with `coefficients`,
#'   `folds`, `weight`, `score`, and `dimensions` components) printed via
#'   the registered `print` method.
#' @export
summary.cdmc_dr_fit <- function(object, ...) {
  coefs <- object$effect$coefficients
  coefficient_table <- data.frame(
    term = names(coefs),
    estimate = unname(as.numeric(coefs)),
    stringsAsFactors = FALSE
  )

  fold_summaries <- object$baseline$fold_summaries %||% list()
  if (length(fold_summaries) > 0L) {
    folds <- data.frame(
      fold = vapply(fold_summaries, function(s) s$fold %||% NA_integer_, integer(1)),
      lambda = vapply(fold_summaries, function(s) s$lambda %||% NA_real_, numeric(1)),
      lambda_method = vapply(fold_summaries, function(s) s$lambda_method %||% NA_character_, character(1)),
      applied_max_weight = vapply(
        fold_summaries,
        function(s) {
          v <- s$applied_max_weight %||% NA_real_
          if (is.null(v) || length(v) == 0L) NA_real_ else as.numeric(v)
        },
        numeric(1)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    folds <- data.frame(
      fold = integer(0),
      lambda = numeric(0),
      lambda_method = character(0),
      applied_max_weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  diagnostics <- object$weight_diagnostics %||% list()
  weight <- list(
    method = object$weight_method,
    observed_ess = diagnostics$overall$observed$ess %||% NA_real_,
    observed_n = diagnostics$overall$observed$n_weights %||% NA_integer_,
    dr_sample_ess = diagnostics$overall$dr_sample$ess %||% NA_real_,
    dr_sample_n = diagnostics$overall$dr_sample$n_weights %||% NA_integer_
  )

  out <- list(
    coefficients = coefficient_table,
    folds = folds,
    weight = weight,
    score = list(
      dr_score = object$dr_score %||% "aipw",
      baseline_weighting = object$baseline_weighting %||% "gps"
    ),
    dimensions = list(
      n_units = object$n_units,
      n_times = object$n_times,
      n_folds = object$n_folds,
      lag_order = object$lag_order,
      n_dr_cells = sum(object$effect$sample_mask, na.rm = TRUE)
    )
  )
  class(out) <- "summary.cdmc_dr_fit"
  out
}

#' @export
print.summary.cdmc_dr_fit <- function(x, ...) {
  cat("causaldosemc cross-fitted DR summary\n")
  cat(sprintf(
    "  units x times: %d x %d   folds: %d   lag order: %d   DR cells: %d\n",
    x$dimensions$n_units, x$dimensions$n_times, x$dimensions$n_folds,
    x$dimensions$lag_order, x$dimensions$n_dr_cells
  ))
  cat(sprintf(
    "  score: %s   train baseline weighting: %s   weight method: %s\n",
    x$score$dr_score, x$score$baseline_weighting, x$weight$method
  ))
  if (is.finite(x$weight$dr_sample_ess) && isTRUE(x$weight$dr_sample_n > 0L)) {
    cat(sprintf(
      "  DR sample weight ESS: %.6g / %d (%.1f%%)\n",
      x$weight$dr_sample_ess, x$weight$dr_sample_n,
      100 * x$weight$dr_sample_ess / x$weight$dr_sample_n
    ))
  }

  if (nrow(x$folds) > 0L) {
    lambda_values <- x$folds$lambda[is.finite(x$folds$lambda)]
    if (length(lambda_values) > 0L) {
      cat(sprintf(
        "  fold lambda: [%.6g, %.6g] (methods: %s)\n",
        min(lambda_values), max(lambda_values),
        paste(unique(stats::na.omit(x$folds$lambda_method)), collapse = ", ")
      ))
    }
    caps <- x$folds$applied_max_weight[is.finite(x$folds$applied_max_weight)]
    if (length(caps) > 0L) {
      cat(sprintf("  fold applied weight cap: [%.6g, %.6g]\n", min(caps), max(caps)))
    }
  }

  if (nrow(x$coefficients) > 0L) {
    cat("  DR effect coefficients:\n")
    print(x$coefficients, row.names = FALSE)
  }

  invisible(x)
}

#' Side-by-side comparison of two DR scores
#'
#' Compares dose-lag coefficients across two [cdmc_dr_fit()] objects fitted
#' on the same panel with different `dr_score` values (typically `"aipw"`
#' and `"plr"`). Both fits should otherwise share data, lag order, weight
#' method, and folds; differing scores let the user assess robustness of
#' the headline estimates to the orthogonal score choice.
#'
#' @param aipw,plr `cdmc_dr_fit` objects to contrast. The names refer to
#'   the typical use; the function does not enforce a particular `dr_score`
#'   value, so it can also be used to compare any two compatible fits.
#'
#' @return A data frame with one row per dose-lag term containing both
#'   estimates, their absolute and relative differences, and a flag for
#'   sign agreement. Sign agreement uses a `0` tolerance equal to
#'   `sqrt(.Machine$double.eps) * max(|estimates|)`.
#' @export
cdmc_dr_score_contrast <- function(aipw, plr) {
  if (!inherits(aipw, "cdmc_dr_fit") || !inherits(plr, "cdmc_dr_fit")) {
    stop("Both arguments must be 'cdmc_dr_fit' objects.", call. = FALSE)
  }
  if (!identical(aipw$lag_order, plr$lag_order)) {
    warning("Comparing fits with different lag_order; aligning on shared term names.", call. = FALSE)
  }

  aipw_coefs <- aipw$effect$coefficients
  plr_coefs <- plr$effect$coefficients
  shared <- intersect(names(aipw_coefs), names(plr_coefs))
  if (length(shared) == 0L) {
    stop("The two fits share no coefficient names.", call. = FALSE)
  }

  a <- as.numeric(aipw_coefs[shared])
  b <- as.numeric(plr_coefs[shared])
  scale_tol <- sqrt(.Machine$double.eps) * max(abs(c(a, b)), 1)
  data.frame(
    term = shared,
    aipw_estimate = a,
    plr_estimate = b,
    abs_diff = b - a,
    rel_diff = ifelse(abs(a) > scale_tol, (b - a) / a, NA_real_),
    sign_agree = (a > scale_tol & b > scale_tol) | (a < -scale_tol & b < -scale_tol) |
      (abs(a) <= scale_tol & abs(b) <= scale_tol),
    stringsAsFactors = FALSE
  )
}
