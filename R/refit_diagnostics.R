cdmc_resolve_refit_tuning_flag <- function(object, rerun_tuning = FALSE) {
  if (!is.logical(rerun_tuning) || length(rerun_tuning) != 1L || is.na(rerun_tuning)) {
    stop("rerun_tuning must be TRUE or FALSE.", call. = FALSE)
  }

  rerun_tuning <- isTRUE(rerun_tuning)
  if (rerun_tuning && identical(object$lambda_tuning$method, "fixed")) {
    stop(
      "rerun_tuning = TRUE requires a source cdmc_fit object with automatically selected lambda.",
      call. = FALSE
    )
  }

  rerun_tuning
}

cdmc_refit_baseline <- function(object, drop_mask, rerun_tuning = FALSE, verbose = FALSE) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  if (!is.matrix(drop_mask) || !identical(dim(drop_mask), dim(object$eligible_mask))) {
    stop("drop_mask must be a logical matrix matching the panel dimensions.", call. = FALSE)
  }

  rerun_tuning <- cdmc_resolve_refit_tuning_flag(object, rerun_tuning = rerun_tuning)
  build_bootstrap_fit_spec <- get("cdmc_build_bootstrap_fit_spec", mode = "function")
  refit_spec <- build_bootstrap_fit_spec(object, rerun_tuning = rerun_tuning)
  fit_control <- refit_spec$fit_spec
  objective <- object$objective %||% fit_control$objective %||% "staged"
  fit_mask <- object$optimization_mask %||% object$eligible_mask
  train_mask <- fit_mask & !drop_mask

  if (identical(objective, "joint")) {
    validate_joint_support <- get("cdmc_validate_joint_support", mode = "function")
    validate_joint_support(train_mask)
  } else {
    cdmc_validate_support(train_mask)
  }

  refit_data <- object$data
  if (".cdmc_observed" %in% names(refit_data)) {
    refit_data <- refit_data[!is.na(refit_data$.cdmc_observed) & refit_data$.cdmc_observed, , drop = FALSE]
  }

  prepared <- cdmc_prepare_panel(
    data = refit_data[, unique(c(object$outcome, object$dose, object$unit, object$time, object$covariates)), drop = FALSE],
    outcome = object$outcome,
    dose = object$dose,
    unit = object$unit,
    time = object$time,
    covariates = object$covariates,
    zero_tolerance = object$zero_tolerance
  )
  build_joint_effect_design <- get("cdmc_build_joint_effect_design", mode = "function")
  fit_model_components <- get("cdmc_fit_model_components", mode = "function")
  select_fit_lambda <- get("cdmc_select_fit_lambda", mode = "function")
  rank_max <- fit_control$rank_max %||% object$rank_max
  joint_design <- if (identical(objective, "joint") && !identical(object$effect_model, "none")) {
    build_joint_effect_design(
      dose_matrix = prepared$dose_matrix,
      fit_mask = train_mask,
      lag_order = object$lag_order,
      model = object$effect_model,
      df = object$effect_df %||% 4L
    )
  } else {
    NULL
  }

  objective_x_matrices <- if (identical(objective, "joint") && !is.null(joint_design)) {
    c(prepared$x_matrices, joint_design$matrices)
  } else {
    prepared$x_matrices
  }

  lambda_selection_result <- select_fit_lambda(
    y_matrix = prepared$y_matrix,
    x_matrices = objective_x_matrices,
    mask = train_mask,
    weight_matrix = object$weight_matrix,
    lambda = fit_control$lambda,
    lambda_fraction = fit_control$lambda_fraction %||% 0.25,
    lambda_selection = fit_control$lambda_selection %||% "heuristic",
    lambda_grid = fit_control$lambda_grid %||% NULL,
    nlambda = fit_control$nlambda %||% 5L,
    lambda_min_ratio = fit_control$lambda_min_ratio %||% 0.05,
    cv_rounds = fit_control$cv_rounds %||% 5L,
    cv_block_size = fit_control$cv_block_size %||% 2L,
    cv_workers = fit_control$cv_workers %||% 1L,
    rank_max = rank_max,
    outer_maxit = fit_control$outer_maxit,
    fe_maxit = fit_control$fe_maxit,
    soft_maxit = fit_control$soft_maxit,
    tol = fit_control$tol,
    fe_tol = fit_control$fe_tol,
    verbose = verbose
  )

  baseline_fit <- fit_model_components(
    prepared = prepared,
    objective = objective,
    fit_mask = train_mask,
    weight_matrix = object$weight_matrix,
    lambda = lambda_selection_result$lambda,
    rank_max = rank_max,
    lag_order = object$lag_order,
    effect_model = object$effect_model,
    effect_df = object$effect_df %||% 4L,
    outer_maxit = fit_control$outer_maxit,
    fe_maxit = fit_control$fe_maxit,
    soft_maxit = fit_control$soft_maxit,
    tol = fit_control$tol,
    fe_tol = fit_control$fe_tol,
    verbose = verbose,
    joint_design = joint_design
  )$baseline

  baseline_fit$lambda <- lambda_selection_result$lambda
  baseline_fit$lambda_tuning <- lambda_selection_result$lambda_tuning
  baseline_fit$rerun_tuning <- rerun_tuning
  baseline_fit
}

cdmc_trim_drop_mask <- function(mask, drop_mask) {
  trimmed_mask <- matrix(FALSE, nrow = nrow(mask), ncol = ncol(mask))
  row_counts <- rowSums(mask)
  col_counts <- colSums(mask)
  remaining_total <- sum(mask)
  target_indices <- which(drop_mask & mask, arr.ind = TRUE)

  for (index in seq_len(nrow(target_indices))) {
    row_index <- target_indices[index, 1L]
    col_index <- target_indices[index, 2L]

    if (row_counts[row_index] <= 1L) {
      next
    }
    if (col_counts[col_index] <= 1L) {
      next
    }
    if ((remaining_total - 1L) <= (nrow(mask) + ncol(mask))) {
      next
    }

    trimmed_mask[row_index, col_index] <- TRUE
    row_counts[row_index] <- row_counts[row_index] - 1L
    col_counts[col_index] <- col_counts[col_index] - 1L
    remaining_total <- remaining_total - 1L
  }

  trimmed_mask
}

cdmc_build_placebo_mask <- function(object, periods) {
  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods > 0L)) {
    stop("placebo periods must contain one or more nonpositive integers.", call. = FALSE)
  }

  n_units <- nrow(object$dose_matrix)
  n_times <- ncol(object$dose_matrix)
  target_mask <- matrix(FALSE, nrow = n_units, ncol = n_times)

  for (unit_index in seq_len(n_units)) {
    treated_times <- which(cdmc_active_dose_mask(object$dose_matrix[unit_index, ], zero_tolerance = object$zero_tolerance))
    if (length(treated_times) == 0L) {
      next
    }

    first_treated <- min(treated_times)
    candidate_times <- first_treated + periods
    candidate_times <- candidate_times[candidate_times >= 1L & candidate_times <= n_times]
    if (length(candidate_times) == 0L) {
      next
    }

    valid_times <- candidate_times[cdmc_zero_dose_mask(object$dose_matrix[unit_index, candidate_times], zero_tolerance = object$zero_tolerance)]
    if (length(valid_times) == 0L) {
      next
    }

    target_mask[unit_index, valid_times] <- TRUE
  }

  target_mask & object$eligible_mask
}

cdmc_summarize_refit_test <- function(object, target_mask, baseline_fit, label, periods) {
  if (!any(target_mask)) {
    stop(sprintf("No eligible cells matched the requested %s periods.", label), call. = FALSE)
  }

  tau_values <- object$y_matrix[target_mask] - baseline_fit$baseline_hat[target_mask]
  sample_size <- length(tau_values)
  mean_tau <- mean(tau_values)
  sd_tau <- if (sample_size > 1L) stats::sd(tau_values) else 0
  standard_error <- if (sample_size > 1L) sd_tau / sqrt(sample_size) else NA_real_
  t_statistic <- if (sample_size > 1L && standard_error > 0) mean_tau / standard_error else NA_real_
  p_value <- if (sample_size > 1L && is.finite(t_statistic)) {
    2 * stats::pt(-abs(t_statistic), df = sample_size - 1L)
  } else {
    NA_real_
  }

  cells <- object$data[, c(object$unit, object$time), drop = FALSE]
  cells$pseudo_tau <- cdmc_flatten_matrix(object$y_matrix - baseline_fit$baseline_hat)
  cells$eligible_control <- cdmc_flatten_matrix(object$eligible_mask)
  cells <- cells[cdmc_flatten_matrix(target_mask), , drop = FALSE]

  list(
    periods = periods,
    n = sample_size,
    mean_tau = mean_tau,
    sd_tau = sd_tau,
    standard_error = standard_error,
    t_statistic = t_statistic,
    p_value = p_value,
    cells = cells,
    baseline_fit = baseline_fit,
    rerun_tuning = baseline_fit$rerun_tuning %||% FALSE,
    refit_lambda = baseline_fit$lambda %||% NA_real_,
    refit_lambda_method = baseline_fit$lambda_tuning$method %||% NA_character_
  )
}

cdmc_placebo_test <- function(object, periods = -2:0, rerun_tuning = FALSE, verbose = FALSE) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  target_mask <- cdmc_trim_drop_mask(
    mask = object$eligible_mask,
    drop_mask = cdmc_build_placebo_mask(object, periods = periods)
  )
  baseline_fit <- cdmc_refit_baseline(
    object,
    drop_mask = target_mask,
    rerun_tuning = rerun_tuning,
    verbose = verbose
  )

  result <- cdmc_summarize_refit_test(
    object = object,
    target_mask = target_mask,
    baseline_fit = baseline_fit,
    label = "placebo",
    periods = sort(unique(as.integer(periods)))
  )
  result$call <- match.call()
  result$fit_object <- object

  class(result) <- "cdmc_placebo_test"
  result
}

cdmc_carryover_refit_test <- function(object, periods = 1L, rerun_tuning = FALSE, verbose = FALSE) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods < 1L)) {
    stop("carryover periods must contain one or more positive integers.", call. = FALSE)
  }

  exit_distance <- cdmc_exit_distances(
    dose_matrix = object$dose_matrix,
    zero_tolerance = object$zero_tolerance
  )
  target_mask <- cdmc_trim_drop_mask(
    mask = object$eligible_mask,
    drop_mask = object$eligible_mask & cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance) & exit_distance %in% periods
  )
  baseline_fit <- cdmc_refit_baseline(
    object,
    drop_mask = target_mask,
    rerun_tuning = rerun_tuning,
    verbose = verbose
  )

  result <- cdmc_summarize_refit_test(
    object = object,
    target_mask = target_mask,
    baseline_fit = baseline_fit,
    label = "carryover",
    periods = periods
  )
  result$call <- match.call()
  result$fit_object <- object

  class(result) <- "cdmc_carryover_refit_test"
  result
}

print.cdmc_placebo_test <- function(x, ...) {
  cat("causaldosemc placebo diagnostic\n")
  cat(sprintf("  periods tested: %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  matched placebo cells: %d\n", x$n))
  cat(sprintf("  rerun tuning in refit: %s\n", if (isTRUE(x$rerun_tuning)) "yes" else "no"))
  if (is.finite(x$refit_lambda)) {
    cat(sprintf("  refit lambda: %.6g\n", x$refit_lambda))
  }
  if (!is.na(x$refit_lambda_method)) {
    cat(sprintf("  refit lambda selection: %s\n", x$refit_lambda_method))
  }
  cat(sprintf("  mean pseudo effect: %.6g\n", x$mean_tau))
  if (is.finite(x$standard_error)) {
    cat(sprintf("  standard error: %.6g\n", x$standard_error))
  }
  if (is.finite(x$t_statistic)) {
    cat(sprintf("  t statistic: %.6g\n", x$t_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  p value: %.6g\n", x$p_value))
  }

  invisible(x)
}

print.cdmc_carryover_refit_test <- function(x, ...) {
  cat("causaldosemc refit carryover diagnostic\n")
  cat(sprintf("  periods tested: %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  matched exit cells: %d\n", x$n))
  cat(sprintf("  rerun tuning in refit: %s\n", if (isTRUE(x$rerun_tuning)) "yes" else "no"))
  if (is.finite(x$refit_lambda)) {
    cat(sprintf("  refit lambda: %.6g\n", x$refit_lambda))
  }
  if (!is.na(x$refit_lambda_method)) {
    cat(sprintf("  refit lambda selection: %s\n", x$refit_lambda_method))
  }
  cat(sprintf("  mean pseudo effect: %.6g\n", x$mean_tau))
  if (is.finite(x$standard_error)) {
    cat(sprintf("  standard error: %.6g\n", x$standard_error))
  }
  if (is.finite(x$t_statistic)) {
    cat(sprintf("  t statistic: %.6g\n", x$t_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  p value: %.6g\n", x$p_value))
  }

  invisible(x)
}
