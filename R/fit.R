cdmc_resolve_fit_objective <- function(objective) {
  match.arg(objective, c("staged", "joint"))
}

cdmc_build_empty_effect_fit <- function(y_matrix, baseline_hat, lag_order, model = "none") {
  list(
    coefficients = numeric(0),
    fitted = matrix(0, nrow = nrow(y_matrix), ncol = ncol(y_matrix)),
    tau = y_matrix - baseline_hat,
    sample_mask = matrix(FALSE, nrow = nrow(y_matrix), ncol = ncol(y_matrix)),
    lag_order = lag_order,
    design_columns = character(0),
    model = model,
    design_info = NULL,
    weights = NULL,
    residuals = numeric(0)
  )
}

cdmc_build_joint_effect_design <- function(
  dose_matrix,
  fit_mask,
  lag_order = 0L,
  model = c("linear", "spline") ,
  df = 4L
) {
  model <- match.arg(model)
  lag_array <- cdmc_build_lagged_doses(dose_matrix, lag_order = lag_order)
  lag_array[is.na(lag_array)] <- 0

  sample_indices <- which(fit_mask, arr.ind = TRUE)
  history <- cdmc_build_sample_history(
    lag_array = lag_array,
    sample_indices = sample_indices,
    tau = numeric(nrow(sample_indices))
  )

  design_info <- cdmc_build_response_design(
    history = history,
    model = model,
    df = df
  )

  design_matrices <- vector("list", length(design_info$column_names))
  names(design_matrices) <- design_info$column_names

  if (identical(model, "linear")) {
    for (index in seq_len(dim(lag_array)[3L])) {
      lag_name <- dimnames(lag_array)[[3L]][index]
      design_matrices[[lag_name]] <- lag_array[, , index]
    }
  } else {
    for (index in seq_along(design_info$lag_names)) {
      lag_name <- design_info$lag_names[[index]]
      lag_values <- as.vector(lag_array[, , lag_name])
      centered <- cdmc_apply_spline_basis(lag_values, design_info$basis_spec[[lag_name]])
      current_names <- paste0(lag_name, "_ns", seq_len(ncol(centered)))

      for (column_index in seq_len(ncol(centered))) {
        design_matrices[[current_names[[column_index]]]] <- matrix(
          centered[, column_index],
          nrow = nrow(dose_matrix),
          ncol = ncol(dose_matrix)
        )
      }
    }
  }

  list(
    matrices = design_matrices[design_info$column_names],
    design_info = design_info,
    sample_mask = fit_mask
  )
}

cdmc_fit_model_components <- function(
  prepared,
  objective = c("staged", "joint"),
  fit_mask,
  weight_matrix = NULL,
  lambda,
  rank_max,
  lag_order = 0L,
  effect_model = c("linear", "spline", "none"),
  effect_df = 4L,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  verbose = FALSE,
  joint_design = NULL
) {
  objective <- match.arg(objective)
  effect_model <- match.arg(effect_model)

  if (identical(objective, "staged")) {
    baseline_fit <- cdmc_fit_baseline(
      y_matrix = prepared$y_matrix,
      x_matrices = prepared$x_matrices,
      mask = fit_mask,
      weight_matrix = weight_matrix,
      lambda = lambda,
      rank_max = rank_max,
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      verbose = verbose
    )

    effect_fit <- if (identical(effect_model, "linear")) {
      cdmc_fit_linear_effect(
        y_matrix = prepared$y_matrix,
        baseline_hat = baseline_fit$baseline_hat,
        dose_matrix = prepared$dose_matrix,
        lag_order = lag_order,
        weight_matrix = weight_matrix
      )
    } else if (identical(effect_model, "spline")) {
      cdmc_fit_spline_effect(
        y_matrix = prepared$y_matrix,
        baseline_hat = baseline_fit$baseline_hat,
        dose_matrix = prepared$dose_matrix,
        lag_order = lag_order,
        df = effect_df,
        weight_matrix = weight_matrix
      )
    } else {
      cdmc_build_empty_effect_fit(
        y_matrix = prepared$y_matrix,
        baseline_hat = baseline_fit$baseline_hat,
        lag_order = lag_order,
        model = effect_model
      )
    }

    baseline_fit$objective <- objective
    baseline_fit$optimization_mask <- fit_mask
    return(list(baseline = baseline_fit, effect = effect_fit))
  }

  if (identical(effect_model, "none")) {
    baseline_fit <- cdmc_fit_baseline(
      y_matrix = prepared$y_matrix,
      x_matrices = prepared$x_matrices,
      mask = fit_mask,
      weight_matrix = weight_matrix,
      lambda = lambda,
      rank_max = rank_max,
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      verbose = verbose
    )
    baseline_fit$objective <- objective
    baseline_fit$optimization_mask <- fit_mask

    return(list(
      baseline = baseline_fit,
      effect = cdmc_build_empty_effect_fit(
        y_matrix = prepared$y_matrix,
        baseline_hat = baseline_fit$baseline_hat,
        lag_order = lag_order,
        model = effect_model
      )
    ))
  }

  joint_design <- joint_design %||% cdmc_build_joint_effect_design(
    dose_matrix = prepared$dose_matrix,
    fit_mask = fit_mask,
    lag_order = lag_order,
    model = effect_model,
    df = effect_df
  )

  baseline_fit <- cdmc_fit_baseline(
    y_matrix = prepared$y_matrix,
    x_matrices = c(prepared$x_matrices, joint_design$matrices),
    mask = fit_mask,
    weight_matrix = weight_matrix,
    lambda = lambda,
    rank_max = rank_max,
    outer_maxit = outer_maxit,
    fe_maxit = fe_maxit,
    soft_maxit = soft_maxit,
    tol = tol,
    fe_tol = fe_tol,
    verbose = verbose
  )

  baseline_gamma_names <- names(prepared$x_matrices)
  treatment_gamma_names <- names(joint_design$matrices)
  all_gamma <- baseline_fit$nuisance$gamma
  baseline_gamma <- if (length(baseline_gamma_names) > 0L) all_gamma[baseline_gamma_names] else numeric(0)
  treatment_gamma <- if (length(treatment_gamma_names) > 0L) all_gamma[treatment_gamma_names] else numeric(0)

  baseline_nuisance_fitted <-
    matrix(baseline_fit$nuisance$alpha, nrow = prepared$n_units, ncol = prepared$n_times) +
    matrix(baseline_fit$nuisance$beta, nrow = prepared$n_units, ncol = prepared$n_times, byrow = TRUE) +
    cdmc_covariate_contribution(
      x_matrices = prepared$x_matrices,
      gamma = baseline_gamma,
      n_units = prepared$n_units,
      n_times = prepared$n_times
    )
  treatment_fitted <- cdmc_covariate_contribution(
    x_matrices = joint_design$matrices,
    gamma = treatment_gamma,
    n_units = prepared$n_units,
    n_times = prepared$n_times
  )

  baseline_fit$nuisance$gamma <- baseline_gamma
  baseline_fit$nuisance$fitted <- baseline_nuisance_fitted
  baseline_fit$baseline_hat <- baseline_nuisance_fitted + baseline_fit$low_rank
  baseline_fit$joint_treatment_coefficients <- treatment_gamma
  baseline_fit$joint_treatment_fitted <- treatment_fitted
  baseline_fit$joint_full_fitted <- baseline_fit$baseline_hat + treatment_fitted
  baseline_fit$objective <- objective
  baseline_fit$optimization_mask <- fit_mask

  effect_fit <- list(
    coefficients = treatment_gamma,
    fitted = treatment_fitted,
    tau = prepared$y_matrix - baseline_fit$baseline_hat,
    sample_mask = fit_mask,
    lag_order = lag_order,
    design_columns = treatment_gamma_names,
    model = effect_model,
    design_info = joint_design$design_info,
    weights = if (is.null(weight_matrix)) NULL else weight_matrix[fit_mask],
    residuals = (prepared$y_matrix - baseline_fit$baseline_hat - treatment_fitted)[fit_mask]
  )

  list(baseline = baseline_fit, effect = effect_fit)
}

cdmc_select_fit_lambda <- function(
  y_matrix,
  x_matrices,
  mask,
  weight_matrix = NULL,
  lambda = NULL,
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
  rank_max,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  verbose = FALSE
) {
  lambda_selection <- match.arg(lambda_selection)

  if (is.null(lambda)) {
    if (identical(lambda_selection, "cv")) {
      lambda_tuning <- cdmc_tune_lambda(
        y_matrix = y_matrix,
        x_matrices = x_matrices,
        mask = mask,
        weight_matrix = weight_matrix,
        rank_max = rank_max,
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
        outer_maxit = outer_maxit,
        fe_maxit = fe_maxit,
        soft_maxit = soft_maxit,
        tol = tol,
        fe_tol = fe_tol,
        verbose = verbose
      )
      lambda <- lambda_tuning$selected_lambda
    } else {
      reference_lambda0 <- cdmc_default_lambda(
        y_matrix = y_matrix,
        x_matrices = x_matrices,
        mask = mask,
        weight_matrix = weight_matrix,
        lambda_fraction = 1,
        fe_maxit = fe_maxit,
        fe_tol = fe_tol
      )
      lambda <- lambda_fraction * reference_lambda0
      lambda_tuning <- list(
        method = "heuristic",
        selected_lambda = lambda,
        reference_lambda0 = reference_lambda0,
        lambda_fraction = lambda_fraction
      )
    }
  } else {
    lambda_tuning <- list(
      method = "fixed",
      selected_lambda = lambda
    )
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("lambda must be a single nonnegative numeric value.", call. = FALSE)
  }

  list(
    lambda = as.numeric(lambda),
    lambda_tuning = lambda_tuning
  )
}

cdmc_fit <- function(
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
  lambda_selection = c("cv", "heuristic"),
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  workers = 1L,
  cv_top_k = NULL,
  cv_coarse_to_fine = FALSE,
  cv_coarse_nlambda = NULL,
  cv_warm_starts = FALSE,
  lambda_control = NULL,
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
) {
  effect_model <- match.arg(effect_model)
  lambda_selection <- match.arg(lambda_selection)
  objective <- cdmc_resolve_fit_objective(objective)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.numeric(lambda_fraction) || length(lambda_fraction) != 1L ||
      lambda_fraction <= 0 || lambda_fraction > 1) {
    stop("lambda_fraction must be a scalar in (0, 1].", call. = FALSE)
  }

  lag_order <- as.integer(lag_order)
  washout <- as.integer(washout)
  outer_maxit <- as.integer(outer_maxit)
  fe_maxit <- as.integer(fe_maxit)
  soft_maxit <- as.integer(soft_maxit)
  nlambda <- as.integer(nlambda)
  cv_rounds <- as.integer(cv_rounds)
  cv_block_size <- as.integer(cv_block_size)
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) || workers < 1 || workers != floor(workers)) {
    stop("workers must be a positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
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

  resolved_lambda_control <- cdmc_resolve_lambda_control(
    lambda_control,
    lambda = lambda,
    fraction = lambda_fraction,
    selection = lambda_selection,
    grid = lambda_grid,
    nlambda = nlambda,
    min_ratio = lambda_min_ratio,
    cv_rounds = cv_rounds,
    cv_block_size = cv_block_size,
    cv_top_k = cv_top_k,
    cv_coarse_to_fine = cv_coarse_to_fine,
    cv_coarse_nlambda = cv_coarse_nlambda,
    cv_warm_starts = cv_warm_starts
  )
  lambda <- resolved_lambda_control$lambda
  lambda_fraction <- resolved_lambda_control$fraction
  lambda_selection <- resolved_lambda_control$selection
  lambda_grid <- resolved_lambda_control$grid
  nlambda <- resolved_lambda_control$nlambda
  lambda_min_ratio <- resolved_lambda_control$min_ratio
  cv_rounds <- resolved_lambda_control$cv_rounds
  cv_block_size <- resolved_lambda_control$cv_block_size
  cv_top_k <- resolved_lambda_control$cv_top_k
  cv_coarse_to_fine <- resolved_lambda_control$cv_coarse_to_fine
  cv_coarse_nlambda <- resolved_lambda_control$cv_coarse_nlambda
  cv_warm_starts <- resolved_lambda_control$cv_warm_starts

  if (!is.numeric(effect_df) || length(effect_df) != 1L || !is.finite(effect_df) || effect_df < 1) {
    stop("effect_df must be a positive integer.", call. = FALSE)
  }
  effect_df <- as.integer(effect_df)

  prepared <- cdmc_prepare_panel(
    data = data,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    zero_tolerance = zero_tolerance
  )

  eligible_mask <- cdmc_build_eligible_mask(
    dose_matrix = prepared$dose_matrix,
    zero_tolerance = zero_tolerance,
    washout = washout
  )
  optimization_mask <- if (identical(objective, "joint")) prepared$observed_mask else eligible_mask
  weight_info <- cdmc_prepare_panel_weights(
    weights = weights,
    data = prepared$data,
    n_units = prepared$n_units,
    n_times = prepared$n_times,
    eligible_mask = optimization_mask
  )
  prepared$data <- weight_info$data
  fit_weight_matrix <- if (weight_info$supplied) weight_info$matrix else NULL
  if (!is.null(fit_weight_matrix)) {
    eligible_mask <- eligible_mask & fit_weight_matrix > 0
    optimization_mask <- optimization_mask & fit_weight_matrix > 0
  }
  if (identical(objective, "joint")) {
    validate_joint_support <- get("cdmc_validate_joint_support", mode = "function")
    validate_joint_support(optimization_mask)
  } else {
    cdmc_validate_support(optimization_mask)
  }

  if (is.null(rank_max)) {
    rank_max <- min(5L, min(prepared$n_units, prepared$n_times) - 1L)
  }

  if (!is.numeric(rank_max) || length(rank_max) != 1L || rank_max < 1L) {
    stop("rank_max must be a positive scalar.", call. = FALSE)
  }
  rank_max <- as.integer(rank_max)

  joint_design <- if (identical(objective, "joint") && !identical(effect_model, "none")) {
    cdmc_build_joint_effect_design(
      dose_matrix = prepared$dose_matrix,
      fit_mask = optimization_mask,
      lag_order = lag_order,
      model = effect_model,
      df = effect_df
    )
  } else {
    NULL
  }
  objective_x_matrices <- if (identical(objective, "joint") && !is.null(joint_design)) {
    c(prepared$x_matrices, joint_design$matrices)
  } else {
    prepared$x_matrices
  }

  lambda_selection_result <- cdmc_select_fit_lambda(
    y_matrix = prepared$y_matrix,
    x_matrices = objective_x_matrices,
    mask = optimization_mask,
    weight_matrix = fit_weight_matrix,
    lambda = lambda,
    lambda_fraction = lambda_fraction,
    lambda_selection = lambda_selection,
    lambda_grid = lambda_grid,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    cv_rounds = cv_rounds,
    cv_block_size = cv_block_size,
    cv_workers = workers,
    cv_top_k = cv_top_k,
    cv_coarse_to_fine = cv_coarse_to_fine,
    cv_coarse_nlambda = cv_coarse_nlambda,
    cv_warm_starts = cv_warm_starts,
    rank_max = rank_max,
    outer_maxit = outer_maxit,
    fe_maxit = fe_maxit,
    soft_maxit = soft_maxit,
    tol = tol,
    fe_tol = fe_tol,
    verbose = verbose
  )
  lambda <- lambda_selection_result$lambda
  lambda_tuning <- lambda_selection_result$lambda_tuning

  fit_components <- cdmc_fit_model_components(
    prepared = prepared,
    objective = objective,
    fit_mask = optimization_mask,
    weight_matrix = fit_weight_matrix,
    lambda = lambda,
    rank_max = rank_max,
    lag_order = lag_order,
    effect_model = effect_model,
    effect_df = effect_df,
    outer_maxit = outer_maxit,
    fe_maxit = fe_maxit,
    soft_maxit = soft_maxit,
    tol = tol,
    fe_tol = fe_tol,
    verbose = verbose,
    joint_design = joint_design
  )
  baseline_fit <- fit_components$baseline
  effect_fit <- fit_components$effect

  fitted_panel <- prepared$data
  fitted_panel$.cdmc_eligible_control <- cdmc_flatten_matrix(eligible_mask)
  fitted_panel$.cdmc_y0_hat <- cdmc_flatten_matrix(baseline_fit$baseline_hat)
  fitted_panel$.cdmc_tau <- cdmc_flatten_matrix(effect_fit$tau)
  fitted_panel$.cdmc_tau_model <- cdmc_flatten_matrix(effect_fit$fitted)
  fitted_panel$.cdmc_tau_linear <- cdmc_flatten_matrix(effect_fit$fitted)
  fitted_panel$.cdmc_tau_model[!fitted_panel$.cdmc_observed] <- NA_real_
  fitted_panel$.cdmc_tau_linear[!fitted_panel$.cdmc_observed] <- NA_real_

  result <- list(
    call = match.call(),
    data = fitted_panel,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    weights = if (weight_info$supplied) weight_info$column else NULL,
    n_units = prepared$n_units,
    n_times = prepared$n_times,
    unit_levels = prepared$unit_levels,
    time_levels = prepared$time_levels,
    observed_mask = prepared$observed_mask,
    y_matrix = prepared$y_matrix,
    dose_matrix = prepared$dose_matrix,
    weight_matrix = fit_weight_matrix,
    weight_supplied = weight_info$supplied,
    weight_scale = weight_info$scale,
    eligible_mask = eligible_mask,
    optimization_mask = optimization_mask,
    baseline = baseline_fit,
    effect = effect_fit,
    fit_control = list(
      weights = if (weight_info$supplied) weight_info$column else NULL,
      lambda = if (identical(lambda_tuning$method, "fixed")) lambda else NULL,
      rank_max = rank_max,
      lambda_fraction = lambda_fraction,
      lambda_selection = lambda_selection,
      lambda_grid = lambda_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      workers = workers,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = cv_coarse_to_fine,
      cv_coarse_nlambda = cv_coarse_nlambda,
      cv_warm_starts = cv_warm_starts,
      washout = washout,
      lag_order = lag_order,
      effect_model = effect_model,
      effect_df = effect_df,
      objective = objective,
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      zero_tolerance = zero_tolerance
    ),
    lambda_tuning = lambda_tuning,
    lambda = lambda,
    rank_max = rank_max,
    washout = washout,
    lag_order = lag_order,
    effect_model = effect_model,
    effect_df = effect_df,
    objective = objective,
    zero_tolerance = zero_tolerance
  )

  class(result) <- "cdmc_fit"
  result
}

cdmc_print_empirical_lambda_note <- function(lambda_method) {
  if (!identical(lambda_method, "heuristic")) {
    return(invisible(NULL))
  }

  cat(
    "  note: empirical workflows should prefer lambda_selection = \"cv\" when blocked cross-validation has enough zero-dose support.\n"
  )

  invisible(NULL)
}

print.cdmc_fit <- function(x, ...) {
  cat("causaldosemc fit\n")
  cat(sprintf("  units x times: %d x %d\n", x$n_units, x$n_times))
  cat(sprintf("  objective: %s\n", x$objective %||% "staged"))
  if (identical(x$objective %||% "staged", "joint")) {
    cat(sprintf("  optimization cells: %d\n", sum(x$optimization_mask %||% x$eligible_mask)))
  }
  cat(sprintf("  eligible zero-dose cells: %d\n", sum(x$eligible_mask)))
  cat(sprintf("  weighted fit: %s\n", if (isTRUE(x$weight_supplied)) "yes" else "no"))
  cat(sprintf("  lambda: %.6g\n", x$lambda))
  cat(sprintf("  lambda selection: %s\n", x$lambda_tuning$method))
  if (identical(x$lambda_tuning$method, "cv")) {
    cat(sprintf(
      "  mean cv score: %.6g\n",
      x$lambda_tuning$mean_scores[[x$lambda_tuning$selected_index]]
    ))
  }
  cdmc_print_empirical_lambda_note(x$lambda_tuning$method)
  cat(sprintf("  effective rank: %d\n", x$baseline$effective_rank))
  cat(sprintf("  baseline solver: %s\n", x$baseline$solver))
  cat(sprintf(
    "  baseline convergence: %s after %d outer iterations\n",
    if (x$baseline$converged) "converged" else "not converged",
    x$baseline$outer_iterations
  ))

  if (length(x$effect$coefficients) > 0L) {
    cat(sprintf("  %s effect coefficients:\n", x$effect_model))
    for (coefficient_name in names(x$effect$coefficients)) {
      cat(sprintf("    %s = %.6g\n", coefficient_name, x$effect$coefficients[[coefficient_name]]))
    }
  }

  invisible(x)
}

summary.cdmc_fit <- function(object, ...) {
  output <- list(
    n_units = object$n_units,
    n_times = object$n_times,
    objective = object$objective %||% "staged",
    optimization_cells = sum(object$optimization_mask %||% object$eligible_mask),
    eligible_controls = sum(object$eligible_mask),
    observed_controls = sum(cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance)),
    weighted_fit = isTRUE(object$weight_supplied),
    lambda = object$lambda,
    lambda_method = object$lambda_tuning$method,
    rank_max = object$rank_max,
    effective_rank = object$baseline$effective_rank,
    baseline_solver = object$baseline$solver,
    baseline_converged = object$baseline$converged,
    outer_iterations = object$baseline$outer_iterations,
    relative_change = object$baseline$relative_change,
    washout = object$washout,
    lag_order = object$lag_order,
    effect_model = object$effect_model,
    effect_df = object$effect_df,
    effect_coefficients = object$effect$coefficients
  )
  class(output) <- "summary.cdmc_fit"
  output
}

print.summary.cdmc_fit <- function(x, ...) {
  cat("Summary of causaldosemc fit\n")
  cat(sprintf("  panel size: %d x %d\n", x$n_units, x$n_times))
  cat(sprintf("  objective: %s\n", x$objective %||% "staged"))
  if (identical(x$objective %||% "staged", "joint")) {
    cat(sprintf("  optimization cells: %d\n", x$optimization_cells))
  }
  cat(sprintf("  observed zero-dose cells: %d\n", x$observed_controls))
  cat(sprintf("  eligible zero-dose cells: %d\n", x$eligible_controls))
  cat(sprintf("  lambda: %.6g\n", x$lambda))
  cat(sprintf("  lambda selection: %s\n", x$lambda_method))
  cdmc_print_empirical_lambda_note(x$lambda_method)
  cat(sprintf("  rank_max: %d\n", x$rank_max))
  cat(sprintf("  effective rank: %d\n", x$effective_rank))
  cat(sprintf(
    "  baseline convergence: %s after %d outer iterations\n",
    if (x$baseline_converged) "converged" else "not converged",
    x$outer_iterations
  ))
  cat(sprintf("  final relative change: %.6g\n", x$relative_change))
  cat(sprintf("  washout: %d\n", x$washout))
  cat(sprintf("  lag order: %d\n", x$lag_order))
  cat(sprintf("  effect model: %s\n", x$effect_model))
  if (identical(x$effect_model, "spline")) {
    cat(sprintf("  effect basis df: %d\n", x$effect_df))
  }

  if (length(x$effect_coefficients) > 0L) {
    cat(sprintf("  %s effect coefficients:\n", x$effect_model))
    for (coefficient_name in names(x$effect_coefficients)) {
      cat(sprintf("    %s = %.6g\n", coefficient_name, x$effect_coefficients[[coefficient_name]]))
    }
  }

  invisible(x)
}

predict.cdmc_fit <- function(object, newdata = NULL, type = c("y0", "tau", "tau_model", "tau_linear"), ...) {
  type <- match.arg(type)

  if (!is.null(newdata)) {
    stop(
      "v0.1 only supports predictions on the training panel. newdata is not implemented yet.",
      call. = FALSE
    )
  }

  values <- switch(
    type,
    y0 = object$baseline$baseline_hat,
    tau = object$effect$tau,
    tau_model = object$effect$fitted,
    tau_linear = object$effect$fitted
  )

  output <- object$data[, c(object$unit, object$time), drop = FALSE]
  output$estimate <- cdmc_flatten_matrix(values)
  if (type != "y0" && ".cdmc_observed" %in% names(object$data)) {
    output$estimate[!object$data$.cdmc_observed] <- NA_real_
  }
  output
}
