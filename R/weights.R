cdmc_resolve_sample_weights <- function(sample_weights, data) {
  if (is.null(sample_weights)) {
    return(NULL)
  }

  resolved <- if (is.character(sample_weights) && length(sample_weights) == 1L) {
    if (!sample_weights %in% names(data)) {
      stop(sprintf("sample_weights column '%s' was not found in data.", sample_weights), call. = FALSE)
    }
    data[[sample_weights]]
  } else if (is.numeric(sample_weights)) {
    sample_weights
  } else {
    stop("sample_weights must be NULL, a column name, or a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(resolved) || any(!is.finite(resolved)) || any(resolved < 0)) {
    stop("Resolved sample weights must be finite and nonnegative.", call. = FALSE)
  }

  if (length(resolved) != nrow(data)) {
    stop("Resolved sample weights must have length equal to nrow(data).", call. = FALSE)
  }

  resolved
}

cdmc_panel_observed_rows <- function(data) {
  if (!is.data.frame(data) || !".cdmc_observed" %in% names(data)) {
    return(rep(TRUE, nrow(data)))
  }

  !is.na(data$.cdmc_observed) & data$.cdmc_observed
}

cdmc_expand_weight_vector <- function(resolved, data, argument_name) {
  observed_rows <- cdmc_panel_observed_rows(data)

  if (length(resolved) == nrow(data)) {
    return(as.numeric(resolved))
  }

  if (".cdmc_input_index" %in% names(data)) {
    input_index <- data$.cdmc_input_index
    valid_index <- !is.na(input_index)
    max_index <- if (any(valid_index)) max(input_index[valid_index]) else 0L

    if (length(resolved) == max_index) {
      expanded <- rep(NA_real_, nrow(data))
      expanded[valid_index] <- resolved[input_index[valid_index]]
      return(expanded)
    }
  }

  if (length(resolved) == sum(observed_rows)) {
    expanded <- rep(NA_real_, nrow(data))
    expanded[observed_rows] <- resolved
    return(expanded)
  }

  stop(
    sprintf("Resolved %s must have length equal to nrow(data) or to the number of observed rows.", argument_name),
    call. = FALSE
  )
}

cdmc_resolve_weight_vector <- function(
  weights,
  data,
  n_units = NULL,
  n_times = NULL,
  argument_name = "weights"
) {
  if (is.null(weights)) {
    return(NULL)
  }

  resolved <- if (inherits(weights, "cdmc_cbps_weights") || inherits(weights, "cdmc_entropy_balance_weights") || inherits(weights, "cdmc_kernel_balance_weights") || inherits(weights, "cdmc_adaptive_balance_weights")) {
    weights$weights
  } else if (is.list(weights) && !is.null(weights$weights)) {
    weights$weights
  } else if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) {
      stop(sprintf("%s column '%s' was not found in data.", argument_name, weights), call. = FALSE)
    }
    data[[weights]]
  } else if (is.matrix(weights)) {
    if (!is.null(n_units) && !is.null(n_times) && !identical(dim(weights), c(n_units, n_times))) {
      stop(
        sprintf("%s matrix must match the panel dimensions.", argument_name),
        call. = FALSE
      )
    }
    cdmc_flatten_matrix(weights)
  } else if (is.numeric(weights)) {
    weights
  } else {
    stop(
      sprintf(
        "%s must be NULL, a column name, a numeric vector, a numeric matrix, or an object with a numeric weights component.",
        argument_name
      ),
      call. = FALSE
    )
  }

  resolved <- cdmc_expand_weight_vector(
    resolved = resolved,
    data = data,
    argument_name = argument_name
  )

  observed_rows <- cdmc_panel_observed_rows(data)
  resolved[!observed_rows & is.na(resolved)] <- 0

  if (!is.numeric(resolved) || any(!is.finite(resolved[observed_rows])) || any(resolved[observed_rows] < 0)) {
    stop(sprintf("Resolved %s must be finite and nonnegative.", argument_name), call. = FALSE)
  }

  if (any(!is.finite(resolved[!observed_rows]))) {
    resolved[!observed_rows] <- 0
  }
  if (any(resolved[!observed_rows] < 0)) {
    stop(sprintf("Resolved %s must be nonnegative on padded rows.", argument_name), call. = FALSE)
  }

  as.numeric(resolved)
}

cdmc_prepare_panel_weights <- function(
  weights,
  data,
  n_units,
  n_times,
  eligible_mask = NULL,
  column_name = ".cdmc_weight"
) {
  resolved <- cdmc_resolve_weight_vector(
    weights = weights,
    data = data,
    n_units = n_units,
    n_times = n_times,
    argument_name = "weights"
  )

  if (is.null(resolved)) {
    return(list(
      data = data,
      supplied = FALSE,
      column = NULL,
      vector = rep(1, nrow(data)),
      matrix = matrix(1, nrow = n_units, ncol = n_times),
      scale = 1
    ))
  }

  positive_weights <- if (is.null(eligible_mask)) {
    resolved[resolved > 0]
  } else {
    resolved[cdmc_flatten_matrix(eligible_mask) & resolved > 0]
  }

  if (length(positive_weights) == 0L) {
    stop(
      "weights must leave at least one positive-weight eligible zero-dose observation.",
      call. = FALSE
    )
  }

  scale <- mean(positive_weights)
  normalized <- resolved / scale
  data[[column_name]] <- normalized

  list(
    data = data,
    supplied = TRUE,
    column = column_name,
    vector = normalized,
    matrix = matrix(normalized, nrow = n_units, ncol = n_times, byrow = TRUE),
    scale = scale
  )
}

cdmc_resolve_entropy_balance_degree <- function(degree) {
  if (!is.numeric(degree) || length(degree) != 1L || !is.finite(degree) || degree < 1) {
    stop("degree must be a positive integer.", call. = FALSE)
  }

  as.integer(degree)
}

cdmc_resolve_entropy_balance_iterations <- function(iterations) {
  if (!is.numeric(iterations) || length(iterations) != 1L || !is.finite(iterations) || iterations < 1) {
    stop("iterations must be a positive integer.", call. = FALSE)
  }

  as.integer(iterations)
}

cdmc_resolve_entropy_balance_reltol <- function(reltol) {
  if (!is.numeric(reltol) || length(reltol) != 1L || !is.finite(reltol) || reltol <= 0) {
    stop("reltol must be a single positive numeric value.", call. = FALSE)
  }

  as.numeric(reltol)
}

cdmc_resolve_kernel_balance_centers <- function(n_centers) {
  if (!is.numeric(n_centers) || length(n_centers) != 1L || !is.finite(n_centers) || n_centers < 1) {
    stop("n_centers must be a positive integer.", call. = FALSE)
  }

  as.integer(n_centers)
}

cdmc_resolve_kernel_balance_bandwidth <- function(bandwidth) {
  if (is.null(bandwidth)) {
    return(NULL)
  }

  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("bandwidth must be NULL or a single positive numeric value.", call. = FALSE)
  }

  as.numeric(bandwidth)
}

cdmc_weighted_mean <- function(x, weights) {
  sum(weights * x) / sum(weights)
}

cdmc_weighted_rms <- function(x, weights) {
  sqrt(sum(weights * x ^ 2) / sum(weights))
}

cdmc_center_weighted_matrix <- function(x, weights) {
  if (ncol(x) < 1L) {
    return(matrix(0, nrow = nrow(x), ncol = 0L))
  }

  centers <- vapply(
    seq_len(ncol(x)),
    function(index) cdmc_weighted_mean(x[, index], weights),
    numeric(1)
  )
  centered <- sweep(x, 2, centers, FUN = "-")
  rms <- vapply(
    seq_len(ncol(centered)),
    function(index) cdmc_weighted_rms(centered[, index], weights),
    numeric(1)
  )
  keep <- is.finite(rms) & rms > sqrt(.Machine$double.eps)
  centered[, keep, drop = FALSE]
}

cdmc_scale_weighted_matrix <- function(x, weights) {
  if (ncol(x) < 1L) {
    return(matrix(0, nrow = nrow(x), ncol = 0L))
  }

  rms <- vapply(
    seq_len(ncol(x)),
    function(index) cdmc_weighted_rms(x[, index], weights),
    numeric(1)
  )
  keep <- is.finite(rms) & rms > sqrt(.Machine$double.eps)
  if (!any(keep)) {
    return(matrix(0, nrow = nrow(x), ncol = 0L))
  }

  sweep(x[, keep, drop = FALSE], 2, rms[keep], FUN = "/")
}

cdmc_select_kernel_balance_centers <- function(kernel_input, n_centers) {
  if (nrow(kernel_input) < 1L || ncol(kernel_input) < 1L) {
    return(matrix(0, nrow = 0L, ncol = ncol(kernel_input)))
  }

  n_centers <- min(as.integer(n_centers), nrow(kernel_input))
  if (n_centers >= nrow(kernel_input)) {
    return(kernel_input)
  }

  scores <- if (ncol(kernel_input) == 1L) {
    kernel_input[, 1L]
  } else {
    tryCatch(
      stats::prcomp(kernel_input, center = FALSE, scale. = FALSE)$x[, 1L],
      error = function(...) rowMeans(kernel_input)
    )
  }
  ordered_rows <- order(scores)
  selected_positions <- unique(round(seq(1, length(ordered_rows), length.out = n_centers)))
  centers <- kernel_input[ordered_rows[selected_positions], , drop = FALSE]
  unique_centers <- unique.data.frame(as.data.frame(centers, stringsAsFactors = FALSE, check.names = FALSE))
  as.matrix(unique_centers)
}

cdmc_default_kernel_balance_bandwidth <- function(kernel_input, centers) {
  if (nrow(centers) > 1L) {
    center_sq <- rowSums(centers ^ 2)
    distance_sq <- outer(center_sq, center_sq, "+") - 2 * tcrossprod(centers)
    distance_sq[distance_sq < 0] <- 0
    distances <- sqrt(distance_sq[upper.tri(distance_sq)])
    distances <- distances[is.finite(distances) & distances > sqrt(.Machine$double.eps)]
    if (length(distances) > 0L) {
      return(as.numeric(stats::median(distances)))
    }
  }

  if (ncol(kernel_input) > 0L) {
    fallback <- sqrt(mean(rowSums(kernel_input ^ 2)))
    if (is.finite(fallback) && fallback > 0) {
      return(as.numeric(fallback))
    }
  }

  1
}

cdmc_compute_rbf_basis <- function(kernel_input, centers, bandwidth) {
  if (nrow(kernel_input) < 1L || ncol(kernel_input) < 1L || nrow(centers) < 1L) {
    return(matrix(0, nrow = nrow(kernel_input), ncol = 0L))
  }

  input_sq <- rowSums(kernel_input ^ 2)
  center_sq <- rowSums(centers ^ 2)
  distance_sq <- outer(input_sq, center_sq, "+") - 2 * tcrossprod(kernel_input, centers)
  distance_sq[distance_sq < 0] <- 0
  basis <- exp(-0.5 * distance_sq / (bandwidth ^ 2))
  colnames(basis) <- paste0("rbf_center", seq_len(ncol(basis)))
  basis
}

cdmc_prepare_entropy_balance_constraints <- function(
  observed_data,
  dose,
  covariates,
  time,
  time_effects,
  model,
  df,
  spline_covariates,
  degree,
  standardize,
  base_weights
) {
  formula <- cdmc_build_gps_formula(
    data = observed_data,
    dose = dose,
    covariates = covariates,
    time = if (isTRUE(time_effects)) time else dose,
    gps_time_effects = time_effects,
    gps_model = model,
    gps_df = df,
    gps_spline_covariates = spline_covariates
  )
  basis_terms <- stats::delete.response(stats::terms(formula))
  basis_matrix <- stats::model.matrix(basis_terms, data = observed_data)

  if (!"(Intercept)" %in% colnames(basis_matrix)) {
    basis_matrix <- cbind("(Intercept)" = 1, basis_matrix)
  }

  main_basis <- basis_matrix[, setdiff(colnames(basis_matrix), "(Intercept)"), drop = FALSE]
  if (ncol(main_basis) > 0L) {
    basis_centers <- vapply(
      seq_len(ncol(main_basis)),
      function(index) cdmc_weighted_mean(main_basis[, index], base_weights),
      numeric(1)
    )
    main_centered <- sweep(main_basis, 2, basis_centers, FUN = "-")
    main_rms <- vapply(
      seq_len(ncol(main_centered)),
      function(index) cdmc_weighted_rms(main_centered[, index], base_weights),
      numeric(1)
    )
    keep_main <- is.finite(main_rms) & main_rms > sqrt(.Machine$double.eps)
    main_centered <- main_centered[, keep_main, drop = FALSE]
  } else {
    main_centered <- matrix(0, nrow = nrow(observed_data), ncol = 0L)
  }

  interaction_basis_main <- if (length(covariates) > 0L) {
    interaction_formula <- cdmc_build_gps_formula(
      data = observed_data,
      dose = dose,
      covariates = covariates,
      time = dose,
      gps_time_effects = FALSE,
      gps_model = model,
      gps_df = df,
      gps_spline_covariates = spline_covariates
    )
    interaction_terms <- stats::delete.response(stats::terms(interaction_formula))
    interaction_matrix <- stats::model.matrix(interaction_terms, data = observed_data)
    interaction_matrix[, setdiff(colnames(interaction_matrix), "(Intercept)"), drop = FALSE]
  } else if (isTRUE(time_effects)) {
    main_basis
  } else {
    matrix(0, nrow = nrow(observed_data), ncol = 0L)
  }

  if (ncol(interaction_basis_main) > 0L) {
    interaction_centers <- vapply(
      seq_len(ncol(interaction_basis_main)),
      function(index) cdmc_weighted_mean(interaction_basis_main[, index], base_weights),
      numeric(1)
    )
    interaction_centered <- sweep(interaction_basis_main, 2, interaction_centers, FUN = "-")
    interaction_rms <- vapply(
      seq_len(ncol(interaction_centered)),
      function(index) cdmc_weighted_rms(interaction_centered[, index], base_weights),
      numeric(1)
    )
    keep_interaction <- is.finite(interaction_rms) & interaction_rms > sqrt(.Machine$double.eps)
    interaction_centered <- interaction_centered[, keep_interaction, drop = FALSE]
  } else {
    interaction_centered <- matrix(0, nrow = nrow(observed_data), ncol = 0L)
  }

  dose_centered <- observed_data[[dose]] - cdmc_weighted_mean(observed_data[[dose]], base_weights)
  interaction_basis <- cbind("(Intercept)" = 1, interaction_centered)
  interaction_constraints <- lapply(seq_len(degree), function(power) {
    interaction_block <- sweep(interaction_basis, 1, dose_centered ^ power, FUN = "*")
    colnames(interaction_block) <- paste0("dose_power", power, ":", colnames(interaction_basis))
    interaction_block
  })

  raw_constraints <- cbind(main_centered, do.call(cbind, interaction_constraints))
  if (ncol(raw_constraints) < 1L) {
    return(list(
      formula = formula,
      raw_constraints = matrix(0, nrow = nrow(observed_data), ncol = 0L),
      constraint_matrix = matrix(0, nrow = nrow(observed_data), ncol = 0L),
      feature_names = character(0)
    ))
  }

  raw_rms <- vapply(
    seq_len(ncol(raw_constraints)),
    function(index) cdmc_weighted_rms(raw_constraints[, index], base_weights),
    numeric(1)
  )
  keep_constraints <- is.finite(raw_rms) & raw_rms > sqrt(.Machine$double.eps)
  raw_constraints <- raw_constraints[, keep_constraints, drop = FALSE]
  raw_rms <- raw_rms[keep_constraints]

  constraint_matrix <- if (isTRUE(standardize) && ncol(raw_constraints) > 0L) {
    sweep(raw_constraints, 2, raw_rms, FUN = "/")
  } else {
    raw_constraints
  }

  if (ncol(constraint_matrix) > 0L) {
    qr_fit <- qr(constraint_matrix)
    if (qr_fit$rank < ncol(constraint_matrix)) {
      keep_rank <- sort(qr_fit$pivot[seq_len(qr_fit$rank)])
      constraint_matrix <- constraint_matrix[, keep_rank, drop = FALSE]
      raw_constraints <- raw_constraints[, keep_rank, drop = FALSE]
    }
  }

  list(
    formula = formula,
    raw_constraints = raw_constraints,
    constraint_matrix = constraint_matrix,
    feature_names = colnames(raw_constraints)
  )
}

cdmc_prepare_kernel_balance_constraints <- function(
  observed_data,
  dose,
  covariates,
  time,
  time_effects,
  model,
  df,
  spline_covariates,
  degree,
  n_centers,
  bandwidth,
  standardize,
  base_weights
) {
  formula <- cdmc_build_gps_formula(
    data = observed_data,
    dose = dose,
    covariates = covariates,
    time = if (isTRUE(time_effects)) time else dose,
    gps_time_effects = time_effects,
    gps_model = model,
    gps_df = df,
    gps_spline_covariates = spline_covariates
  )
  basis_terms <- stats::delete.response(stats::terms(formula))
  basis_matrix <- stats::model.matrix(basis_terms, data = observed_data)

  if (!"(Intercept)" %in% colnames(basis_matrix)) {
    basis_matrix <- cbind("(Intercept)" = 1, basis_matrix)
  }

  main_basis <- basis_matrix[, setdiff(colnames(basis_matrix), "(Intercept)"), drop = FALSE]
  main_centered <- cdmc_center_weighted_matrix(main_basis, base_weights)

  kernel_formula <- cdmc_build_gps_formula(
    data = observed_data,
    dose = dose,
    covariates = covariates,
    time = dose,
    gps_time_effects = FALSE,
    gps_model = model,
    gps_df = df,
    gps_spline_covariates = spline_covariates
  )
  kernel_terms <- stats::delete.response(stats::terms(kernel_formula))
  kernel_matrix <- stats::model.matrix(kernel_terms, data = observed_data)
  kernel_input <- kernel_matrix[, setdiff(colnames(kernel_matrix), "(Intercept)"), drop = FALSE]
  kernel_input <- cdmc_scale_weighted_matrix(kernel_input, base_weights)

  centers <- cdmc_select_kernel_balance_centers(kernel_input, n_centers = n_centers)
  applied_bandwidth <- cdmc_resolve_kernel_balance_bandwidth(bandwidth)
  if (is.null(applied_bandwidth)) {
    applied_bandwidth <- cdmc_default_kernel_balance_bandwidth(kernel_input, centers)
  }
  rbf_basis <- cdmc_compute_rbf_basis(kernel_input, centers, bandwidth = applied_bandwidth)
  rbf_centered <- cdmc_center_weighted_matrix(rbf_basis, base_weights)

  dose_centered <- observed_data[[dose]] - cdmc_weighted_mean(observed_data[[dose]], base_weights)
  interaction_basis <- cbind("(Intercept)" = 1, rbf_centered)
  interaction_constraints <- lapply(seq_len(degree), function(power) {
    interaction_block <- sweep(interaction_basis, 1, dose_centered ^ power, FUN = "*")
    colnames(interaction_block) <- paste0("dose_power", power, ":", colnames(interaction_basis))
    interaction_block
  })

  raw_constraints <- cbind(main_centered, do.call(cbind, interaction_constraints))
  if (ncol(raw_constraints) < 1L) {
    return(list(
      formula = formula,
      raw_constraints = matrix(0, nrow = nrow(observed_data), ncol = 0L),
      constraint_matrix = matrix(0, nrow = nrow(observed_data), ncol = 0L),
      feature_names = character(0),
      centers = centers,
      applied_bandwidth = applied_bandwidth
    ))
  }

  raw_rms <- vapply(
    seq_len(ncol(raw_constraints)),
    function(index) cdmc_weighted_rms(raw_constraints[, index], base_weights),
    numeric(1)
  )
  keep_constraints <- is.finite(raw_rms) & raw_rms > sqrt(.Machine$double.eps)
  raw_constraints <- raw_constraints[, keep_constraints, drop = FALSE]
  raw_rms <- raw_rms[keep_constraints]

  constraint_matrix <- if (isTRUE(standardize) && ncol(raw_constraints) > 0L) {
    sweep(raw_constraints, 2, raw_rms, FUN = "/")
  } else {
    raw_constraints
  }

  if (ncol(constraint_matrix) > 0L) {
    qr_fit <- qr(constraint_matrix)
    if (qr_fit$rank < ncol(constraint_matrix)) {
      keep_rank <- sort(qr_fit$pivot[seq_len(qr_fit$rank)])
      constraint_matrix <- constraint_matrix[, keep_rank, drop = FALSE]
      raw_constraints <- raw_constraints[, keep_rank, drop = FALSE]
    }
  }

  list(
    formula = formula,
    raw_constraints = raw_constraints,
    constraint_matrix = constraint_matrix,
    feature_names = colnames(raw_constraints),
    centers = centers,
    applied_bandwidth = applied_bandwidth
  )
}

cdmc_entropy_balance_weight_map <- function(constraint_matrix, base_weights, theta) {
  if (ncol(constraint_matrix) < 1L) {
    return(base_weights)
  }

  eta <- as.vector(constraint_matrix %*% theta)
  max_eta <- max(eta)
  weighted_exp <- base_weights * exp(eta - max_eta)
  # Equivalent to sum(base_weights) * weighted_exp / sum(weighted_exp); use a
  # scalar pre-scaling factor to avoid an extra n-vector temporary.
  weighted_exp * (sum(base_weights) / sum(weighted_exp))
}

cdmc_fit_entropy_balance_solver <- function(constraint_matrix, base_weights, iterations, reltol) {
  if (ncol(constraint_matrix) < 1L) {
    return(list(
      par = numeric(0),
      value = 0,
      convergence = 0L,
      counts = c(fn = 0L, gr = 0L),
      weights = base_weights,
      max_abs_standardized_balance = 0
    ))
  }

  objective <- function(theta) {
    eta <- as.vector(constraint_matrix %*% theta)
    max_eta <- max(eta)
    log(sum(base_weights * exp(eta - max_eta))) + max_eta
  }

  gradient_hessian <- function(theta) {
    weights <- cdmc_entropy_balance_weight_map(constraint_matrix, base_weights, theta)
    total_weight <- sum(weights)
    # Hot path: original code allocated four n x p temporaries (a sweep
    # and two copies of `centered * sqrt(w/W)`). Use the FWL identity
    #   H = C' diag(w/W) C - g g'
    # which collapses to one n x p allocation (`weighted_C`) plus a single
    # p x p crossprod and outer product.
    weighted_C <- constraint_matrix * weights
    gradient <- colSums(weighted_C) / total_weight
    hessian <- crossprod(constraint_matrix, weighted_C) / total_weight -
      tcrossprod(gradient)

    list(
      weights = weights,
      gradient = gradient,
      hessian = hessian
    )
  }

  theta <- rep(0, ncol(constraint_matrix))
  value <- objective(theta)
  fn_count <- 1L
  gr_count <- 0L
  tolerance <- max(1e-8, sqrt(reltol))
  converged <- FALSE

  for (iteration in seq_len(iterations)) {
    state <- gradient_hessian(theta)
    gr_count <- gr_count + 1L
    max_gradient <- if (length(state$gradient) > 0L) max(abs(state$gradient)) else 0
    if (max_gradient <= tolerance) {
      converged <- TRUE
      break
    }

    step <- tryCatch(
      qr.solve(state$hessian + diag(1e-8, ncol(state$hessian)), state$gradient),
      error = function(...) state$gradient
    )
    if (any(!is.finite(step))) {
      step <- state$gradient
    }

    step_scale <- 1
    improved <- FALSE
    while (step_scale > 1e-8) {
      candidate_theta <- theta - step_scale * step
      candidate_value <- objective(candidate_theta)
      fn_count <- fn_count + 1L
      if (is.finite(candidate_value) && candidate_value <= value + 1e-12) {
        theta <- candidate_theta
        value <- candidate_value
        improved <- TRUE
        break
      }
      step_scale <- step_scale / 2
    }

    if (!improved) {
      theta <- theta - 0.1 * state$gradient
      value <- objective(theta)
      fn_count <- fn_count + 1L
    }
  }

  final_state <- gradient_hessian(theta)
  gr_count <- gr_count + 1L

  list(
    par = theta,
    value = value,
    convergence = if (converged || max(abs(final_state$gradient)) <= tolerance) 0L else 1L,
    counts = c(fn = fn_count, gr = gr_count),
    weights = final_state$weights,
    max_abs_standardized_balance = if (length(final_state$gradient) > 0L) max(abs(final_state$gradient)) else 0
  )
}

cdmc_entropy_balance_summary <- function(raw_constraints, base_weights, weights) {
  if (ncol(raw_constraints) < 1L) {
    return(data.frame())
  }

  base_total <- sum(base_weights)
  weight_total <- sum(weights)
  data.frame(
    feature = colnames(raw_constraints),
    base_moment = colSums(raw_constraints * base_weights) / base_total,
    weighted_moment = colSums(raw_constraints * weights) / weight_total,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cdmc_entropy_balance_weights <- function(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  degree = 1L,
  standardize = TRUE,
  iterations = 1000L,
  reltol = 1e-8,
  sample_weights = NULL,
  max_weight = "adaptive"
) {
  data <- as.data.frame(data)
  model <- match.arg(model)

  if (!is.character(dose) || length(dose) != 1L || !dose %in% names(data)) {
    stop("dose must be a single column name present in data.", call. = FALSE)
  }
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }
  covariates <- covariates %||% character(0)
  if (!is.logical(time_effects) || length(time_effects) != 1L || is.na(time_effects)) {
    stop("time_effects must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(time_effects) && (!is.character(time) || length(time) != 1L || !time %in% names(data))) {
    stop("time must be a single column name present in data when time_effects = TRUE.", call. = FALSE)
  }
  if (length(covariates) == 0L && !isTRUE(time_effects)) {
    stop("At least one balancing covariate or time_effects = TRUE is required for entropy balancing.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("df must be a positive integer.", call. = FALSE)
  }

  df <- as.integer(df)
  degree <- cdmc_resolve_entropy_balance_degree(degree)
  iterations <- cdmc_resolve_entropy_balance_iterations(iterations)
  reltol <- cdmc_resolve_entropy_balance_reltol(reltol)
  max_weight <- cdmc_resolve_max_weight_spec(max_weight)

  observed_rows <- cdmc_panel_observed_rows(data)
  if (!any(observed_rows)) {
    stop("Entropy balancing requires at least one observed row.", call. = FALSE)
  }

  observed_data <- data[observed_rows, , drop = FALSE]
  missing_covariates <- setdiff(covariates, names(observed_data))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("Missing covariate columns: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  spline_covariates <- cdmc_resolve_gps_spline_covariates(covariates, spline_covariates)
  sample_weights_full <- cdmc_resolve_sample_weights(sample_weights, data)
  base_weights <- if (is.null(sample_weights_full)) rep(1, nrow(observed_data)) else sample_weights_full[observed_rows]
  if (sum(base_weights) <= sqrt(.Machine$double.eps)) {
    stop("Entropy balancing requires positive total sample weight on observed rows.", call. = FALSE)
  }

  prepared <- cdmc_prepare_entropy_balance_constraints(
    observed_data = observed_data,
    dose = dose,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    degree = degree,
    standardize = standardize,
    base_weights = base_weights
  )
  solver <- cdmc_fit_entropy_balance_solver(
    constraint_matrix = prepared$constraint_matrix,
    base_weights = base_weights,
    iterations = iterations,
    reltol = reltol
  )

  raw_weights_observed <- solver$weights
  applied_max_weight <- cdmc_resolve_internal_max_weight(raw_weights_observed, max_weight = max_weight)
  clipped_weights_observed <- cdmc_cap_internal_weights(raw_weights_observed, max_weight = applied_max_weight)
  balance_summary <- cdmc_entropy_balance_summary(
    raw_constraints = prepared$raw_constraints,
    base_weights = base_weights,
    weights = clipped_weights_observed
  )
  max_abs_balance <- if (nrow(balance_summary) > 0L) max(abs(balance_summary$weighted_moment)) else 0
  if (!(solver$convergence == 0L || solver$max_abs_standardized_balance <= max(1e-6, sqrt(reltol)))) {
    stop(
      sprintf(
        "Entropy balancing did not converge (code %d, max standardized imbalance %.6g).",
        solver$convergence,
        solver$max_abs_standardized_balance
      ),
      call. = FALSE
    )
  }

  expanded_weights <- numeric(nrow(data))
  expanded_raw_weights <- numeric(nrow(data))
  expanded_weights[observed_rows] <- clipped_weights_observed
  expanded_raw_weights[observed_rows] <- raw_weights_observed

  result <- list(
    call = match.call(),
    weight_method = "entropy_balance",
    formula = prepared$formula,
    weights = expanded_weights,
    raw_weights = expanded_raw_weights,
    gps_model = model,
    gps_df = df,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    degree = degree,
    standardize = standardize,
    iterations = iterations,
    reltol = reltol,
    requested_max_weight = max_weight,
    applied_max_weight = applied_max_weight,
    sample_weights = sample_weights_full,
    converged = solver$convergence == 0L || solver$max_abs_standardized_balance <= max(1e-6, sqrt(reltol)),
    optimization = solver,
    balance_summary = balance_summary,
    max_abs_balance = max_abs_balance,
    max_abs_standardized_balance = solver$max_abs_standardized_balance
  )

  class(result) <- "cdmc_entropy_balance_weights"
  result
}

print.cdmc_entropy_balance_weights <- function(x, ...) {
  cat("causaldosemc entropy-balance weights\n")
  cat(sprintf("  formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  observations: %d\n", length(x$weights)))
  cat(sprintf("  model: %s\n", x$model))
  cat(sprintf("  degree: %d\n", x$degree))
  cat(sprintf("  standardize constraints: %s\n", if (isTRUE(x$standardize)) "yes" else "no"))
  if (identical(x$requested_max_weight, "adaptive")) {
    if (is.null(x$applied_max_weight)) {
      cat("  weight cap: adaptive\n")
    } else {
      cat(sprintf("  weight cap: adaptive (%.6g)\n", x$applied_max_weight))
    }
  } else if (is.null(x$requested_max_weight)) {
    cat("  weight cap: none\n")
  } else {
    cat(sprintf("  weight cap: %.6g\n", x$requested_max_weight))
  }
  cat(sprintf("  max abs balance: %.6g\n", x$max_abs_balance))
  invisible(x)
}

cdmc_kernel_balance_weights <- function(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  degree = 1L,
  n_centers = 25L,
  bandwidth = NULL,
  standardize = TRUE,
  iterations = 1000L,
  reltol = 1e-8,
  sample_weights = NULL,
  max_weight = "adaptive"
) {
  data <- as.data.frame(data)
  model <- match.arg(model)

  if (!is.character(dose) || length(dose) != 1L || !dose %in% names(data)) {
    stop("dose must be a single column name present in data.", call. = FALSE)
  }
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }
  covariates <- covariates %||% character(0)
  if (!is.logical(time_effects) || length(time_effects) != 1L || is.na(time_effects)) {
    stop("time_effects must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(time_effects) && (!is.character(time) || length(time) != 1L || !time %in% names(data))) {
    stop("time must be a single column name present in data when time_effects = TRUE.", call. = FALSE)
  }
  if (length(covariates) == 0L && !isTRUE(time_effects)) {
    stop("At least one balancing covariate or time_effects = TRUE is required for kernel balancing.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("df must be a positive integer.", call. = FALSE)
  }
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize)) {
    stop("standardize must be TRUE or FALSE.", call. = FALSE)
  }

  df <- as.integer(df)
  degree <- cdmc_resolve_entropy_balance_degree(degree)
  n_centers <- cdmc_resolve_kernel_balance_centers(n_centers)
  bandwidth <- cdmc_resolve_kernel_balance_bandwidth(bandwidth)
  iterations <- cdmc_resolve_entropy_balance_iterations(iterations)
  reltol <- cdmc_resolve_entropy_balance_reltol(reltol)
  max_weight <- cdmc_resolve_max_weight_spec(max_weight)

  observed_rows <- cdmc_panel_observed_rows(data)
  if (!any(observed_rows)) {
    stop("Kernel balancing requires at least one observed row.", call. = FALSE)
  }

  observed_data <- data[observed_rows, , drop = FALSE]
  missing_covariates <- setdiff(covariates, names(observed_data))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("Missing covariate columns: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  spline_covariates <- cdmc_resolve_gps_spline_covariates(covariates, spline_covariates)
  sample_weights_full <- cdmc_resolve_sample_weights(sample_weights, data)
  base_weights <- if (is.null(sample_weights_full)) rep(1, nrow(observed_data)) else sample_weights_full[observed_rows]
  if (sum(base_weights) <= sqrt(.Machine$double.eps)) {
    stop("Kernel balancing requires positive total sample weight on observed rows.", call. = FALSE)
  }

  fit_kernel_balance <- function(current_bandwidth, current_n_centers) {
    prepared_fit <- cdmc_prepare_kernel_balance_constraints(
      observed_data = observed_data,
      dose = dose,
      covariates = covariates,
      time = time,
      time_effects = time_effects,
      model = model,
      df = df,
      spline_covariates = spline_covariates,
      degree = degree,
      n_centers = current_n_centers,
      bandwidth = current_bandwidth,
      standardize = standardize,
      base_weights = base_weights
    )
    solver_fit <- cdmc_fit_entropy_balance_solver(
      constraint_matrix = prepared_fit$constraint_matrix,
      base_weights = base_weights,
      iterations = iterations,
      reltol = reltol
    )

    list(prepared = prepared_fit, solver = solver_fit)
  }

  kernel_fit <- fit_kernel_balance(bandwidth, n_centers)
  prepared <- kernel_fit$prepared
  solver <- kernel_fit$solver
  balance_tolerance <- max(5e-2, sqrt(reltol))

  if (!(solver$convergence == 0L || solver$max_abs_standardized_balance <= balance_tolerance)) {
    center_grid <- unique(pmax(1L, c(n_centers, floor(n_centers / 2), floor(n_centers / 4), 1L)))
    bandwidth_grid <- if (is.null(bandwidth)) {
      unique(prepared$applied_bandwidth * c(1, 2 ^ seq_len(6)))
    } else {
      prepared$applied_bandwidth
    }
    bandwidth_grid <- bandwidth_grid[is.finite(bandwidth_grid) & bandwidth_grid > 0]

    best_score <- solver$max_abs_standardized_balance
    for (candidate_n_centers in center_grid) {
      for (candidate_bandwidth in bandwidth_grid) {
        if (candidate_n_centers == n_centers && candidate_bandwidth == prepared$applied_bandwidth) {
          next
        }
        candidate_fit <- fit_kernel_balance(candidate_bandwidth, candidate_n_centers)
        candidate_score <- candidate_fit$solver$max_abs_standardized_balance
        if (is.finite(candidate_score) && (!is.finite(best_score) || candidate_score < best_score)) {
          prepared <- candidate_fit$prepared
          solver <- candidate_fit$solver
          best_score <- candidate_score
        }
        if (solver$convergence == 0L || solver$max_abs_standardized_balance <= balance_tolerance) {
          break
        }
      }
      if (solver$convergence == 0L || solver$max_abs_standardized_balance <= balance_tolerance) {
        break
      }
    }
  }

  raw_weights_observed <- solver$weights
  applied_max_weight <- cdmc_resolve_internal_max_weight(raw_weights_observed, max_weight = max_weight)
  clipped_weights_observed <- cdmc_cap_internal_weights(raw_weights_observed, max_weight = applied_max_weight)
  balance_summary <- cdmc_entropy_balance_summary(
    raw_constraints = prepared$raw_constraints,
    base_weights = base_weights,
    weights = clipped_weights_observed
  )
  max_abs_balance <- if (nrow(balance_summary) > 0L) max(abs(balance_summary$weighted_moment)) else 0
  if (!(solver$convergence == 0L || solver$max_abs_standardized_balance <= balance_tolerance)) {
    stop(
      sprintf(
        "Kernel balancing did not converge (code %d, max standardized imbalance %.6g).",
        solver$convergence,
        solver$max_abs_standardized_balance
      ),
      call. = FALSE
    )
  }

  expanded_weights <- numeric(nrow(data))
  expanded_raw_weights <- numeric(nrow(data))
  expanded_weights[observed_rows] <- clipped_weights_observed
  expanded_raw_weights[observed_rows] <- raw_weights_observed

  result <- list(
    call = match.call(),
    weight_method = "kernel_balance",
    formula = prepared$formula,
    weights = expanded_weights,
    raw_weights = expanded_raw_weights,
    gps_model = model,
    gps_df = df,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    degree = degree,
    n_centers = n_centers,
    actual_n_centers = nrow(prepared$centers),
    bandwidth = bandwidth,
    applied_bandwidth = prepared$applied_bandwidth,
    standardize = standardize,
    iterations = iterations,
    reltol = reltol,
    balance_tolerance = balance_tolerance,
    requested_max_weight = max_weight,
    applied_max_weight = applied_max_weight,
    sample_weights = sample_weights_full,
    converged = solver$convergence == 0L || solver$max_abs_standardized_balance <= balance_tolerance,
    optimization = solver,
    balance_summary = balance_summary,
    max_abs_balance = max_abs_balance,
    max_abs_standardized_balance = solver$max_abs_standardized_balance
  )

  class(result) <- "cdmc_kernel_balance_weights"
  result
}

print.cdmc_kernel_balance_weights <- function(x, ...) {
  cat("causaldosemc kernel-balance weights\n")
  cat(sprintf("  formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  observations: %d\n", length(x$weights)))
  cat(sprintf("  model: %s\n", x$model))
  cat(sprintf("  degree: %d\n", x$degree))
  cat(sprintf("  centers: %d\n", x$actual_n_centers))
  cat(sprintf("  bandwidth: %.6g\n", x$applied_bandwidth))
  cat(sprintf("  standardize constraints: %s\n", if (isTRUE(x$standardize)) "yes" else "no"))
  if (identical(x$requested_max_weight, "adaptive")) {
    if (is.null(x$applied_max_weight)) {
      cat("  weight cap: adaptive\n")
    } else {
      cat(sprintf("  weight cap: adaptive (%.6g)\n", x$applied_max_weight))
    }
  } else if (is.null(x$requested_max_weight)) {
    cat("  weight cap: none\n")
  } else {
    cat(sprintf("  weight cap: %.6g\n", x$requested_max_weight))
  }
  cat(sprintf("  max abs balance: %.6g\n", x$max_abs_balance))
  invisible(x)
}

cdmc_resolve_adaptive_balance_methods <- function(methods = NULL) {
  allowed_methods <- c("cbps", "entropy_balance", "kernel_balance")

  if (is.null(methods)) {
    resolved <- c("entropy_balance", "kernel_balance")
    if (requireNamespace("CBPS", quietly = TRUE)) {
      resolved <- c("cbps", resolved)
    }
    return(resolved)
  }

  if (!is.character(methods) || length(methods) < 1L) {
    stop("methods must be NULL or a nonempty character vector.", call. = FALSE)
  }

  resolved <- unique(as.character(methods))
  unknown_methods <- setdiff(resolved, allowed_methods)
  if (length(unknown_methods) > 0L) {
    stop(
      sprintf(
        "Unknown adaptive balance methods: %s.",
        paste(unknown_methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if ("cbps" %in% resolved && !requireNamespace("CBPS", quietly = TRUE)) {
    stop("The CBPS package must be installed when methods includes 'cbps'.", call. = FALSE)
  }

  resolved
}

cdmc_prepare_cbps_internal_fit <- function(
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
  sample_weights = NULL,
  max_weight = "adaptive"
) {
  cbps_fit <- cdmc_cbps_weights(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    time_effects = gps_time_effects,
    model = gps_model,
    df = gps_df,
    spline_covariates = gps_spline_covariates,
    standardize = cbps_standardize,
    method = cbps_method,
    iterations = cbps_iterations,
    twostep = cbps_twostep,
    sample_weights = sample_weights
  )

  cbps_fit$raw_weights <- cbps_fit$weights
  cbps_fit$requested_max_weight <- cdmc_resolve_max_weight_spec(max_weight)
  cbps_fit$applied_max_weight <- cdmc_resolve_internal_max_weight(cbps_fit$raw_weights, max_weight = cbps_fit$requested_max_weight)
  cbps_fit$weights <- cdmc_cap_internal_weights(cbps_fit$raw_weights, max_weight = cbps_fit$applied_max_weight)
  cbps_fit$weight_method <- "cbps"
  cbps_fit$gps_time_effects <- gps_time_effects
  cbps_fit$gps_model <- gps_model
  cbps_fit$gps_df <- gps_df
  cbps_fit$gps_spline_covariates <- gps_spline_covariates
  cbps_fit$gps_bandwidth <- NULL
  cbps_fit
}

cdmc_prepare_adaptive_balance_score_basis <- function(
  methods,
  observed_data,
  dose,
  covariates,
  time,
  time_effects,
  model,
  df,
  spline_covariates,
  entropy_balance_degree,
  kernel_balance_degree,
  kernel_balance_centers,
  kernel_balance_bandwidth,
  base_weights
) {
  if ("kernel_balance" %in% methods) {
    return(cdmc_prepare_kernel_balance_constraints(
      observed_data = observed_data,
      dose = dose,
      covariates = covariates,
      time = time,
      time_effects = time_effects,
      model = model,
      df = df,
      spline_covariates = spline_covariates,
      degree = kernel_balance_degree,
      n_centers = kernel_balance_centers,
      bandwidth = kernel_balance_bandwidth,
      standardize = TRUE,
      base_weights = base_weights
    ))
  }

  cdmc_prepare_entropy_balance_constraints(
    observed_data = observed_data,
    dose = dose,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    degree = entropy_balance_degree,
    standardize = TRUE,
    base_weights = base_weights
  )
}

cdmc_score_balance_weights <- function(weights, constraint_matrix) {
  weights <- as.numeric(weights)
  valid <- is.finite(weights) & weights >= 0
  if (!any(valid)) {
    return(Inf)
  }

  weights <- weights[valid]
  if (ncol(constraint_matrix) < 1L) {
    return(0)
  }
  score_matrix <- constraint_matrix[valid, , drop = FALSE]
  total_weight <- sum(weights)
  if (!is.finite(total_weight) || total_weight <= 0) {
    return(Inf)
  }

  imbalance <- colSums(score_matrix * weights) / total_weight
  if (length(imbalance) < 1L) {
    return(0)
  }

  max(abs(imbalance))
}

cdmc_adaptive_balance_score_row <- function(method, fit, observed_rows, score_basis) {
  observed_weights <- fit$weights[observed_rows]
  ess <- cdmc_effective_sample_size(observed_weights)
  data.frame(
    method = method,
    fit_success = TRUE,
    fit_error = NA_character_,
    max_abs_standardized_balance = cdmc_score_balance_weights(
      weights = observed_weights,
      constraint_matrix = score_basis$constraint_matrix
    ),
    ess = ess,
    ess_fraction = if (is.finite(ess)) ess / sum(observed_rows) else NA_real_,
    applied_max_weight = fit$applied_max_weight %||% NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

cdmc_select_adaptive_balance_index <- function(candidate_scores) {
  valid_rows <- which(!is.na(candidate_scores$fit_success) & candidate_scores$fit_success)
  if (length(valid_rows) < 1L) {
    return(integer(0))
  }

  balance_rank <- candidate_scores$max_abs_standardized_balance[valid_rows]
  balance_rank[!is.finite(balance_rank)] <- Inf
  ess_rank <- candidate_scores$ess_fraction[valid_rows]
  ess_rank[!is.finite(ess_rank)] <- -Inf
  max_weight_rank <- candidate_scores$applied_max_weight[valid_rows]
  max_weight_rank[!is.finite(max_weight_rank)] <- Inf

  valid_rows[order(balance_rank, -ess_rank, max_weight_rank, candidate_scores$method[valid_rows])][1L]
}

cdmc_fit_adaptive_balance_candidate <- function(
  method,
  data,
  dose,
  covariates,
  time,
  time_effects,
  model,
  df,
  spline_covariates,
  cbps_standardize,
  cbps_method,
  cbps_iterations,
  cbps_twostep,
  entropy_balance_degree,
  entropy_balance_standardize,
  entropy_balance_iterations,
  entropy_balance_reltol,
  kernel_balance_degree,
  kernel_balance_centers,
  kernel_balance_bandwidth,
  kernel_balance_standardize,
  kernel_balance_iterations,
  kernel_balance_reltol,
  sample_weights,
  max_weight
) {
  switch(
    method,
    cbps = cdmc_prepare_cbps_internal_fit(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = time_effects,
      gps_model = model,
      gps_df = df,
      gps_spline_covariates = spline_covariates,
      cbps_standardize = cbps_standardize,
      cbps_method = cbps_method,
      cbps_iterations = cbps_iterations,
      cbps_twostep = cbps_twostep,
      sample_weights = sample_weights,
      max_weight = max_weight
    ),
    entropy_balance = cdmc_entropy_balance_weights(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      time_effects = time_effects,
      model = model,
      df = df,
      spline_covariates = spline_covariates,
      degree = entropy_balance_degree,
      standardize = entropy_balance_standardize,
      iterations = entropy_balance_iterations,
      reltol = entropy_balance_reltol,
      sample_weights = sample_weights,
      max_weight = max_weight
    ),
    kernel_balance = cdmc_kernel_balance_weights(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      time_effects = time_effects,
      model = model,
      df = df,
      spline_covariates = spline_covariates,
      degree = kernel_balance_degree,
      n_centers = kernel_balance_centers,
      bandwidth = kernel_balance_bandwidth,
      standardize = kernel_balance_standardize,
      iterations = kernel_balance_iterations,
      reltol = kernel_balance_reltol,
      sample_weights = sample_weights,
      max_weight = max_weight
    ),
    stop(sprintf("Unknown adaptive balance candidate method '%s'.", method), call. = FALSE)
  )
}

cdmc_adaptive_balance_weights <- function(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  methods = NULL,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
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
  sample_weights = NULL,
  max_weight = "adaptive"
) {
  data <- as.data.frame(data)
  model <- match.arg(model)
  cbps_method <- match.arg(cbps_method)

  if (!is.character(dose) || length(dose) != 1L || !dose %in% names(data)) {
    stop("dose must be a single column name present in data.", call. = FALSE)
  }
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }
  covariates <- covariates %||% character(0)
  if (!is.logical(time_effects) || length(time_effects) != 1L || is.na(time_effects)) {
    stop("time_effects must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(time_effects) && (!is.character(time) || length(time) != 1L || !time %in% names(data))) {
    stop("time must be a single column name present in data when time_effects = TRUE.", call. = FALSE)
  }
  if (length(covariates) == 0L && !isTRUE(time_effects)) {
    stop("At least one balancing covariate or time_effects = TRUE is required for adaptive balancing.", call. = FALSE)
  }
  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("df must be a positive integer.", call. = FALSE)
  }
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

  methods <- cdmc_resolve_adaptive_balance_methods(methods)
  df <- as.integer(df)
  spline_covariates <- cdmc_resolve_gps_spline_covariates(covariates, spline_covariates)
  entropy_balance_degree <- cdmc_resolve_entropy_balance_degree(entropy_balance_degree)
  entropy_balance_iterations <- cdmc_resolve_entropy_balance_iterations(entropy_balance_iterations)
  entropy_balance_reltol <- cdmc_resolve_entropy_balance_reltol(entropy_balance_reltol)
  kernel_balance_degree <- cdmc_resolve_entropy_balance_degree(kernel_balance_degree)
  kernel_balance_centers <- cdmc_resolve_kernel_balance_centers(kernel_balance_centers)
  kernel_balance_bandwidth <- cdmc_resolve_kernel_balance_bandwidth(kernel_balance_bandwidth)
  kernel_balance_iterations <- cdmc_resolve_entropy_balance_iterations(kernel_balance_iterations)
  kernel_balance_reltol <- cdmc_resolve_entropy_balance_reltol(kernel_balance_reltol)
  cbps_iterations <- cdmc_resolve_cbps_iterations(cbps_iterations)
  max_weight <- cdmc_resolve_max_weight_spec(max_weight)

  observed_rows <- cdmc_panel_observed_rows(data)
  if (!any(observed_rows)) {
    stop("Adaptive balancing requires at least one observed row.", call. = FALSE)
  }

  observed_data <- data[observed_rows, , drop = FALSE]
  missing_covariates <- setdiff(covariates, names(observed_data))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("Missing covariate columns: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  sample_weights_full <- cdmc_resolve_sample_weights(sample_weights, data)
  base_weights <- if (is.null(sample_weights_full)) rep(1, nrow(observed_data)) else sample_weights_full[observed_rows]
  if (sum(base_weights) <= sqrt(.Machine$double.eps)) {
    stop("Adaptive balancing requires positive total sample weight on observed rows.", call. = FALSE)
  }

  score_basis <- cdmc_prepare_adaptive_balance_score_basis(
    methods = methods,
    observed_data = observed_data,
    dose = dose,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    entropy_balance_degree = entropy_balance_degree,
    kernel_balance_degree = kernel_balance_degree,
    kernel_balance_centers = kernel_balance_centers,
    kernel_balance_bandwidth = kernel_balance_bandwidth,
    base_weights = base_weights
  )

  candidate_fits <- vector("list", length(methods))
  names(candidate_fits) <- methods
  candidate_rows <- vector("list", length(methods))

  for (method_index in seq_along(methods)) {
    current_method <- methods[method_index]
    current_fit <- tryCatch(
      cdmc_fit_adaptive_balance_candidate(
        method = current_method,
        data = data,
        dose = dose,
        covariates = covariates,
        time = time,
        time_effects = time_effects,
        model = model,
        df = df,
        spline_covariates = spline_covariates,
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
        sample_weights = sample_weights_full,
        max_weight = max_weight
      ),
      error = function(error) error
    )

    if (inherits(current_fit, "error")) {
      candidate_rows[[method_index]] <- data.frame(
        method = current_method,
        fit_success = FALSE,
        fit_error = conditionMessage(current_fit),
        max_abs_standardized_balance = NA_real_,
        ess = NA_real_,
        ess_fraction = NA_real_,
        applied_max_weight = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      next
    }

    candidate_fits[[method_index]] <- current_fit
    candidate_rows[[method_index]] <- cdmc_adaptive_balance_score_row(
      method = current_method,
      fit = current_fit,
      observed_rows = observed_rows,
      score_basis = score_basis
    )
  }

  candidate_scores <- do.call(rbind, candidate_rows)
  selected_index <- cdmc_select_adaptive_balance_index(candidate_scores)
  if (length(selected_index) < 1L) {
    stop(
      sprintf(
        "Adaptive balancing failed for all candidate methods (%s).",
        paste(sprintf("%s: %s", candidate_scores$method, candidate_scores$fit_error), collapse = "; ")
      ),
      call. = FALSE
    )
  }

  candidate_scores$selected <- FALSE
  candidate_scores$selected[selected_index] <- TRUE
  selected_method <- candidate_scores$method[selected_index]
  selected_fit <- candidate_fits[[selected_method]]

  result <- list(
    call = match.call(),
    weight_method = "adaptive_balance",
    selected_method = selected_method,
    methods = methods,
    formula = selected_fit$formula,
    weights = selected_fit$weights,
    raw_weights = selected_fit$raw_weights %||% selected_fit$weights,
    covariates = covariates,
    time = time,
    time_effects = time_effects,
    model = model,
    df = df,
    spline_covariates = spline_covariates,
    candidate_scores = candidate_scores,
    selected_fit = selected_fit,
    requested_max_weight = selected_fit$requested_max_weight %||% max_weight,
    applied_max_weight = selected_fit$applied_max_weight %||% NULL,
    sample_weights = sample_weights_full,
    converged = TRUE,
    max_abs_standardized_balance = candidate_scores$max_abs_standardized_balance[selected_index],
    max_abs_balance = selected_fit$max_abs_balance %||% NULL
  )

  class(result) <- "cdmc_adaptive_balance_weights"
  result
}

print.cdmc_adaptive_balance_weights <- function(x, ...) {
  cat("causaldosemc adaptive balance weights\n")
  cat(sprintf("  selected method: %s\n", x$selected_method))
  cat(sprintf("  candidate methods: %s\n", paste(x$methods, collapse = ", ")))
  if (!is.null(x$applied_max_weight) && is.finite(x$applied_max_weight)) {
    cat(sprintf("  selected weight cap: %.6g\n", x$applied_max_weight))
  }
  if (!is.null(x$max_abs_standardized_balance) && is.finite(x$max_abs_standardized_balance)) {
    cat(sprintf("  selected max abs standardized balance: %.6g\n", x$max_abs_standardized_balance))
  }
  if (is.data.frame(x$candidate_scores) && nrow(x$candidate_scores) > 0L) {
    for (row_index in seq_len(nrow(x$candidate_scores))) {
      row <- x$candidate_scores[row_index, , drop = FALSE]
      if (!isTRUE(row$fit_success)) {
        cat(sprintf("  %s: failed (%s)\n", row$method, row$fit_error))
      } else {
        cat(sprintf(
          "  %s: imbalance %.6g, ESS fraction %.3f%s\n",
          row$method,
          row$max_abs_standardized_balance,
          row$ess_fraction,
          if (isTRUE(row$selected)) " [selected]" else ""
        ))
      }
    }
  }
  invisible(x)
}

cdmc_call_cbps <- function(...) {
  withCallingHandlers(
    CBPS::CBPS(...),
    warning = function(condition) {
      if (grepl("Treatment vector is numeric", conditionMessage(condition), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

cdmc_cbps_weights <- function(
  data,
  dose,
  covariates = NULL,
  time = NULL,
  time_effects = FALSE,
  model = c("linear", "spline"),
  df = 4L,
  spline_covariates = covariates,
  standardize = TRUE,
  method = c("over", "exact"),
  iterations = 1000L,
  twostep = TRUE,
  sample_weights = NULL,
  ...
) {
  cdmc_assert_installed("CBPS")

  data <- as.data.frame(data)
  method <- match.arg(method)
  model <- match.arg(model)

  if (!is.character(dose) || length(dose) != 1L || !dose %in% names(data)) {
    stop("dose must be a single column name present in data.", call. = FALSE)
  }

  if (!is.null(covariates) && !is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }
  covariates <- covariates %||% character(0)

  if (!is.logical(time_effects) || length(time_effects) != 1L || is.na(time_effects)) {
    stop("time_effects must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(time_effects) && (!is.character(time) || length(time) != 1L || !time %in% names(data))) {
    stop("time must be a single column name present in data when time_effects = TRUE.", call. = FALSE)
  }

  if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df < 1) {
    stop("df must be a positive integer.", call. = FALSE)
  }
  df <- as.integer(df)

  if (length(covariates) == 0L && !isTRUE(time_effects)) {
    stop("At least one balancing covariate or time_effects = TRUE is required for CBPS weighting.", call. = FALSE)
  }

  missing_covariates <- setdiff(covariates, names(data))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("Missing covariate columns: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  observed_rows <- cdmc_panel_observed_rows(data)
  if (!any(observed_rows)) {
    stop("CBPS weighting requires at least one observed row.", call. = FALSE)
  }

  observed_data <- data[observed_rows, , drop = FALSE]
  spline_covariates <- cdmc_resolve_gps_spline_covariates(covariates, spline_covariates)
  sample_weights_full <- cdmc_resolve_sample_weights(sample_weights, data)
  sample_weights_fit <- if (is.null(sample_weights_full)) NULL else sample_weights_full[observed_rows]
  formula <- cdmc_build_gps_formula(
    data = observed_data,
    dose = dose,
    covariates = covariates,
    time = if (isTRUE(time_effects)) time else dose,
    gps_time_effects = time_effects,
    gps_model = model,
    gps_df = df,
    gps_spline_covariates = spline_covariates
  )
  requested_method <- method
  fit <- tryCatch(
    cdmc_call_cbps(
      formula = formula,
      data = observed_data,
      ATT = 0,
      iterations = as.integer(iterations),
      standardize = standardize,
      method = method,
      twostep = twostep,
      sample.weights = sample_weights_fit,
      ...
    ),
    error = function(condition) {
      fallback_needed <- identical(requested_method, "over") && grepl(
        "infinite value in the weighting matrix|just-identified version",
        conditionMessage(condition),
        ignore.case = TRUE
      )

      if (!fallback_needed) {
        stop(condition)
      }

      method <<- "exact"
      cdmc_call_cbps(
        formula = formula,
        data = observed_data,
        ATT = 0,
        iterations = as.integer(iterations),
        standardize = standardize,
        method = method,
        twostep = twostep,
        sample.weights = sample_weights_fit,
        ...
      )
    }
  )

  weights <- fit$weights
  if (!is.numeric(weights) || length(weights) != nrow(observed_data) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("CBPS did not return a valid nonnegative weight vector for the supplied data.", call. = FALSE)
  }

  expanded_weights <- numeric(nrow(data))
  expanded_weights[observed_rows] <- as.numeric(weights)

  result <- list(
    call = match.call(),
    formula = formula,
    weights = expanded_weights,
    requested_method = requested_method,
    method = method,
    standardize = standardize,
    sample_weights = sample_weights_full,
    converged = fit$converged,
    fit = fit
  )

  class(result) <- "cdmc_cbps_weights"
  result
}

print.cdmc_cbps_weights <- function(x, ...) {
  cat("causaldosemc CBPS weights\n")
  cat(sprintf("  formula: %s\n", paste(deparse(x$formula), collapse = " ")))
  cat(sprintf("  observations: %d\n", length(x$weights)))
  cat(sprintf("  method: %s\n", x$method))
  if (!identical(x$requested_method, x$method)) {
    cat(sprintf("  requested method: %s\n", x$requested_method))
  }
  cat(sprintf("  standardized: %s\n", if (x$standardize) "yes" else "no"))
  if (!is.null(x$converged)) {
    cat(sprintf("  convergence code: %s\n", paste(x$converged, collapse = ", ")))
  }

  invisible(x)
}