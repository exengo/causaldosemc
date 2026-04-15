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

  resolved <- if (inherits(weights, "cdmc_cbps_weights")) {
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
    CBPS::CBPS(
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
      CBPS::CBPS(
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