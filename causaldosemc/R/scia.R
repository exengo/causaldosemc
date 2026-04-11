cdmc_prepare_scia_sample <- function(
  object,
  lags,
  outcome_proxy = c("tau", "y"),
  covariates = object$covariates,
  include_current_covariates = TRUE,
  include_covariate_lags = TRUE,
  include_unit_effects = TRUE,
  include_time_effects = TRUE
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  lags <- as.integer(lags)
  if (!is.numeric(lags) || length(lags) != 1L || lags < 1L) {
    stop("lags must be an integer greater than or equal to 1.", call. = FALSE)
  }

  outcome_proxy <- match.arg(outcome_proxy)
  covariates <- covariates %||% character(0)
  if (!is.character(covariates)) {
    stop("covariates must be NULL or a character vector of column names.", call. = FALSE)
  }

  missing_covariates <- setdiff(covariates, names(object$data))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("Missing covariate columns: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  dose_lag_array <- cdmc_build_lagged_doses(object$dose_matrix, lag_order = lags)
  outcome_matrix <- switch(
    outcome_proxy,
    tau = object$effect$tau,
    y = object$y_matrix
  )
  outcome_lag_array <- cdmc_build_lagged_doses(outcome_matrix, lag_order = lags)

  valid_history <- apply(!is.na(dose_lag_array), c(1, 2), all) &
    apply(!is.na(outcome_lag_array), c(1, 2), all)

  covariate_matrices <- list()
  covariate_lag_arrays <- list()
  if (length(covariates) > 0L) {
    for (covariate in covariates) {
      covariate_matrix <- matrix(
        object$data[[covariate]],
        nrow = object$n_units,
        ncol = object$n_times,
        byrow = TRUE
      )
      covariate_matrices[[covariate]] <- covariate_matrix

      if (isTRUE(include_covariate_lags)) {
        covariate_lag_arrays[[covariate]] <- cdmc_build_lagged_doses(covariate_matrix, lag_order = lags)
        valid_history <- valid_history & apply(!is.na(covariate_lag_arrays[[covariate]]), c(1, 2), all)
      }
    }
  }

  sample_indices <- which(valid_history, arr.ind = TRUE)
  if (nrow(sample_indices) == 0L) {
    stop("No observations retain the requested lag history for SCIA screening.", call. = FALSE)
  }

  linear_indices <- (sample_indices[, 1L] - 1L) * object$n_times + sample_indices[, 2L]
  sample <- object$data[linear_indices, c(object$unit, object$time), drop = FALSE]
  sample$.cdmc_unit_index <- sample_indices[, 1L]
  sample$.cdmc_time_index <- sample_indices[, 2L]
  sample$.cdmc_current_dose <- object$dose_matrix[sample_indices]
  sample$.cdmc_current_active <- as.integer(
    cdmc_active_dose_mask(object$dose_matrix[sample_indices], zero_tolerance = object$zero_tolerance)
  )

  dose_history_terms <- character(lags)
  outcome_history_terms <- character(lags)
  for (lag_index in seq_len(lags)) {
    dose_name <- paste0("dose_lag", lag_index)
    proxy_name <- paste0(outcome_proxy, "_lag", lag_index)
    sample[[dose_name]] <- dose_lag_array[, , lag_index + 1L][sample_indices]
    sample[[proxy_name]] <- outcome_lag_array[, , lag_index + 1L][sample_indices]
    dose_history_terms[[lag_index]] <- dose_name
    outcome_history_terms[[lag_index]] <- proxy_name
  }

  current_covariate_terms <- character(0)
  lagged_covariate_terms <- character(0)
  if (length(covariates) > 0L) {
    if (isTRUE(include_current_covariates)) {
      for (covariate in covariates) {
        sample[[covariate]] <- covariate_matrices[[covariate]][sample_indices]
        current_covariate_terms <- c(current_covariate_terms, covariate)
      }
    }

    if (isTRUE(include_covariate_lags)) {
      for (covariate in covariates) {
        for (lag_index in seq_len(lags)) {
          lagged_name <- paste0(covariate, "_lag", lag_index)
          sample[[lagged_name]] <- covariate_lag_arrays[[covariate]][, , lag_index + 1L][sample_indices]
          lagged_covariate_terms <- c(lagged_covariate_terms, lagged_name)
        }
      }
    }
  }

  restricted_terms <- c(dose_history_terms, current_covariate_terms, lagged_covariate_terms)
  if (isTRUE(include_unit_effects)) {
    restricted_terms <- c(restricted_terms, "factor(.cdmc_unit_index)")
  }
  if (isTRUE(include_time_effects)) {
    restricted_terms <- c(restricted_terms, "factor(.cdmc_time_index)")
  }

  restricted_formula <- stats::reformulate(
    termlabels = restricted_terms,
    response = ".cdmc_current_dose"
  )
  augmented_formula <- stats::reformulate(
    termlabels = c(restricted_terms, outcome_history_terms),
    response = ".cdmc_current_dose"
  )

  list(
    sample = sample,
    restricted_formula = restricted_formula,
    augmented_formula = augmented_formula,
    dose_history_terms = dose_history_terms,
    outcome_history_terms = outcome_history_terms,
    current_covariate_terms = current_covariate_terms,
    lagged_covariate_terms = lagged_covariate_terms,
    lags = lags,
    outcome_proxy = outcome_proxy,
    covariates = covariates,
    include_current_covariates = include_current_covariates,
    include_covariate_lags = include_covariate_lags,
    include_unit_effects = include_unit_effects,
    include_time_effects = include_time_effects
  )
}

cdmc_scia_screen_table <- function(augmented_fit, outcome_history_terms) {
  coefficient_table <- summary(augmented_fit)$coefficients
  available_terms <- intersect(outcome_history_terms, rownames(coefficient_table))

  if (length(available_terms) == 0L) {
    return(data.frame(
      term = character(0),
      estimate = numeric(0),
      std_error = numeric(0),
      statistic = numeric(0),
      p_value = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    term = available_terms,
    estimate = unname(coefficient_table[available_terms, 1L]),
    std_error = unname(coefficient_table[available_terms, 2L]),
    statistic = unname(coefficient_table[available_terms, 3L]),
    p_value = unname(coefficient_table[available_terms, 4L]),
    stringsAsFactors = FALSE
  )
}

cdmc_scia_test <- function(
  object,
  lags = 1L,
  outcome_proxy = c("tau", "y"),
  covariates = object$covariates,
  include_current_covariates = TRUE,
  include_covariate_lags = TRUE,
  include_unit_effects = TRUE,
  include_time_effects = TRUE,
  alpha = 0.05
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).", call. = FALSE)
  }

  prepared <- cdmc_prepare_scia_sample(
    object = object,
    lags = lags,
    outcome_proxy = outcome_proxy,
    covariates = covariates,
    include_current_covariates = include_current_covariates,
    include_covariate_lags = include_covariate_lags,
    include_unit_effects = include_unit_effects,
    include_time_effects = include_time_effects
  )

  restricted_fit <- stats::lm(formula = prepared$restricted_formula, data = prepared$sample)
  augmented_fit <- stats::lm(formula = prepared$augmented_formula, data = prepared$sample)
  comparison <- stats::anova(restricted_fit, augmented_fit)

  f_statistic <- if (nrow(comparison) >= 2L) comparison$F[2L] else NA_real_
  p_value <- if (nrow(comparison) >= 2L) comparison$`Pr(>F)`[2L] else NA_real_
  delta_r_squared <- summary(augmented_fit)$adj.r.squared - summary(restricted_fit)$adj.r.squared

  result <- list(
    call = match.call(),
    lags = prepared$lags,
    outcome_proxy = prepared$outcome_proxy,
    covariates = prepared$covariates,
    include_current_covariates = prepared$include_current_covariates,
    include_covariate_lags = prepared$include_covariate_lags,
    include_unit_effects = prepared$include_unit_effects,
    include_time_effects = prepared$include_time_effects,
    alpha = alpha,
    sample_size = nrow(prepared$sample),
    restricted_formula = prepared$restricted_formula,
    augmented_formula = prepared$augmented_formula,
    restricted_fit = restricted_fit,
    augmented_fit = augmented_fit,
    comparison = comparison,
    f_statistic = f_statistic,
    p_value = p_value,
    delta_r_squared = delta_r_squared,
    passed = is.finite(p_value) && p_value >= alpha,
    screen_table = cdmc_scia_screen_table(augmented_fit, prepared$outcome_history_terms),
    sample = prepared$sample,
    fit_object = object
  )

  class(result) <- "cdmc_scia_test"
  result
}

print.cdmc_scia_test <- function(x, ...) {
  cat("causaldosemc SCIA screening diagnostic\n")
  cat(sprintf("  lag order screened: %d\n", x$lags))
  cat(sprintf("  outcome proxy: %s\n", x$outcome_proxy))
  cat(sprintf("  sample size: %d\n", x$sample_size))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  if (is.finite(x$f_statistic)) {
    cat(sprintf("  incremental F statistic: %.6g\n", x$f_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  p value: %.6g\n", x$p_value))
    cat(sprintf("  passed: %s\n", if (x$passed) "yes" else "no"))
  }
  cat(sprintf("  adjusted R-squared increment: %.6g\n", x$delta_r_squared))

  if (nrow(x$screen_table) > 0L) {
    cat("  lagged outcome terms:\n")
    print(x$screen_table, row.names = FALSE)
  }

  invisible(x)
}