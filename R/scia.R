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
    restricted_terms = restricted_terms,
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

cdmc_format_scia_block_name <- function(lag_indices, outcome_history_terms) {
  if (length(lag_indices) == 1L) {
    return(outcome_history_terms[[lag_indices]])
  }

  paste0("lags_", paste(lag_indices, collapse = "_"))
}

cdmc_resolve_scia_restriction_blocks <- function(lags, outcome_history_terms, restriction_blocks = NULL) {
  if (is.null(restriction_blocks)) {
    resolved <- as.list(seq_len(lags))
    names(resolved) <- vapply(seq_len(lags), function(index) {
      cdmc_format_scia_block_name(index, outcome_history_terms)
    }, character(1))
    return(resolved)
  }

  if (is.numeric(restriction_blocks)) {
    restriction_blocks <- list(restriction_blocks)
  }

  if (!is.list(restriction_blocks) || length(restriction_blocks) < 1L) {
    stop("restriction_blocks must be NULL, a numeric vector, or a nonempty list of lag indices.", call. = FALSE)
  }

  resolved <- lapply(restriction_blocks, function(block) {
    if (!is.numeric(block) || length(block) < 1L || any(!is.finite(block))) {
      stop("Each SCIA restriction block must contain one or more finite lag indices.", call. = FALSE)
    }

    block <- sort(unique(as.integer(block)))
    if (any(block < 1L) || any(block > lags)) {
      stop(
        sprintf("SCIA restriction lag indices must lie between 1 and %d.", lags),
        call. = FALSE
      )
    }

    block
  })

  block_names <- names(restriction_blocks)
  if (is.null(block_names)) {
    block_names <- rep("", length(resolved))
  }
  block_names[!nzchar(block_names)] <- vapply(resolved[!nzchar(block_names)], function(block) {
    cdmc_format_scia_block_name(block, outcome_history_terms)
  }, character(1))

  if (anyDuplicated(block_names)) {
    stop("SCIA restriction block names must be unique.", call. = FALSE)
  }

  names(resolved) <- block_names
  resolved
}

cdmc_scia_conditioning_table <- function(prepared) {
  rows <- list(
    dose_history = prepared$dose_history_terms,
    current_covariates = prepared$current_covariate_terms,
    lagged_covariates = prepared$lagged_covariate_terms,
    unit_effects = if (isTRUE(prepared$include_unit_effects)) "factor(.cdmc_unit_index)" else character(0),
    time_effects = if (isTRUE(prepared$include_time_effects)) "factor(.cdmc_time_index)" else character(0)
  )

  rows <- rows[vapply(rows, length, integer(1)) > 0L]
  if (length(rows) < 1L) {
    return(data.frame(
      block = character(0),
      n_terms = integer(0),
      terms = character(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    block = names(rows),
    n_terms = vapply(rows, length, integer(1)),
    terms = vapply(rows, function(row) paste(row, collapse = ", "), character(1)),
    stringsAsFactors = FALSE
  )
}

cdmc_cluster_robust_vcov <- function(model, cluster_ids) {
  # CR1-style cluster-robust covariance for stats::lm objects.
  # V = (X'X)^{-1} * sum_g X_g' e_g e_g' X_g * (X'X)^{-1} * G/(G-1) * (n-1)/(n-k)
  X <- stats::model.matrix(model)
  e <- stats::residuals(model)
  if (length(cluster_ids) != nrow(X)) {
    return(NULL)
  }
  keep <- !is.na(e)
  X <- X[keep, , drop = FALSE]
  e <- e[keep]
  cluster_ids <- cluster_ids[keep]
  groups <- split(seq_along(e), cluster_ids)
  G <- length(groups)
  n <- nrow(X)
  k <- ncol(X)
  if (G < 2L || n <= k) {
    return(NULL)
  }
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(err) NULL)
  if (is.null(XtX_inv)) return(NULL)
  meat <- matrix(0, k, k)
  for (idx in groups) {
    Xg <- X[idx, , drop = FALSE]
    eg <- e[idx]
    score <- crossprod(Xg, eg)
    meat <- meat + tcrossprod(score)
  }
  finite_correction <- (G / (G - 1L)) * ((n - 1L) / (n - k))
  V <- XtX_inv %*% meat %*% XtX_inv * finite_correction
  rownames(V) <- colnames(X)
  colnames(V) <- colnames(X)
  V
}

cdmc_cluster_wald_test <- function(model, terms, cluster_ids) {
  # Wald F-style test on a subset of coefficients with cluster-robust V.
  V <- cdmc_cluster_robust_vcov(model, cluster_ids)
  if (is.null(V)) {
    return(list(f_statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_))
  }
  beta <- stats::coef(model)
  available <- intersect(terms, names(beta))
  if (length(available) == 0L) {
    return(list(f_statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_))
  }
  beta_k <- beta[available]
  V_kk <- V[available, available, drop = FALSE]
  V_kk_inv <- tryCatch(solve(V_kk), error = function(e) NULL)
  if (is.null(V_kk_inv)) {
    return(list(f_statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_))
  }
  q <- length(available)
  G <- length(unique(cluster_ids[!is.na(stats::residuals(model))]))
  wald <- as.numeric(t(beta_k) %*% V_kk_inv %*% beta_k)
  f_stat <- wald / q
  df1 <- q
  df2 <- max(G - 1L, 1L)
  p_value <- stats::pf(f_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  list(f_statistic = f_stat, df1 = df1, df2 = df2, p_value = p_value)
}

cdmc_scia_compare_restrictions <- function(prepared, restricted_fit, restriction_blocks, alpha = 0.05, p_adjust_method = c("holm", "none"), cluster_ids = NULL) {
  p_adjust_method <- match.arg(p_adjust_method)

  rows <- lapply(names(restriction_blocks), function(block_name) {
    block_lags <- restriction_blocks[[block_name]]
    block_terms <- prepared$outcome_history_terms[block_lags]
    block_formula <- stats::reformulate(
      termlabels = c(prepared$restricted_terms, block_terms),
      response = ".cdmc_current_dose"
    )
    block_fit <- stats::lm(formula = block_formula, data = prepared$sample)
    cluster_test <- if (!is.null(cluster_ids)) {
      cdmc_cluster_wald_test(block_fit, terms = block_terms, cluster_ids = cluster_ids)
    } else {
      list(f_statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, p_value = NA_real_)
    }
    iid_comparison <- stats::anova(restricted_fit, block_fit)

    data.frame(
      restriction_name = block_name,
      lags = paste(block_lags, collapse = ","),
      terms = paste(block_terms, collapse = ", "),
      n_terms = length(block_terms),
      f_statistic = if (is.finite(cluster_test$f_statistic)) cluster_test$f_statistic else if (nrow(iid_comparison) >= 2L) iid_comparison$F[2L] else NA_real_,
      p_value = if (is.finite(cluster_test$p_value)) cluster_test$p_value else if (nrow(iid_comparison) >= 2L) iid_comparison$`Pr(>F)`[2L] else NA_real_,
      iid_p_value = if (nrow(iid_comparison) >= 2L) iid_comparison$`Pr(>F)`[2L] else NA_real_,
      cluster_p_value = cluster_test$p_value,
      delta_r_squared = summary(block_fit)$adj.r.squared - summary(restricted_fit)$adj.r.squared,
      stringsAsFactors = FALSE
    )
  })

  table <- do.call(rbind, rows)
  if (nrow(table) < 1L) {
    table$adjusted_p_value <- numeric(0)
    table$passed <- logical(0)
    return(table)
  }

  adjusted <- if (identical(p_adjust_method, "none")) {
    table$p_value
  } else {
    stats::p.adjust(table$p_value, method = p_adjust_method)
  }
  table$adjusted_p_value <- adjusted
  table$passed <- is.na(adjusted) | adjusted >= alpha
  table
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
  restriction_blocks = NULL,
  p_adjust_method = c("holm", "none"),
  alpha = 0.05
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).", call. = FALSE)
  }

  p_adjust_method <- match.arg(p_adjust_method)

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
  cluster_ids <- prepared$sample$.cdmc_unit_index
  cluster_global_test <- cdmc_cluster_wald_test(
    augmented_fit,
    terms = prepared$outcome_history_terms,
    cluster_ids = cluster_ids
  )
  resolved_restriction_blocks <- cdmc_resolve_scia_restriction_blocks(
    lags = prepared$lags,
    outcome_history_terms = prepared$outcome_history_terms,
    restriction_blocks = restriction_blocks
  )
  restriction_table <- cdmc_scia_compare_restrictions(
    prepared = prepared,
    restricted_fit = restricted_fit,
    restriction_blocks = resolved_restriction_blocks,
    alpha = alpha,
    p_adjust_method = p_adjust_method,
    cluster_ids = cluster_ids
  )

  iid_f <- if (nrow(comparison) >= 2L) comparison$F[2L] else NA_real_
  iid_p <- if (nrow(comparison) >= 2L) comparison$`Pr(>F)`[2L] else NA_real_
  f_statistic <- if (is.finite(cluster_global_test$f_statistic)) cluster_global_test$f_statistic else iid_f
  p_value <- if (is.finite(cluster_global_test$p_value)) cluster_global_test$p_value else iid_p
  delta_r_squared <- summary(augmented_fit)$adj.r.squared - summary(restricted_fit)$adj.r.squared
  global_passed <- is.finite(p_value) && p_value >= alpha
  failed_restrictions <- if (nrow(restriction_table) > 0L) {
    restriction_table$restriction_name[!restriction_table$passed]
  } else {
    character(0)
  }

  result <- list(
    call = match.call(),
    lags = prepared$lags,
    outcome_proxy = prepared$outcome_proxy,
    covariates = prepared$covariates,
    include_current_covariates = prepared$include_current_covariates,
    include_covariate_lags = prepared$include_covariate_lags,
    include_unit_effects = prepared$include_unit_effects,
    include_time_effects = prepared$include_time_effects,
    conditioning_table = cdmc_scia_conditioning_table(prepared),
    restriction_blocks = resolved_restriction_blocks,
    p_adjust_method = p_adjust_method,
    alpha = alpha,
    sample_size = nrow(prepared$sample),
    restricted_formula = prepared$restricted_formula,
    augmented_formula = prepared$augmented_formula,
    restricted_fit = restricted_fit,
    augmented_fit = augmented_fit,
    comparison = comparison,
    f_statistic = f_statistic,
    p_value = p_value,
    iid_f_statistic = iid_f,
    iid_p_value = iid_p,
    cluster_f_statistic = cluster_global_test$f_statistic,
    cluster_p_value = cluster_global_test$p_value,
    cluster_df1 = cluster_global_test$df1,
    cluster_df2 = cluster_global_test$df2,
    delta_r_squared = delta_r_squared,
    global_passed = global_passed,
    restriction_table = restriction_table,
    failed_restrictions = failed_restrictions,
    passed = global_passed && length(failed_restrictions) == 0L,
    screen_table = cdmc_scia_screen_table(augmented_fit, prepared$outcome_history_terms),
    sample = prepared$sample,
    fit_object = object
  )

  class(result) <- "cdmc_scia_test"
  result
}

print.cdmc_scia_test <- function(x, ...) {
  cat("causaldosemc SCIA restriction diagnostic\n")
  cat(sprintf("  lag order screened: %d\n", x$lags))
  cat(sprintf("  outcome proxy: %s\n", x$outcome_proxy))
  cat(sprintf("  sample size: %d\n", x$sample_size))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  if (is.finite(x$f_statistic)) {
    cat(sprintf("  global incremental F statistic: %.6g\n", x$f_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  global p value: %.6g\n", x$p_value))
    cat(sprintf("  global passed: %s\n", if (isTRUE(x$global_passed)) "yes" else "no"))
  }
  cat(sprintf("  adjusted R-squared increment: %.6g\n", x$delta_r_squared))
  cat(sprintf("  restriction p adjustment: %s\n", x$p_adjust_method))
  cat(sprintf("  overall passed: %s\n", if (x$passed) "yes" else "no"))

  if (!is.null(x$conditioning_table) && nrow(x$conditioning_table) > 0L) {
    cat("  conditioning blocks:\n")
    print(x$conditioning_table, row.names = FALSE)
  }

  if (!is.null(x$restriction_table) && nrow(x$restriction_table) > 0L) {
    cat("  restriction checks:\n")
    print(x$restriction_table[, c("restriction_name", "lags", "p_value", "adjusted_p_value", "passed"), drop = FALSE], row.names = FALSE)
    if (length(x$failed_restrictions) > 0L) {
      cat(sprintf("  failed restrictions: %s\n", paste(x$failed_restrictions, collapse = ", ")))
    }
  }

  if (nrow(x$screen_table) > 0L) {
    cat("  lagged outcome terms:\n")
    print(x$screen_table, row.names = FALSE)
  }

  invisible(x)
}