cdmc_exit_distances <- function(dose_matrix, zero_tolerance) {
  n_units <- nrow(dose_matrix)
  n_times <- ncol(dose_matrix)
  distances <- matrix(NA_integer_, nrow = n_units, ncol = n_times)

  for (unit_index in seq_len(n_units)) {
    last_positive <- NA_integer_

    for (time_index in seq_len(n_times)) {
      if (cdmc_active_dose_mask(dose_matrix[unit_index, time_index], zero_tolerance = zero_tolerance)) {
        last_positive <- time_index
        next
      }

      if (is.na(last_positive)) {
        next
      }

      distances[unit_index, time_index] <- time_index - last_positive
    }
  }

  distances
}

cdmc_carryover_test <- function(object, periods = 1L) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods < 1L)) {
    stop("periods must contain one or more positive integers.", call. = FALSE)
  }

  exit_distance <- cdmc_exit_distances(
    dose_matrix = object$dose_matrix,
    zero_tolerance = object$zero_tolerance
  )

  target_mask <- cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance) & exit_distance %in% periods
  if (!any(target_mask)) {
    stop(
      "No post-exit zero-dose observations match the requested carryover periods.",
      call. = FALSE
    )
  }

  tau_values <- object$effect$tau[target_mask]
  sample_size <- length(tau_values)
  mean_tau <- mean(tau_values)
  sd_tau <- if (sample_size > 1L) stats::sd(tau_values) else 0
  naive_se <- if (sample_size > 1L) sd_tau / sqrt(sample_size) else NA_real_

  # Cluster-robust SE for the sample mean, clustering by unit. Residuals from
  # the same unit are correlated through the shared low-rank baseline, so the
  # iid-mean SE is anti-conservative. Use CR1-style: SE^2 = (G/(G-1)) *
  # sum_g (sum_{i in g} (tau_i - mean))^2 / n^2.
  target_rows <- object$data[, c(object$unit, object$time), drop = FALSE]
  target_rows$exit_distance <- cdmc_flatten_matrix(exit_distance)
  target_rows$tau <- cdmc_flatten_matrix(object$effect$tau)
  target_rows$eligible_control <- cdmc_flatten_matrix(object$eligible_mask)
  target_rows <- target_rows[cdmc_flatten_matrix(target_mask), , drop = FALSE]

  cluster_ids <- target_rows[[object$unit]]
  n_clusters <- length(unique(cluster_ids))
  if (sample_size > 1L && n_clusters > 1L) {
    centred <- tau_values - mean_tau
    cluster_sums <- tapply(centred, cluster_ids, sum)
    cluster_se <- sqrt((n_clusters / (n_clusters - 1L)) * sum(cluster_sums^2) / sample_size^2)
  } else {
    cluster_se <- NA_real_
  }

  # Prefer the cluster-robust SE when available; fall back to naive SE so the
  # test still runs on degenerate one-cluster cases.
  standard_error <- if (is.finite(cluster_se) && cluster_se > 0) cluster_se else naive_se
  effective_df <- if (is.finite(cluster_se)) max(n_clusters - 1L, 1L) else max(sample_size - 1L, 1L)
  t_statistic <- if (is.finite(standard_error) && standard_error > 0) mean_tau / standard_error else NA_real_
  p_value <- if (is.finite(t_statistic)) {
    2 * stats::pt(-abs(t_statistic), df = effective_df)
  } else {
    NA_real_
  }

  result <- list(
    call = match.call(),
    periods = periods,
    n = sample_size,
    n_clusters = n_clusters,
    mean_tau = mean_tau,
    sd_tau = sd_tau,
    standard_error = standard_error,
    naive_standard_error = naive_se,
    cluster_standard_error = cluster_se,
    t_statistic = t_statistic,
    p_value = p_value,
    cells = target_rows,
    fit_object = object
  )

  class(result) <- "cdmc_carryover_test"
  result
}

print.cdmc_carryover_test <- function(x, ...) {
  cat("causaldosemc carryover diagnostic\n")
  cat(sprintf("  periods tested: %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  matched zero-dose exit cells: %d\n", x$n))
  if (!is.null(x$n_clusters)) {
    cat(sprintf("  unit clusters: %d\n", x$n_clusters))
  }
  cat(sprintf("  mean residual effect: %.6g\n", x$mean_tau))
  if (is.finite(x$standard_error)) {
    se_label <- if (!is.null(x$cluster_standard_error) && is.finite(x$cluster_standard_error)) {
      "cluster-robust standard error"
    } else {
      "standard error"
    }
    cat(sprintf("  %s: %.6g\n", se_label, x$standard_error))
  }
  if (is.finite(x$t_statistic)) {
    cat(sprintf("  t statistic: %.6g\n", x$t_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  p value: %.6g\n", x$p_value))
  }

  invisible(x)
}

cdmc_entry_distances <- function(dose_matrix, zero_tolerance) {
  # For each cell, distance (in periods) to the nearest *future* active dose.
  # NA if no subsequent activation. Mirrors cdmc_exit_distances backward.
  n_units <- nrow(dose_matrix)
  n_times <- ncol(dose_matrix)
  distances <- matrix(NA_integer_, nrow = n_units, ncol = n_times)

  for (unit_index in seq_len(n_units)) {
    next_positive <- NA_integer_

    for (time_index in rev(seq_len(n_times))) {
      if (cdmc_active_dose_mask(dose_matrix[unit_index, time_index], zero_tolerance = zero_tolerance)) {
        next_positive <- time_index
        next
      }

      if (is.na(next_positive)) {
        next
      }

      distances[unit_index, time_index] <- next_positive - time_index
    }
  }

  distances
}

cdmc_anticipation_test <- function(object, periods = 1L) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods < 1L)) {
    stop("periods must contain one or more positive integers.", call. = FALSE)
  }

  entry_distance <- cdmc_entry_distances(
    dose_matrix = object$dose_matrix,
    zero_tolerance = object$zero_tolerance
  )

  target_mask <- cdmc_zero_dose_mask(object$dose_matrix, zero_tolerance = object$zero_tolerance) &
    entry_distance %in% periods
  if (!any(target_mask)) {
    stop(
      "No pre-entry zero-dose observations match the requested anticipation periods.",
      call. = FALSE
    )
  }

  tau_values <- object$effect$tau[target_mask]
  sample_size <- length(tau_values)
  mean_tau <- mean(tau_values)
  sd_tau <- if (sample_size > 1L) stats::sd(tau_values) else 0
  naive_se <- if (sample_size > 1L) sd_tau / sqrt(sample_size) else NA_real_

  target_rows <- object$data[, c(object$unit, object$time), drop = FALSE]
  target_rows$entry_distance <- cdmc_flatten_matrix(entry_distance)
  target_rows$tau <- cdmc_flatten_matrix(object$effect$tau)
  target_rows$eligible_control <- cdmc_flatten_matrix(object$eligible_mask)
  target_rows <- target_rows[cdmc_flatten_matrix(target_mask), , drop = FALSE]

  cluster_ids <- target_rows[[object$unit]]
  n_clusters <- length(unique(cluster_ids))
  if (sample_size > 1L && n_clusters > 1L) {
    centred <- tau_values - mean_tau
    cluster_sums <- tapply(centred, cluster_ids, sum)
    cluster_se <- sqrt((n_clusters / (n_clusters - 1L)) * sum(cluster_sums^2) / sample_size^2)
  } else {
    cluster_se <- NA_real_
  }

  standard_error <- if (is.finite(cluster_se) && cluster_se > 0) cluster_se else naive_se
  effective_df <- if (is.finite(cluster_se)) max(n_clusters - 1L, 1L) else max(sample_size - 1L, 1L)
  t_statistic <- if (is.finite(standard_error) && standard_error > 0) mean_tau / standard_error else NA_real_
  p_value <- if (is.finite(t_statistic)) {
    2 * stats::pt(-abs(t_statistic), df = effective_df)
  } else {
    NA_real_
  }

  result <- list(
    call = match.call(),
    periods = periods,
    n = sample_size,
    n_clusters = n_clusters,
    mean_tau = mean_tau,
    sd_tau = sd_tau,
    standard_error = standard_error,
    naive_standard_error = naive_se,
    cluster_standard_error = cluster_se,
    t_statistic = t_statistic,
    p_value = p_value,
    cells = target_rows,
    fit_object = object
  )

  class(result) <- "cdmc_anticipation_test"
  result
}

print.cdmc_anticipation_test <- function(x, ...) {
  cat("causaldosemc anticipation (leads) diagnostic\n")
  cat(sprintf("  periods tested (lead lengths): %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  matched pre-entry zero-dose cells: %d\n", x$n))
  if (!is.null(x$n_clusters)) {
    cat(sprintf("  unit clusters: %d\n", x$n_clusters))
  }
  cat(sprintf("  mean residual effect: %.6g\n", x$mean_tau))
  if (is.finite(x$standard_error)) {
    se_label <- if (!is.null(x$cluster_standard_error) && is.finite(x$cluster_standard_error)) {
      "cluster-robust standard error"
    } else {
      "standard error"
    }
    cat(sprintf("  %s: %.6g\n", se_label, x$standard_error))
  }
  if (is.finite(x$t_statistic)) {
    cat(sprintf("  t statistic: %.6g\n", x$t_statistic))
  }
  if (is.finite(x$p_value)) {
    cat(sprintf("  p value: %.6g\n", x$p_value))
  }

  invisible(x)
}
