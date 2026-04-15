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
  standard_error <- if (sample_size > 1L) sd_tau / sqrt(sample_size) else NA_real_
  t_statistic <- if (sample_size > 1L && standard_error > 0) mean_tau / standard_error else NA_real_
  p_value <- if (sample_size > 1L && is.finite(t_statistic)) {
    2 * stats::pt(-abs(t_statistic), df = sample_size - 1L)
  } else {
    NA_real_
  }

  target_rows <- object$data[, c(object$unit, object$time), drop = FALSE]
  target_rows$exit_distance <- cdmc_flatten_matrix(exit_distance)
  target_rows$tau <- cdmc_flatten_matrix(object$effect$tau)
  target_rows$eligible_control <- cdmc_flatten_matrix(object$eligible_mask)
  target_rows <- target_rows[cdmc_flatten_matrix(target_mask), , drop = FALSE]

  result <- list(
    call = match.call(),
    periods = periods,
    n = sample_size,
    mean_tau = mean_tau,
    sd_tau = sd_tau,
    standard_error = standard_error,
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
  cat(sprintf("  mean residual effect: %.6g\n", x$mean_tau))
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
