cdmc_extract_diagnostic_tau <- function(object) {
  if (inherits(object, "cdmc_placebo_test") || inherits(object, "cdmc_carryover_refit_test")) {
    return(object$cells$pseudo_tau)
  }

  if (inherits(object, "cdmc_carryover_test")) {
    return(object$cells$tau)
  }

  stop(
    "object must inherit from 'cdmc_placebo_test', 'cdmc_carryover_refit_test', or 'cdmc_carryover_test'.",
    call. = FALSE
  )
}

cdmc_equivalence_test <- function(object, margin, alpha = 0.05) {
  if (!is.numeric(margin) || length(margin) != 1L || !is.finite(margin) || margin <= 0) {
    stop("margin must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).", call. = FALSE)
  }

  tau_values <- cdmc_extract_diagnostic_tau(object)
  sample_size <- length(tau_values)
  mean_tau <- mean(tau_values)
  sd_tau <- if (sample_size > 1L) stats::sd(tau_values) else 0
  standard_error <- if (sample_size > 1L) sd_tau / sqrt(sample_size) else NA_real_
  degrees_freedom <- max(sample_size - 1L, 0L)

  if (sample_size < 2L || !is.finite(standard_error) || standard_error <= 0) {
    lower_p_value <- if (mean_tau > -margin) 0 else 1
    upper_p_value <- if (mean_tau < margin) 0 else 1
    lower_t_statistic <- NA_real_
    upper_t_statistic <- NA_real_
  } else {
    lower_t_statistic <- (mean_tau + margin) / standard_error
    upper_t_statistic <- (mean_tau - margin) / standard_error
    lower_p_value <- 1 - stats::pt(lower_t_statistic, df = degrees_freedom)
    upper_p_value <- stats::pt(upper_t_statistic, df = degrees_freedom)
  }

  result <- list(
    call = match.call(),
    diagnostic_class = class(object)[1L],
    periods = object$periods %||% NULL,
    n = sample_size,
    mean_tau = mean_tau,
    sd_tau = sd_tau,
    standard_error = standard_error,
    margin = margin,
    alpha = alpha,
    lower_t_statistic = lower_t_statistic,
    upper_t_statistic = upper_t_statistic,
    lower_p_value = lower_p_value,
    upper_p_value = upper_p_value,
    p_value = max(lower_p_value, upper_p_value),
    equivalent = max(lower_p_value, upper_p_value) < alpha,
    diagnostic_object = object,
    fit_object = object$fit_object %||% NULL
  )

  class(result) <- "cdmc_equivalence_test"
  result
}

print.cdmc_equivalence_test <- function(x, ...) {
  cat("causaldosemc equivalence diagnostic\n")
  cat(sprintf("  source diagnostic: %s\n", x$diagnostic_class))
  if (!is.null(x$periods)) {
    cat(sprintf("  periods tested: %s\n", paste(x$periods, collapse = ", ")))
  }
  cat(sprintf("  observations: %d\n", x$n))
  cat(sprintf("  mean effect: %.6g\n", x$mean_tau))
  cat(sprintf("  equivalence margin: %.6g\n", x$margin))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  cat(sprintf("  lower one-sided p value: %.6g\n", x$lower_p_value))
  cat(sprintf("  upper one-sided p value: %.6g\n", x$upper_p_value))
  cat(sprintf("  equivalence p value: %.6g\n", x$p_value))
  cat(sprintf("  equivalent: %s\n", if (x$equivalent) "yes" else "no"))

  invisible(x)
}