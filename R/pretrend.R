cdmc_joint_placebo_test <- function(
  object,
  periods = -3:-1,
  alpha = 0.05,
  equivalence_margin = NULL,
  verbose = FALSE
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods >= 0L)) {
    stop("periods must contain one or more strictly negative integers.", call. = FALSE)
  }

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).", call. = FALSE)
  }

  placebo_results <- lapply(periods, function(period) {
    placebo <- cdmc_placebo_test(object, periods = period, verbose = verbose)
    placebo$fit_object <- NULL
    placebo
  })

  summary_table <- data.frame(
    period = periods,
    n = vapply(placebo_results, `[[`, integer(1), "n"),
    mean_tau = vapply(placebo_results, `[[`, numeric(1), "mean_tau"),
    p_value = vapply(placebo_results, `[[`, numeric(1), "p_value"),
    stringsAsFactors = FALSE
  )

  if (is.null(equivalence_margin)) {
    fisher_statistic <- -2 * sum(log(summary_table$p_value))
    joint_p_value <- stats::pchisq(fisher_statistic, df = 2L * nrow(summary_table), lower.tail = FALSE)

    result <- list(
      call = match.call(),
      periods = periods,
      alpha = alpha,
      equivalence_margin = NULL,
      tests = summary_table,
      method = "fisher",
      statistic = fisher_statistic,
      joint_p_value = joint_p_value,
      passed = joint_p_value >= alpha,
      components = placebo_results,
      fit_object = object
    )
  } else {
    equivalence_results <- lapply(placebo_results, function(placebo) {
      cdmc_equivalence_test(placebo, margin = equivalence_margin, alpha = alpha)
    })
    equivalence_p_values <- vapply(equivalence_results, `[[`, numeric(1), "p_value")

    summary_table$equivalence_p_value <- equivalence_p_values
    summary_table$equivalent <- vapply(equivalence_results, `[[`, logical(1), "equivalent")

    joint_p_value <- min(1, length(equivalence_p_values) * max(equivalence_p_values))

    result <- list(
      call = match.call(),
      periods = periods,
      alpha = alpha,
      equivalence_margin = equivalence_margin,
      tests = summary_table,
      method = "bonferroni-equivalence",
      statistic = NA_real_,
      joint_p_value = joint_p_value,
      passed = joint_p_value < alpha,
      components = placebo_results,
      equivalence_components = equivalence_results,
      fit_object = object
    )
  }

  class(result) <- "cdmc_joint_placebo_test"
  result
}

print.cdmc_joint_placebo_test <- function(x, ...) {
  cat("causaldosemc joint placebo diagnostic\n")
  cat(sprintf("  periods tested: %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  if (!is.null(x$equivalence_margin)) {
    cat(sprintf("  equivalence margin: %.6g\n", x$equivalence_margin))
  }
  cat(sprintf("  joint p value: %.6g\n", x$joint_p_value))
  cat(sprintf("  passed: %s\n", if (x$passed) "yes" else "no"))

  invisible(x)
}