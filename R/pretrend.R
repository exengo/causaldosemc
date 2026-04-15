cdmc_joint_placebo_period_label <- function(period) {
  if (period < 0L) {
    return(sprintf("m%d", abs(as.integer(period))))
  }
  if (period > 0L) {
    return(sprintf("p%d", as.integer(period)))
  }

  "0"
}

cdmc_validate_joint_placebo_period_block <- function(periods) {
  periods <- sort(unique(as.integer(periods)))
  if (length(periods) == 0L || any(!is.finite(periods)) || any(periods >= 0L)) {
    stop("Placebo windows must contain one or more strictly negative integers.", call. = FALSE)
  }

  periods
}

cdmc_resolve_joint_placebo_windows <- function(periods = -3:-1, placebo_windows = NULL) {
  if (is.null(placebo_windows)) {
    periods <- cdmc_validate_joint_placebo_period_block(periods)
    resolved <- as.list(periods)
    names(resolved) <- paste0(
      "period_",
      vapply(periods, cdmc_joint_placebo_period_label, character(1))
    )
    return(resolved)
  }

  if (!is.list(placebo_windows) || length(placebo_windows) < 1L) {
    stop("placebo_windows must be NULL or a nonempty named list of negative placebo periods.", call. = FALSE)
  }

  resolved <- lapply(placebo_windows, cdmc_validate_joint_placebo_period_block)
  window_names <- names(placebo_windows)
  if (is.null(window_names)) {
    window_names <- rep("", length(resolved))
  }
  window_names[!nzchar(window_names)] <- vapply(resolved[!nzchar(window_names)], function(window_periods) {
    if (length(window_periods) == 1L) {
      return(paste0("period_", cdmc_joint_placebo_period_label(window_periods)))
    }

    paste0("window_", paste(vapply(window_periods, cdmc_joint_placebo_period_label, character(1)), collapse = "_"))
  }, character(1))

  if (anyDuplicated(window_names)) {
    stop("placebo_windows names must be unique.", call. = FALSE)
  }

  names(resolved) <- window_names
  resolved
}

cdmc_joint_placebo_summary_table <- function(placebo_windows, placebo_results, alpha, equivalence_results = NULL) {
  summary_table <- data.frame(
    window_name = names(placebo_windows),
    period = vapply(placebo_windows, function(window_periods) {
      if (length(window_periods) == 1L) window_periods else NA_integer_
    }, integer(1)),
    window_periods = vapply(placebo_windows, function(window_periods) {
      paste(window_periods, collapse = ", ")
    }, character(1)),
    n_periods = vapply(placebo_windows, length, integer(1)),
    n = vapply(placebo_results, `[[`, integer(1), "n"),
    mean_tau = vapply(placebo_results, `[[`, numeric(1), "mean_tau"),
    p_value = vapply(placebo_results, `[[`, numeric(1), "p_value"),
    stringsAsFactors = FALSE
  )

  if (is.null(equivalence_results)) {
    summary_table$passed_window <- is.na(summary_table$p_value) | summary_table$p_value >= alpha
    return(summary_table)
  }

  summary_table$equivalence_p_value <- vapply(equivalence_results, `[[`, numeric(1), "p_value")
  summary_table$equivalent <- vapply(equivalence_results, `[[`, logical(1), "equivalent")
  summary_table$passed_window <- summary_table$equivalent
  summary_table
}

cdmc_joint_placebo_test <- function(
  object,
  periods = -3:-1,
  placebo_windows = NULL,
  alpha = 0.05,
  equivalence_margin = NULL,
  rerun_tuning = FALSE,
  verbose = FALSE
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  placebo_windows <- cdmc_resolve_joint_placebo_windows(
    periods = periods,
    placebo_windows = placebo_windows
  )
  periods <- sort(unique(unlist(placebo_windows, use.names = FALSE)))

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).", call. = FALSE)
  }

  placebo_results <- lapply(placebo_windows, function(window_periods) {
    placebo <- cdmc_placebo_test(
      object,
      periods = window_periods,
      rerun_tuning = rerun_tuning,
      verbose = verbose
    )
    placebo$fit_object <- NULL
    placebo
  })
  names(placebo_results) <- names(placebo_windows)

  summary_table <- cdmc_joint_placebo_summary_table(
    placebo_windows = placebo_windows,
    placebo_results = placebo_results,
    alpha = alpha
  )

  if (is.null(equivalence_margin)) {
    fisher_statistic <- -2 * sum(log(summary_table$p_value))
    joint_p_value <- stats::pchisq(fisher_statistic, df = 2L * nrow(summary_table), lower.tail = FALSE)

    result <- list(
      call = match.call(),
      periods = periods,
      placebo_windows = placebo_windows,
      alpha = alpha,
      equivalence_margin = NULL,
      window_table = summary_table,
      tests = summary_table,
      method = "fisher",
      statistic = fisher_statistic,
      joint_p_value = joint_p_value,
      passed = joint_p_value >= alpha,
      rerun_tuning = isTRUE(rerun_tuning),
      components = placebo_results,
      fit_object = object
    )
  } else {
    equivalence_results <- lapply(placebo_results, function(placebo) {
      cdmc_equivalence_test(placebo, margin = equivalence_margin, alpha = alpha)
    })
    names(equivalence_results) <- names(placebo_windows)
    summary_table <- cdmc_joint_placebo_summary_table(
      placebo_windows = placebo_windows,
      placebo_results = placebo_results,
      alpha = alpha,
      equivalence_results = equivalence_results
    )
    equivalence_p_values <- summary_table$equivalence_p_value

    joint_p_value <- min(1, length(equivalence_p_values) * max(equivalence_p_values))

    result <- list(
      call = match.call(),
      periods = periods,
      placebo_windows = placebo_windows,
      alpha = alpha,
      equivalence_margin = equivalence_margin,
      window_table = summary_table,
      tests = summary_table,
      method = "bonferroni-equivalence",
      statistic = NA_real_,
      joint_p_value = joint_p_value,
      passed = joint_p_value < alpha,
      rerun_tuning = isTRUE(rerun_tuning),
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
  cat(sprintf("  pooled periods: %s\n", paste(x$periods, collapse = ", ")))
  cat(sprintf("  placebo windows: %d\n", nrow(x$tests)))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  cat(sprintf("  rerun tuning in refits: %s\n", if (isTRUE(x$rerun_tuning)) "yes" else "no"))
  if (!is.null(x$equivalence_margin)) {
    cat(sprintf("  equivalence margin: %.6g\n", x$equivalence_margin))
  }
  cat(sprintf("  joint p value: %.6g\n", x$joint_p_value))
  cat(sprintf("  passed: %s\n", if (x$passed) "yes" else "no"))
  if (nrow(x$tests) > 0L) {
    display_columns <- c("window_name", "window_periods", "n", "mean_tau", "p_value", "passed_window")
    if ("equivalence_p_value" %in% names(x$tests)) {
      display_columns <- c("window_name", "window_periods", "n", "mean_tau", "equivalence_p_value", "passed_window")
    }
    cat("  window summary:\n")
    print(x$tests[, display_columns, drop = FALSE], row.names = FALSE)
  }

  invisible(x)
}