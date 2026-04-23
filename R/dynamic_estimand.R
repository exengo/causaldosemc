cdmc_dynamic_lag_names <- function(object) {
  if (inherits(object, "cdmc_dose_response")) {
    return(object$design_info$lag_names)
  }

  if (inherits(object, "cdmc_fit")) {
    if (identical(object$effect_model, "none")) {
      stop(
        "cdmc_dynamic_estimand() requires a cdmc_fit object with effect_model = 'linear' or 'spline', or a cdmc_dose_response() object.",
        call. = FALSE
      )
    }

    return(object$effect$design_info$lag_names %||% paste0("dose_lag", seq.int(0L, object$lag_order)))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    return(cdmc_dr_design_columns(object))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
    call. = FALSE
  )
}

cdmc_dynamic_history_template <- function(lag_names, n_rows = 1L) {
  cdmc_named_numeric_data_frame(lag_names, n_rows = n_rows)
}

cdmc_resolve_dynamic_history <- function(history, lag_names, argument_name = "history") {
  if (is.null(history)) {
    stop(sprintf("%s must be supplied.", argument_name), call. = FALSE)
  }

  if (is.matrix(history)) {
    history <- as.data.frame(history, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (is.data.frame(history)) {
    if (nrow(history) < 1L) {
      stop(sprintf("%s must contain at least one row.", argument_name), call. = FALSE)
    }

    missing_columns <- setdiff(lag_names, names(history))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf("%s is missing required lag columns: %s.", argument_name, paste(missing_columns, collapse = ", ")),
        call. = FALSE
      )
    }

    return(history[, lag_names, drop = FALSE])
  }

  if (!is.numeric(history)) {
    stop(
      sprintf("%s must be numeric, a matrix, or a data frame.", argument_name),
      call. = FALSE
    )
  }

  resolved <- cdmc_dynamic_history_template(lag_names)

  if (length(history) == 1L && is.null(names(history))) {
    resolved[[lag_names[[1L]]]] <- as.numeric(history)
    return(resolved)
  }

  if (length(history) == length(lag_names) && is.null(names(history))) {
    resolved[1L, lag_names] <- as.numeric(history)
    return(resolved)
  }

  if (is.null(names(history))) {
    stop(
      sprintf("Unnamed numeric %s inputs must have length 1 or %d.", argument_name, length(lag_names)),
      call. = FALSE
    )
  }

  unknown_names <- setdiff(names(history), lag_names)
  if (length(unknown_names) > 0L) {
    stop(
      sprintf("Unknown lag names in %s: %s.", argument_name, paste(unknown_names, collapse = ", ")),
      call. = FALSE
    )
  }

  resolved[1L, names(history)] <- as.numeric(history)
  resolved
}

cdmc_align_dynamic_histories <- function(history, reference_history) {
  n_paths <- max(nrow(history), nrow(reference_history))

  if (nrow(history) == 1L && n_paths > 1L) {
    history <- history[rep(1L, n_paths), , drop = FALSE]
  }
  if (nrow(reference_history) == 1L && n_paths > 1L) {
    reference_history <- reference_history[rep(1L, n_paths), , drop = FALSE]
  }

  if (nrow(history) != nrow(reference_history)) {
    stop(
      "history and reference_history must have the same number of rows, unless one of them has exactly one row.",
      call. = FALSE
    )
  }

  list(history = history, reference_history = reference_history)
}

cdmc_resolve_dynamic_path_weights <- function(path_weights, n_paths) {
  if (is.null(path_weights)) {
    return(rep(1, n_paths))
  }

  if (!is.numeric(path_weights) || any(!is.finite(path_weights)) || any(path_weights <= 0)) {
    stop("path_weights must contain finite positive values when supplied.", call. = FALSE)
  }

  if (length(path_weights) == 1L) {
    return(rep(as.numeric(path_weights), n_paths))
  }

  if (length(path_weights) != n_paths) {
    stop("path_weights must have length 1 or match the number of requested paths.", call. = FALSE)
  }

  as.numeric(path_weights)
}

cdmc_resolve_dynamic_labels <- function(history, labels = NULL) {
  if (is.null(labels)) {
    labels <- rownames(history)
    if (!is.null(labels) && identical(labels, as.character(seq_len(nrow(history))))) {
      labels <- NULL
    }
  }

  if (is.null(labels) || length(labels) != nrow(history)) {
    labels <- paste0("path", seq_len(nrow(history)))
  }

  labels <- as.character(labels)
  labels[is.na(labels) | labels == ""] <- paste0("path", which(is.na(labels) | labels == ""))
  make.unique(labels)
}

cdmc_prepare_dynamic_design_history <- function(history) {
  transform(
    history,
    row = seq_len(nrow(history)),
    col = seq_len(nrow(history)),
    tau = 0
  )
}

cdmc_predict_fit_dynamic_response <- function(object, history) {
  if (identical(object$effect_model, "none")) {
    stop(
      "cdmc_dynamic_estimand() is not available for cdmc_fit objects with effect_model = 'none'.",
      call. = FALSE
    )
  }

  if (identical(object$effect_model, "spline") && is.null(object$effect$design_info)) {
    stop(
      "Spline dynamic prediction is unavailable because the fitted cdmc_fit object did not retain a spline effect design.",
      call. = FALSE
    )
  }

  design_info <- cdmc_build_response_design(
    history = cdmc_prepare_dynamic_design_history(history),
    model = object$effect_model,
    df = object$effect$design_info$df %||% object$effect_df %||% 4L,
    basis_spec = object$effect$design_info$basis_spec
  )
  coefficients <- numeric(length(design_info$column_names))
  names(coefficients) <- design_info$column_names

  if (length(object$effect$coefficients) > 0L) {
    common_columns <- intersect(names(object$effect$coefficients), names(coefficients))
    coefficients[common_columns] <- object$effect$coefficients[common_columns]
  }

  as.vector(design_info$design %*% coefficients)
}

cdmc_predict_dr_dynamic_response <- function(object, history) {
  coefficient_vector <- cdmc_dr_coefficient_vector(object)
  as.vector(as.matrix(history[, names(coefficient_vector), drop = FALSE]) %*% coefficient_vector)
}

cdmc_predict_dynamic_response <- function(object, history) {
  if (inherits(object, "cdmc_dose_response")) {
    return(stats::predict(object, history = history, type = "response")$estimate)
  }

  if (inherits(object, "cdmc_fit")) {
    return(cdmc_predict_fit_dynamic_response(object, history = history))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    return(cdmc_predict_dr_dynamic_response(object, history = history))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
    call. = FALSE
  )
}

cdmc_predict_fit_dynamic_slope <- function(object, history, eps) {
  if (identical(object$effect_model, "linear")) {
    coefficient <- unname(object$effect$coefficients[["dose_lag0"]] %||% 0)
    return(rep(coefficient, nrow(history)))
  }

  history_up <- history
  history_down <- history
  history_up$dose_lag0 <- history_up$dose_lag0 + eps
  history_down$dose_lag0 <- history_down$dose_lag0 - eps

  response_up <- cdmc_predict_fit_dynamic_response(object, history = history_up)
  response_down <- cdmc_predict_fit_dynamic_response(object, history = history_down)
  (response_up - response_down) / (2 * eps)
}

cdmc_predict_dynamic_slope <- function(object, history, eps) {
  if (inherits(object, "cdmc_dose_response")) {
    return(stats::predict(object, history = history, type = "slope", eps = eps)$estimate)
  }

  if (inherits(object, "cdmc_fit")) {
    return(cdmc_predict_fit_dynamic_slope(object, history = history, eps = eps))
  }

  if (inherits(object, "cdmc_dr_fit")) {
    coefficient <- unname(cdmc_dr_coefficient_vector(object)[["dose_lag0"]] %||% 0)
    return(rep(coefficient, nrow(history)))
  }

  stop(
    "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
    call. = FALSE
  )
}

cdmc_dynamic_stat_prefix <- function(type) {
  switch(
    type,
    response = "dynamic_response",
    contrast = "dynamic_contrast",
    slope = "dynamic_slope"
  )
}

cdmc_dynamic_estimate_vector <- function(estimate_table, type, aggregate, path_weights) {
  prefix <- cdmc_dynamic_stat_prefix(type)
  path_names <- paste0(prefix, "_", make.names(estimate_table$label, unique = TRUE))
  individual <- estimate_table$estimate
  names(individual) <- path_names

  if (identical(aggregate, "individual")) {
    return(individual)
  }

  mean_estimate <- stats::weighted.mean(estimate_table$estimate, w = path_weights)
  mean_name <- paste0("mean_", prefix)

  if (identical(aggregate, "mean")) {
    return(stats::setNames(mean_estimate, mean_name))
  }

  c(individual, stats::setNames(mean_estimate, mean_name))
}

cdmc_dynamic_estimand <- function(
  object,
  history,
  reference_history = NULL,
  type = c("response", "contrast", "slope"),
  aggregate = c("individual", "mean", "both"),
  path_weights = NULL,
  eps = 1e-4,
  labels = NULL
) {
  if (!inherits(object, "cdmc_fit") &&
      !inherits(object, "cdmc_dr_fit") &&
      !inherits(object, "cdmc_dose_response")) {
    stop(
      "object must inherit from 'cdmc_fit', 'cdmc_dr_fit', or 'cdmc_dose_response'.",
      call. = FALSE
    )
  }

  type <- match.arg(type)
  aggregate <- match.arg(aggregate)

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("eps must be a single positive numeric value.", call. = FALSE)
  }

  lag_names <- cdmc_dynamic_lag_names(object)
  history <- cdmc_resolve_dynamic_history(history, lag_names = lag_names, argument_name = "history")
  path_labels <- cdmc_resolve_dynamic_labels(history, labels = labels)

  estimate_table <- data.frame(
    label = path_labels,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(type, "contrast")) {
    reference_history <- cdmc_resolve_dynamic_history(
      reference_history %||% cdmc_dynamic_history_template(lag_names, n_rows = 1L),
      lag_names = lag_names,
      argument_name = "reference_history"
    )
    aligned <- cdmc_align_dynamic_histories(history, reference_history)
    history <- aligned$history
    reference_history <- aligned$reference_history
    estimate_table$effect <- cdmc_predict_dynamic_response(object, history = history)
    estimate_table$reference_effect <- cdmc_predict_dynamic_response(object, history = reference_history)
    estimate_table$estimate <- estimate_table$effect - estimate_table$reference_effect
  } else if (identical(type, "response")) {
    if (!is.null(reference_history)) {
      stop("reference_history is only supported when type = 'contrast'.", call. = FALSE)
    }
    estimate_table$estimate <- cdmc_predict_dynamic_response(object, history = history)
  } else {
    if (!is.null(reference_history)) {
      stop("reference_history is only supported when type = 'contrast'.", call. = FALSE)
    }
    estimate_table$estimate <- cdmc_predict_dynamic_slope(object, history = history, eps = eps)
  }

  path_weights <- cdmc_resolve_dynamic_path_weights(path_weights, n_paths = nrow(history))
  estimate_table$path_weight <- path_weights
  estimates <- cdmc_dynamic_estimate_vector(
    estimate_table = estimate_table,
    type = type,
    aggregate = aggregate,
    path_weights = path_weights
  )

  result <- list(
    call = match.call(),
    source_object = object,
    source_class = class(object)[1L],
    type = type,
    aggregate = aggregate,
    history = history,
    reference_history = reference_history,
    path_weights = path_weights,
    labels = path_labels,
    lag_names = lag_names,
    eps = eps,
    estimate_table = estimate_table,
    estimates = estimates
  )

  class(result) <- "cdmc_dynamic_estimand"
  result
}

print.cdmc_dynamic_estimand <- function(x, ...) {
  cat("causaldosemc dynamic estimand\n")
  cat(sprintf("  source: %s\n", x$source_class))
  cat(sprintf("  type: %s\n", x$type))
  cat(sprintf("  paths: %d\n", nrow(x$history)))
  cat(sprintf("  aggregate: %s\n", x$aggregate))

  print(x$estimate_table, row.names = FALSE)

  invisible(x)
}

summary.cdmc_dynamic_estimand <- function(object, ...) {
  object$estimate_table
}