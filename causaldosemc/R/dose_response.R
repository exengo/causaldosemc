cdmc_prepare_effect_sample <- function(
  object,
  lag_order = object$lag_order,
  include_zero_dose = TRUE
) {
  lag_array <- cdmc_build_lagged_doses(object$dose_matrix, lag_order = lag_order)
  valid_history <- apply(!is.na(lag_array), c(1, 2), all) & !is.na(object$effect$tau)
  exposure_magnitude <- cdmc_lag_exposure_magnitude(lag_array)

  sample_mask <- if (include_zero_dose) {
    valid_history
  } else {
    valid_history & exposure_magnitude > 0
  }

  sample_indices <- which(sample_mask, arr.ind = TRUE)
  history <- data.frame(
    row = sample_indices[, 1L],
    col = sample_indices[, 2L],
    tau = object$effect$tau[sample_indices],
    stringsAsFactors = FALSE
  )

  for (index in seq_len(dim(lag_array)[3L])) {
    lag_name <- dimnames(lag_array)[[3L]][index]
    history[[lag_name]] <- lag_array[, , index][sample_indices]
  }

  list(
    lag_order = lag_order,
    lag_names = dimnames(lag_array)[[3L]],
    sample_mask = sample_mask,
    history = history
  )
}

cdmc_resolve_effect_weights <- function(weights, object, sample_mask) {
  if (is.null(weights) && !is.null(object$fit_control$weights)) {
    weights <- object$fit_control$weights
  }

  if (is.null(weights)) {
    return(NULL)
  }

  resolved <- cdmc_resolve_weight_vector(
    weights = weights,
    data = object$data,
    n_units = nrow(object$dose_matrix),
    n_times = ncol(object$dose_matrix),
    argument_name = "weights"
  )

  resolved[cdmc_flatten_matrix(sample_mask)]
}

cdmc_build_spline_basis_spec <- function(x, df) {
  boundary_knots <- range(x)
  interior_candidates <- sort(unique(x[x > boundary_knots[1L] & x < boundary_knots[2L]]))
  n_knots <- min(max(df - 1L, 0L), length(interior_candidates))

  if (n_knots < 1L) {
    return(list(
      type = "linear",
      zero_basis = 0,
      ncol = 1L
    ))
  }

  knot_probs <- seq(0, 1, length.out = n_knots + 2L)[seq_len(n_knots) + 1L]
  knots <- as.numeric(stats::quantile(interior_candidates, probs = knot_probs, names = FALSE, type = 1L))
  knots <- unique(knots[knots > boundary_knots[1L] & knots < boundary_knots[2L]])

  if (length(knots) < 1L) {
    return(list(
      type = "linear",
      zero_basis = 0,
      ncol = 1L
    ))
  }

  zero_basis <- as.numeric(splines::ns(0, knots = knots, Boundary.knots = boundary_knots))

  list(
    type = "ns",
    knots = knots,
    boundary_knots = boundary_knots,
    zero_basis = zero_basis,
    ncol = length(zero_basis)
  )
}

cdmc_apply_spline_basis <- function(x, basis_spec) {
  if (identical(basis_spec$type, "linear")) {
    transformed <- matrix(x, ncol = 1L)
  } else {
    transformed <- splines::ns(
      x,
      knots = basis_spec$knots,
      Boundary.knots = basis_spec$boundary_knots
    )
  }

  sweep(transformed, 2, basis_spec$zero_basis, FUN = "-")
}

cdmc_build_response_design <- function(history, model = c("linear", "spline"), df = 4L, basis_spec = NULL) {
  model <- match.arg(model)
  lag_names <- setdiff(names(history), c("row", "col", "tau"))

  if (model == "linear") {
    design <- as.matrix(history[, lag_names, drop = FALSE])
    return(list(
      design = design,
      basis_spec = NULL,
      column_names = colnames(design),
      model = model,
      lag_names = lag_names
    ))
  }

  df <- as.integer(df)
  if (df < 1L) {
    stop("df must be a positive integer for spline dose-response models.", call. = FALSE)
  }

  basis_spec <- basis_spec %||% vector("list", length(lag_names))
  names(basis_spec) <- lag_names
  design_parts <- vector("list", length(lag_names))
  column_names <- character(0)

  for (index in seq_along(lag_names)) {
    lag_name <- lag_names[[index]]
    x <- history[[lag_name]]

    if (is.null(basis_spec[[lag_name]])) {
      basis_spec[[lag_name]] <- cdmc_build_spline_basis_spec(x, df = df)
    }

    centered <- cdmc_apply_spline_basis(x, basis_spec[[lag_name]])
    current_names <- paste0(lag_name, "_ns", seq_len(ncol(centered)))
    colnames(centered) <- current_names

    design_parts[[index]] <- centered
    column_names <- c(column_names, current_names)
  }

  design <- do.call(cbind, design_parts)

  list(
    design = design,
    basis_spec = basis_spec,
    column_names = column_names,
    model = model,
    lag_names = lag_names,
    df = df
  )
}

cdmc_fit_weighted_regression <- function(design, response, weights = NULL) {
  if (is.null(weights)) {
    fit <- stats::lm.fit(x = design, y = response)
  } else {
    fit <- stats::lm.wfit(x = design, y = response, w = weights)
  }

  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- colnames(design)
  fitted_values <- as.vector(design %*% coefficients)

  list(
    coefficients = coefficients,
    fitted_values = fitted_values,
    residuals = response - fitted_values,
    weights = weights
  )
}

cdmc_dose_response <- function(
  object,
  model = c("linear", "spline"),
  lag_order = object$lag_order,
  df = 4L,
  include_zero_dose = TRUE,
  weights = NULL
) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  model <- match.arg(model)
  prepared <- cdmc_prepare_effect_sample(
    object = object,
    lag_order = lag_order,
    include_zero_dose = include_zero_dose
  )
  observation_weights <- cdmc_resolve_effect_weights(
    weights = weights,
    object = object,
    sample_mask = prepared$sample_mask
  )
  design_info <- cdmc_build_response_design(
    history = prepared$history,
    model = model,
    df = df
  )
  fit <- cdmc_fit_weighted_regression(
    design = design_info$design,
    response = prepared$history$tau,
    weights = observation_weights
  )

  result <- list(
    call = match.call(),
    model = model,
    lag_order = lag_order,
    include_zero_dose = include_zero_dose,
    weights = weights,
    fit_object = object,
    history = prepared$history,
    design_info = design_info,
    coefficients = fit$coefficients,
    fitted_values = fit$fitted_values,
    residuals = fit$residuals
  )

  class(result) <- "cdmc_dose_response"
  result
}

cdmc_prepare_prediction_history <- function(object, dose = NULL, history = NULL) {
  lag_names <- object$design_info$lag_names

  if (!is.null(history)) {
    history <- as.data.frame(history)
    missing_columns <- setdiff(lag_names, names(history))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf("history is missing required lag columns: %s.", paste(missing_columns, collapse = ", ")),
        call. = FALSE
      )
    }
    return(history[, lag_names, drop = FALSE])
  }

  if (is.null(dose)) {
    stop("Provide either dose or history when predicting from a dose-response fit.", call. = FALSE)
  }

  prediction_history <- as.data.frame(matrix(0, nrow = length(dose), ncol = length(lag_names)))
  names(prediction_history) <- lag_names
  prediction_history[["dose_lag0"]] <- dose
  prediction_history
}

predict.cdmc_dose_response <- function(object, dose = NULL, history = NULL, type = c("response", "slope"), eps = 1e-4, ...) {
  type <- match.arg(type)
  prediction_history <- cdmc_prepare_prediction_history(object, dose = dose, history = history)
  design_info <- cdmc_build_response_design(
    history = transform(prediction_history, row = seq_len(nrow(prediction_history)), col = seq_len(nrow(prediction_history)), tau = 0),
    model = object$model,
    df = object$design_info$df %||% 1L,
    basis_spec = object$design_info$basis_spec
  )

  response <- as.vector(design_info$design %*% object$coefficients)
  if (type == "response") {
    output <- prediction_history
    output$estimate <- response
    return(output)
  }

  history_up <- prediction_history
  history_down <- prediction_history
  history_up$dose_lag0 <- history_up$dose_lag0 + eps
  history_down$dose_lag0 <- history_down$dose_lag0 - eps

  design_up <- cdmc_build_response_design(
    history = transform(history_up, row = seq_len(nrow(history_up)), col = seq_len(nrow(history_up)), tau = 0),
    model = object$model,
    df = object$design_info$df %||% 1L,
    basis_spec = object$design_info$basis_spec
  )
  design_down <- cdmc_build_response_design(
    history = transform(history_down, row = seq_len(nrow(history_down)), col = seq_len(nrow(history_down)), tau = 0),
    model = object$model,
    df = object$design_info$df %||% 1L,
    basis_spec = object$design_info$basis_spec
  )

  estimate_up <- as.vector(design_up$design %*% object$coefficients)
  estimate_down <- as.vector(design_down$design %*% object$coefficients)
  slope <- (estimate_up - estimate_down) / (history_up$dose_lag0 - history_down$dose_lag0)

  output <- prediction_history
  output$estimate <- slope
  output
}

print.cdmc_dose_response <- function(x, ...) {
  cat("causaldosemc dose-response fit\n")
  cat(sprintf("  model: %s\n", x$model))
  cat(sprintf("  lag order: %d\n", x$lag_order))
  cat(sprintf("  sample size: %d\n", nrow(x$history)))
  cat(sprintf("  include zero-dose sample: %s\n", if (x$include_zero_dose) "yes" else "no"))
  if (length(x$coefficients) > 0L) {
    cat(sprintf("  coefficient count: %d\n", length(x$coefficients)))
  }

  invisible(x)
}
