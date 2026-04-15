cdmc_prepare_effect_history <- function(
  y_matrix,
  baseline_hat,
  dose_matrix,
  lag_order = 0L,
  weight_matrix = NULL
) {
  tau_matrix <- y_matrix - baseline_hat
  lag_array <- cdmc_build_lagged_doses(dose_matrix, lag_order = lag_order)

  valid_history <- apply(!is.na(lag_array), c(1, 2), all)
  exposure_magnitude <- cdmc_lag_exposure_magnitude(lag_array)

  effect_sample_mask <- valid_history & !is.na(tau_matrix) & exposure_magnitude > 0
  if (!is.null(weight_matrix)) {
    effect_sample_mask <- effect_sample_mask & weight_matrix > 0
  }

  if (!any(effect_sample_mask)) {
    return(list(
      tau_matrix = tau_matrix,
      lag_array = lag_array,
      sample_mask = effect_sample_mask,
      sample_indices = NULL,
      history = NULL,
      weights = NULL
    ))
  }

  sample_indices <- which(effect_sample_mask, arr.ind = TRUE)
  history <- data.frame(
    row = sample_indices[, 1L],
    col = sample_indices[, 2L],
    tau = tau_matrix[sample_indices],
    stringsAsFactors = FALSE
  )

  for (index in seq_len(dim(lag_array)[3L])) {
    lag_name <- dimnames(lag_array)[[3L]][index]
    history[[lag_name]] <- lag_array[, , index][sample_indices]
  }

  list(
    tau_matrix = tau_matrix,
    lag_array = lag_array,
    sample_mask = effect_sample_mask,
    sample_indices = sample_indices,
    history = history,
    weights = if (is.null(weight_matrix)) NULL else weight_matrix[sample_indices]
  )
}

cdmc_fit_linear_effect <- function(
  y_matrix,
  baseline_hat,
  dose_matrix,
  lag_order = 0L,
  weight_matrix = NULL
) {
  prepared <- cdmc_prepare_effect_history(
    y_matrix = y_matrix,
    baseline_hat = baseline_hat,
    dose_matrix = dose_matrix,
    lag_order = lag_order,
    weight_matrix = weight_matrix
  )
  fitted_matrix <- matrix(0, nrow = nrow(y_matrix), ncol = ncol(y_matrix))

  if (!any(prepared$sample_mask)) {
    return(list(
      coefficients = numeric(0),
      fitted = fitted_matrix,
      tau = prepared$tau_matrix,
      sample_mask = prepared$sample_mask,
      lag_order = lag_order,
      design_columns = character(0),
      model = "linear",
      design_info = NULL,
      weights = NULL,
      residuals = numeric(0)
    ))
  }

  design <- as.matrix(prepared$history[, dimnames(prepared$lag_array)[[3L]], drop = FALSE])
  colnames(design) <- dimnames(prepared$lag_array)[[3L]]

  fit <- if (is.null(prepared$weights)) {
    stats::lm.fit(x = design, y = prepared$history$tau)
  } else {
    stats::lm.wfit(x = design, y = prepared$history$tau, w = prepared$weights)
  }
  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- colnames(design)

  fitted_values <- as.vector(design %*% coefficients)
  fitted_matrix[prepared$sample_indices] <- fitted_values

  list(
    coefficients = coefficients,
    fitted = fitted_matrix,
    tau = prepared$tau_matrix,
    sample_mask = prepared$sample_mask,
    lag_order = lag_order,
    design_columns = colnames(design),
    model = "linear",
    design_info = list(
      design = design,
      basis_spec = NULL,
      column_names = colnames(design),
      model = "linear",
      lag_names = colnames(design),
      df = NULL
    ),
    weights = prepared$weights,
    residuals = prepared$history$tau - fitted_values
  )
}

cdmc_fit_spline_effect <- function(
  y_matrix,
  baseline_hat,
  dose_matrix,
  lag_order = 0L,
  df = 4L,
  weight_matrix = NULL
) {
  prepared <- cdmc_prepare_effect_history(
    y_matrix = y_matrix,
    baseline_hat = baseline_hat,
    dose_matrix = dose_matrix,
    lag_order = lag_order,
    weight_matrix = weight_matrix
  )
  fitted_matrix <- matrix(0, nrow = nrow(y_matrix), ncol = ncol(y_matrix))

  if (!any(prepared$sample_mask)) {
    return(list(
      coefficients = numeric(0),
      fitted = fitted_matrix,
      tau = prepared$tau_matrix,
      sample_mask = prepared$sample_mask,
      lag_order = lag_order,
      design_columns = character(0),
      model = "spline",
      design_info = NULL,
      weights = NULL,
      residuals = numeric(0)
    ))
  }

  design_info <- cdmc_build_response_design(
    history = prepared$history,
    model = "spline",
    df = df
  )
  fit <- cdmc_fit_weighted_regression(
    design = design_info$design,
    response = prepared$history$tau,
    weights = prepared$weights
  )
  fitted_matrix[prepared$sample_indices] <- fit$fitted_values

  list(
    coefficients = fit$coefficients,
    fitted = fitted_matrix,
    tau = prepared$tau_matrix,
    sample_mask = prepared$sample_mask,
    lag_order = lag_order,
    design_columns = design_info$column_names,
    model = "spline",
    design_info = design_info,
    weights = prepared$weights,
    residuals = fit$residuals
  )
}
