cdmc_validate_panel_data <- function(
  data,
  outcome,
  dose,
  unit,
  time,
  covariates,
  zero_tolerance
) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.", call. = FALSE)
  }

  if (!".cdmc_input_index" %in% names(data)) {
    data$.cdmc_input_index <- seq_len(nrow(data))
  }

  required_columns <- c(outcome, dose, unit, time, covariates)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "Missing required columns: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(data[[outcome]])) {
    stop("outcome must be numeric.", call. = FALSE)
  }

  if (!is.numeric(data[[dose]])) {
    stop("dose must be numeric.", call. = FALSE)
  }

  observed_rows <- if (".cdmc_observed" %in% names(data)) {
    if (!is.logical(data$.cdmc_observed)) {
      stop(".cdmc_observed must be logical when supplied.", call. = FALSE)
    }

    !is.na(data$.cdmc_observed) & data$.cdmc_observed
  } else {
    rep(TRUE, nrow(data))
  }

  if (any(is.na(data[observed_rows, required_columns, drop = FALSE]))) {
    stop(
      "Observed rows must have complete outcome, dose, unit, time, and covariate data.",
      call. = FALSE
    )
  }

  unit_time_key <- paste(data[[unit]], data[[time]], sep = "\r")
  if (anyDuplicated(unit_time_key) > 0L) {
    stop(
      "Each unit-time combination must appear exactly once.",
      call. = FALSE
    )
  }

  unit_levels <- unique(data[[unit]][order(data[[unit]])])
  time_levels <- unique(data[[time]][order(data[[time]])])

  ordered_data <- data[
    order(match(data[[unit]], unit_levels), match(data[[time]], time_levels)),
    ,
    drop = FALSE
  ]

  grid_columns <- list(
    rep(unit_levels, each = length(time_levels)),
    rep(time_levels, times = length(unit_levels))
  )
  names(grid_columns) <- c(unit, time)
  grid_data <- as.data.frame(grid_columns, stringsAsFactors = FALSE, check.names = FALSE)

  ordered_key <- paste(ordered_data[[unit]], ordered_data[[time]], sep = "\r")
  grid_key <- paste(grid_data[[unit]], grid_data[[time]], sep = "\r")
  matched_index <- match(grid_key, ordered_key)
  matched_columns <- ordered_data[matched_index, setdiff(names(ordered_data), c(unit, time)), drop = FALSE]
  full_data <- cbind(grid_data, matched_columns)
  full_data$.cdmc_observed <- !is.na(full_data$.cdmc_input_index)

  observed_full_rows <- full_data$.cdmc_observed
  full_data[[dose]][observed_full_rows & abs(full_data[[dose]]) <= zero_tolerance] <- 0

  list(
    data = full_data,
    unit_levels = unit_levels,
    time_levels = time_levels
  )
}

cdmc_zero_dose_mask <- function(dose_matrix, zero_tolerance) {
  !is.na(dose_matrix) & abs(dose_matrix) <= zero_tolerance
}

cdmc_active_dose_mask <- function(dose_matrix, zero_tolerance) {
  !is.na(dose_matrix) & !cdmc_zero_dose_mask(dose_matrix, zero_tolerance = zero_tolerance)
}

cdmc_lag_exposure_magnitude <- function(lag_array) {
  n_units <- dim(lag_array)[1L]
  n_times <- dim(lag_array)[2L]
  lag_values <- lag_array
  lag_values[is.na(lag_values)] <- 0

  exposure_magnitude <- rowSums(abs(matrix(
    lag_values,
    nrow = n_units * n_times,
    ncol = dim(lag_array)[3L]
  )))

  matrix(
    exposure_magnitude,
    nrow = n_units,
    ncol = n_times,
    dimnames = dimnames(lag_array)[1:2]
  )
}

cdmc_prepare_panel <- function(
  data,
  outcome,
  dose,
  unit,
  time,
  covariates,
  zero_tolerance
) {
  validated <- cdmc_validate_panel_data(
    data = data,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    zero_tolerance = zero_tolerance
  )

  ordered_data <- validated$data
  n_units <- length(validated$unit_levels)
  n_times <- length(validated$time_levels)

  ordered_data$.cdmc_unit_index <- rep(seq_len(n_units), each = n_times)
  ordered_data$.cdmc_time_index <- rep(seq_len(n_times), times = n_units)

  y_matrix <- matrix(
    ordered_data[[outcome]],
    nrow = n_units,
    ncol = n_times,
    byrow = TRUE
  )

  dose_matrix <- matrix(
    ordered_data[[dose]],
    nrow = n_units,
    ncol = n_times,
    byrow = TRUE
  )

  x_matrices <- list()
  if (length(covariates) > 0L) {
    for (covariate in covariates) {
      x_matrices[[covariate]] <- matrix(
        ordered_data[[covariate]],
        nrow = n_units,
        ncol = n_times,
        byrow = TRUE
      )
    }
  }

  observed_mask <- matrix(
    ordered_data$.cdmc_observed,
    nrow = n_units,
    ncol = n_times,
    byrow = TRUE
  )

  list(
    data = ordered_data,
    y_matrix = y_matrix,
    dose_matrix = dose_matrix,
    x_matrices = x_matrices,
    observed_mask = observed_mask,
    n_units = n_units,
    n_times = n_times,
    unit_levels = validated$unit_levels,
    time_levels = validated$time_levels
  )
}

cdmc_build_eligible_mask <- function(dose_matrix, zero_tolerance, washout = 0L) {
  eligible_mask <- cdmc_zero_dose_mask(dose_matrix, zero_tolerance = zero_tolerance)

  if (washout <= 0L) {
    return(eligible_mask)
  }

  n_times <- ncol(dose_matrix)
  for (unit_index in seq_len(nrow(dose_matrix))) {
    treated_times <- which(cdmc_active_dose_mask(dose_matrix[unit_index, , drop = FALSE], zero_tolerance = zero_tolerance))
    if (length(treated_times) == 0L) {
      next
    }

    for (treated_time in treated_times) {
      upper <- min(n_times, treated_time + washout)
      if (treated_time < upper) {
        eligible_mask[unit_index, seq.int(treated_time + 1L, upper)] <- FALSE
      }
    }
  }

  eligible_mask & cdmc_zero_dose_mask(dose_matrix, zero_tolerance = zero_tolerance)
}

cdmc_validate_mask_support <- function(mask, unit_message, time_message, size_message) {
  if (any(rowSums(mask) == 0L)) {
    stop(unit_message, call. = FALSE)
  }

  if (any(colSums(mask) == 0L)) {
    stop(time_message, call. = FALSE)
  }

  if (sum(mask) <= (nrow(mask) + ncol(mask))) {
    stop(size_message, call. = FALSE)
  }
}

cdmc_validate_support <- function(eligible_mask) {
  cdmc_validate_mask_support(
    mask = eligible_mask,
    unit_message = "Each unit must retain at least one eligible zero-dose observation after washout masking.",
    time_message = "Each time period must retain at least one eligible zero-dose observation after washout masking.",
    size_message = "The eligible zero-dose sample is too small to support two-way effects plus low-rank estimation."
  )
}

cdmc_validate_joint_support <- function(fit_mask) {
  cdmc_validate_mask_support(
    mask = fit_mask,
    unit_message = "Each unit must retain at least one observed estimation cell under the joint objective.",
    time_message = "Each time period must retain at least one observed estimation cell under the joint objective.",
    size_message = "The observed estimation sample is too small to support two-way effects plus low-rank estimation under the joint objective."
  )
}

cdmc_build_lagged_doses <- function(dose_matrix, lag_order = 0L) {
  lag_order <- as.integer(lag_order)
  n_units <- nrow(dose_matrix)
  n_times <- ncol(dose_matrix)
  lag_names <- paste0("dose_lag", seq.int(0L, lag_order))

  lag_array <- array(
    NA_real_,
    dim = c(n_units, n_times, lag_order + 1L),
    dimnames = list(NULL, NULL, lag_names)
  )

  lag_array[, , 1L] <- dose_matrix
  if (lag_order == 0L) {
    return(lag_array)
  }

  for (lag_index in seq_len(lag_order)) {
    lag_array[, seq.int(lag_index + 1L, n_times), lag_index + 1L] <-
      dose_matrix[, seq_len(n_times - lag_index), drop = FALSE]
  }

  lag_array
}
