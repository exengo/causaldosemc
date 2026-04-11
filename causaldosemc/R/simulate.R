simulate_cdmc_data <- function(
  n_units = 20L,
  n_times = 12L,
  rank = 2L,
  beta = 1,
  lag_beta = NULL,
  n_covariates = 1L,
  noise_sd = 0.1,
  switch_on_prob = 0.2,
  switch_off_prob = 0.35,
  signed_dose = FALSE,
  negative_dose_prob = 0.5,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!is.logical(signed_dose) || length(signed_dose) != 1L) {
    stop("signed_dose must be a single logical value.", call. = FALSE)
  }

  if (!is.numeric(negative_dose_prob) || length(negative_dose_prob) != 1L ||
      !is.finite(negative_dose_prob) || negative_dose_prob < 0 || negative_dose_prob > 1) {
    stop("negative_dose_prob must be a scalar in [0, 1].", call. = FALSE)
  }

  lag_beta <- if (is.null(lag_beta)) numeric(0) else lag_beta
  units <- paste0("u", seq_len(n_units))
  times <- seq_len(n_times)

  loadings <- matrix(stats::rnorm(n_units * rank), nrow = n_units, ncol = rank)
  factors <- matrix(stats::rnorm(n_times * rank), nrow = n_times, ncol = rank)
  unit_fe <- stats::rnorm(n_units, sd = 0.5)
  time_fe <- stats::rnorm(n_times, sd = 0.5)

  x_matrices <- list()
  gamma <- numeric(0)
  if (n_covariates > 0L) {
    gamma <- seq(0.4, 0.4 * n_covariates, length.out = n_covariates)
    for (index in seq_len(n_covariates)) {
      x_matrices[[paste0("x", index)]] <- matrix(
        stats::rnorm(n_units * n_times),
        nrow = n_units,
        ncol = n_times
      )
    }
  }

  dose_matrix <- matrix(0, nrow = n_units, ncol = n_times)
  for (unit_index in seq_len(n_units)) {
    exposed <- FALSE
    for (time_index in seq_len(n_times)) {
      if (!exposed && stats::runif(1) < switch_on_prob) {
        exposed <- TRUE
      } else if (exposed && stats::runif(1) < switch_off_prob) {
        exposed <- FALSE
      }

      if (exposed) {
        dose_value <- stats::runif(1, min = 0.25, max = 2)
        if (signed_dose) {
          dose_value <- dose_value * if (stats::runif(1) < negative_dose_prob) -1 else 1
        }
        dose_matrix[unit_index, time_index] <- dose_value
      }
    }
  }

  baseline_matrix <- matrix(unit_fe, nrow = n_units, ncol = n_times) +
    matrix(time_fe, nrow = n_units, ncol = n_times, byrow = TRUE) +
    loadings %*% t(factors)

  if (length(x_matrices) > 0L) {
    for (index in seq_along(x_matrices)) {
      baseline_matrix <- baseline_matrix + gamma[index] * x_matrices[[index]]
    }
  }

  lag_coefficients <- c(beta, lag_beta)
  lag_array <- cdmc_build_lagged_doses(dose_matrix, lag_order = length(lag_beta))
  effect_matrix <- matrix(0, nrow = n_units, ncol = n_times)
  for (index in seq_along(lag_coefficients)) {
    lag_slice <- lag_array[, , index]
    lag_slice[is.na(lag_slice)] <- 0
    effect_matrix <- effect_matrix + lag_coefficients[index] * lag_slice
  }

  outcome_matrix <- baseline_matrix + effect_matrix +
    matrix(stats::rnorm(n_units * n_times, sd = noise_sd), nrow = n_units, ncol = n_times)

  output <- data.frame(
    unit = rep(units, each = n_times),
    time = rep(times, times = n_units),
    y = cdmc_flatten_matrix(outcome_matrix),
    dose = cdmc_flatten_matrix(dose_matrix)
  )

  if (length(x_matrices) > 0L) {
    for (name in names(x_matrices)) {
      output[[name]] <- cdmc_flatten_matrix(x_matrices[[name]])
    }
  }

  output
}
