test_that("SCIA screen returns a diagnostic object on standard fits", {
  panel <- simulate_cdmc_data(
    n_units = 16,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 7171
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 7171
  )

  diagnostic <- cdmc_scia_test(fit, lags = 1, outcome_proxy = "tau")

  expect_s3_class(diagnostic, "cdmc_scia_test")
  expect_equal(diagnostic$lags, 1)
  expect_true(diagnostic$sample_size > 0)
  expect_true("tau_lag1" %in% diagnostic$screen_table$term)
})

test_that("SCIA screen detects lagged outcome driven treatment assignment", {
  panel <- simulate_cdmc_data(
    n_units = 18,
    n_times = 12,
    rank = 2,
    beta = 0,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0,
    switch_off_prob = 1,
    seed = 7272
  )

  baseline_y <- matrix(panel$y, nrow = 18, ncol = 12, byrow = TRUE)
  dose_matrix <- matrix(0, nrow = 18, ncol = 12)
  for (unit_index in seq_len(nrow(dose_matrix))) {
    for (time_index in 2:ncol(dose_matrix)) {
      latent_dose <- 0.55 * baseline_y[unit_index, time_index - 1L] + stats::rnorm(1, sd = 0.08)
      if (abs(latent_dose) >= 0.5 && stats::runif(1) < 0.55) {
        dose_matrix[unit_index, time_index] <- max(min(latent_dose, 1.8), -1.8)
      }
    }
  }

  for (time_index in seq_len(ncol(dose_matrix))) {
    if (all(abs(dose_matrix[, time_index]) > 0)) {
      dose_matrix[1L, time_index] <- 0
    }
  }

  panel$dose <- as.vector(t(dose_matrix))
  panel$y <- panel$y + 0.8 * panel$dose

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 7272
  )

  diagnostic <- cdmc_scia_test(fit, lags = 1, outcome_proxy = "y")

  expect_s3_class(diagnostic, "cdmc_scia_test")
  expect_true(is.finite(diagnostic$p_value))
  expect_true(diagnostic$p_value < 0.05)
})