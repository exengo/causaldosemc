test_that("CBPS weights can be estimated and reused in dose-response fitting", {
  testthat::skip_if_not_installed("CBPS")

  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 707
  )

  cbps_weights <- cdmc_cbps_weights(
    data = panel,
    dose = "dose",
    covariates = "x1",
    iterations = 250L
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = cbps_weights,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 707
  )

  dr <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1
  )

  expect_s3_class(cbps_weights, "cdmc_cbps_weights")
  expect_length(cbps_weights$weights, nrow(panel))
  expect_true(all(is.finite(cbps_weights$weights)))
  expect_true(fit$weight_supplied)
  expect_s3_class(dr, "cdmc_dose_response")
})

test_that("entropy-balance weights can be estimated and reused in dose-response fitting", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 717
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.1 + panel$x1[active]^2 + ifelse(panel$time[active] > stats::median(panel$time), 0.6, -0.1)
  )

  entropy_weights <- cdmc_entropy_balance_weights(
    data = panel,
    dose = "dose",
    covariates = "x1",
    time = "time",
    time_effects = TRUE,
    iterations = 500L,
    reltol = 1e-8
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = entropy_weights,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 717
  )

  dr <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1
  )

  expect_s3_class(entropy_weights, "cdmc_entropy_balance_weights")
  expect_length(entropy_weights$weights, nrow(panel))
  expect_true(all(is.finite(entropy_weights$weights)))
  expect_true(is.finite(entropy_weights$max_abs_balance))
  expect_true(
    entropy_weights$max_abs_balance <= max(abs(entropy_weights$balance_summary$base_moment)) + 1e-8
  )
  expect_true(fit$weight_supplied)
  expect_s3_class(dr, "cdmc_dose_response")
})

test_that("kernel-balance weights can be estimated and reused in dose-response fitting", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 727
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.2 + sin(panel$x1[active]) + ifelse(panel$time[active] > stats::median(panel$time), 0.5, -0.1)
  )

  kernel_weights <- cdmc_kernel_balance_weights(
    data = panel,
    dose = "dose",
    covariates = "x1",
    time = "time",
    time_effects = TRUE,
    n_centers = 6L,
    iterations = 400L,
    reltol = 1e-8
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = kernel_weights,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 727
  )

  dr <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1
  )

  expect_s3_class(kernel_weights, "cdmc_kernel_balance_weights")
  expect_length(kernel_weights$weights, nrow(panel))
  expect_true(all(is.finite(kernel_weights$weights)))
  expect_true(is.finite(kernel_weights$max_abs_balance))
  expect_true(kernel_weights$actual_n_centers >= 1L)
  expect_true(is.finite(kernel_weights$applied_bandwidth))
  expect_true(fit$weight_supplied)
  expect_s3_class(dr, "cdmc_dose_response")
})

test_that("adaptive-balance weights can be estimated and reused in dose-response fitting", {
  panel <- simulate_cdmc_data(
    n_units = 20,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.12,
    switch_off_prob = 0.5,
    seed = 829
  )
  original_dose <- panel$dose
  active <- abs(original_dose) > 0
  panel$dose[active] <- sign(original_dose[active]) * (
    0.25 + sin(panel$x1[active]) + ifelse(panel$time[active] > stats::median(panel$time), 0.4, -0.1)
  )

  adaptive_weights <- cdmc_adaptive_balance_weights(
    data = panel,
    dose = "dose",
    covariates = "x1",
    time = "time",
    time_effects = TRUE,
    methods = c("entropy_balance", "kernel_balance"),
    kernel_balance_centers = 6L,
    entropy_balance_iterations = 400L,
    kernel_balance_iterations = 400L
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = adaptive_weights,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 829
  )

  dr <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1
  )

  expect_s3_class(adaptive_weights, "cdmc_adaptive_balance_weights")
  expect_true(adaptive_weights$selected_method %in% c("entropy_balance", "kernel_balance"))
  expect_length(adaptive_weights$weights, nrow(panel))
  expect_true(all(is.finite(adaptive_weights$weights)))
  expect_true(is.data.frame(adaptive_weights$candidate_scores))
  expect_true(all(c("entropy_balance", "kernel_balance") %in% adaptive_weights$candidate_scores$method))
  expect_true(any(adaptive_weights$candidate_scores$selected))
  expect_true(fit$weight_supplied)
  expect_s3_class(dr, "cdmc_dose_response")
})