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