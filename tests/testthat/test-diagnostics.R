test_that("carryover diagnostic detects positive post-exit residual effects", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.5,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 7
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
    seed = 7
  )

  diagnostic <- cdmc_carryover_test(fit, periods = 1)

  expect_s3_class(diagnostic, "cdmc_carryover_test")
  expect_true(diagnostic$n > 0)
  expect_true(diagnostic$mean_tau > 0)
  expect_true(is.finite(diagnostic$p_value))
  expect_true(diagnostic$p_value < 0.05)
})
