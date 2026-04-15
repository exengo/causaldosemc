test_that("cdmc_fit can tune lambda with blocked cross-validation", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 99
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 99
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(fit$lambda %in% c(0.4, 0.2, 0.1))
  expect_equal(fit$lambda, fit$lambda_tuning$selected_lambda)
  expect_equal(length(fit$lambda_tuning$lambda_grid), 3)
  expect_true(all(is.finite(fit$lambda_tuning$mean_scores)))
  expect_true(all(fit$lambda_tuning$holdout_counts > 0))
})

test_that("cdmc_fit can tune lambda with observation weights", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 9,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 199
  )
  panel$w <- 0.5 + abs(panel$x1)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = "w",
    lambda = NULL,
    lambda_selection = "cv",
    lambda_grid = c(0.4, 0.2, 0.1),
    cv_rounds = 2,
    cv_block_size = 1,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 199
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(fit$weight_supplied)
  expect_equal(fit$lambda_tuning$method, "cv")
  expect_true(all(is.finite(fit$lambda_tuning$mean_scores)))
})
