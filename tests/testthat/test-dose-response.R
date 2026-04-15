test_that("spline dose-response predicts larger effects at larger doses", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 505
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
    washout = 1,
    lag_order = 1,
    seed = 505
  )

  dr <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 1,
    df = 3,
    include_zero_dose = TRUE
  )

  pred <- predict(dr, dose = c(0, 1))

  expect_s3_class(dr, "cdmc_dose_response")
  expect_true(nrow(pred) == 2)
  expect_true(pred$estimate[2] > pred$estimate[1])
})

test_that("weighted linear dose-response can be fit from residual effects", {
  panel <- simulate_cdmc_data(
    n_units = 8,
    n_times = 8,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.45,
    seed = 606
  )
  panel$w <- rep(seq(1, 2, length.out = 8), each = 8)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 606
  )

  dr <- cdmc_dose_response(
    fit,
    model = "linear",
    lag_order = 1,
    weights = "w"
  )

  expect_s3_class(dr, "cdmc_dose_response")
  expect_true(all(c("dose_lag0", "dose_lag1") %in% names(dr$coefficients)))
})

test_that("slope prediction uses a symmetric finite difference under signed doses", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 909
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

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
    lag_order = 0,
    seed = 909
  )
  dr <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 0,
    df = 3,
    include_zero_dose = TRUE
  )

  eps <- 1e-4
  slope <- predict(dr, dose = -0.5, type = "slope", eps = eps)
  response_up <- predict(dr, dose = -0.5 + eps, type = "response")
  response_down <- predict(dr, dose = -0.5 - eps, type = "response")
  manual_slope <- (response_up$estimate - response_down$estimate) / (2 * eps)

  expect_s3_class(dr, "cdmc_dose_response")
  expect_equal(slope$estimate, manual_slope, tolerance = 1e-6)
})

test_that("gam dose-response supports anchored nonlinear prediction", {
  testthat::skip_if_not_installed("mgcv")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 9393
  )
  panel$y <- panel$y + 0.8 * panel$dose^2

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
    lag_order = 0,
    seed = 9393
  )
  dr <- cdmc_dose_response(
    fit,
    model = "gam",
    lag_order = 0,
    df = 4,
    include_zero_dose = TRUE
  )

  zero_pred <- predict(dr, dose = 0)
  positive_pred <- predict(dr, dose = 0.5)
  slope <- predict(dr, dose = 0.5, type = "slope")

  expect_s3_class(dr, "cdmc_dose_response")
  expect_identical(dr$model, "gam")
  expect_equal(zero_pred$estimate, 0, tolerance = 1e-6)
  expect_true(positive_pred$estimate > zero_pred$estimate)
  expect_true(is.finite(slope$estimate))
})

test_that("tree dose-response supports anchored nonlinear prediction", {
  testthat::skip_if_not_installed("rpart")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 9494
  )
  panel$y <- panel$y + ifelse(panel$dose > 0.25, 0.8, 0)

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
    lag_order = 0,
    seed = 9494
  )
  dr <- cdmc_dose_response(
    fit,
    model = "tree",
    lag_order = 0,
    include_zero_dose = TRUE
  )

  zero_pred <- predict(dr, dose = 0)
  positive_pred <- predict(dr, dose = 0.8)
  slope <- predict(dr, dose = 0.8, type = "slope")

  expect_s3_class(dr, "cdmc_dose_response")
  expect_identical(dr$model, "tree")
  expect_true("dose_lag0" %in% names(dr$coefficients))
  expect_equal(zero_pred$estimate, 0, tolerance = 1e-6)
  expect_true(positive_pred$estimate >= zero_pred$estimate)
  expect_true(is.finite(slope$estimate))
})

test_that("forest dose-response supports anchored nonlinear prediction", {
  testthat::skip_if_not_installed("ranger")

  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 9595
  )
  panel$y <- panel$y + ifelse(panel$dose > 0.25, 0.8, 0)

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
    lag_order = 0,
    seed = 9595
  )
  dr <- cdmc_dose_response(
    fit,
    model = "forest",
    lag_order = 0,
    forest_trees = 50,
    include_zero_dose = TRUE
  )

  zero_pred <- predict(dr, dose = 0)
  positive_pred <- predict(dr, dose = 0.8)
  slope <- predict(dr, dose = 0.8, type = "slope")

  expect_s3_class(dr, "cdmc_dose_response")
  expect_identical(dr$model, "forest")
  expect_identical(dr$forest_trees, 50L)
  expect_true("dose_lag0" %in% names(dr$coefficients))
  expect_equal(zero_pred$estimate, 0, tolerance = 1e-6)
  expect_true(positive_pred$estimate >= zero_pred$estimate)
  expect_true(is.finite(slope$estimate))
})
