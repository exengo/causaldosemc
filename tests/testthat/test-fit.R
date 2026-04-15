test_that("core estimator recovers simple distributed lag effects", {
  panel <- simulate_cdmc_data(
    n_units = 14,
    n_times = 12,
    rank = 2,
    beta = 1.1,
    lag_beta = 0.35,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    seed = 42
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
    seed = 42
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(length(fit$effect$coefficients) == 2)
  expect_equal(
    unname(fit$effect$coefficients["dose_lag0"]),
    1.1,
    tolerance = 0.45
  )
  expect_equal(
    unname(fit$effect$coefficients["dose_lag1"]),
    0.35,
    tolerance = 0.45
  )
})

test_that("core estimator supports the optional joint objective", {
  panel <- simulate_cdmc_data(
    n_units = 14,
    n_times = 12,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    seed = 424
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
    objective = "joint",
    seed = 424
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_identical(fit$objective, "joint")
  expect_identical(fit$fit_control$objective, "joint")
  expect_true(all(fit$optimization_mask == fit$observed_mask))
  expect_equal(unname(fit$effect$coefficients["dose_lag0"]), 0.9, tolerance = 0.55)
  expect_equal(unname(fit$effect$coefficients["dose_lag1"]), 0.25, tolerance = 0.55)
  expect_true(all(is.finite(fit$data$.cdmc_y0_hat[fit$data$.cdmc_observed])))
})

test_that("core estimator accepts observation weights in the main fit", {
  panel <- simulate_cdmc_data(
    n_units = 14,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    seed = 420
  )
  panel$w <- rep(seq(0.9, 1.1, length.out = 12), times = 14)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = "w",
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 420
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(fit$weight_supplied)
  expect_equal(fit$fit_control$weights, ".cdmc_weight")
  expect_true(all(fit$data$.cdmc_weight > 0))
  expect_equal(unname(fit$effect$coefficients["dose_lag0"]), 1, tolerance = 0.5)
  expect_equal(unname(fit$effect$coefficients["dose_lag1"]), 0.3, tolerance = 0.5)
})

test_that("core estimator aligns numeric weight vectors on unbalanced panels", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.18,
    switch_off_prob = 0.45,
    seed = 430
  )
  drop_index <- head(which(panel$dose != 0), 2L)
  panel <- panel[-drop_index, ]
  weights <- 0.5 + abs(panel$x1)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weights = weights,
    lambda = 0.2,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 430
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(fit$weight_supplied)
  expect_true(all(fit$data$.cdmc_weight[fit$data$.cdmc_observed] > 0))
  expect_true(all(fit$data$.cdmc_weight[!fit$data$.cdmc_observed] == 0))
})

test_that("core estimator recovers linear effects under signed doses", {
  panel <- simulate_cdmc_data(
    n_units = 14,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 4242
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
    seed = 4242
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(any(panel$dose < 0))
  expect_equal(unname(fit$effect$coefficients["dose_lag0"]), 1, tolerance = 0.45)
  expect_equal(unname(fit$effect$coefficients["dose_lag1"]), 0.25, tolerance = 0.45)
})

test_that("core estimator supports spline second-stage effects", {
  panel <- simulate_cdmc_data(
    n_units = 14,
    n_times = 12,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    seed = 5252
  )
  panel$y <- panel$y + 0.9 * panel$dose^2

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
    lag_order = 0,
    effect_model = "spline",
    effect_df = 3,
    seed = 5252
  )

  fitted_effect <- predict(fit, type = "tau_model")
  active <- cbind(fitted_effect, dose = panel$dose)
  active <- active[active$dose > 0, , drop = FALSE]
  dose_split <- stats::median(active$dose)

  expect_s3_class(fit, "cdmc_fit")
  expect_identical(fit$effect_model, "spline")
  expect_identical(fit$effect_df, 3L)
  expect_true(length(fit$effect$coefficients) > 1L)
  expect_true(".cdmc_tau_model" %in% names(fit$data))
  expect_equal(fitted_effect$estimate, fit$data$.cdmc_tau_model)
  expect_true(mean(active$estimate[active$dose > dose_split]) > mean(active$estimate[active$dose <= dose_split]))
})