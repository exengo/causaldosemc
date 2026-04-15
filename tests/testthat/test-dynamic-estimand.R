test_that("dynamic estimands reproduce built-in linear path contrasts", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 8181
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
    soft_maxit = 200,
    seed = 8181
  )

  history <- data.frame(
    dose_lag0 = c(1, 0.5),
    dose_lag1 = c(0.5, 1),
    row.names = c("frontloaded", "backloaded")
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = history,
    type = "contrast",
    aggregate = "both"
  )
  coefficient_vector <- fit$effect$coefficients[c("dose_lag0", "dose_lag1")]
  manual <- as.vector(as.matrix(history) %*% coefficient_vector)

  expect_s3_class(dynamic, "cdmc_dynamic_estimand")
  expect_equal(dynamic$estimate_table$estimate, manual)
  expect_equal(unname(dynamic$estimates[["mean_dynamic_contrast"]]), mean(manual))
})

test_that("dynamic estimands reproduce spline dose-response slopes", {
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
    seed = 8282
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
    seed = 8282
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 0,
    df = 3,
    include_zero_dose = TRUE
  )

  history <- data.frame(dose_lag0 = c(-0.5, 0.5))
  dynamic <- cdmc_dynamic_estimand(
    dose_response,
    history = history,
    type = "slope",
    aggregate = "both",
    path_weights = c(1, 2)
  )
  manual <- predict(dose_response, history = history, type = "slope")$estimate

  expect_s3_class(dynamic, "cdmc_dynamic_estimand")
  expect_equal(dynamic$estimate_table$estimate, manual)
  expect_equal(
    unname(dynamic$estimates[["mean_dynamic_slope"]]),
    stats::weighted.mean(manual, w = c(1, 2))
  )
})

test_that("dynamic estimands reproduce DR path contrasts", {
  panel <- simulate_cdmc_data(
    n_units = 18,
    n_times = 10,
    rank = 2,
    beta = 0.9,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.18,
    switch_off_prob = 0.42,
    seed = 8383
  )

  fit <- cdmc_dr_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    weight_method = "gaussian_gps",
    n_folds = 2,
    lambda = 0.2,
    rank_max = 3,
    washout = 0,
    lag_order = 1,
    seed = 8383
  )

  history <- data.frame(
    dose_lag0 = c(1, 0.25),
    dose_lag1 = c(0.5, 1)
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = history,
    type = "contrast",
    aggregate = "individual"
  )
  coefficient_vector <- fit$effect$coefficients[c("dose_lag0", "dose_lag1")]
  manual <- as.vector(as.matrix(history) %*% coefficient_vector)

  expect_s3_class(dynamic, "cdmc_dynamic_estimand")
  expect_equal(dynamic$estimate_table$estimate, manual)
})

test_that("bootstrap supports dynamic estimands from fitted models", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 1,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.05,
    switch_on_prob = 0.15,
    switch_off_prob = 0.45,
    seed = 8484
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
    seed = 8484
  )
  dynamic <- cdmc_dynamic_estimand(
    fit,
    history = data.frame(dose_lag0 = c(1, 0.5), dose_lag1 = c(0.5, 1)),
    type = "contrast",
    aggregate = "both"
  )

  boot <- cdmc_bootstrap(dynamic, n_boot = 2, statistics = "estimate", seed = 8484)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true(all(c(
    "dynamic_contrast_path1",
    "dynamic_contrast_path2",
    "mean_dynamic_contrast"
  ) %in% colnames(boot$draws)))
  expect_equal(boot$simultaneous$method, "max-t")
  expect_true(all(c(
    "simultaneous_conf_low",
    "simultaneous_conf_high"
  ) %in% names(boot$summary)))
  expect_true(all(boot$summary$simultaneous_conf_low <= boot$summary$simultaneous_conf_high))
  expect_true(boot$simultaneous$n_success <= boot$n_boot)
})

test_that("bootstrap supports dynamic estimands built from dose-response objects", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 8,
    rank = 2,
    beta = 0.2,
    lag_beta = NULL,
    n_covariates = 1,
    noise_sd = 0.04,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    signed_dose = TRUE,
    negative_dose_prob = 0.5,
    seed = 8585
  )
  panel$y <- panel$y + 0.7 * panel$dose^2

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
    seed = 8585
  )
  dose_response <- cdmc_dose_response(
    fit,
    model = "spline",
    lag_order = 0,
    df = 3,
    include_zero_dose = TRUE
  )
  dynamic <- cdmc_dynamic_estimand(
    dose_response,
    history = data.frame(dose_lag0 = c(-0.5, 0.5)),
    type = "slope",
    aggregate = "mean"
  )

  boot <- cdmc_bootstrap(dynamic, n_boot = 2, statistics = "estimate", seed = 8585)

  expect_s3_class(boot, "cdmc_bootstrap")
  expect_true("mean_dynamic_slope" %in% colnames(boot$draws))
})