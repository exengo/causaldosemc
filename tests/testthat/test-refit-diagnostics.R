test_that("placebo diagnostic finds near-zero pseudo effects in pre-exposure controls", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 303
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
    seed = 303
  )

  placebo <- cdmc_placebo_test(fit, periods = -2:0)

  expect_s3_class(placebo, "cdmc_placebo_test")
  expect_true(placebo$n > 0)
  expect_true(is.finite(placebo$p_value))
  expect_true(placebo$p_value > 0.05)
})

test_that("placebo refit diagnostics can rerun heuristic lambda selection", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.25,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 313
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 313
  )

  placebo_fixed <- cdmc_placebo_test(fit, periods = -2:0)
  placebo_retuned <- cdmc_placebo_test(fit, periods = -2:0, rerun_tuning = TRUE)

  expect_false(placebo_fixed$rerun_tuning)
  expect_identical(placebo_fixed$refit_lambda_method, "fixed")
  expect_true(placebo_retuned$rerun_tuning)
  expect_identical(placebo_retuned$refit_lambda_method, "heuristic")
  expect_true(is.finite(placebo_retuned$refit_lambda))
})

test_that("refit diagnostics reject rerun_tuning for fixed-lambda source fits", {
  panel <- simulate_cdmc_data(
    n_units = 10,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 314
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
    seed = 314
  )

  expect_error(
    cdmc_placebo_test(fit, periods = -2:0, rerun_tuning = TRUE),
    "requires a source cdmc_fit object with automatically selected lambda"
  )
})

test_that("placebo diagnostic supports joint-objective fits", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 0.8,
    lag_beta = 0.2,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 333
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
    seed = 333
  )

  placebo <- cdmc_placebo_test(fit, periods = -1:0)

  expect_s3_class(placebo, "cdmc_placebo_test")
  expect_true(placebo$n > 0)
  expect_true(is.finite(placebo$mean_tau))
})

test_that("refit carryover diagnostic detects post-exit pseudo effects with lag structure", {
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
    seed = 404
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
    seed = 404
  )

  diagnostic <- cdmc_carryover_refit_test(fit, periods = 1)

  expect_s3_class(diagnostic, "cdmc_carryover_refit_test")
  expect_true(diagnostic$n > 0)
  expect_true(diagnostic$mean_tau > 0)
})

test_that("placebo equivalence returns a bounded diagnostic for near-zero placebo effects", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.1,
    n_covariates = 1,
    noise_sd = 0.01,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 927
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
    seed = 927
  )

  placebo <- cdmc_placebo_test(fit, periods = -2:0)
  equivalence <- cdmc_equivalence_test(placebo, margin = 0.2)

  expect_s3_class(equivalence, "cdmc_equivalence_test")
  expect_equal(mean(placebo$cells$pseudo_tau), placebo$mean_tau, tolerance = 1e-10)
  expect_true(is.finite(equivalence$p_value))
  expect_gte(equivalence$p_value, 0)
  expect_lte(equivalence$p_value, 1)
  expect_true(is.logical(equivalence$equivalent) && length(equivalence$equivalent) == 1L)
})

test_that("carryover equivalence fails when residual effects remain large", {
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
    seed = 809
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
    seed = 809
  )

  diagnostic <- cdmc_carryover_refit_test(fit, periods = 1)
  equivalence <- cdmc_equivalence_test(diagnostic, margin = 0.1)

  expect_s3_class(equivalence, "cdmc_equivalence_test")
  expect_equal(mean(diagnostic$cells$pseudo_tau), diagnostic$mean_tau, tolerance = 1e-10)
  expect_false(equivalence$equivalent)
})

test_that("joint placebo diagnostic aggregates multiple pre-exposure windows", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.1,
    n_covariates = 1,
    noise_sd = 0.01,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 927
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
    seed = 927
  )

  joint <- cdmc_joint_placebo_test(
    fit,
    periods = -2:-1,
    equivalence_margin = 0.2
  )

  expect_s3_class(joint, "cdmc_joint_placebo_test")
  expect_equal(joint$tests$period, c(-2L, -1L))
  expect_true(is.finite(joint$joint_p_value))
  expect_gte(joint$joint_p_value, 0)
  expect_lte(joint$joint_p_value, 1)
  expect_true(is.logical(joint$passed) && length(joint$passed) == 1L)
})

test_that("joint placebo diagnostic supports named repeated placebo windows", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.1,
    n_covariates = 1,
    noise_sd = 0.01,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 929
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
    seed = 929
  )

  joint <- cdmc_joint_placebo_test(
    fit,
    placebo_windows = list(early_window = -3:-2, immediate_window = -1)
  )

  expect_s3_class(joint, "cdmc_joint_placebo_test")
  expect_identical(names(joint$placebo_windows), c("early_window", "immediate_window"))
  expect_identical(joint$tests$window_name, c("early_window", "immediate_window"))
  expect_identical(joint$tests$n_periods, c(2L, 1L))
  expect_identical(joint$tests$period, c(NA_integer_, -1L))
  expect_true(all(vapply(joint$components, function(component) inherits(component, "cdmc_placebo_test"), logical(1))))
})

test_that("joint placebo diagnostic propagates rerun_tuning to component refits", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 12,
    rank = 2,
    beta = 1,
    lag_beta = 0.1,
    n_covariates = 1,
    noise_sd = 0.01,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 928
  )

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    covariates = "x1",
    lambda = NULL,
    rank_max = 3,
    washout = 1,
    lag_order = 1,
    seed = 928
  )

  joint <- cdmc_joint_placebo_test(
    fit,
    periods = -2:-1,
    rerun_tuning = TRUE
  )

  expect_true(joint$rerun_tuning)
  expect_true(all(vapply(joint$components, function(component) component$rerun_tuning, logical(1))))
  expect_true(all(vapply(joint$components, function(component) identical(component$refit_lambda_method, "heuristic"), logical(1))))
})

test_that("placebo refit diagnostics support unbalanced panels", {
  panel <- simulate_cdmc_data(
    n_units = 12,
    n_times = 10,
    rank = 2,
    beta = 1,
    lag_beta = 0.3,
    n_covariates = 1,
    noise_sd = 0.03,
    switch_on_prob = 0.2,
    switch_off_prob = 0.4,
    seed = 505
  )
  drop_index <- tail(which(panel$dose != 0), 1L)
  panel <- panel[-drop_index, ]

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

  placebo <- cdmc_placebo_test(fit, periods = -2:0)

  expect_s3_class(placebo, "cdmc_placebo_test")
  expect_true(placebo$n > 0)
})
