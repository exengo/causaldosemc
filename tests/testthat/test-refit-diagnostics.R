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

test_that("placebo equivalence passes for near-zero placebo effects", {
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
  expect_true(equivalence$equivalent)
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
  expect_true(joint$passed)
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
