test_that("cdmc_fit pads unbalanced panels onto the full unit-time grid", {
  panel <- simulate_cdmc_data(n_units = 4, n_times = 5, seed = 1)
  drop_index <- which(panel$dose != 0)[1L]
  panel <- panel[-drop_index, ]

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    lambda = 0.2,
    rank_max = 2,
    lag_order = 1,
    seed = 1
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_equal(nrow(fit$data), 20)
  expect_equal(sum(fit$data$.cdmc_observed), nrow(panel))
  expect_false(any(fit$eligible_mask[!fit$observed_mask]))
  expect_true(all(is.na(fit$data$y[!fit$data$.cdmc_observed])))
  expect_true(all(is.na(predict(fit, type = "tau_model")$estimate[!fit$data$.cdmc_observed])))
})

test_that("cdmc_fit supports signed doses with zero as the control state", {
  panel <- simulate_cdmc_data(n_units = 4, n_times = 5, signed_dose = TRUE, seed = 2)

  fit <- cdmc_fit(
    data = panel,
    outcome = "y",
    dose = "dose",
    unit = "unit",
    time = "time",
    lambda = 0.2,
    rank_max = 2,
    lag_order = 1,
    seed = 2
  )

  expect_s3_class(fit, "cdmc_fit")
  expect_true(any(panel$dose < 0))
  expect_true(any(abs(fit$dose_matrix) > fit$zero_tolerance))
})
