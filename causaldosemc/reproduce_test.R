pkgload::load_all('.')
panel <- simulate_cdmc_data(
  n_units = 20,
  n_times = 10,
  rank = 2,
  beta = 0.65,
  lag_beta = NULL,
  n_covariates = 1,
  noise_sd = 0.04,
  switch_on_prob = 0.18,
  switch_off_prob = 0.42,
  seed = 1768
)
original_dose <- panel$dose
active <- abs(original_dose) > 0
panel$dose[active] <- sign(original_dose[active]) * (
  0.2 + sin(2 * panel$x1[active]) + (panel$time[active] / max(panel$time))^2
)
panel$y <- panel$y + 0.6 * (panel$dose - original_dose)

fit <- cdmc_dr_fit(
  data = panel,
  outcome = 'y',
  dose = 'dose',
  unit = 'unit',
  time = 'time',
  covariates = 'x1',
  weight_method = 'gaussian_gps',
  gps_model = 'linear',
  n_folds = 2,
  lambda = 0.2,
  rank_max = 3,
  washout = 0,
  lag_order = 0,
  seed = 1768
)

scan <- cdmc_sensitivity_scan(
  fit,
  washout_grid = fit$washout,
  zero_tolerance_grid = fit$zero_tolerance,
  weight_specs = list(gam = list(weight_method = 'gaussian_gps', gps_model = 'gam', gps_df = 4))
)

print(scan$results[, c('weight_scenario', 'fit_ok', 'fit_error', 'weight_mode')])
