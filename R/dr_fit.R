cdmc_assign_unit_folds <- function(unit_levels, n_folds) {
  n_folds <- as.integer(n_folds)
  if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
    stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
  }
  if (n_folds > length(unit_levels)) {
    stop("n_folds cannot exceed the number of units.", call. = FALSE)
  }

  shuffled_units <- sample(unit_levels)
  fold_sizes <- rep(length(unit_levels) %/% n_folds, n_folds)
  fold_sizes[seq_len(length(unit_levels) %% n_folds)] <- fold_sizes[seq_len(length(unit_levels) %% n_folds)] + 1L

  fold_id <- integer(length(unit_levels))
  names(fold_id) <- unit_levels
  start_index <- 1L

  for (fold_index in seq_len(n_folds)) {
    end_index <- start_index + fold_sizes[fold_index] - 1L
    current_units <- shuffled_units[seq.int(start_index, end_index)]
    fold_id[current_units] <- fold_index
    start_index <- end_index + 1L
  }

  data.frame(
    unit = unit_levels,
    fold = unname(fold_id[unit_levels]),
    stringsAsFactors = FALSE
  )
}

cdmc_resolve_unit_folds <- function(unit_levels, n_folds, fold_assignments = NULL) {
  if (is.null(fold_assignments)) {
    return(cdmc_assign_unit_folds(unit_levels, n_folds = n_folds))
  }

  if (!is.data.frame(fold_assignments)) {
    stop("fold_assignments must be NULL or a data frame with columns 'unit' and 'fold'.", call. = FALSE)
  }

  required_columns <- c("unit", "fold")
  missing_columns <- setdiff(required_columns, names(fold_assignments))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "fold_assignments is missing required columns: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  provided_units <- as.character(fold_assignments$unit)
  if (anyNA(provided_units) || anyDuplicated(provided_units)) {
    stop("fold_assignments must contain each unit exactly once with no missing values.", call. = FALSE)
  }

  unit_keys <- as.character(unit_levels)
  missing_units <- setdiff(unit_keys, provided_units)
  extra_units <- setdiff(provided_units, unit_keys)
  if (length(missing_units) > 0L || length(extra_units) > 0L) {
    message_parts <- character(0)
    if (length(missing_units) > 0L) {
      message_parts <- c(message_parts, sprintf("missing units: %s", paste(missing_units, collapse = ", ")))
    }
    if (length(extra_units) > 0L) {
      message_parts <- c(message_parts, sprintf("unknown units: %s", paste(extra_units, collapse = ", ")))
    }

    stop(
      sprintf("fold_assignments must match the current unit set exactly (%s).", paste(message_parts, collapse = "; ")),
      call. = FALSE
    )
  }

  unit_index <- match(unit_keys, provided_units)
  fold_values_raw <- fold_assignments$fold[unit_index]
  if (anyNA(fold_values_raw) || any(!is.finite(fold_values_raw))) {
    stop("fold_assignments$fold must be finite for every unit.", call. = FALSE)
  }
  if (any(abs(fold_values_raw - round(fold_values_raw)) > sqrt(.Machine$double.eps))) {
    stop("fold_assignments$fold must contain integer-valued fold labels.", call. = FALSE)
  }

  fold_values <- as.integer(round(fold_values_raw))
  if (any(fold_values < 1L)) {
    stop("fold_assignments$fold must contain positive fold labels.", call. = FALSE)
  }

  unique_folds <- sort(unique(fold_values))
  if (!identical(unique_folds, seq_along(unique_folds))) {
    stop("fold_assignments$fold must use consecutive fold labels 1, 2, ..., K.", call. = FALSE)
  }
  if (length(unique_folds) < 2L) {
    stop("fold_assignments must define at least two folds.", call. = FALSE)
  }

  if (!is.null(n_folds)) {
    n_folds <- as.integer(n_folds)
    if (!is.numeric(n_folds) || length(n_folds) != 1L || n_folds < 2L) {
      stop("n_folds must be an integer greater than or equal to 2.", call. = FALSE)
    }
    if (length(unique_folds) != n_folds) {
      stop(
        sprintf(
          "Provided fold_assignments define %d folds, but n_folds = %d.",
          length(unique_folds),
          n_folds
        ),
        call. = FALSE
      )
    }
  }

  data.frame(
    unit = unit_levels,
    fold = fold_values,
    stringsAsFactors = FALSE
  )
}

cdmc_resolve_dr_weight_method <- function(weights, weight_method) {
  if (is.null(weight_method)) {
    return(if (is.null(weights)) "gaussian_gps" else "external")
  }

  resolved <- match.arg(weight_method, c("external", "gaussian_gps", "kernel_gps", "cbps"))
  if (identical(resolved, "external") && is.null(weights)) {
    stop("weights must be supplied when weight_method = 'external'.", call. = FALSE)
  }
  if (!identical(resolved, "external") && !is.null(weights)) {
    stop("Do not supply weights when weight_method uses an internal weighting nuisance model.", call. = FALSE)
  }

  resolved
}

cdmc_resolve_gps_bandwidth <- function(gps_bandwidth) {
  if (is.null(gps_bandwidth)) {
    return(NULL)
  }

  if (!is.numeric(gps_bandwidth) || length(gps_bandwidth) != 1L || !is.finite(gps_bandwidth) || gps_bandwidth <= 0) {
    stop("gps_bandwidth must be NULL or a single positive numeric value.", call. = FALSE)
  }

  as.numeric(gps_bandwidth)
}

cdmc_resolve_cbps_iterations <- function(cbps_iterations) {
  if (!is.numeric(cbps_iterations) || length(cbps_iterations) != 1L || !is.finite(cbps_iterations) || cbps_iterations < 1) {
    stop("cbps_iterations must be a positive integer.", call. = FALSE)
  }

  as.integer(cbps_iterations)
}

cdmc_resolve_gps_forest_trees <- function(gps_forest_trees) {
  if (!is.numeric(gps_forest_trees) || length(gps_forest_trees) != 1L || !is.finite(gps_forest_trees) || gps_forest_trees < 1) {
    stop("gps_forest_trees must be a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_trees)
}

cdmc_resolve_gps_forest_mtry <- function(gps_forest_mtry) {
  if (is.null(gps_forest_mtry)) {
    return(NULL)
  }

  if (!is.numeric(gps_forest_mtry) || length(gps_forest_mtry) != 1L || !is.finite(gps_forest_mtry) || gps_forest_mtry < 1) {
    stop("gps_forest_mtry must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_mtry)
}

cdmc_resolve_gps_forest_min_node_size <- function(gps_forest_min_node_size) {
  if (is.null(gps_forest_min_node_size)) {
    return(NULL)
  }

  if (!is.numeric(gps_forest_min_node_size) || length(gps_forest_min_node_size) != 1L || !is.finite(gps_forest_min_node_size) || gps_forest_min_node_size < 1) {
    stop("gps_forest_min_node_size must be NULL or a positive integer.", call. = FALSE)
  }

  as.integer(gps_forest_min_node_size)
}

cdmc_cap_internal_weights <- function(weights, max_weight = NULL) {
  if (is.null(max_weight)) {
    return(as.numeric(weights))
  }

  if (!is.numeric(max_weight) || length(max_weight) != 1L || !is.finite(max_weight) || max_weight <= 0) {
    stop("max_weight must be a single positive numeric value when supplied.", call. = FALSE)
  }

  as.numeric(pmin(weights, max_weight))
}

cdmc_resolve_gps_spline_covariates <- function(covariates, gps_spline_covariates) {
  if (is.null(gps_spline_covariates)) {
    return(character(0))
  }

  if (!is.character(gps_spline_covariates)) {
    stop("gps_spline_covariates must be NULL or a character vector of column names.", call. = FALSE)
  }

  missing_covariates <- setdiff(gps_spline_covariates, covariates %||% character(0))
  if (length(missing_covariates) > 0L) {
    stop(
      sprintf("gps_spline_covariates must be a subset of weight_covariates. Unknown names: %s.", paste(missing_covariates, collapse = ", ")),
      call. = FALSE
    )
  }

  gps_spline_covariates
}

cdmc_is_gps_spline_candidate <- function(x, gps_df) {
  if (!is.numeric(x)) {
    return(FALSE)
  }

  length(unique(stats::na.omit(x))) > gps_df
}

cdmc_build_gam_gps_formula <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects,
  gps_df = 4L,
  gps_spline_covariates = covariates
) {
  cdmc_assert_installed("mgcv")

  smooth_covariates <- cdmc_resolve_gps_spline_covariates(covariates, gps_spline_covariates)
  termlabels <- vapply(covariates %||% character(0), function(covariate) {
    if (covariate %in% smooth_covariates &&
        cdmc_is_gps_spline_candidate(data[[covariate]], gps_df = gps_df)) {
      sprintf("s(%s, k = %d)", covariate, gps_df)
    } else {
      covariate
    }
  }, character(1))

  if (isTRUE(gps_time_effects)) {
    termlabels <- c(termlabels, sprintf("factor(%s)", time))
  }

  formula_text <- sprintf(
    "%s ~ %s",
    dose,
    if (length(termlabels) == 0L) "1" else paste(termlabels, collapse = " + ")
  )
  formula_environment <- new.env(parent = baseenv())
  formula_environment$s <- mgcv::s
  stats::as.formula(formula_text, env = formula_environment)
}

cdmc_build_gps_formula <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = covariates
) {
  gps_model <- match.arg(gps_model)
  if (identical(gps_model, "gam")) {
    return(cdmc_build_gam_gps_formula(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = gps_time_effects,
      gps_df = gps_df,
      gps_spline_covariates = gps_spline_covariates
    ))
  }

  spline_covariates <- if (identical(gps_model, "spline")) {
    cdmc_resolve_gps_spline_covariates(covariates, gps_spline_covariates)
  } else {
    character(0)
  }

  termlabels <- vapply(covariates %||% character(0), function(covariate) {
    if (identical(gps_model, "spline") &&
        covariate %in% spline_covariates &&
        cdmc_is_gps_spline_candidate(data[[covariate]], gps_df = gps_df)) {
      sprintf("splines::ns(%s, df = %d)", covariate, gps_df)
    } else {
      covariate
    }
  }, character(1))

  if (isTRUE(gps_time_effects)) {
    termlabels <- c(termlabels, sprintf("factor(%s)", time))
  }

  stats::reformulate(termlabels = termlabels, response = dose)
}

cdmc_build_gps_forest_frame <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects
) {
  predictor_names <- covariates %||% character(0)
  forest_data <- data.frame(.cdmc_gps_dose = data[[dose]], check.names = FALSE)

  for (covariate in predictor_names) {
    forest_data[[covariate]] <- data[[covariate]]
  }

  if (isTRUE(gps_time_effects)) {
    forest_data$.cdmc_gps_time_factor <- factor(data[[time]])
    predictor_names <- c(predictor_names, ".cdmc_gps_time_factor")
  }

  list(
    data = forest_data,
    predictor_names = predictor_names,
    formula = stats::reformulate(termlabels = predictor_names, response = ".cdmc_gps_dose")
  )
}

cdmc_predict_gps_mean_model <- function(object, newdata) {
  if (inherits(object, "cdmc_gps_constant_mean")) {
    return(rep(object$mean, nrow(newdata)))
  }
  if (inherits(object, "ranger")) {
    return(as.numeric(predict(object, data = newdata)$predictions))
  }

  as.numeric(stats::predict(object, newdata = newdata))
}

cdmc_fit_gps_mean_model <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL
) {
  gps_model <- match.arg(gps_model)
  if (identical(gps_model, "forest")) {
    cdmc_assert_installed("ranger")
    forest_frame <- cdmc_build_gps_forest_frame(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = gps_time_effects
    )

    if (length(forest_frame$predictor_names) == 0L) {
      fit_object <- structure(list(mean = mean(data[[dose]], na.rm = TRUE)), class = "cdmc_gps_constant_mean")
    } else {
      fit_object <- ranger::ranger(
        dependent.variable.name = ".cdmc_gps_dose",
        data = forest_frame$data,
        num.trees = gps_forest_trees,
        mtry = gps_forest_mtry,
        min.node.size = gps_forest_min_node_size,
        respect.unordered.factors = "order",
        seed = 1L,
        num.threads = 1L
      )
    }

    return(list(
      formula = forest_frame$formula,
      fit = fit_object,
      prepare_newdata = function(newdata) {
        forest_newdata <- cdmc_build_gps_forest_frame(
          data = newdata,
          dose = dose,
          covariates = covariates,
          time = time,
          gps_time_effects = gps_time_effects
        )$data
        forest_newdata[, forest_frame$predictor_names, drop = FALSE]
      }
    ))
  }

  formula <- cdmc_build_gps_formula(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates
  )

  fit_object <- if (identical(gps_model, "gam")) {
    mgcv::gam(
      formula = formula,
      data = data,
      method = "REML"
    )
  } else if (identical(gps_model, "tree")) {
    cdmc_assert_installed("rpart")
    rpart::rpart(
      formula = formula,
      data = data,
      method = "anova"
    )
  } else {
    stats::lm(formula = formula, data = data)
  }

  list(
    formula = formula,
    fit = fit_object,
    prepare_newdata = function(newdata) newdata
  )
}

cdmc_fit_kernel_density <- function(x, bandwidth = NULL, n = 512L) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) == 0L) {
    stop("Kernel density estimation requires at least one finite observation.", call. = FALSE)
  }

  resolved_bandwidth <- bandwidth %||% suppressWarnings(stats::bw.nrd0(x))
  if (!is.finite(resolved_bandwidth) || resolved_bandwidth <= sqrt(.Machine$double.eps)) {
    fallback_scale <- stats::sd(x)
    if (!is.finite(fallback_scale) || fallback_scale <= sqrt(.Machine$double.eps)) {
      fallback_scale <- 1
    }
    resolved_bandwidth <- fallback_scale * length(x)^(-1 / 5)
  }

  density_fit <- stats::density(
    x = x,
    bw = resolved_bandwidth,
    n = max(256L, as.integer(n)),
    from = min(x) - 4 * resolved_bandwidth,
    to = max(x) + 4 * resolved_bandwidth
  )

  list(
    bandwidth = resolved_bandwidth,
    x = density_fit$x,
    y = density_fit$y
  )
}

cdmc_predict_kernel_density <- function(object, values) {
  approximated <- stats::approx(
    x = object$x,
    y = object$y,
    xout = values,
    yleft = 0,
    yright = 0,
    rule = 1L
  )$y

  pmax(as.numeric(approximated), .Machine$double.xmin)
}

cdmc_fit_gaussian_gps <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL
) {
  mean_model <- cdmc_fit_gps_mean_model(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size
  )
  mean_fit <- mean_model$fit
  fitted_mean <- cdmc_predict_gps_mean_model(mean_fit, newdata = mean_model$prepare_newdata(data))
  residuals <- data[[dose]] - fitted_mean
  sigma_cond <- sqrt(mean(residuals ^ 2))
  if (!is.finite(sigma_cond) || sigma_cond <= sqrt(.Machine$double.eps)) {
    sigma_cond <- sqrt(.Machine$double.eps)
  }

  observed_dose <- stats::na.omit(data[[dose]])
  marginal_mean <- mean(observed_dose)
  marginal_sd <- stats::sd(observed_dose)
  if (!is.finite(marginal_sd) || marginal_sd <= sqrt(.Machine$double.eps)) {
    marginal_sd <- sigma_cond
  }

  list(
    weight_method = "gaussian_gps",
    formula = mean_model$formula,
    lm_fit = mean_fit,
    sigma_cond = sigma_cond,
    marginal_mean = marginal_mean,
    marginal_sd = marginal_sd,
    covariates = covariates,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_bandwidth = NULL,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_prepare_newdata = mean_model$prepare_newdata
  )
}

cdmc_fit_kernel_gps <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL
) {
  mean_model <- cdmc_fit_gps_mean_model(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size
  )
  mean_fit <- mean_model$fit
  fitted_mean <- cdmc_predict_gps_mean_model(mean_fit, newdata = mean_model$prepare_newdata(data))
  residuals <- data[[dose]] - fitted_mean
  residual_density <- cdmc_fit_kernel_density(residuals, bandwidth = gps_bandwidth)
  marginal_density <- cdmc_fit_kernel_density(data[[dose]], bandwidth = gps_bandwidth)

  list(
    weight_method = "kernel_gps",
    formula = mean_model$formula,
    lm_fit = mean_fit,
    residual_density = residual_density,
    marginal_density = marginal_density,
    covariates = covariates,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_bandwidth = residual_density$bandwidth,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    gps_prepare_newdata = mean_model$prepare_newdata
  )
}

cdmc_predict_gaussian_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  conditional_mean <- cdmc_predict_gps_mean_model(object$lm_fit, newdata = object$gps_prepare_newdata(newdata))
  conditional_density <- stats::dnorm(
    newdata[[dose]],
    mean = conditional_mean,
    sd = object$sigma_cond
  )
  numerator <- if (isTRUE(stabilize)) {
    stats::dnorm(
      newdata[[dose]],
      mean = object$marginal_mean,
      sd = object$marginal_sd
    )
  } else {
    rep(1, nrow(newdata))
  }

  weights <- numerator / pmax(conditional_density, .Machine$double.xmin)
  cdmc_cap_internal_weights(weights, max_weight = max_weight)
}

cdmc_predict_kernel_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  conditional_mean <- cdmc_predict_gps_mean_model(object$lm_fit, newdata = object$gps_prepare_newdata(newdata))
  residuals <- newdata[[dose]] - conditional_mean
  conditional_density <- cdmc_predict_kernel_density(object$residual_density, residuals)
  numerator <- if (isTRUE(stabilize)) {
    cdmc_predict_kernel_density(object$marginal_density, newdata[[dose]])
  } else {
    rep(1, nrow(newdata))
  }

  weights <- numerator / pmax(conditional_density, .Machine$double.xmin)
  cdmc_cap_internal_weights(weights, max_weight = max_weight)
}

cdmc_fit_cbps_weights <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  max_weight = NULL
) {
  cbps_fit <- cdmc_cbps_weights(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    time_effects = gps_time_effects,
    model = gps_model,
    df = gps_df,
    spline_covariates = gps_spline_covariates,
    standardize = cbps_standardize,
    method = cbps_method,
    iterations = cbps_iterations,
    twostep = cbps_twostep
  )

  cbps_fit$weights <- cdmc_cap_internal_weights(cbps_fit$weights, max_weight = max_weight)
  cbps_fit$weight_method <- "cbps"
  cbps_fit$gps_time_effects <- gps_time_effects
  cbps_fit$gps_model <- gps_model
  cbps_fit$gps_df <- gps_df
  cbps_fit$gps_spline_covariates <- gps_spline_covariates
  cbps_fit$gps_bandwidth <- NULL
  cbps_fit
}

cdmc_fit_internal_gps <- function(
  weight_method,
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL
) {
  weight_method <- match.arg(weight_method, c("gaussian_gps", "kernel_gps"))

  if (identical(weight_method, "gaussian_gps")) {
    return(cdmc_fit_gaussian_gps(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = gps_time_effects,
      gps_model = gps_model,
      gps_df = gps_df,
      gps_spline_covariates = gps_spline_covariates,
      gps_forest_trees = gps_forest_trees,
      gps_forest_mtry = gps_forest_mtry,
      gps_forest_min_node_size = gps_forest_min_node_size
    ))
  }

  cdmc_fit_kernel_gps(
    data = data,
    dose = dose,
    covariates = covariates,
    time = time,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_bandwidth = gps_bandwidth,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size
  )
}

cdmc_predict_internal_gps_weights <- function(
  object,
  newdata,
  dose,
  stabilize = TRUE,
  max_weight = NULL
) {
  if (identical(object$weight_method, "gaussian_gps")) {
    return(cdmc_predict_gaussian_gps_weights(
      object = object,
      newdata = newdata,
      dose = dose,
      stabilize = stabilize,
      max_weight = max_weight
    ))
  }

  if (identical(object$weight_method, "kernel_gps")) {
    return(cdmc_predict_kernel_gps_weights(
      object = object,
      newdata = newdata,
      dose = dose,
      stabilize = stabilize,
      max_weight = max_weight
    ))
  }

  stop("Internal GPS prediction requires a fitted gaussian_gps or kernel_gps nuisance model.", call. = FALSE)
}

cdmc_validate_holdout_support <- function(eligible_mask) {
  if (any(rowSums(eligible_mask) == 0L)) {
    stop(
      "Each holdout unit must retain at least one eligible zero-dose observation for cross-fitted baseline prediction.",
      call. = FALSE
    )
  }
}

cdmc_extract_time_basis <- function(baseline_fit, n_times) {
  soft_fit <- baseline_fit$soft_fit
  if (is.null(soft_fit) || is.null(soft_fit$v) || length(soft_fit$d) == 0L) {
    return(matrix(0, nrow = n_times, ncol = 0L))
  }

  sweep(soft_fit$v, 2, soft_fit$d, FUN = "*")
}

cdmc_predict_holdout_baseline <- function(object, prepared, eligible_mask, weight_matrix = NULL) {
  if (!inherits(object, "cdmc_fit")) {
    stop("object must inherit from 'cdmc_fit'.", call. = FALSE)
  }

  cdmc_validate_holdout_support(eligible_mask)

  n_units <- prepared$n_units
  n_times <- prepared$n_times
  nuisance_fit <- object$baseline$nuisance
  basis <- cdmc_extract_time_basis(object$baseline, n_times = n_times)
  x_beta <- cdmc_covariate_contribution(
    x_matrices = prepared$x_matrices,
    gamma = nuisance_fit$gamma,
    n_units = n_units,
    n_times = n_times
  )

  baseline_hat <- matrix(0, nrow = n_units, ncol = n_times)
  low_rank_hat <- matrix(0, nrow = n_units, ncol = n_times)
  alpha_hat <- numeric(n_units)

  for (unit_index in seq_len(n_units)) {
    row_mask <- eligible_mask[unit_index, ]
    response <- prepared$y_matrix[unit_index, row_mask] - nuisance_fit$beta[row_mask] - x_beta[unit_index, row_mask]
    design <- cbind(1, basis[row_mask, , drop = FALSE])
    row_weights <- if (is.null(weight_matrix)) NULL else weight_matrix[unit_index, row_mask]

    regression_fit <- if (is.null(row_weights)) {
      stats::lm.fit(x = design, y = response)
    } else {
      stats::lm.wfit(x = design, y = response, w = row_weights)
    }

    coefficients <- regression_fit$coefficients
    coefficients[is.na(coefficients)] <- 0
    alpha_hat[unit_index] <- coefficients[[1L]]

    if (length(coefficients) > 1L) {
      low_rank_hat[unit_index, ] <- as.vector(basis %*% coefficients[-1L])
    }
  }

  baseline_hat <- matrix(alpha_hat, nrow = n_units, ncol = n_times) +
    matrix(nuisance_fit$beta, nrow = n_units, ncol = n_times, byrow = TRUE) +
    x_beta +
    low_rank_hat

  list(
    baseline_hat = baseline_hat,
    alpha = alpha_hat,
    low_rank = low_rank_hat,
    time_basis = basis
  )
}

cdmc_prepare_dr_fold_sample <- function(object, prepared, baseline_hat, weight_matrix = NULL, fold_id = NA_integer_) {
  tau_matrix <- prepared$y_matrix - baseline_hat
  lag_array <- cdmc_build_lagged_doses(prepared$dose_matrix, lag_order = object$lag_order)
  valid_history <- apply(!is.na(lag_array), c(1, 2), all)
  exposure_magnitude <- cdmc_lag_exposure_magnitude(lag_array)

  sample_mask <- valid_history & !is.na(tau_matrix) & exposure_magnitude > 0
  if (!is.null(weight_matrix)) {
    sample_mask <- sample_mask & weight_matrix > 0
  }

  sample_indices <- which(sample_mask, arr.ind = TRUE)
  if (nrow(sample_indices) == 0L) {
    return(list(
      sample_mask = sample_mask,
      design = matrix(0, nrow = 0L, ncol = object$lag_order + 1L),
      pseudo_tau = numeric(0),
      tau_observed = numeric(0),
      tau_model = numeric(0),
      sample_table = data.frame()
    ))
  }

  design <- do.call(
    cbind,
    lapply(seq_len(dim(lag_array)[3L]), function(index) lag_array[, , index][sample_indices])
  )
  colnames(design) <- dimnames(lag_array)[[3L]]

  coefficient_vector <- numeric(ncol(design))
  names(coefficient_vector) <- colnames(design)
  if (length(object$effect$coefficients) > 0L) {
    common_columns <- intersect(names(object$effect$coefficients), names(coefficient_vector))
    coefficient_vector[common_columns] <- object$effect$coefficients[common_columns]
  }

  tau_observed <- tau_matrix[sample_indices]
  tau_model <- as.vector(design %*% coefficient_vector)
  sample_weights <- if (is.null(weight_matrix)) rep(1, length(tau_observed)) else weight_matrix[sample_indices]
  pseudo_tau <- tau_model + sample_weights * (tau_observed - tau_model)

  linear_indices <- (sample_indices[, 1L] - 1L) * prepared$n_times + sample_indices[, 2L]
  sample_table <- prepared$data[linear_indices, c(names(prepared$data)), drop = FALSE]
  sample_table[, colnames(design)] <- design
  sample_table$.cdmc_fold_id <- fold_id
  sample_table$.cdmc_tau_observed <- tau_observed
  sample_table$.cdmc_tau_model <- tau_model
  sample_table$.cdmc_tau_dr <- pseudo_tau
  sample_table$.cdmc_weight_used <- sample_weights

  list(
    sample_mask = sample_mask,
    design = design,
    pseudo_tau = pseudo_tau,
    tau_observed = tau_observed,
    tau_model = tau_model,
    sample_weights = sample_weights,
    sample_indices = sample_indices,
    sample_table = sample_table
  )
}

cdmc_dr_fit <- function(
  data,
  outcome,
  dose,
  unit,
  time,
  covariates = NULL,
  weights = NULL,
  weight_method = NULL,
  weight_covariates = covariates,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest"),
  gps_df = 4L,
  gps_spline_covariates = weight_covariates,
  gps_bandwidth = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  cbps_standardize = TRUE,
  cbps_method = c("over", "exact"),
  cbps_iterations = 1000L,
  cbps_twostep = TRUE,
  stabilize_weights = TRUE,
  max_weight = NULL,
  n_folds = 2L,
  fold_assignments = NULL,
  lambda = NULL,
  rank_max = NULL,
  lambda_fraction = 0.25,
  lambda_selection = c("heuristic", "cv"),
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  washout = 0L,
  lag_order = 0L,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  zero_tolerance = 1e-8,
  seed = NULL,
  verbose = FALSE
) {
  n_folds_requested <- if (missing(n_folds)) NULL else n_folds
  lambda_selection <- match.arg(lambda_selection)
  gps_model <- match.arg(gps_model)
  cbps_method <- match.arg(cbps_method)
  weight_method <- cdmc_resolve_dr_weight_method(weights = weights, weight_method = weight_method)

  if (!is.null(weight_covariates) && !is.character(weight_covariates)) {
    stop("weight_covariates must be NULL or a character vector of column names.", call. = FALSE)
  }

  if (!is.numeric(gps_df) || length(gps_df) != 1L || !is.finite(gps_df) || gps_df < 1) {
    stop("gps_df must be a positive integer.", call. = FALSE)
  }
  gps_df <- as.integer(gps_df)
  if (identical(gps_model, "gam") && gps_df < 3L) {
    stop("gps_df must be at least 3 when gps_model = 'gam'.", call. = FALSE)
  }
  gps_spline_covariates <- cdmc_resolve_gps_spline_covariates(weight_covariates, gps_spline_covariates)
  gps_bandwidth <- cdmc_resolve_gps_bandwidth(gps_bandwidth)
  gps_forest_trees <- cdmc_resolve_gps_forest_trees(gps_forest_trees)
  gps_forest_mtry <- cdmc_resolve_gps_forest_mtry(gps_forest_mtry)
  gps_forest_min_node_size <- cdmc_resolve_gps_forest_min_node_size(gps_forest_min_node_size)
  cbps_iterations <- cdmc_resolve_cbps_iterations(cbps_iterations)

  if (!is.logical(cbps_standardize) || length(cbps_standardize) != 1L || is.na(cbps_standardize)) {
    stop("cbps_standardize must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(cbps_twostep) || length(cbps_twostep) != 1L || is.na(cbps_twostep)) {
    stop("cbps_twostep must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  full_prepared <- cdmc_prepare_panel(
    data = data,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    zero_tolerance = zero_tolerance
  )
  full_data <- full_prepared$data
  full_data$.cdmc_global_unit_index <- full_data$.cdmc_unit_index
  full_data$.cdmc_global_time_index <- full_data$.cdmc_time_index

  full_eligible_mask <- cdmc_build_eligible_mask(
    dose_matrix = full_prepared$dose_matrix,
    zero_tolerance = zero_tolerance,
    washout = washout
  )
  if (identical(weight_method, "external")) {
    weight_seed_info <- cdmc_prepare_panel_weights(
      weights = weights,
      data = full_data,
      n_units = full_prepared$n_units,
      n_times = full_prepared$n_times,
      eligible_mask = full_eligible_mask,
      column_name = ".cdmc_dr_weight_raw"
    )
    full_data <- weight_seed_info$data
    full_eligible_mask <- full_eligible_mask & weight_seed_info$matrix > 0
  } else {
    missing_weight_covariates <- setdiff(weight_covariates, names(full_data))
    if (length(missing_weight_covariates) > 0L) {
      stop(
        sprintf("Missing weight_covariates columns: %s.", paste(missing_weight_covariates, collapse = ", ")),
        call. = FALSE
      )
    }
    if (identical(weight_method, "cbps") && !gps_model %in% c("linear", "spline")) {
      stop(
        "gps_model is not available when weight_method = 'cbps'. Use 'linear' or 'spline', or switch to gaussian_gps or kernel_gps.",
        call. = FALSE
      )
    }
    if (identical(weight_method, "cbps") &&
        length(weight_covariates %||% character(0)) == 0L &&
        !isTRUE(gps_time_effects)) {
      stop(
        "Internal CBPS weighting requires at least one weight_covariate or gps_time_effects = TRUE.",
        call. = FALSE
      )
    }
  }
  cdmc_validate_support(full_eligible_mask)

  fold_assignments <- cdmc_resolve_unit_folds(
    unit_levels = full_prepared$unit_levels,
    n_folds = if (is.null(fold_assignments)) n_folds else n_folds_requested,
    fold_assignments = fold_assignments
  )
  n_folds <- length(unique(fold_assignments$fold))
  fold_lookup <- setNames(fold_assignments$fold, fold_assignments$unit)

  baseline_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_model_oof <- matrix(0, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  tau_dr_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  fitted_dr_oof <- matrix(0, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  sample_mask_oof <- matrix(FALSE, nrow = full_prepared$n_units, ncol = full_prepared$n_times)
  weight_oof <- matrix(NA_real_, nrow = full_prepared$n_units, ncol = full_prepared$n_times)

  design_parts <- vector("list", n_folds)
  pseudo_parts <- vector("list", n_folds)
  fold_samples <- vector("list", n_folds)
  fold_summaries <- vector("list", n_folds)
  sample_counter <- 0L

  for (fold_index in seq_len(n_folds)) {
    holdout_units <- fold_assignments$unit[fold_assignments$fold == fold_index]
    train_data <- full_data[!full_data[[unit]] %in% holdout_units, , drop = FALSE]
    holdout_data <- full_data[full_data[[unit]] %in% holdout_units, , drop = FALSE]

    if (verbose) {
      message(sprintf("cross-fit fold %d/%d: %d holdout units", fold_index, n_folds, length(holdout_units)))
    }

    gps_fit <- NULL
    holdout_gps_fit <- NULL
    if (identical(weight_method, "cbps")) {
      gps_fit <- cdmc_fit_cbps_weights(
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates,
        cbps_standardize = cbps_standardize,
        cbps_method = cbps_method,
        cbps_iterations = cbps_iterations,
        cbps_twostep = cbps_twostep,
        max_weight = max_weight
      )
      holdout_gps_fit <- cdmc_fit_cbps_weights(
        data = holdout_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates,
        cbps_standardize = cbps_standardize,
        cbps_method = cbps_method,
        cbps_iterations = cbps_iterations,
        cbps_twostep = cbps_twostep,
        max_weight = max_weight
      )
      train_data$.cdmc_dr_weight_raw <- gps_fit$weights
      holdout_data$.cdmc_dr_weight_raw <- holdout_gps_fit$weights
    } else if (!identical(weight_method, "external")) {
      gps_fit <- cdmc_fit_internal_gps(
        weight_method = weight_method,
        data = train_data,
        dose = dose,
        covariates = weight_covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = gps_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates,
        gps_bandwidth = gps_bandwidth,
        gps_forest_trees = gps_forest_trees,
        gps_forest_mtry = gps_forest_mtry,
        gps_forest_min_node_size = gps_forest_min_node_size
      )
      train_data$.cdmc_dr_weight_raw <- cdmc_predict_internal_gps_weights(
        object = gps_fit,
        newdata = train_data,
        dose = dose,
        stabilize = stabilize_weights,
        max_weight = max_weight
      )
      holdout_data$.cdmc_dr_weight_raw <- cdmc_predict_internal_gps_weights(
        object = gps_fit,
        newdata = holdout_data,
        dose = dose,
        stabilize = stabilize_weights,
        max_weight = max_weight
      )
    }

    fold_fit <- cdmc_fit(
      data = train_data,
      outcome = outcome,
      dose = dose,
      unit = unit,
      time = time,
      covariates = covariates,
      weights = ".cdmc_dr_weight_raw",
      lambda = lambda,
      rank_max = rank_max,
      lambda_fraction = lambda_fraction,
      lambda_selection = lambda_selection,
      lambda_grid = lambda_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      washout = washout,
      lag_order = lag_order,
      effect_model = "linear",
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      zero_tolerance = zero_tolerance,
      seed = NULL,
      verbose = FALSE
    )

    holdout_prepared <- cdmc_prepare_panel(
      data = holdout_data,
      outcome = outcome,
      dose = dose,
      unit = unit,
      time = time,
      covariates = covariates,
      zero_tolerance = zero_tolerance
    )
    holdout_eligible_mask <- cdmc_build_eligible_mask(
      dose_matrix = holdout_prepared$dose_matrix,
      zero_tolerance = zero_tolerance,
      washout = washout
    )
    holdout_weight_info <- cdmc_prepare_panel_weights(
      weights = ".cdmc_dr_weight_raw",
      data = holdout_prepared$data,
      n_units = holdout_prepared$n_units,
      n_times = holdout_prepared$n_times,
      eligible_mask = holdout_eligible_mask,
      column_name = ".cdmc_weight"
    )
    holdout_prepared$data <- holdout_weight_info$data
    holdout_weight_matrix <- holdout_weight_info$matrix
    holdout_eligible_mask <- holdout_eligible_mask & holdout_weight_matrix > 0

    holdout_baseline <- cdmc_predict_holdout_baseline(
      object = fold_fit,
      prepared = holdout_prepared,
      eligible_mask = holdout_eligible_mask,
      weight_matrix = holdout_weight_matrix
    )
    fold_sample <- cdmc_prepare_dr_fold_sample(
      object = fold_fit,
      prepared = holdout_prepared,
      baseline_hat = holdout_baseline$baseline_hat,
      weight_matrix = holdout_weight_matrix,
      fold_id = fold_index
    )

    global_indices <- cbind(
      holdout_prepared$data$.cdmc_global_unit_index,
      holdout_prepared$data$.cdmc_global_time_index
    )
    baseline_oof[global_indices] <- cdmc_flatten_matrix(holdout_baseline$baseline_hat)
    tau_oof[global_indices] <- cdmc_flatten_matrix(holdout_prepared$y_matrix - holdout_baseline$baseline_hat)
    weight_oof[global_indices] <- holdout_prepared$data$.cdmc_weight

    if (nrow(fold_sample$design) > 0L) {
      sample_counter <- sample_counter + 1L
      design_parts[[sample_counter]] <- fold_sample$design
      pseudo_parts[[sample_counter]] <- fold_sample$pseudo_tau
      fold_samples[[sample_counter]] <- fold_sample$sample_table

      sample_linear_indices <- (fold_sample$sample_indices[, 1L] - 1L) * holdout_prepared$n_times + fold_sample$sample_indices[, 2L]
      sample_global_indices <- cbind(
        holdout_prepared$data$.cdmc_global_unit_index[sample_linear_indices],
        holdout_prepared$data$.cdmc_global_time_index[sample_linear_indices]
      )
      sample_mask_oof[sample_global_indices] <- TRUE
      tau_model_oof[sample_global_indices] <- fold_sample$tau_model
      tau_dr_oof[sample_global_indices] <- fold_sample$pseudo_tau
    }

    fold_summaries[[fold_index]] <- list(
      fold = fold_index,
      holdout_units = holdout_units,
      weight_method = weight_method,
      lambda = fold_fit$lambda,
      baseline_solver = fold_fit$baseline$solver,
      gps_formula = if (is.null(gps_fit)) NULL else gps_fit$formula,
      gps_weight_method = if (is.null(gps_fit)) NULL else gps_fit$weight_method,
      gps_model = if (is.null(gps_fit)) NULL else gps_fit$gps_model,
      gps_df = if (is.null(gps_fit)) NULL else gps_fit$gps_df,
      gps_bandwidth = if (is.null(gps_fit)) NULL else gps_fit$gps_bandwidth,
      gps_forest_trees = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_trees,
      gps_forest_mtry = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_mtry,
      gps_forest_min_node_size = if (is.null(gps_fit)) NULL else gps_fit$gps_forest_min_node_size,
      gps_sigma_cond = if (is.null(gps_fit)) NULL else gps_fit$sigma_cond %||% NULL,
      cbps_requested_method = if (is.null(gps_fit) || !identical(weight_method, "cbps")) NULL else gps_fit$requested_method,
      cbps_train_method = if (is.null(gps_fit) || !identical(weight_method, "cbps")) NULL else gps_fit$method,
      cbps_holdout_method = if (is.null(holdout_gps_fit) || !identical(weight_method, "cbps")) NULL else holdout_gps_fit$method,
      cbps_standardize = if (identical(weight_method, "cbps")) cbps_standardize else NULL,
      cbps_iterations = if (identical(weight_method, "cbps")) cbps_iterations else NULL,
      cbps_twostep = if (identical(weight_method, "cbps")) cbps_twostep else NULL,
      sample_size = nrow(fold_sample$design)
    )
  }

  design_parts <- design_parts[seq_len(sample_counter)]
  pseudo_parts <- pseudo_parts[seq_len(sample_counter)]
  fold_samples <- fold_samples[seq_len(sample_counter)]

  if (sample_counter == 0L) {
    stop("No nonzero-dose histories were available for the doubly robust stage.", call. = FALSE)
  }

  final_design <- do.call(rbind, design_parts)
  final_pseudo_tau <- unlist(pseudo_parts, use.names = FALSE)
  final_effect_fit <- cdmc_fit_weighted_regression(
    design = final_design,
    response = final_pseudo_tau,
    weights = NULL
  )

  for (sample_table in fold_samples) {
    design <- as.matrix(sample_table[, paste0("dose_lag", seq.int(0L, lag_order)), drop = FALSE])
    fitted_values <- as.vector(design %*% final_effect_fit$coefficients)
    sample_global_indices <- cbind(
      sample_table$.cdmc_global_unit_index,
      sample_table$.cdmc_global_time_index
    )
    fitted_dr_oof[sample_global_indices] <- fitted_values
  }

  fitted_panel <- full_data[, setdiff(names(full_data), c(".cdmc_global_unit_index", ".cdmc_global_time_index")), drop = FALSE]
  fitted_panel$.cdmc_fold_id <- unname(fold_lookup[fitted_panel[[unit]]])
  fitted_panel$.cdmc_y0_oof <- cdmc_flatten_matrix(baseline_oof)
  fitted_panel$.cdmc_tau_oof <- cdmc_flatten_matrix(tau_oof)
  fitted_panel$.cdmc_tau_model_oof <- cdmc_flatten_matrix(tau_model_oof)
  fitted_panel$.cdmc_tau_dr <- cdmc_flatten_matrix(tau_dr_oof)
  fitted_panel$.cdmc_tau_dr_linear <- cdmc_flatten_matrix(fitted_dr_oof)
  fitted_panel$.cdmc_dr_sample <- cdmc_flatten_matrix(sample_mask_oof)
  fitted_panel$.cdmc_weight <- cdmc_flatten_matrix(weight_oof)
  fitted_panel$.cdmc_tau_model_oof[!fitted_panel$.cdmc_observed] <- NA_real_
  fitted_panel$.cdmc_tau_dr_linear[!fitted_panel$.cdmc_observed] <- NA_real_

  result <- list(
    call = match.call(),
    data = fitted_panel,
    outcome = outcome,
    dose = dose,
    unit = unit,
    time = time,
    covariates = covariates,
    weights = ".cdmc_weight",
    n_units = full_prepared$n_units,
    n_times = full_prepared$n_times,
    unit_levels = full_prepared$unit_levels,
    time_levels = full_prepared$time_levels,
    observed_mask = full_prepared$observed_mask,
    y_matrix = full_prepared$y_matrix,
    dose_matrix = full_prepared$dose_matrix,
    weight_matrix = weight_oof,
    fold_assignments = fold_assignments,
    n_folds = n_folds,
    weight_method = weight_method,
    weight_covariates = weight_covariates,
    gps_time_effects = gps_time_effects,
    gps_model = gps_model,
    gps_df = gps_df,
    gps_spline_covariates = gps_spline_covariates,
    gps_bandwidth = gps_bandwidth,
    gps_forest_trees = gps_forest_trees,
    gps_forest_mtry = gps_forest_mtry,
    gps_forest_min_node_size = gps_forest_min_node_size,
    cbps_standardize = cbps_standardize,
    cbps_method = cbps_method,
    cbps_iterations = cbps_iterations,
    cbps_twostep = cbps_twostep,
    stabilize_weights = stabilize_weights,
    max_weight = max_weight,
    rank_max = rank_max,
    washout = as.integer(washout),
    lag_order = as.integer(lag_order),
    lambda = lambda,
    lambda_selection = lambda_selection,
    zero_tolerance = zero_tolerance,
    fit_control = list(
      weights = if (identical(weight_method, "external")) ".cdmc_weight" else NULL,
      weight_method = weight_method,
      weight_covariates = weight_covariates,
      gps_time_effects = gps_time_effects,
      gps_model = gps_model,
      gps_df = gps_df,
      gps_spline_covariates = gps_spline_covariates,
      gps_bandwidth = gps_bandwidth,
      gps_forest_trees = gps_forest_trees,
      gps_forest_mtry = gps_forest_mtry,
      gps_forest_min_node_size = gps_forest_min_node_size,
      cbps_standardize = cbps_standardize,
      cbps_method = cbps_method,
      cbps_iterations = cbps_iterations,
      cbps_twostep = cbps_twostep,
      stabilize_weights = stabilize_weights,
      max_weight = max_weight,
      n_folds = n_folds,
      lambda = lambda,
      rank_max = rank_max,
      lambda_fraction = lambda_fraction,
      lambda_selection = lambda_selection,
      lambda_grid = lambda_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      washout = as.integer(washout),
      lag_order = as.integer(lag_order),
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      zero_tolerance = zero_tolerance
    ),
    baseline = list(
      baseline_hat = baseline_oof,
      fold_summaries = fold_summaries
    ),
    effect = list(
      coefficients = final_effect_fit$coefficients,
      fitted = fitted_dr_oof,
      tau = tau_oof,
      tau_model = tau_model_oof,
      tau_dr = tau_dr_oof,
      sample_mask = sample_mask_oof,
      design_columns = colnames(final_design),
      residuals = final_effect_fit$residuals
    )
  )

  class(result) <- "cdmc_dr_fit"
  result
}

print.cdmc_dr_fit <- function(x, ...) {
  cat("causaldosemc cross-fitted DR fit\n")
  cat(sprintf("  units x times: %d x %d\n", x$n_units, x$n_times))
  cat(sprintf("  folds: %d\n", x$n_folds))
  cat(sprintf("  weight method: %s\n", x$weight_method))
  if (identical(x$weight_method, "cbps")) {
    cat(sprintf("  cbps method: %s\n", x$cbps_method))
    cat(sprintf("  cbps standardized: %s\n", if (isTRUE(x$cbps_standardize)) "yes" else "no"))
  }
  if (identical(x$weight_method, "gaussian_gps") || identical(x$weight_method, "kernel_gps")) {
    cat(sprintf("  gps model: %s\n", x$gps_model))
    if (identical(x$gps_model, "spline")) {
      cat(sprintf("  gps spline df: %d\n", x$gps_df))
    }
    if (identical(x$gps_model, "gam")) {
      cat(sprintf("  gps smooth basis k: %d\n", x$gps_df))
    }
    if (identical(x$gps_model, "forest")) {
      cat(sprintf("  gps forest trees: %d\n", x$gps_forest_trees))
    }
    if (identical(x$weight_method, "kernel_gps")) {
      cat(sprintf("  gps kernel bandwidth: %.6g\n", x$gps_bandwidth))
    }
  }
  cat(sprintf("  DR sample cells: %d\n", sum(x$effect$sample_mask, na.rm = TRUE)))

  if (length(x$effect$coefficients) > 0L) {
    cat("  DR effect coefficients:\n")
    for (coefficient_name in names(x$effect$coefficients)) {
      cat(sprintf("    %s = %.6g\n", coefficient_name, x$effect$coefficients[[coefficient_name]]))
    }
  }

  invisible(x)
}