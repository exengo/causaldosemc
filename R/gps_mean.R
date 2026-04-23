# GPS conditional-mean models used by the DR pipeline.
# Extracted from dr_fit.R for code locality with PLR (Robinson) score.

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
  gps_model = c("linear", "spline", "gam", "tree", "forest", "boost"),
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

cdmc_available_gps_stack_models <- function() {
  models <- c("linear", "spline")

  if (requireNamespace("mgcv", quietly = TRUE)) {
    models <- c(models, "gam")
  }
  if (requireNamespace("rpart", quietly = TRUE)) {
    models <- c(models, "tree")
  }
  if (requireNamespace("ranger", quietly = TRUE)) {
    models <- c(models, "forest")
  }
  if (requireNamespace("gbm", quietly = TRUE)) {
    models <- c(models, "boost")
  }

  models
}

cdmc_resolve_gps_stack_models <- function(gps_stack_models = NULL) {
  allowed_models <- c("linear", "spline", "gam", "tree", "forest", "boost")

  if (is.null(gps_stack_models)) {
    return(cdmc_available_gps_stack_models())
  }

  if (!is.character(gps_stack_models) || length(gps_stack_models) < 1L) {
    stop("gps_stack_models must be NULL or a nonempty character vector of base learner names.", call. = FALSE)
  }

  gps_stack_models <- unique(gps_stack_models)
  unknown_models <- setdiff(gps_stack_models, allowed_models)
  if (length(unknown_models) > 0L) {
    stop(
      sprintf("Unknown gps_stack_models entries: %s.", paste(unknown_models, collapse = ", ")),
      call. = FALSE
    )
  }

  if ("gam" %in% gps_stack_models) {
    cdmc_assert_installed("mgcv")
  }
  if ("tree" %in% gps_stack_models) {
    cdmc_assert_installed("rpart")
  }
  if ("forest" %in% gps_stack_models) {
    cdmc_assert_installed("ranger")
  }
  if ("boost" %in% gps_stack_models) {
    cdmc_assert_installed("gbm")
  }

  gps_stack_models
}

cdmc_fit_gps_stack_weights <- function(prediction_matrix, response) {
  prediction_matrix <- as.matrix(prediction_matrix)
  if (ncol(prediction_matrix) == 1L) {
    weights <- 1
    names(weights) <- colnames(prediction_matrix)
    return(weights)
  }

  fit <- stats::lm.fit(x = prediction_matrix, y = response)
  weights <- as.numeric(fit$coefficients)
  weights[!is.finite(weights)] <- 0
  weights <- pmax(weights, 0)

  if (sum(weights) <= sqrt(.Machine$double.eps)) {
    weights <- rep(1 / ncol(prediction_matrix), ncol(prediction_matrix))
  } else {
    weights <- weights / sum(weights)
  }

  names(weights) <- colnames(prediction_matrix)
  weights
}

cdmc_build_stack_gps_formula <- function(dose, stack_models) {
  stats::as.formula(
    sprintf("%s ~ %s", dose, paste(paste0("stack_", stack_models), collapse = " + "))
  )
}

cdmc_predict_gps_mean_model <- function(object, newdata) {
  if (inherits(object, "cdmc_gps_constant_mean")) {
    return(rep(object$mean, nrow(newdata)))
  }
  if (inherits(object, "cdmc_gps_stack_mean")) {
    prediction_matrix <- do.call(
      cbind,
      lapply(object$components, function(component) {
        cdmc_predict_gps_mean_model(
          component$fit,
          newdata = component$prepare_newdata(newdata)
        )
      })
    )
    prediction_matrix <- as.matrix(prediction_matrix)
    colnames(prediction_matrix) <- object$model_names
    return(as.numeric(prediction_matrix %*% object$weights))
  }
  if (inherits(object, "gbm")) {
    return(as.numeric(stats::predict(object, newdata = newdata, n.trees = object$n.trees, type = "response")))
  }
  if (inherits(object, "ranger")) {
    return(as.numeric(stats::predict(object, data = newdata)$predictions))
  }

  as.numeric(stats::predict(object, newdata = newdata))
}

cdmc_fit_gps_mean_model <- function(
  data,
  dose,
  covariates,
  time,
  gps_time_effects = TRUE,
  gps_model = c("linear", "spline", "gam", "tree", "forest", "stack", "boost"),
  gps_df = 4L,
  gps_spline_covariates = covariates,
  gps_stack_models = NULL,
  gps_forest_trees = 200L,
  gps_forest_mtry = NULL,
  gps_forest_min_node_size = NULL,
  gps_boost_trees = 200L,
  gps_boost_depth = 2L,
  gps_boost_shrinkage = 0.05,
  gps_boost_min_obs_node = 10L
) {
  gps_model <- match.arg(gps_model)
  if (identical(gps_model, "stack")) {
    stack_models <- cdmc_resolve_gps_stack_models(gps_stack_models)
    component_models <- lapply(stack_models, function(stack_model) {
      cdmc_fit_gps_mean_model(
        data = data,
        dose = dose,
        covariates = covariates,
        time = time,
        gps_time_effects = gps_time_effects,
        gps_model = stack_model,
        gps_df = gps_df,
        gps_spline_covariates = gps_spline_covariates,
        gps_stack_models = NULL,
        gps_forest_trees = gps_forest_trees,
        gps_forest_mtry = gps_forest_mtry,
        gps_forest_min_node_size = gps_forest_min_node_size,
        gps_boost_trees = gps_boost_trees,
        gps_boost_depth = gps_boost_depth,
        gps_boost_shrinkage = gps_boost_shrinkage,
        gps_boost_min_obs_node = gps_boost_min_obs_node
      )
    })

    prediction_matrix <- do.call(
      cbind,
      lapply(component_models, function(component_model) {
        cdmc_predict_gps_mean_model(
          component_model$fit,
          newdata = component_model$prepare_newdata(data)
        )
      })
    )
    prediction_matrix <- as.matrix(prediction_matrix)
    colnames(prediction_matrix) <- stack_models
    fit_object <- structure(
      list(
        model_names = stack_models,
        weights = cdmc_fit_gps_stack_weights(prediction_matrix, response = data[[dose]]),
        components = lapply(component_models, function(component_model) {
          list(
            fit = component_model$fit,
            prepare_newdata = component_model$prepare_newdata
          )
        })
      ),
      class = "cdmc_gps_stack_mean"
    )

    return(list(
      formula = cdmc_build_stack_gps_formula(dose, stack_models),
      fit = fit_object,
      prepare_newdata = function(newdata) newdata,
      gps_stack_models = stack_models
    ))
  }

  if (identical(gps_model, "boost")) {
    cdmc_assert_installed("gbm")
    boost_frame <- cdmc_build_gps_forest_frame(
      data = data,
      dose = dose,
      covariates = covariates,
      time = time,
      gps_time_effects = gps_time_effects
    )

    if (length(boost_frame$predictor_names) == 0L) {
      fit_object <- structure(list(mean = mean(data[[dose]], na.rm = TRUE)), class = "cdmc_gps_constant_mean")
    } else {
      fit_object <- gbm::gbm(
        formula = boost_frame$formula,
        data = boost_frame$data,
        distribution = "gaussian",
        n.trees = gps_boost_trees,
        interaction.depth = gps_boost_depth,
        shrinkage = gps_boost_shrinkage,
        n.minobsinnode = gps_boost_min_obs_node,
        bag.fraction = 1,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE,
        n.cores = 1L
      )
    }

    return(list(
      formula = boost_frame$formula,
      fit = fit_object,
      prepare_newdata = function(newdata) {
        boost_newdata <- cdmc_build_gps_forest_frame(
          data = newdata,
          dose = dose,
          covariates = covariates,
          time = time,
          gps_time_effects = gps_time_effects
        )$data
        boost_newdata[, boost_frame$predictor_names, drop = FALSE]
      },
      gps_stack_models = NULL
    ))
  }

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
      },
      gps_stack_models = NULL
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
    prepare_newdata = function(newdata) newdata,
    gps_stack_models = NULL
  )
}
