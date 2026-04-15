cdmc_masked_row_means <- function(x, mask, weights = NULL) {
  weighted_mask <- if (is.null(weights)) mask else weights * mask
  x[!mask] <- 0
  counts <- rowSums(weighted_mask)
  if (any(counts <= sqrt(.Machine$double.eps))) {
    stop(
      "Cannot estimate unit fixed effects because at least one unit has no eligible observations.",
      call. = FALSE
    )
  }

  rowSums(x * weighted_mask) / counts
}

cdmc_masked_col_means <- function(x, mask, weights = NULL) {
  weighted_mask <- if (is.null(weights)) mask else weights * mask
  x[!mask] <- 0
  counts <- colSums(weighted_mask)
  if (any(counts <= sqrt(.Machine$double.eps))) {
    stop(
      "Cannot estimate time fixed effects because at least one time period has no eligible observations.",
      call. = FALSE
    )
  }

  colSums(x * weighted_mask) / counts
}

cdmc_covariate_contribution <- function(x_matrices, gamma, n_units, n_times) {
  contribution <- matrix(0, nrow = n_units, ncol = n_times)
  if (length(x_matrices) == 0L) {
    return(contribution)
  }

  for (index in seq_along(x_matrices)) {
    contribution <- contribution + gamma[index] * x_matrices[[index]]
  }

  contribution
}

cdmc_estimate_gamma <- function(y_matrix, l_matrix, x_matrices, alpha, beta, mask, weights = NULL) {
  if (length(x_matrices) == 0L) {
    return(numeric(0))
  }

  observed_indices <- which(mask, arr.ind = TRUE)
  response <- y_matrix[observed_indices] -
    l_matrix[observed_indices] -
    alpha[observed_indices[, 1L]] -
    beta[observed_indices[, 2L]]

  design <- do.call(
    cbind,
    lapply(x_matrices, function(x_matrix) x_matrix[observed_indices])
  )

  fit <- if (is.null(weights)) {
    stats::lm.fit(x = design, y = response)
  } else {
    stats::lm.wfit(x = design, y = response, w = weights[observed_indices])
  }
  coefficients <- fit$coefficients
  coefficients[is.na(coefficients)] <- 0
  names(coefficients) <- names(x_matrices)
  coefficients
}

cdmc_fit_nuisance <- function(
  y_matrix,
  l_matrix,
  x_matrices,
  mask,
  weights = NULL,
  start = NULL,
  maxit = 200L,
  tol = 1e-8
) {
  n_units <- nrow(y_matrix)
  n_times <- ncol(y_matrix)
  n_covariates <- length(x_matrices)

  alpha <- if (!is.null(start$alpha)) start$alpha else numeric(n_units)
  beta <- if (!is.null(start$beta)) start$beta else numeric(n_times)
  gamma <- if (!is.null(start$gamma)) start$gamma else rep(0, n_covariates)
  if (n_covariates > 0L) {
    names(gamma) <- names(x_matrices)
  } else {
    gamma <- numeric(0)
  }

  x_beta <- cdmc_covariate_contribution(x_matrices, gamma, n_units, n_times)
  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(maxit)) {
    adjusted <- y_matrix - l_matrix - x_beta

    alpha_new <- cdmc_masked_row_means(
      adjusted - matrix(beta, nrow = n_units, ncol = n_times, byrow = TRUE),
      mask,
      weights = weights
    )

    beta_new <- cdmc_masked_col_means(
      adjusted - matrix(alpha_new, nrow = n_units, ncol = n_times),
      mask,
      weights = weights
    )

    shift <- mean(alpha_new)
    alpha_new <- alpha_new - shift
    beta_new <- beta_new + shift

    gamma_new <- cdmc_estimate_gamma(
      y_matrix = y_matrix,
      l_matrix = l_matrix,
      x_matrices = x_matrices,
      alpha = alpha_new,
      beta = beta_new,
      mask = mask,
      weights = weights
    )

    x_beta_new <- cdmc_covariate_contribution(
      x_matrices = x_matrices,
      gamma = gamma_new,
      n_units = n_units,
      n_times = n_times
    )

    max_delta <- max(
      max(abs(alpha_new - alpha)),
      max(abs(beta_new - beta)),
      if (length(gamma_new) > 0L) max(abs(gamma_new - gamma)) else 0
    )

    alpha <- alpha_new
    beta <- beta_new
    gamma <- gamma_new
    x_beta <- x_beta_new

    if (max_delta < tol) {
      converged <- TRUE
      break
    }
  }

  fitted <- matrix(alpha, nrow = n_units, ncol = n_times) +
    matrix(beta, nrow = n_units, ncol = n_times, byrow = TRUE) +
    x_beta

  list(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    fitted = fitted,
    iterations = iteration,
    converged = converged
  )
}
