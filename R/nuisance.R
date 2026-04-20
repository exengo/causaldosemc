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

  # --- Loop-invariant pre-computations ---
  # (1) Weighted mask, its transpose, and per-axis counts.
  # Using matrix-vector products instead of broadcasting alpha/beta into full
  # n×T matrices saves ~2 O(n×T) allocations per inner iteration.
  weighted_mask <- if (is.null(weights)) mask else weights * mask
  weighted_mask[!mask] <- 0
  weighted_mask_t <- t(weighted_mask)
  counts_row <- rowSums(weighted_mask)
  counts_col <- colSums(weighted_mask)

  if (any(counts_row <= sqrt(.Machine$double.eps))) {
    stop(
      "Cannot estimate unit fixed effects because at least one unit has no eligible observations.",
      call. = FALSE
    )
  }
  if (any(counts_col <= sqrt(.Machine$double.eps))) {
    stop(
      "Cannot estimate time fixed effects because at least one time period has no eligible observations.",
      call. = FALSE
    )
  }

  # (2) Gamma estimation pre-computations.
  # mask, x_matrices, and l_matrix are all fixed for the duration of this call.
  # Pre-computing the design matrix and its QR factorization avoids repeating
  # which(mask, arr.ind=TRUE), cbind/lapply extractions, and lm.fit QR on every
  # inner iteration — together these account for ~20% of nuisance iteration time.
  if (n_covariates > 0L) {
    obs_idx <- which(mask, arr.ind = TRUE)
    obs_row <- obs_idx[, 1L]
    obs_col <- obs_idx[, 2L]
    x_at_obs <- do.call(cbind, lapply(x_matrices, function(x) x[obs_idx]))
    yl_at_obs <- y_matrix[obs_idx] - l_matrix[obs_idx]
    weights_at_obs <- if (!is.null(weights)) weights[obs_idx] else NULL
    sw <- if (!is.null(weights_at_obs)) sqrt(pmax(weights_at_obs, 0)) else NULL
    gamma_qr <- qr(if (is.null(sw)) x_at_obs else x_at_obs * sw)
  } else {
    obs_idx <- obs_row <- obs_col <- yl_at_obs <- weights_at_obs <- sw <- gamma_qr <- NULL
  }

  converged <- FALSE
  iteration <- 0L

  for (iteration in seq_len(maxit)) {
    adjusted <- y_matrix - l_matrix - x_beta
    adjusted[!mask] <- 0

    # Compute adjusted * weighted_mask once and reuse for both row and col sums.
    adj_weighted <- adjusted * weighted_mask
    adj_row_sums <- rowSums(adj_weighted)
    adj_col_sums <- colSums(adj_weighted)

    # Row means of (adjusted - beta_broadcast): use matrix-vector product.
    alpha_new <- (adj_row_sums - as.vector(weighted_mask %*% beta)) / counts_row

    # Col means of (adjusted - alpha_new_broadcast): same approach, pre-shift alpha.
    beta_new <- (adj_col_sums - as.vector(weighted_mask_t %*% alpha_new)) / counts_col

    shift <- mean(alpha_new)
    alpha_new <- alpha_new - shift
    beta_new  <- beta_new  + shift

    if (n_covariates > 0L) {
      response <- yl_at_obs - alpha_new[obs_row] - beta_new[obs_col]
      rhs <- if (is.null(sw)) response else response * sw
      gamma_new <- qr.coef(gamma_qr, rhs)
      gamma_new[is.na(gamma_new)] <- 0
      names(gamma_new) <- names(x_matrices)
    } else {
      gamma_new <- numeric(0)
    }

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
    beta  <- beta_new
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
