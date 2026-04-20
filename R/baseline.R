cdmc_spectral_norm <- function(x) {
  if (!is.matrix(x)) {
    stop("x must be a matrix.", call. = FALSE)
  }

  if (all(abs(x) <= sqrt(.Machine$double.eps))) {
    return(0)
  }

  base::svd(x, nu = 0L, nv = 0L)$d[[1L]]
}

cdmc_shrink_singular_values <- function(x, lambda, rank_max) {
  # Use truncated SVD (irlba) when the matrix is large and rank_max is small
  # relative to matrix dimensions — avoids computing all singular values/vectors.
  n_min <- min(dim(x))
  use_truncated <- n_min > 50L &&
                   rank_max < n_min %/% 3L &&
                   requireNamespace("irlba", quietly = TRUE)

  svd_fit <- if (use_truncated) {
    k <- min(rank_max + 3L, n_min - 1L)
    tryCatch(
      irlba::irlba(x, nv = k, nu = k, maxit = 300L),
      error = function(e) base::svd(x)
    )
  } else {
    base::svd(x)
  }

  shrunk <- pmax(svd_fit$d - lambda, 0)
  keep <- which(shrunk > sqrt(.Machine$double.eps))
  if (length(keep) == 0L) {
    return(list(
      matrix = matrix(0, nrow = nrow(x), ncol = ncol(x)),
      u = NULL,
      d = numeric(0),
      v = NULL
    ))
  }

  keep <- keep[seq_len(min(length(keep), rank_max))]
  u <- svd_fit$u[, keep, drop = FALSE]
  v <- svd_fit$v[, keep, drop = FALSE]
  d <- shrunk[keep]

  list(
    matrix = sweep(u, 2, d, FUN = "*") %*% t(v),
    u = u,
    d = d,
    v = v
  )
}

cdmc_fit_weighted_low_rank <- function(
  residual_matrix,
  weight_matrix,
  lambda,
  rank_max,
  start = NULL,
  maxit = 100L,
  tol = 1e-5,
  verbose = FALSE
) {
  if (!is.matrix(weight_matrix) || !identical(dim(weight_matrix), dim(residual_matrix))) {
    stop("weight_matrix must be a numeric matrix matching residual_matrix.", call. = FALSE)
  }

  positive_weights <- weight_matrix[weight_matrix > 0]
  if (length(positive_weights) == 0L) {
    stop("weight_matrix must contain at least one positive entry.", call. = FALSE)
  }

  l_matrix <- start %||% matrix(0, nrow = nrow(residual_matrix), ncol = ncol(residual_matrix))
  step_size <- 1 / max(positive_weights)
  converged <- FALSE
  relative_change <- Inf
  iteration <- 0L
  svt_fit <- list(u = NULL, d = numeric(0), v = NULL)

  for (iteration in seq_len(maxit)) {
    update_matrix <- l_matrix - step_size * (weight_matrix * (l_matrix - residual_matrix))
    svt_fit <- cdmc_shrink_singular_values(
      x = update_matrix,
      lambda = step_size * lambda,
      rank_max = rank_max
    )

    relative_change <- cdmc_relative_change(svt_fit$matrix, l_matrix)
    if (verbose) {
      message(sprintf(
        "    weighted svt iteration %d: relative low-rank change = %.6g",
        iteration,
        relative_change
      ))
    }

    l_matrix <- svt_fit$matrix
    if (relative_change < tol) {
      converged <- TRUE
      break
    }
  }

  c(
    svt_fit,
    list(
      matrix = l_matrix,
      iterations = iteration,
      converged = converged,
      relative_change = relative_change
    )
  )
}

cdmc_default_lambda <- function(
  y_matrix,
  x_matrices,
  mask,
  weight_matrix = NULL,
  lambda_fraction = 0.25,
  fe_maxit = 200L,
  fe_tol = 1e-8
) {
  cdmc_assert_installed("softImpute")

  nuisance_zero <- cdmc_fit_nuisance(
    y_matrix = y_matrix,
    l_matrix = matrix(0, nrow = nrow(y_matrix), ncol = ncol(y_matrix)),
    x_matrices = x_matrices,
    mask = mask,
    weights = weight_matrix,
    maxit = fe_maxit,
    tol = fe_tol
  )

  residual_matrix <- y_matrix - nuisance_zero$fitted
  if (is.null(weight_matrix)) {
    residual_matrix[!mask] <- NA_real_
    return(lambda_fraction * softImpute::lambda0(residual_matrix))
  }

  residual_matrix[!mask] <- 0
  lambda_fraction * cdmc_spectral_norm((weight_matrix * mask) * residual_matrix)
}

cdmc_fit_baseline <- function(
  y_matrix,
  x_matrices,
  mask,
  weight_matrix = NULL,
  lambda,
  rank_max,
  fit_start = NULL,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  verbose = FALSE
) {
  cdmc_assert_installed("softImpute")

  n_units <- nrow(y_matrix)
  n_times <- ncol(y_matrix)
  optimization_weights <- if (is.null(weight_matrix)) NULL else weight_matrix * mask

  fit_start <- if (is.list(fit_start)) fit_start else NULL
  l_matrix <- fit_start$l_matrix %||% matrix(0, nrow = n_units, ncol = n_times)
  nuisance_start <- fit_start$nuisance %||% NULL
  warm_start <- fit_start$warm %||% NULL
  converged <- FALSE
  relative_change <- Inf
  soft_fit <- NULL
  iteration <- 0L

  pb <- if (verbose) cdmc_progress_bar("Baseline fitting", total = outer_maxit) else NULL

  for (iteration in seq_len(outer_maxit)) {
    nuisance_fit <- cdmc_fit_nuisance(
      y_matrix = y_matrix,
      l_matrix = l_matrix,
      x_matrices = x_matrices,
      mask = mask,
      weights = weight_matrix,
      start = nuisance_start,
      maxit = fe_maxit,
      tol = fe_tol
    )

    residual_matrix <- y_matrix - nuisance_fit$fitted
    if (is.null(optimization_weights)) {
      residual_matrix[!mask] <- NA_real_

      soft_fit <- softImpute::softImpute(
        x = residual_matrix,
        rank.max = rank_max,
        lambda = lambda,
        type = "svd",
        thresh = tol,
        maxit = soft_maxit,
        trace.it = verbose,
        warm.start = warm_start
      )

      l_new <- cdmc_reconstruct_low_rank(
        fit = soft_fit,
        n_units = n_units,
        n_times = n_times
      )
    } else {
      residual_matrix[!mask] <- 0
      soft_fit <- cdmc_fit_weighted_low_rank(
        residual_matrix = residual_matrix,
        weight_matrix = optimization_weights,
        lambda = lambda,
        rank_max = rank_max,
        start = if (!is.null(warm_start$matrix)) warm_start$matrix else l_matrix,
        maxit = soft_maxit,
        tol = tol,
        verbose = verbose
      )
      l_new <- soft_fit$matrix
    }

    relative_change <- cdmc_relative_change(l_new, l_matrix)

    l_matrix <- l_new
    nuisance_start <- nuisance_fit
    warm_start <- soft_fit

    if (verbose) {
      if (!is.null(pb)) {
        cdmc_progress_update(pb)
      } else {
        message(sprintf(
          "outer iteration %d: relative low-rank change = %.6g",
          iteration,
          relative_change
        ))
      }
    }

    if (relative_change < tol) {
      converged <- TRUE
      break
    }
  }

  cdmc_progress_done(pb)

  final_nuisance <- cdmc_fit_nuisance(
    y_matrix = y_matrix,
    l_matrix = l_matrix,
    x_matrices = x_matrices,
    mask = mask,
    weights = weight_matrix,
    start = nuisance_start,
    maxit = fe_maxit,
    tol = fe_tol
  )

  list(
    low_rank = l_matrix,
    nuisance = final_nuisance,
    baseline_hat = final_nuisance$fitted + l_matrix,
    lambda = lambda,
    rank_max = rank_max,
    outer_iterations = iteration,
    converged = converged,
    relative_change = relative_change,
    soft_fit = soft_fit,
    fit_start = list(
      l_matrix = l_matrix,
      nuisance = final_nuisance,
      warm = soft_fit
    ),
    solver = if (is.null(optimization_weights)) "softImpute" else "weighted_svt",
    effective_rank = if (is.null(soft_fit)) 0L else sum(soft_fit$d > sqrt(.Machine$double.eps))
  )
}
