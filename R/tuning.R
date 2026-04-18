cdmc_resolve_lambda_grid <- function(
  lambda_grid,
  y_matrix,
  x_matrices,
  mask,
  weight_matrix = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  fe_maxit = 200L,
  fe_tol = 1e-8
) {
  if (!is.null(lambda_grid)) {
    if (!is.numeric(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) {
      stop("lambda_grid must contain finite nonnegative numeric values.", call. = FALSE)
    }

    resolved <- sort(unique(as.numeric(lambda_grid)), decreasing = TRUE)
    if (length(resolved) == 0L) {
      stop("lambda_grid must contain at least one value.", call. = FALSE)
    }

    return(resolved)
  }

  nlambda <- as.integer(nlambda)
  if (!is.numeric(nlambda) || length(nlambda) != 1L || nlambda < 1L) {
    stop("nlambda must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(lambda_min_ratio) || length(lambda_min_ratio) != 1L ||
      lambda_min_ratio <= 0 || lambda_min_ratio > 1) {
    stop("lambda_min_ratio must be a scalar in (0, 1].", call. = FALSE)
  }

  lambda_max <- cdmc_default_lambda(
    y_matrix = y_matrix,
    x_matrices = x_matrices,
    mask = mask,
    weight_matrix = weight_matrix,
    lambda_fraction = 1,
    fe_maxit = fe_maxit,
    fe_tol = fe_tol
  )

  if (!is.finite(lambda_max) || lambda_max < 0) {
    stop("Unable to compute a valid reference lambda for tuning.", call. = FALSE)
  }

  if (lambda_max <= sqrt(.Machine$double.eps) || nlambda == 1L) {
    return(lambda_max)
  }

  lambda_path <- lambda_max * exp(seq(0, log(lambda_min_ratio), length.out = nlambda))
  sort(unique(lambda_path), decreasing = TRUE)
}

cdmc_select_holdout_mask <- function(mask, block_size = 2L) {
  block_size <- as.integer(block_size)
  if (block_size < 1L) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }

  n_units <- nrow(mask)
  n_times <- ncol(mask)
  block_size <- min(block_size, n_times)

  holdout_mask <- matrix(FALSE, nrow = n_units, ncol = n_times)
  row_counts <- rowSums(mask)
  col_counts <- colSums(mask)

  for (unit_index in sample(seq_len(n_units))) {
    if (row_counts[unit_index] <= block_size) {
      next
    }

    max_start <- n_times - block_size + 1L
    if (max_start < 1L) {
      next
    }

    candidate_starts <- seq_len(max_start)
    candidate_starts <- candidate_starts[vapply(
      candidate_starts,
      function(start_index) {
        all(mask[unit_index, seq.int(start_index, start_index + block_size - 1L)])
      },
      logical(1)
    )]

    if (length(candidate_starts) == 0L) {
      next
    }

    for (start_index in sample(candidate_starts)) {
      time_indices <- seq.int(start_index, start_index + block_size - 1L)
      if (row_counts[unit_index] - block_size < 1L) {
        next
      }

      if (any(col_counts[time_indices] <= 1L)) {
        next
      }

      holdout_mask[unit_index, time_indices] <- TRUE
      row_counts[unit_index] <- row_counts[unit_index] - block_size
      col_counts[time_indices] <- col_counts[time_indices] - 1L
      break
    }
  }

  remaining_total <- sum(mask) - sum(holdout_mask)
  if ((!any(holdout_mask) || remaining_total <= (n_units + n_times)) && block_size > 1L) {
    return(cdmc_select_holdout_mask(mask = mask, block_size = block_size - 1L))
  }

  attr(holdout_mask, "block_size_used") <- if (any(holdout_mask)) block_size else 0L
  holdout_mask
}

cdmc_tune_lambda <- function(
  y_matrix,
  x_matrices,
  mask,
  weight_matrix = NULL,
  rank_max,
  lambda_grid = NULL,
  nlambda = 5L,
  lambda_min_ratio = 0.05,
  cv_rounds = 5L,
  cv_block_size = 2L,
  cv_workers = 1L,
  cv_top_k = NULL,
  cv_coarse_to_fine = FALSE,
  cv_coarse_nlambda = NULL,
  cv_warm_starts = FALSE,
  outer_maxit = 20L,
  fe_maxit = 200L,
  soft_maxit = 100L,
  tol = 1e-5,
  fe_tol = 1e-8,
  verbose = FALSE
) {
  lambda_grid <- cdmc_resolve_lambda_grid(
    lambda_grid = lambda_grid,
    y_matrix = y_matrix,
    x_matrices = x_matrices,
    mask = mask,
    weight_matrix = weight_matrix,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio,
    fe_maxit = fe_maxit,
    fe_tol = fe_tol
  )

  cv_rounds <- as.integer(cv_rounds)
  cv_block_size <- as.integer(cv_block_size)
  if (cv_rounds < 1L) {
    stop("cv_rounds must be a positive integer.", call. = FALSE)
  }
  if (cv_block_size < 1L) {
    stop("cv_block_size must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(cv_workers) || length(cv_workers) != 1L || !is.finite(cv_workers) || cv_workers < 1 || cv_workers != floor(cv_workers)) {
    stop("cv_workers must be a positive integer.", call. = FALSE)
  }
  cv_workers <- as.integer(cv_workers)
  if (!is.null(cv_top_k)) {
    if (!is.numeric(cv_top_k) || length(cv_top_k) != 1L || !is.finite(cv_top_k) || cv_top_k < 1 || cv_top_k != floor(cv_top_k)) {
      stop("cv_top_k must be NULL or a positive integer.", call. = FALSE)
    }
    cv_top_k <- as.integer(cv_top_k)
  }
  if (!is.logical(cv_coarse_to_fine) || length(cv_coarse_to_fine) != 1L || is.na(cv_coarse_to_fine)) {
    stop("cv_coarse_to_fine must be TRUE or FALSE.", call. = FALSE)
  }
  cv_coarse_to_fine <- isTRUE(cv_coarse_to_fine)
  if (!is.null(cv_coarse_nlambda)) {
    if (!is.numeric(cv_coarse_nlambda) || length(cv_coarse_nlambda) != 1L || !is.finite(cv_coarse_nlambda) || cv_coarse_nlambda < 3 || cv_coarse_nlambda != floor(cv_coarse_nlambda)) {
      stop("cv_coarse_nlambda must be NULL or an integer >= 3.", call. = FALSE)
    }
    cv_coarse_nlambda <- as.integer(cv_coarse_nlambda)
  }
  if (!is.logical(cv_warm_starts) || length(cv_warm_starts) != 1L || is.na(cv_warm_starts)) {
    stop("cv_warm_starts must be TRUE or FALSE.", call. = FALSE)
  }
  cv_warm_starts <- isTRUE(cv_warm_starts)

  if (cv_coarse_to_fine && length(lambda_grid) >= 5L) {
    coarse_n <- cv_coarse_nlambda %||% min(length(lambda_grid), max(3L, ceiling(sqrt(length(lambda_grid)))))
    coarse_n <- min(coarse_n, length(lambda_grid))

    coarse_indices <- sort(unique(round(seq(1, length(lambda_grid), length.out = coarse_n))))
    coarse_grid <- lambda_grid[coarse_indices]

    coarse_result <- cdmc_tune_lambda(
      y_matrix = y_matrix,
      x_matrices = x_matrices,
      mask = mask,
      weight_matrix = weight_matrix,
      rank_max = rank_max,
      lambda_grid = coarse_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_workers = cv_workers,
      cv_top_k = cv_top_k,
      cv_coarse_to_fine = FALSE,
      cv_coarse_nlambda = NULL,
      cv_warm_starts = cv_warm_starts,
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      verbose = verbose
    )

    coarse_selected_index <- which.min(abs(lambda_grid - coarse_result$selected_lambda))
    fine_indices <- seq.int(max(1L, coarse_selected_index - 1L), min(length(lambda_grid), coarse_selected_index + 1L))
    fine_grid <- lambda_grid[fine_indices]

    fine_result <- cdmc_tune_lambda(
      y_matrix = y_matrix,
      x_matrices = x_matrices,
      mask = mask,
      weight_matrix = weight_matrix,
      rank_max = rank_max,
      lambda_grid = fine_grid,
      nlambda = nlambda,
      lambda_min_ratio = lambda_min_ratio,
      cv_rounds = cv_rounds,
      cv_block_size = cv_block_size,
      cv_workers = cv_workers,
      cv_top_k = NULL,
      cv_coarse_to_fine = FALSE,
      cv_coarse_nlambda = NULL,
      cv_warm_starts = cv_warm_starts,
      outer_maxit = outer_maxit,
      fe_maxit = fe_maxit,
      soft_maxit = soft_maxit,
      tol = tol,
      fe_tol = fe_tol,
      verbose = verbose
    )

    fine_result$cv_coarse_to_fine <- TRUE
    fine_result$cv_coarse_nlambda <- coarse_n
    fine_result$cv_coarse_grid <- coarse_grid
    fine_result$cv_coarse_selected_lambda <- coarse_result$selected_lambda
    fine_result$cv_fine_grid <- fine_grid
    return(fine_result)
  }
  if (sum(mask) <= (nrow(mask) + ncol(mask) + 1L)) {
    stop(
      "Too few eligible zero-dose observations for blocked cross-validation.",
      call. = FALSE
    )
  }

  use_parallel <- cv_workers > 1L
  if (use_parallel && identical(.Platform$OS.type, "windows")) {
    warning(
      "Parallel CV tuning currently uses multicore execution and is not available on Windows. Falling back to sequential execution.",
      call. = FALSE
    )
    use_parallel <- FALSE
    cv_workers <- 1L
  }
  worker_count <- if (use_parallel) min(cv_workers, length(lambda_grid)) else 1L
  use_warm_starts <- cv_warm_starts && !use_parallel
  use_screening <- !is.null(cv_top_k) && cv_top_k < length(lambda_grid)
  active_indices <- seq_along(lambda_grid)

  scores <- matrix(
    NA_real_,
    nrow = length(lambda_grid),
    ncol = cv_rounds,
    dimnames = list(paste0("lambda_", seq_along(lambda_grid)), paste0("round", seq_len(cv_rounds)))
  )
  holdout_counts <- integer(cv_rounds)
  block_sizes_used <- integer(cv_rounds)

  for (round_index in seq_len(cv_rounds)) {
    holdout_mask <- cdmc_select_holdout_mask(mask = mask, block_size = cv_block_size)
    if (!any(holdout_mask)) {
      stop(
        "Unable to construct a non-empty blocked holdout set. Reduce cv_block_size or use heuristic lambda selection.",
        call. = FALSE
      )
    }

    train_mask <- mask & !holdout_mask
    cdmc_validate_support(train_mask)

    holdout_counts[round_index] <- sum(holdout_mask)
    block_sizes_used[round_index] <- attr(holdout_mask, "block_size_used")

    if (verbose) {
      message(sprintf(
        "cv round %d/%d: held out %d control cells",
        round_index,
        cv_rounds,
        holdout_counts[round_index]
      ))
    }

    evaluate_lambda <- function(lambda_index) {
      if (verbose && !use_parallel) {
        message(sprintf(
          "  evaluating lambda %d/%d = %.6g",
          lambda_index,
          length(lambda_grid),
          lambda_grid[lambda_index]
        ))
      }

      baseline_fit <- cdmc_fit_baseline(
        y_matrix = y_matrix,
        x_matrices = x_matrices,
        mask = train_mask,
        weight_matrix = weight_matrix,
        lambda = lambda_grid[lambda_index],
        rank_max = rank_max,
        outer_maxit = outer_maxit,
        fe_maxit = fe_maxit,
        soft_maxit = soft_maxit,
        tol = tol,
        fe_tol = fe_tol,
        verbose = FALSE
      )

      holdout_errors <- (y_matrix[holdout_mask] - baseline_fit$baseline_hat[holdout_mask]) ^ 2
      if (is.null(weight_matrix)) {
        mean(holdout_errors)
      } else {
        stats::weighted.mean(holdout_errors, w = weight_matrix[holdout_mask])
      }
    }

    lambda_indices <- if (use_screening) active_indices else seq_along(lambda_grid)

    lambda_scores <- if (use_parallel) {
      unlist(
        parallel::mclapply(
          lambda_indices,
          evaluate_lambda,
          mc.cores = worker_count,
          mc.set.seed = TRUE
        ),
        use.names = FALSE
      )
    } else {
      lambda_scores_seq <- numeric(length(lambda_indices))
      lambda_fit_start <- NULL
      for (position in seq_along(lambda_indices)) {
        lambda_index <- lambda_indices[[position]]
        if (verbose) {
          message(sprintf(
            "  evaluating lambda %d/%d = %.6g",
            lambda_index,
            length(lambda_grid),
            lambda_grid[lambda_index]
          ))
        }

        baseline_fit <- cdmc_fit_baseline(
          y_matrix = y_matrix,
          x_matrices = x_matrices,
          mask = train_mask,
          weight_matrix = weight_matrix,
          lambda = lambda_grid[lambda_index],
          rank_max = rank_max,
          fit_start = if (use_warm_starts) lambda_fit_start else NULL,
          outer_maxit = outer_maxit,
          fe_maxit = fe_maxit,
          soft_maxit = soft_maxit,
          tol = tol,
          fe_tol = fe_tol,
          verbose = FALSE
        )

        holdout_errors <- (y_matrix[holdout_mask] - baseline_fit$baseline_hat[holdout_mask]) ^ 2
        lambda_scores_seq[position] <- if (is.null(weight_matrix)) {
          mean(holdout_errors)
        } else {
          stats::weighted.mean(holdout_errors, w = weight_matrix[holdout_mask])
        }
        if (use_warm_starts) {
          lambda_fit_start <- baseline_fit$fit_start
        }
      }

      lambda_scores_seq
    }

    scores[lambda_indices, round_index] <- lambda_scores

    if (use_screening && round_index < cv_rounds && length(active_indices) > cv_top_k) {
      partial_means <- rowMeans(scores[active_indices, seq_len(round_index), drop = FALSE], na.rm = TRUE)
      keep_candidates <- active_indices[order(partial_means, decreasing = FALSE)[seq_len(cv_top_k)]]
      active_indices <- sort(keep_candidates)

      if (verbose) {
        message(sprintf(
          "  retaining top %d lambda candidates for remaining rounds",
          length(active_indices)
        ))
      }
    }
  }

  rounds_evaluated <- rowSums(!is.na(scores))
  mean_scores <- rowMeans(scores, na.rm = TRUE)
  if (use_screening) {
    mean_scores[rounds_evaluated < cv_rounds] <- Inf
  }
  best_index <- which.min(mean_scores)

  list(
    method = "cv",
    lambda_grid = lambda_grid,
    scores = scores,
    mean_scores = mean_scores,
    selected_lambda = lambda_grid[best_index],
    selected_index = best_index,
    cv_rounds = cv_rounds,
    cv_block_size_requested = cv_block_size,
    cv_workers = worker_count,
    cv_parallel = use_parallel,
    cv_top_k = cv_top_k,
    cv_screening = use_screening,
    cv_rounds_evaluated = rounds_evaluated,
    cv_coarse_to_fine = FALSE,
    cv_coarse_nlambda = NULL,
    cv_coarse_grid = NULL,
    cv_coarse_selected_lambda = NULL,
    cv_fine_grid = lambda_grid,
    cv_warm_starts = use_warm_starts,
    cv_block_sizes_used = block_sizes_used,
    holdout_counts = holdout_counts
  )
}
