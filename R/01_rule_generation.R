
make_cutoff_candidates <- function(x, parent,
                                   n_total,
                                   min_support_prop = 0.05,
                                   min_span_prop = 0.20) {
  active <- which(parent == 1)
  n_active <- length(active)

  n_min <- ceiling(n_total * min_support_prop)
  span <- ceiling(n_total * min_span_prop)

  if (n_active < 2L * n_min) return(numeric(0))

  x_active <- sort(x[active])
  if (anyNA(x_active)) {
    stop("Missing predictor values are not supported in cutoff generation.")
  }

  runs <- rle(x_active)
  values <- runs$values
  freq <- runs$lengths

  if (length(values) <= 2L) return(numeric(0))

  left <- cumsum(freq)
  right <- n_active - left
  eligible <- which(left >= n_min & right >= n_min)
  if (length(eligible) == 0L) return(numeric(0))

  lower_rank <- n_min
  upper_rank <- n_active - n_min
  rank_width <- upper_rank - lower_rank
  n_intervals <- floor(rank_width / span)

  if (n_intervals == 0L) {
    targets <- (lower_rank + upper_rank) / 2
  } else {
    covered_width <- n_intervals * span
    grid_start <- lower_rank + (rank_width - covered_width) / 2
    targets <- grid_start + seq.int(0L, n_intervals) * span
  }

  idx <- vapply(
    targets,
    FUN.VALUE = integer(1),
    FUN = function(target) {
      candidates <- eligible[left[eligible] >= target]
      if (length(candidates) == 0L) return(NA_integer_)
      candidates[1L]
    }
  )
  idx <- unique(idx[!is.na(idx)])
  values[idx]
}

new_condition <- function(var, direction, cutoff) {
  list(
    var = as.integer(var),
    direction = as.character(direction),
    cutoff = as.numeric(cutoff)
  )
}

canonical_conditions <- function(conditions) {
  if (length(conditions) == 0L) return(conditions)
  ord <- order(
    vapply(conditions, `[[`, integer(1), "var"),
    vapply(conditions, `[[`, character(1), "direction"),
    vapply(conditions, `[[`, numeric(1), "cutoff")
  )
  conditions[ord]
}

make_rule_key <- function(conditions) {
  conditions <- canonical_conditions(conditions)
  paste(
    vapply(
      conditions,
      FUN.VALUE = character(1),
      FUN = function(z) {
        paste(
          z$var,
          z$direction,
          sprintf("%.17g", z$cutoff),
          sep = ":"
        )
      }
    ),
    collapse = "&"
  )
}

format_rule <- function(conditions, var_names, digits = 3L) {
  conditions <- canonical_conditions(conditions)
  paste(
    vapply(
      conditions,
      FUN.VALUE = character(1),
      FUN = function(z) {
        vn <- var_names[z$var]
        if (z$direction == "eq") {
          paste0(vn, " = ", formatC(z$cutoff, digits = digits, format = "fg"))
        } else {
          paste0(
            vn, " ", z$direction, " ",
            formatC(z$cutoff, digits = digits, format = "fg")
          )
        }
      }
    ),
    collapse = " AND "
  )
}

evaluate_conditions <- function(X, conditions) {
  out <- rep(1, nrow(X))
  for (z in conditions) {
    if (z$direction == "<=") out <- out * as.integer(X[, z$var] <= z$cutoff)
    if (z$direction == ">")  out <- out * as.integer(X[, z$var] > z$cutoff)
    if (z$direction == "eq") out <- out * as.integer(X[, z$var] == z$cutoff)
  }
  as.numeric(out)
}

project_rules <- function(rule_pool, X) {
  X <- as_numeric_matrix(X)
  if (length(rule_pool$conditions) == 0L) {
    return(matrix(numeric(0), nrow = nrow(X), ncol = 0L))
  }

  Z <- vapply(
    rule_pool$conditions,
    FUN.VALUE = numeric(nrow(X)),
    FUN = function(cond) evaluate_conditions(X, cond)
  )
  if (is.vector(Z)) Z <- matrix(Z, ncol = 1L)
  colnames(Z) <- rule_pool$key
  Z
}

fit_loss <- function(Z, outcome_type, y = NULL, time = NULL, status = NULL) {
  Z <- as.matrix(Z)
  if (ncol(Z) == 0L) {
    if (outcome_type == "continuous") {
      return(sum((y - mean(y))^2))
    }
    if (outcome_type == "binary") {
      p <- min(max(mean(y), 1e-8), 1 - 1e-8)
      return(-2 * sum(y * log(p) + (1 - y) * log(1 - p)))
    }
    if (outcome_type == "survival") {
      fit <- survival::coxph(
        survival::Surv(time, status) ~ 1,
        ties = "efron",
        singular.ok = TRUE
      )
      return(-2 * as.numeric(logLik(fit)))
    }
  }

  if (outcome_type == "continuous") {
    fit <- stats::lm.fit(cbind(1, Z), y)
    return(sum(fit$residuals^2))
  }

  if (outcome_type == "binary") {
    fit <- suppressWarnings(stats::glm.fit(
      x = cbind(1, Z), y = y, family = stats::binomial()
    ))
    if (!is.finite(fit$deviance)) return(Inf)
    return(fit$deviance)
  }

  dat <- data.frame(time = time, status = status, Z, check.names = FALSE)
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time, status) ~ .,
      data = dat,
      ties = "efron",
      singular.ok = TRUE,
      x = FALSE,
      y = FALSE,
      model = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(Inf)
  val <- -2 * as.numeric(logLik(fit))
  if (!is.finite(val)) Inf else val
}

gable_generate_basis <- function(
    X,
    outcome_type = c("binary", "continuous", "survival"),
    y = NULL,
    time = NULL,
    status = NULL,
    Mmax = 10L,
    degree = 2L,
    min_support_prop = 0.05,
    min_span_prop = 0.20,
    improvement_tol = 1e-8
) {
  outcome_type <- match.arg(outcome_type)
  X <- as_numeric_matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  var_names <- colnames(X)

  is_binary_x <- vapply(seq_len(p), function(j) {
    length(unique(X[, j])) == 2L
  }, logical(1))

  basis <- matrix(numeric(0), nrow = n, ncol = 0L)
  rules <- list()
  keys <- character(0)
  current_loss <- fit_loss(
    basis, outcome_type,
    y = y, time = time, status = status
  )

  while (length(rules) < Mmax) {
    best <- NULL
    best_loss <- Inf

    parent_ids <- c(0L, seq_along(rules))

    for (parent_id in parent_ids) {
      if (parent_id == 0L) {
        parent_vec <- rep(1, n)
        parent_conditions <- list()
        used_vars <- integer(0)
      } else {
        parent_vec <- basis[, parent_id]
        parent_conditions <- rules[[parent_id]]$conditions
        used_vars <- vapply(
          parent_conditions, `[[`, integer(1), "var"
        )
      }

      parent_degree <- length(parent_conditions)
      if (parent_degree >= degree) next

      candidate_vars <- setdiff(seq_len(p), used_vars)

      for (j in candidate_vars) {
        if (is_binary_x[j]) {
          vals <- sort(unique(X[, j]))
          upper <- vals[2L]
          z1 <- parent_vec * as.integer(X[, j] == upper)
          z0 <- parent_vec * as.integer(X[, j] != upper)

          n_min <- ceiling(n * min_support_prop)
          if (sum(z1) < n_min || sum(z0) < n_min) next

          cond <- c(parent_conditions, list(new_condition(j, "eq", upper)))
          key <- make_rule_key(cond)
          if (key %in% keys) next

          loss <- fit_loss(
            cbind(basis, z1), outcome_type,
            y = y, time = time, status = status
          )
          if (loss < best_loss) {
            best_loss <- loss
            best <- list(
              columns = list(z1),
              conditions = list(cond),
              keys = key
            )
          }
        } else {
          cutoffs <- make_cutoff_candidates(
            x = X[, j],
            parent = parent_vec,
            n_total = n,
            min_support_prop = min_support_prop,
            min_span_prop = min_span_prop
          )
          if (length(cutoffs) == 0L) next

          for (cut in cutoffs) {
            z_le <- parent_vec * as.integer(X[, j] <= cut)
            z_gt <- parent_vec * as.integer(X[, j] > cut)

            cond_le <- c(parent_conditions, list(new_condition(j, "<=", cut)))
            cond_gt <- c(parent_conditions, list(new_condition(j, ">", cut)))
            key_le <- make_rule_key(cond_le)
            key_gt <- make_rule_key(cond_gt)

            if (key_le %in% keys || key_gt %in% keys) next

            loss <- fit_loss(
              cbind(basis, z_le, z_gt), outcome_type,
              y = y, time = time, status = status
            )
            if (loss < best_loss) {
              best_loss <- loss
              best <- list(
                columns = list(z_le, z_gt),
                conditions = list(cond_le, cond_gt),
                keys = c(key_le, key_gt)
              )
            }
          }
        }
      }
    }

    if (is.null(best) ||
        !is.finite(best_loss) ||
        best_loss >= current_loss - improvement_tol) {
      break
    }

    n_add <- length(best$columns)
    if (length(rules) + n_add > Mmax) break

    for (k in seq_len(n_add)) {
      basis <- cbind(basis, best$columns[[k]])
      rules[[length(rules) + 1L]] <- list(
        conditions = canonical_conditions(best$conditions[[k]]),
        key = best$keys[k],
        degree = length(best$conditions[[k]])
      )
      keys <- c(keys, best$keys[k])
    }
    current_loss <- best_loss
  }

  if (ncol(basis) > 0L) colnames(basis) <- keys

  list(
    Basis = basis,
    rules = rules,
    loss = current_loss,
    outcome_type = outcome_type,
    var_names = var_names
  )
}

generate_bootstrap_rule_pool <- function(
    X,
    outcome_type,
    y = NULL,
    time = NULL,
    status = NULL,
    B = 50L,
    Mmax = 10L,
    degree = 2L,
    min_support_prop = 0.05,
    min_span_prop = 0.20,
    improvement_tol = 1e-8,
    seed = NULL
) {
  X <- as_numeric_matrix(X)
  n <- nrow(X)
  var_names <- colnames(X)
  if (!is.null(seed)) set.seed(seed)

  all_rules <- list()
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    fit <- gable_generate_basis(
      X = X[idx, , drop = FALSE],
      outcome_type = outcome_type,
      y = if (!is.null(y)) y[idx] else NULL,
      time = if (!is.null(time)) time[idx] else NULL,
      status = if (!is.null(status)) status[idx] else NULL,
      Mmax = Mmax,
      degree = degree,
      min_support_prop = min_support_prop,
      min_span_prop = min_span_prop,
      improvement_tol = improvement_tol
    )

    if (length(fit$rules) > 0L) {
      for (r in fit$rules) {
        all_rules[[length(all_rules) + 1L]] <- list(
          conditions = r$conditions,
          key = r$key,
          degree = r$degree,
          bootstrap_id = b
        )
      }
    }
  }

  if (length(all_rules) == 0L) {
    return(list(
      key = character(0),
      conditions = list(),
      degree = integer(0),
      count = integer(0),
      count_bootstrap = integer(0),
      rule = character(0),
      var_names = var_names
    ))
  }

  keys <- vapply(all_rules, `[[`, character(1), "key")
  unique_keys <- unique(keys)

  conditions <- lapply(unique_keys, function(k) {
    all_rules[[match(k, keys)]]$conditions
  })
  degree <- vapply(conditions, length, integer(1))
  count <- vapply(unique_keys, function(k) sum(keys == k), integer(1))
  count_bootstrap <- vapply(unique_keys, function(k) {
    ids <- vapply(all_rules[keys == k], `[[`, integer(1), "bootstrap_id")
    length(unique(ids))
  }, integer(1))
  rule <- vapply(
    conditions, format_rule, character(1),
    var_names = var_names, digits = 3L
  )

  ord <- order(-count_bootstrap, -count, unique_keys)
  list(
    key = unique_keys[ord],
    conditions = conditions[ord],
    degree = degree[ord],
    count = count[ord],
    count_bootstrap = count_bootstrap[ord],
    rule = rule[ord],
    var_names = var_names
  )
}
