
clip_brp_linear_matrix <- function(X, lower, upper) {
  X <- as_numeric_matrix(X)
  if (ncol(X) != length(lower) || ncol(X) != length(upper)) {
    stop("The new predictor matrix does not match the fitted linear transform.")
  }

  out <- X
  for (j in seq_len(ncol(out))) {
    out[, j] <- pmin(pmax(out[, j], lower[j]), upper[j])
  }
  out
}

make_brp_feature_matrix <- function(bundle, X_new) {
  required <- bundle$preprocess$predictor_names
  X_new <- as.data.frame(X_new, check.names = FALSE)
  missing <- setdiff(required, names(X_new))
  if (length(missing) > 0L) {
    stop("X_new is missing predictors: ", paste(missing, collapse = ", "))
  }
  X_new <- X_new[, required, drop = FALSE]
  X_new <- as_numeric_matrix(X_new)

  pool <- bundle$preprocess$rule_pool
  Z_raw <- project_rules(pool, X_new)
  if (ncol(Z_raw) > 0L) {
    Z <- sweep(Z_raw, 2, bundle$preprocess$rule_sd, "/")
    colnames(Z) <- paste0("rule__", pool$key)
  } else {
    Z <- matrix(numeric(0), nrow(X_new), 0L)
  }

  L_raw <- clip_brp_linear_matrix(
    X_new,
    lower = bundle$preprocess$lower,
    upper = bundle$preprocess$upper
  )
  L <- sweep(L_raw, 2, bundle$preprocess$linear_sd, "/")
  colnames(L) <- paste0("linear__", required)

  out <- cbind(Z, L)
  expected <- bundle$model$feature_names
  missing_features <- setdiff(expected, colnames(out))
  extra_features <- setdiff(colnames(out), expected)
  if (length(missing_features) > 0L || length(extra_features) > 0L) {
    stop(
      "The reconstructed design matrix does not match the fitted model. ",
      "Missing: ", paste(missing_features, collapse = ", "),
      "; extra: ", paste(extra_features, collapse = ", ")
    )
  }
  out[, expected, drop = FALSE]
}

predict_brp_from_raw_X <- function(bundle, X_new, type = "link") {
  x_new <- make_brp_feature_matrix(bundle, X_new)
  as.numeric(stats::predict(
    bundle$model$cv_fit,
    newx = x_new,
    s = bundle$model$lambda_choice,
    type = type
  ))
}

build_full_term_table <- function(beta, design, predictor_names, B) {
  pool <- design$rule_pool
  feature_names <- c(
    paste0("rule__", pool$key),
    paste0("linear__", predictor_names)
  )

  beta <- as.numeric(beta[feature_names])
  if (length(beta) != length(feature_names) || anyNA(beta)) {
    stop("Unable to align lasso coefficients with the BRP-GABLE design.")
  }

  n_rule <- length(pool$key)
  is_rule <- seq_along(feature_names) <= n_rule
  rule_idx <- ifelse(is_rule, seq_along(feature_names), NA_integer_)
  linear_idx <- ifelse(!is_rule, seq_along(feature_names) - n_rule, NA_integer_)

  scale_sd <- c(design$rule_sd, design$linear_sd)
  coef_original <- beta / scale_sd
  selected <- beta != 0

  term <- character(length(feature_names))
  term[is_rule] <- gsub(
    "[[:space:]]+", " ", trimws(pool$rule)
  )
  term[!is_rule] <- paste0(predictor_names[linear_idx[!is_rule]], " (linear)")

  rule_key <- rep(NA_character_, length(feature_names))
  rule_key[is_rule] <- pool$key
  degree <- rep(1L, length(feature_names))
  degree[is_rule] <- pool$degree
  support <- rep(NA_real_, length(feature_names))
  support[is_rule] <- design$rule_support
  count <- rep(NA_integer_, length(feature_names))
  count[is_rule] <- pool$count
  count_bootstrap <- rep(NA_integer_, length(feature_names))
  count_bootstrap[is_rule] <- pool$count_bootstrap

  rule_variables <- rep(NA_character_, length(feature_names))
  if (n_rule > 0L) {
    rule_variables[is_rule] <- vapply(
      pool$conditions,
      FUN.VALUE = character(1),
      FUN = function(cond) {
        ids <- unique(vapply(cond, `[[`, integer(1), "var"))
        paste(predictor_names[ids], collapse = "; ")
      }
    )
  }
  rule_variables[!is_rule] <- predictor_names[linear_idx[!is_rule]]

  out <- data.frame(
    feature_name = feature_names,
    term = term,
    term_type = ifelse(is_rule, "rule", "linear"),
    rule_key = rule_key,
    variables = rule_variables,
    degree = degree,
    support = support,
    count = count,
    count_bootstrap = count_bootstrap,
    generation_frequency = count_bootstrap / B,
    coefficient_standardized = beta,
    scale_sd = scale_sd,
    coefficient_original = coef_original,
    hazard_ratio = exp(coef_original),
    importance = abs(beta),
    selected = selected,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  selected_max <- max(out$importance[out$selected], na.rm = TRUE)
  if (!is.finite(selected_max) || selected_max <= 0) selected_max <- 1
  out$importance_relative <- 100 * out$importance / selected_max

  selected_order <- which(out$selected)[
    order(-out$importance[out$selected], out$term[out$selected])
  ]
  out$term_id <- NA_character_
  out$term_id[selected_order] <- paste0("T", seq_along(selected_order))

  selected_rule_order <- selected_order[out$term_type[selected_order] == "rule"]
  out$rule_id <- NA_character_
  out$rule_id[selected_rule_order] <- paste0(
    "R", seq_along(selected_rule_order)
  )
  out
}

make_variable_contribution_table <- function(term_table, rule_pool,
                                             predictor_names) {
  selected <- term_table[term_table$selected, , drop = FALSE]
  contribution <- setNames(rep(0, length(predictor_names)), predictor_names)

  for (i in seq_len(nrow(selected))) {
    if (selected$term_type[i] == "linear") {
      vars <- selected$variables[i]
    } else {
      idx <- match(selected$rule_key[i], rule_pool$key)
      if (is.na(idx)) next
      var_ids <- unique(vapply(
        rule_pool$conditions[[idx]], `[[`, integer(1), "var"
      ))
      vars <- predictor_names[var_ids]
    }
    vars <- unique(vars[vars %in% predictor_names])
    if (!length(vars)) next
    contribution[vars] <- contribution[vars] + selected$importance[i] / length(vars)
  }

  out <- data.frame(
    variable = names(contribution),
    contribution = as.numeric(contribution),
    stringsAsFactors = FALSE
  )
  out[order(-out$contribution, out$variable), , drop = FALSE]
}

compute_ale_curve <- function(bundle, variable, n_bins = 20L) {
  X <- bundle$data$X
  if (!variable %in% names(X)) stop("Unknown ALE variable: ", variable)
  x <- as.numeric(X[[variable]])
  if (length(unique(x)) <= 2L) {
    stop(variable, " is binary; use compute_binary_effect() instead.")
  }

  n_bins <- max(2L, as.integer(n_bins))
  breaks <- unique(as.numeric(stats::quantile(
    x,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )))
  if (length(breaks) < 3L) {
    stop("Too few unique quantile boundaries for ALE variable ", variable, ".")
  }

  K <- length(breaks) - 1L
  bin <- findInterval(
    x, breaks,
    all.inside = TRUE,
    rightmost.closed = TRUE
  )
  bin <- pmin(pmax(bin, 1L), K)

  local_effect <- numeric(K)
  n_in_bin <- integer(K)
  for (k in seq_len(K)) {
    ids <- which(bin == k)
    n_in_bin[k] <- length(ids)
    if (!length(ids)) next

    X_low <- X[ids, , drop = FALSE]
    X_high <- X[ids, , drop = FALSE]
    X_low[[variable]] <- breaks[k]
    X_high[[variable]] <- breaks[k + 1L]
    local_effect[k] <- mean(
      predict_brp_from_raw_X(bundle, X_high, type = "link") -
        predict_brp_from_raw_X(bundle, X_low, type = "link")
    )
  }

  accumulated <- cumsum(local_effect)
  center <- stats::weighted.mean(accumulated, n_in_bin)
  data.frame(
    variable = variable,
    x = breaks,
    ale = c(-center, accumulated - center),
    bin_lower = c(NA_real_, breaks[-length(breaks)]),
    bin_upper = c(NA_real_, breaks[-1L]),
    n_in_bin = c(NA_integer_, n_in_bin),
    stringsAsFactors = FALSE
  )
}

compute_binary_effect <- function(bundle, variable) {
  X <- bundle$data$X
  if (!variable %in% names(X)) stop("Unknown binary variable: ", variable)
  values <- sort(unique(X[[variable]]))
  if (length(values) != 2L) stop(variable, " does not have exactly two values.")

  predictions <- lapply(values, function(value) {
    X_new <- X
    X_new[[variable]] <- value
    predict_brp_from_raw_X(bundle, X_new, type = "link")
  })
  means <- vapply(predictions, mean, numeric(1))
  observed_weight <- vapply(values, function(v) mean(X[[variable]] == v), numeric(1))
  center <- sum(means * observed_weight)
  contrast <- means[2L] - means[1L]

  data.frame(
    variable = variable,
    level = values,
    centered_effect = means - center,
    contrast_high_vs_low = contrast,
    relative_hazard_high_vs_low = exp(contrast),
    observed_proportion = observed_weight,
    stringsAsFactors = FALSE
  )
}

compute_jaccard_similarity <- function(activation) {
  activation <- as.matrix(activation)
  if (ncol(activation) == 0L) {
    return(matrix(numeric(0), 0L, 0L))
  }
  storage.mode(activation) <- "double"
  intersection <- crossprod(activation)
  support_count <- colSums(activation)
  union <- outer(support_count, support_count, "+") - intersection
  out <- intersection / union
  out[union == 0] <- NA_real_
  diag(out) <- 1
  dimnames(out) <- list(colnames(activation), colnames(activation))
  out
}

select_representative_rules <- function(selected_rules, n_rules = 3L,
                                        min_support = 0.10,
                                        max_support = 0.90) {
  if (!nrow(selected_rules)) return(selected_rules)
  n_rules <- min(as.integer(n_rules), nrow(selected_rules))
  eligible <- selected_rules[
    is.finite(selected_rules$support) &
      selected_rules$support >= min_support &
      selected_rules$support <= max_support,
    , drop = FALSE
  ]
  if (nrow(eligible) < n_rules) eligible <- selected_rules

  chosen <- integer(0)
  high <- which(eligible$hazard_ratio > 1)
  if (length(high)) chosen <- c(chosen, high[which.max(eligible$importance[high])])
  low <- which(eligible$hazard_ratio < 1)
  if (length(low)) chosen <- c(chosen, low[which.max(eligible$importance[low])])

  remaining <- order(-eligible$importance, eligible$rule_id)
  chosen <- unique(c(chosen, remaining))
  eligible[chosen[seq_len(min(n_rules, length(chosen)))], , drop = FALSE]
}

build_km_summaries <- function(bundle, representative_rules,
                               n_risk_times = 5L) {
  if (!nrow(representative_rules)) {
    return(list(curves = data.frame(), risk_table = data.frame()))
  }

  curve_rows <- list()
  risk_rows <- list()
  for (i in seq_len(nrow(representative_rules))) {
    rid <- representative_rules$rule_id[i]
    z <- bundle$selected_rule_activation[, rid]
    group <- factor(
      ifelse(z == 1, "Satisfied", "Not satisfied"),
      levels = c("Not satisfied", "Satisfied")
    )
    fit <- survival::survfit(
      survival::Surv(bundle$data$time, bundle$data$status) ~ group
    )
    sm <- summary(fit)
    strata <- sub("^group=", "", as.character(sm$strata))
    curves <- data.frame(
      rule_id = rid,
      rule = representative_rules$term[i],
      group = strata,
      time = sm$time,
      survival = sm$surv,
      lower = sm$lower,
      upper = sm$upper,
      n_risk = sm$n.risk,
      n_event = sm$n.event,
      stringsAsFactors = FALSE
    )
    starts <- data.frame(
      rule_id = rid,
      rule = representative_rules$term[i],
      group = levels(group),
      time = 0,
      survival = 1,
      lower = 1,
      upper = 1,
      n_risk = as.integer(table(group)[levels(group)]),
      n_event = 0,
      stringsAsFactors = FALSE
    )
    curve_rows[[i]] <- rbind(starts, curves)

    risk_times <- pretty(
      c(0, max(bundle$data$time, na.rm = TRUE)),
      n = max(2L, as.integer(n_risk_times))
    )
    risk_times <- risk_times[
      risk_times >= 0 & risk_times <= max(bundle$data$time, na.rm = TRUE)
    ]
    risk_sm <- summary(fit, times = risk_times, extend = TRUE)
    risk_rows[[i]] <- data.frame(
      rule_id = rid,
      rule = representative_rules$term[i],
      group = sub("^group=", "", as.character(risk_sm$strata)),
      time = risk_sm$time,
      n_risk = risk_sm$n.risk,
      stringsAsFactors = FALSE
    )
  }
  list(
    curves = do.call(rbind, curve_rows),
    risk_table = do.call(rbind, risk_rows)
  )
}

fit_full_brp_gable_survival <- function(
    real_data = prepare_mgus2_data(),
    config = default_config(),
    output_dir = file.path("results", "real_data", "experiment2"),
    dataset_name = real_data$dataset_name %||% "real_survival",
    full_fit_id = 0L,
    n_ale_bins = 20L,
    n_ale_variables = 3L,
    n_km_rules = 3L
) {
  assert_packages(c("survival", "glmnet"))
  real_data <- validate_real_survival_data(real_data)
  X <- real_data$X
  predictor_names <- names(X)
  seeds <- make_seed_set(config$base_seed, as.integer(full_fit_id))

  pool_start <- proc.time()[["elapsed"]]
  pool_raw <- generate_bootstrap_rule_pool(
    X = X,
    outcome_type = "survival",
    time = real_data$time,
    status = real_data$status,
    B = config$B,
    Mmax = config$Mmax,
    degree = config$degree,
    min_support_prop = config$min_support_prop,
    min_span_prop = config$min_span_prop,
    improvement_tol = config$improvement_tol,
    seed = seeds$brp
  )
  pool_elapsed <- proc.time()[["elapsed"]] - pool_start

  design <- build_brp_design(
    pool_raw, X, X,
    linear_trim = config$linear_trim
  )
  x_full <- cbind(design$rules_train, design$linear_train)
  foldid <- make_foldid(
    n = nrow(X),
    nfolds = config$nfolds,
    seed = seeds$cv,
    strata = NULL
  )

  fit_start <- proc.time()[["elapsed"]]
  cv_fit <- glmnet::cv.glmnet(
    x = x_full,
    y = survival::Surv(real_data$time, real_data$status),
    family = "cox",
    alpha = 1,
    foldid = foldid,
    standardize = FALSE
  )
  fit_elapsed <- proc.time()[["elapsed"]] - fit_start

  coef_matrix <- as.matrix(stats::coef(cv_fit, s = config$lambda_choice))
  beta <- coef_matrix[, 1L]
  names(beta) <- rownames(coef_matrix)

  term_table <- build_full_term_table(
    beta = beta,
    design = design,
    predictor_names = predictor_names,
    B = config$B
  )
  selected_terms <- term_table[term_table$selected, , drop = FALSE]
  selected_terms <- selected_terms[
    order(-selected_terms$importance, selected_terms$term),
    , drop = FALSE
  ]
  selected_rules <- selected_terms[
    selected_terms$term_type == "rule",
    , drop = FALSE
  ]

  activation_all <- project_rules(design$rule_pool, X)
  selected_rule_activation <- activation_all[
    , match(selected_rules$rule_key, design$rule_pool$key), drop = FALSE
  ]
  colnames(selected_rule_activation) <- selected_rules$rule_id

  risk_score <- as.numeric(stats::predict(
    cv_fit,
    newx = x_full,
    s = config$lambda_choice,
    type = "link"
  ))

  bundle <- list(
    dataset_name = dataset_name,
    data = list(
      X = X,
      time = real_data$time,
      status = real_data$status,
      patient_id = seq_len(nrow(X))
    ),
    config = config,
    seeds = list(brp = seeds$brp, cv = seeds$cv),
    runtime = list(
      pool_seconds = pool_elapsed,
      lasso_seconds = fit_elapsed,
      total_seconds = pool_elapsed + fit_elapsed
    ),
    pool_raw = pool_raw,
    preprocess = list(
      predictor_names = predictor_names,
      rule_pool = design$rule_pool,
      rule_support = design$rule_support,
      rule_sd = design$rule_sd,
      linear_sd = design$linear_sd,
      lower = design$lower,
      upper = design$upper
    ),
    model = list(
      cv_fit = cv_fit,
      lambda_choice = config$lambda_choice,
      lambda_value = unname(cv_fit[[config$lambda_choice]]),
      feature_names = colnames(x_full),
      foldid = foldid,
      risk_score = risk_score
    ),
    term_table = term_table,
    selected_terms = selected_terms,
    selected_rules = selected_rules,
    selected_rule_activation = selected_rule_activation
  )

  bundle$variable_contribution <- make_variable_contribution_table(
    term_table,
    design$rule_pool,
    predictor_names
  )
  bundle$jaccard_similarity <- compute_jaccard_similarity(
    selected_rule_activation
  )

  continuous <- predictor_names[
    vapply(X, function(z) length(unique(z)) > 2L, logical(1))
  ]
  ranked_continuous <- bundle$variable_contribution$variable[
    bundle$variable_contribution$variable %in% continuous &
      bundle$variable_contribution$contribution > 0
  ]
  ranked_continuous <- head(ranked_continuous, as.integer(n_ale_variables))
  ale_curves <- lapply(
    ranked_continuous,
    function(v) compute_ale_curve(bundle, v, n_bins = n_ale_bins)
  )
  bundle$ale_curves <- if (length(ale_curves)) {
    do.call(rbind, ale_curves)
  } else {
    data.frame()
  }

  binary <- predictor_names[
    vapply(X, function(z) length(unique(z)) == 2L, logical(1))
  ]
  ranked_binary <- bundle$variable_contribution$variable[
    bundle$variable_contribution$variable %in% binary &
      bundle$variable_contribution$contribution > 0
  ]
  binary_effects <- lapply(
    ranked_binary,
    function(v) compute_binary_effect(bundle, v)
  )
  bundle$binary_effects <- if (length(binary_effects)) {
    do.call(rbind, binary_effects)
  } else {
    data.frame()
  }

  bundle$representative_rules <- select_representative_rules(
    selected_rules,
    n_rules = n_km_rules
  )
  km <- build_km_summaries(bundle, bundle$representative_rules)
  bundle$km_curves <- km$curves
  bundle$km_risk_table <- km$risk_table

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dataset_tag <- gsub("[^A-Za-z0-9_-]+", "_", tolower(dataset_name))
  tag <- paste0(
    dataset_tag,
    "_full_n", nrow(X),
    "_B", config$B,
    "_M", config$Mmax,
    "_deg", config$degree,
    "_support", gsub("\\.", "_", config$min_support_prop),
    "_span", gsub("\\.", "_", config$min_span_prop)
  )

  rds_path <- file.path(output_dir, paste0("brp_gable_", tag, ".rds"))
  saveRDS(bundle, rds_path)
  write.csv(
    selected_terms,
    file.path(output_dir, paste0("selected_terms_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    term_table[term_table$term_type == "rule", , drop = FALSE],
    file.path(output_dir, paste0("candidate_rules_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    bundle$variable_contribution,
    file.path(output_dir, paste0("variable_contribution_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  if (nrow(bundle$ale_curves)) {
    write.csv(
      bundle$ale_curves,
      file.path(output_dir, paste0("ale_curves_", tag, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  if (nrow(bundle$binary_effects)) {
    write.csv(
      bundle$binary_effects,
      file.path(output_dir, paste0("binary_effects_", tag, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  if (length(bundle$jaccard_similarity)) {
    write.csv(
      bundle$jaccard_similarity,
      file.path(output_dir, paste0("jaccard_similarity_", tag, ".csv")),
      row.names = TRUE,
      fileEncoding = "UTF-8"
    )
  }
  if (nrow(bundle$km_risk_table)) {
    write.csv(
      bundle$km_risk_table,
      file.path(output_dir, paste0("km_risk_table_", tag, ".csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }

  attr(bundle, "output_rds") <- rds_path
  invisible(bundle)
}
