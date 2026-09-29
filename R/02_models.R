
build_brp_design <- function(rule_pool, X_train, X_test,
                             linear_trim = 0.025) {
  X_train <- as_numeric_matrix(X_train)
  X_test <- as_numeric_matrix(X_test)

  Z_train_raw <- project_rules(rule_pool, X_train)
  Z_test_raw <- project_rules(rule_pool, X_test)

  if (ncol(Z_train_raw) > 0L) {
    support <- colMeans(Z_train_raw)
    informative <- support > 0 & support < 1
    Z_train_raw <- Z_train_raw[, informative, drop = FALSE]
    Z_test_raw <- Z_test_raw[, informative, drop = FALSE]

    pool2 <- rule_pool
    pool2$key <- rule_pool$key[informative]
    pool2$conditions <- rule_pool$conditions[informative]
    pool2$degree <- rule_pool$degree[informative]
    pool2$count <- rule_pool$count[informative]
    pool2$count_bootstrap <- rule_pool$count_bootstrap[informative]
    pool2$rule <- rule_pool$rule[informative]

    support <- support[informative]
    rule_sd <- sqrt(support * (1 - support))
    Z_train <- sweep(Z_train_raw, 2, rule_sd, "/")
    Z_test <- sweep(Z_test_raw, 2, rule_sd, "/")
    colnames(Z_train) <- paste0("rule__", pool2$key)
    colnames(Z_test) <- colnames(Z_train)
  } else {
    pool2 <- rule_pool
    support <- numeric(0)
    rule_sd <- numeric(0)
    Z_train <- matrix(numeric(0), nrow(X_train), 0L)
    Z_test <- matrix(numeric(0), nrow(X_test), 0L)
  }

  lower <- apply(X_train, 2, stats::quantile, probs = linear_trim,
                 na.rm = TRUE, names = FALSE)
  upper <- apply(X_train, 2, stats::quantile, probs = 1 - linear_trim,
                 na.rm = TRUE, names = FALSE)

  clip_matrix <- function(X) {
    out <- X
    for (j in seq_len(ncol(out))) {
      out[, j] <- pmin(pmax(out[, j], lower[j]), upper[j])
    }
    out
  }

  L_train_raw <- clip_matrix(X_train)
  L_test_raw <- clip_matrix(X_test)
  linear_sd <- apply(L_train_raw, 2, stats::sd)
  linear_sd[!is.finite(linear_sd) | linear_sd == 0] <- 1

  L_train <- sweep(L_train_raw, 2, linear_sd, "/")
  L_test <- sweep(L_test_raw, 2, linear_sd, "/")
  colnames(L_train) <- paste0("linear__", colnames(X_train))
  colnames(L_test) <- colnames(L_train)

  list(
    rules_train = Z_train,
    rules_test = Z_test,
    linear_train = L_train,
    linear_test = L_test,
    rule_support = support,
    rule_sd = rule_sd,
    linear_sd = linear_sd,
    lower = lower,
    upper = upper,
    rule_pool = pool2
  )
}

fit_brp_lasso <- function(design, outcome_type,
                          y_train = NULL,
                          time_train = NULL,
                          status_train = NULL,
                          y_test = NULL,
                          time_test = NULL,
                          status_test = NULL,
                          include_linear = TRUE,
                          foldid,
                          lambda_choice = "lambda.1se") {
  x_train <- if (include_linear) {
    cbind(design$rules_train, design$linear_train)
  } else {
    design$rules_train
  }
  x_test <- if (include_linear) {
    cbind(design$rules_test, design$linear_test)
  } else {
    design$rules_test
  }

  if (ncol(x_train) == 0L) stop("No usable columns for BRP-GABLE lasso.")

  family <- switch(
    outcome_type,
    binary = "binomial",
    continuous = "gaussian",
    survival = "cox"
  )
  yy <- if (outcome_type == "survival") {
    survival::Surv(time_train, status_train)
  } else {
    y_train
  }

  fit <- glmnet::cv.glmnet(
    x = x_train,
    y = yy,
    family = family,
    alpha = 1,
    foldid = foldid,
    standardize = FALSE
  )

  pred_type <- if (outcome_type == "binary") "response" else "link"
  if (outcome_type == "continuous") pred_type <- "response"

  pred <- as.numeric(predict(
    fit, newx = x_test, s = lambda_choice, type = pred_type
  ))

  metric <- switch(
    outcome_type,
    binary = auc_value(y_test, pred),
    continuous = mse_value(y_test, pred),
    survival = cindex_value(time_test, status_test, pred)
  )

  beta <- as.matrix(coef(fit, s = lambda_choice))
  selected_names <- rownames(beta)[as.numeric(beta) != 0]
  selected_names <- setdiff(selected_names, "(Intercept)")

  selected_rule_keys <- sub(
    "^rule__", "",
    selected_names[grepl("^rule__", selected_names)]
  )
  selected_rule_idx <- match(
    selected_rule_keys, design$rule_pool$key,
    nomatch = 0L
  )
  degrees <- design$rule_pool$degree[selected_rule_idx[selected_rule_idx > 0L]]
  n_linear <- sum(grepl("^linear__", selected_names))

  list(
    metric_value = metric,
    complexity = summarize_rule_complexity(degrees, n_linear),
    fit = fit,
    prediction = pred,
    selected_names = selected_names
  )
}

fit_gable_baseline <- function(X_train, outcome_type,
                               y_train = NULL,
                               time_train = NULL,
                               status_train = NULL,
                               X_test,
                               y_test = NULL,
                               time_test = NULL,
                               status_test = NULL,
                               Mmax = 10L,
                               degree = 2L,
                               min_support_prop = 0.05,
                               min_span_prop = 0.20,
                               improvement_tol = 1e-8,
                               step_k = 6) {
  fit0 <- gable_generate_basis(
    X_train,
    outcome_type = outcome_type,
    y = y_train,
    time = time_train,
    status = status_train,
    Mmax = Mmax,
    degree = degree,
    min_support_prop = min_support_prop,
    min_span_prop = min_span_prop,
    improvement_tol = improvement_tol
  )

  Z_train <- fit0$Basis
  if (ncol(Z_train) == 0L) {
    stop("GABLE generated no rules.")
  }
  Z_test <- vapply(
    fit0$rules,
    FUN.VALUE = numeric(nrow(X_test)),
    FUN = function(r) evaluate_conditions(as_numeric_matrix(X_test), r$conditions)
  )
  if (is.vector(Z_test)) Z_test <- matrix(Z_test, ncol = 1L)
  colnames(Z_train) <- paste0("R", seq_len(ncol(Z_train)))
  colnames(Z_test) <- colnames(Z_train)

  if (outcome_type == "binary") {
    dat <- data.frame(y = y_train, Z_train, check.names = FALSE)
    full <- stats::glm(y ~ ., data = dat, family = stats::binomial())
    model <- tryCatch(
      MASS::stepAIC(full, direction = "both", k = step_k, trace = FALSE),
      error = function(e) full
    )
    pred <- as.numeric(predict(model, newdata = as.data.frame(Z_test),
                               type = "response"))
    metric <- auc_value(y_test, pred)
  } else if (outcome_type == "continuous") {
    dat <- data.frame(y = y_train, Z_train, check.names = FALSE)
    full <- stats::lm(y ~ ., data = dat)
    model <- tryCatch(
      MASS::stepAIC(full, direction = "both", k = step_k, trace = FALSE),
      error = function(e) full
    )
    pred <- as.numeric(predict(model, newdata = as.data.frame(Z_test)))
    metric <- mse_value(y_test, pred)
  } else {
    dat <- data.frame(time = time_train, status = status_train,
                      Z_train, check.names = FALSE)
    full <- survival::coxph(
      survival::Surv(time, status) ~ ., data = dat,
      singular.ok = TRUE
    )
    model <- tryCatch(
      MASS::stepAIC(full, direction = "both", k = step_k, trace = FALSE),
      error = function(e) full
    )
    pred <- as.numeric(predict(model, newdata = as.data.frame(Z_test),
                               type = "lp"))
    metric <- cindex_value(time_test, status_test, pred)
  }

  coef_names <- setdiff(names(stats::coef(model)), "(Intercept)")
  idx <- suppressWarnings(as.integer(sub("^R", "", coef_names)))
  idx <- idx[is.finite(idx) & idx >= 1L & idx <= length(fit0$rules)]
  degrees <- vapply(fit0$rules[idx], function(r) r$degree, integer(1))

  list(
    metric_value = metric,
    complexity = summarize_rule_complexity(degrees, 0L),
    fit = model,
    prediction = pred
  )
}

fit_linear_baseline <- function(X_train, outcome_type,
                                y_train = NULL,
                                time_train = NULL,
                                status_train = NULL,
                                X_test,
                                y_test = NULL,
                                time_test = NULL,
                                status_test = NULL) {
  X_train <- as.data.frame(X_train)
  X_test <- as.data.frame(X_test)

  if (outcome_type == "binary") {
    dat <- data.frame(y = y_train, X_train)
    fit <- stats::glm(y ~ ., data = dat, family = stats::binomial())
    pred <- as.numeric(predict(fit, newdata = X_test, type = "response"))
    metric <- auc_value(y_test, pred)
  } else if (outcome_type == "continuous") {
    dat <- data.frame(y = y_train, X_train)
    fit <- stats::lm(y ~ ., data = dat)
    pred <- as.numeric(predict(fit, newdata = X_test))
    metric <- mse_value(y_test, pred)
  } else {
    dat <- data.frame(time = time_train, status = status_train, X_train)
    fit <- survival::coxph(
      survival::Surv(time, status) ~ ., data = dat,
      singular.ok = TRUE
    )
    pred <- as.numeric(predict(fit, newdata = X_test, type = "lp"))
    metric <- cindex_value(time_test, status_test, pred)
  }

  list(
    metric_value = metric,
    complexity = summarize_rule_complexity(
      degrees = integer(0),
      n_linear_selected = sum(
        is.finite(stats::coef(fit)) &
          stats::coef(fit) != 0 &
          names(stats::coef(fit)) != "(Intercept)"
      )
    ),
    fit = fit,
    prediction = pred
  )
}

fit_rulefit_default <- function(X_train, outcome_type,
                                 y_train = NULL,
                                 time_train = NULL,
                                 status_train = NULL,
                                 X_test,
                                 y_test = NULL,
                                 time_test = NULL,
                                 status_test = NULL,
                                 seed = NULL,
                                 removecomplements = NULL,
                                 removeduplicates = NULL,
                                 maxdepth = NULL) {
  if (!is.null(seed)) set.seed(seed)

  fit_pre <- function(formula, data, family) {
    pre_args <- list(
      formula = formula,
      data = data,
      family = family
    )
    if (!is.null(removecomplements)) {
      pre_args$removecomplements <- removecomplements
    }
    if (!is.null(removeduplicates)) {
      pre_args$removeduplicates <- removeduplicates
    }
    if (!is.null(maxdepth)) {
      pre_args$maxdepth <- as.integer(maxdepth)
    }
    do.call(pre::pre, pre_args)
  }

  if (outcome_type == "binary") {
    dat <- data.frame(y = as.factor(y_train), X_train)
    fit <- fit_pre(y ~ ., dat, "binomial")
    pred <- as.numeric(predict(
      fit, newdata = X_test, type = "response",
      penalty.par.val = "lambda.1se"
    ))
    metric <- auc_value(y_test, pred)
  } else if (outcome_type == "continuous") {
    dat <- data.frame(y = y_train, X_train)
    fit <- fit_pre(y ~ ., dat, "gaussian")
    pred <- as.numeric(predict(
      fit, newdata = X_test, type = "response",
      penalty.par.val = "lambda.1se"
    ))
    metric <- mse_value(y_test, pred)
  } else {
    dat <- data.frame(time = time_train, status = status_train, X_train)
    fit <- fit_pre(survival::Surv(time, status) ~ ., dat, "cox")
    pred <- as.numeric(predict(
      fit, newdata = X_test, type = "link",
      penalty.par.val = "lambda.1se"
    ))
    metric <- cindex_value(time_test, status_test, pred)
  }

  cm <- as.data.frame(
    coef(fit, penalty.par.val = "lambda.1se"),
    stringsAsFactors = FALSE
  )
  if ("coefficient" %in% names(cm)) names(cm)[names(cm) == "coefficient"] <- "coef"
  if (!"coef" %in% names(cm)) stop("Unable to find RuleFit coefficient column.")

  cm <- cm[cm$coef != 0, , drop = FALSE]
  if ("rule" %in% names(cm)) {
    cm <- cm[!(cm$rule %in% c("(Intercept)", "Intercept")), , drop = FALSE]
  }

  is_rule <- if ("description" %in% names(cm)) !is.na(cm$description) else rep(FALSE, nrow(cm))
  descriptions <- if ("description" %in% names(cm)) cm$description[is_rule] else character(0)
  degrees <- vapply(
    descriptions,
    FUN.VALUE = integer(1),
    FUN = function(s) {
      if (is.na(s) || !nzchar(s)) return(NA_integer_)
      length(strsplit(s, "&", fixed = TRUE)[[1]])
    }
  )
  degrees <- degrees[!is.na(degrees)]
  n_linear <- sum(!is_rule)

  list(
    metric_value = metric,
    complexity = summarize_rule_complexity(degrees, n_linear),
    fit = fit,
    prediction = pred
  )
}

fit_gbm_baseline <- function(X_train, outcome_type,
                             y_train = NULL,
                             time_train = NULL,
                             status_train = NULL,
                             X_test,
                             y_test = NULL,
                             time_test = NULL,
                             status_test = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X_train <- as.data.frame(X_train)
  X_test <- as.data.frame(X_test)

  if (outcome_type == "binary") {
    dat <- data.frame(y = y_train, X_train)
    fit <- gbm::gbm(
      y ~ ., data = dat, distribution = "bernoulli",
      verbose = FALSE
    )
    pred <- as.numeric(predict(
      fit, X_test, n.trees = fit$n.trees, type = "response"
    ))
    metric <- auc_value(y_test, pred)
  } else if (outcome_type == "continuous") {
    dat <- data.frame(y = y_train, X_train)
    fit <- gbm::gbm(
      y ~ ., data = dat, distribution = "gaussian",
      verbose = FALSE
    )
    pred <- as.numeric(predict(
      fit, X_test, n.trees = fit$n.trees, type = "response"
    ))
    metric <- mse_value(y_test, pred)
  } else {
    dat <- data.frame(time = time_train, status = status_train, X_train)
    fit <- gbm::gbm(
      survival::Surv(time, status) ~ .,
      data = dat, distribution = "coxph",
      verbose = FALSE
    )
    pred <- as.numeric(predict(
      fit, X_test, n.trees = fit$n.trees, type = "link"
    ))
    metric <- cindex_value(time_test, status_test, pred)
  }

  list(
    metric_value = metric,
    complexity = empty_complexity(),
    fit = fit,
    prediction = pred
  )
}

fit_rf_baseline <- function(X_train, outcome_type,
                            y_train = NULL,
                            time_train = NULL,
                            status_train = NULL,
                            X_test,
                            y_test = NULL,
                            time_test = NULL,
                            status_test = NULL,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (outcome_type == "binary") {
    fit <- randomForest::randomForest(
      x = as.data.frame(X_train),
      y = as.factor(y_train)
    )
    pred <- as.numeric(predict(fit, as.data.frame(X_test), type = "prob")[, 2L])
    metric <- auc_value(y_test, pred)
  } else if (outcome_type == "continuous") {
    fit <- randomForest::randomForest(
      x = as.data.frame(X_train),
      y = y_train
    )
    pred <- as.numeric(predict(fit, as.data.frame(X_test)))
    metric <- mse_value(y_test, pred)
  } else {
    dat <- data.frame(
      time = time_train,
      status = status_train,
      X_train,
      check.names = FALSE
    )

    Surv <- survival::Surv

    fit <- randomForestSRC::rfsrc(
      Surv(time, status) ~ .,
      data = dat
    )
    pr <- predict(
      fit,
      newdata = as.data.frame(X_test, check.names = FALSE)
    )
    pred <- as.numeric(pr$predicted)
    metric <- cindex_value(time_test, status_test, pred)
  }

  list(
    metric_value = metric,
    complexity = empty_complexity(),
    fit = fit,
    prediction = pred
  )
}
