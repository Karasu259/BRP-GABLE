
assert_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      "\nInstall them before running the simulation."
    )
  }
  invisible(TRUE)
}

as_numeric_matrix <- function(X) {
  X <- as.data.frame(X, check.names = FALSE)
  bad <- !vapply(X, is.numeric, logical(1))
  if (any(bad)) {
    stop("All predictors must be numeric. Non-numeric columns: ",
         paste(names(X)[bad], collapse = ", "))
  }
  out <- as.matrix(X)
  storage.mode(out) <- "double"
  if (is.null(colnames(out))) colnames(out) <- paste0("X", seq_len(ncol(out)))
  out
}

make_seed_set <- function(base_seed, sim_id) {
  stopifnot(length(base_seed) == 1L, length(sim_id) == 1L)
  list(
    data     = as.integer(base_seed + sim_id),
    split    = as.integer(base_seed + 100000L + sim_id),
    brp      = as.integer(base_seed + 200000L + sim_id),
    cv       = as.integer(base_seed + 300000L + sim_id),
    gable    = as.integer(base_seed + 400000L + sim_id),
    rulefit  = as.integer(base_seed + 500000L + sim_id),
    gbm      = as.integer(base_seed + 600000L + sim_id),
    rf       = as.integer(base_seed + 700000L + sim_id)
  )
}

make_train_test_split <- function(n, train_fraction = 0.5, seed = NULL,
                                  strata = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_train <- floor(n * train_fraction)

  if (is.null(strata)) {
    train <- sample.int(n, n_train, replace = FALSE)
  } else {
    strata <- as.factor(strata)
    train <- integer(0)
    levs <- levels(strata)
    for (lev in levs) {
      ids <- which(strata == lev)
      take <- max(1L, floor(length(ids) * train_fraction))
      take <- min(take, length(ids))
      train <- c(train, sample(ids, take, replace = FALSE))
    }
    if (length(train) > n_train) train <- sample(train, n_train)
    if (length(train) < n_train) {
      rest <- setdiff(seq_len(n), train)
      train <- c(train, sample(rest, n_train - length(train)))
    }
  }

  list(
    train = sort(unique(train)),
    test = setdiff(seq_len(n), train)
  )
}

make_foldid <- function(n, nfolds = 10L, seed = NULL, strata = NULL) {
  if (!is.null(seed)) set.seed(seed)
  nfolds <- min(as.integer(nfolds), n)
  foldid <- integer(n)

  if (is.null(strata)) {
    foldid[sample.int(n)] <- rep(seq_len(nfolds), length.out = n)
  } else {
    strata <- as.factor(strata)
    for (lev in levels(strata)) {
      ids <- sample(which(strata == lev))
      foldid[ids] <- rep(seq_len(nfolds), length.out = length(ids))
    }
  }
  foldid
}

empty_complexity <- function() {
  list(
    n_rules_total = NA_integer_,
    n_rule_terms_selected = NA_integer_,
    n_linear_selected = NA_integer_,
    n_rules_deg1 = NA_integer_,
    n_rules_deg2 = NA_integer_,
    n_rules_deg3 = NA_integer_,
    n_rules_deg4 = NA_integer_,
    n_rules_deg5plus = NA_integer_,
    prop_deg1 = NA_real_,
    prop_deg2 = NA_real_,
    prop_deg_le2 = NA_real_,
    avg_rule_degree = NA_real_
  )
}

summarize_rule_complexity <- function(degrees, n_linear_selected = 0L) {
  degrees <- as.integer(degrees)
  degrees <- degrees[is.finite(degrees) & degrees >= 1L]
  n_rule_terms <- length(degrees)
  n_linear_selected <- as.integer(n_linear_selected)[1L]
  if (!is.finite(n_linear_selected) || n_linear_selected < 0L) {
    n_linear_selected <- 0L
  }

  term_degrees <- c(degrees, rep.int(1L, n_linear_selected))
  n <- length(term_degrees)

  if (n == 0L) {
    return(list(
      n_rules_total = 0L,
      n_rule_terms_selected = 0L,
      n_linear_selected = n_linear_selected,
      n_rules_deg1 = 0L,
      n_rules_deg2 = 0L,
      n_rules_deg3 = 0L,
      n_rules_deg4 = 0L,
      n_rules_deg5plus = 0L,
      prop_deg1 = NA_real_,
      prop_deg2 = NA_real_,
      prop_deg_le2 = NA_real_,
      avg_rule_degree = NA_real_
    ))
  }

  list(
    n_rules_total = n,
    n_rule_terms_selected = n_rule_terms,
    n_linear_selected = n_linear_selected,
    n_rules_deg1 = sum(term_degrees == 1L),
    n_rules_deg2 = sum(term_degrees == 2L),
    n_rules_deg3 = sum(term_degrees == 3L),
    n_rules_deg4 = sum(term_degrees == 4L),
    n_rules_deg5plus = sum(term_degrees >= 5L),
    prop_deg1 = mean(term_degrees == 1L),
    prop_deg2 = mean(term_degrees == 2L),
    prop_deg_le2 = mean(term_degrees <= 2L),
    avg_rule_degree = mean(term_degrees)
  )
}

failed_model_result <- function(method, metric_name, error_message,
                                elapsed_time = NA_real_) {
  c(
    list(
      method = method,
      fit_ok = FALSE,
      metric_name = metric_name,
      metric_value = NA_real_,
      elapsed_time = elapsed_time,
      error_message = as.character(error_message)
    ),
    empty_complexity()
  )
}

successful_model_result <- function(method, metric_name, metric_value,
                                    complexity = empty_complexity(),
                                    elapsed_time = NA_real_) {
  c(
    list(
      method = method,
      fit_ok = TRUE,
      metric_name = metric_name,
      metric_value = as.numeric(metric_value),
      elapsed_time = as.numeric(elapsed_time),
      error_message = NA_character_
    ),
    complexity
  )
}

safe_model <- function(method, metric_name, expr) {
  start <- proc.time()[["elapsed"]]
  tryCatch(
    {
      ans <- force(expr)
      elapsed <- proc.time()[["elapsed"]] - start
      successful_model_result(
        method = method,
        metric_name = metric_name,
        metric_value = ans$metric_value,
        complexity = ans$complexity %||% empty_complexity(),
        elapsed_time = elapsed
      )
    },
    error = function(e) {
      elapsed <- proc.time()[["elapsed"]] - start
      failed_model_result(method, metric_name, conditionMessage(e), elapsed)
    }
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

result_list_to_df <- function(results, sim_id, outcome,
                              scenario, scenario_name) {
  rows <- lapply(results, function(z) {
    as.data.frame(c(
      list(
        Sim = sim_id,
        Outcome = outcome,
        Scenario = scenario,
        Scenario_name = scenario_name
      ),
      z
    ), stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- do.call(rbind, rows)

  numeric_cols <- c(
    "Sim", "Scenario", "metric_value", "elapsed_time",
    "n_rules_total", "n_rule_terms_selected", "n_linear_selected",
    "n_rules_deg1", "n_rules_deg2", "n_rules_deg3",
    "n_rules_deg4", "n_rules_deg5plus",
    "prop_deg1", "prop_deg2", "prop_deg_le2",
    "avg_rule_degree"
  )
  for (nm in intersect(numeric_cols, names(out))) {
    out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  }
  out$fit_ok <- as.logical(out$fit_ok)
  out
}

auc_value <- function(y, score) {
  as.numeric(pROC::auc(pROC::roc(y, score, quiet = TRUE, direction = "<")))
}

mse_value <- function(y, pred) mean((y - pred)^2)

cindex_value <- function(time, status, risk_score) {
  as.numeric(
    survival::concordance(
      survival::Surv(time, status) ~ I(-as.numeric(risk_score))
    )$concordance
  )
}
