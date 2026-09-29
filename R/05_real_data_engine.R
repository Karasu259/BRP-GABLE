
prepare_mgus2_data <- function() {
  dat <- survival::mgus2[, c(
    "futime",
    "death",
    "age",
    "sex",
    "dxyr",
    "hgb",
    "creat",
    "mspike"
  )]

  dat <- stats::na.omit(dat)
  dat$sex <- as.integer(dat$sex == "M")

  names(dat)[names(dat) == "futime"] <- "time"
  names(dat)[names(dat) == "death"] <- "status"

  list(
    X = dat[, setdiff(names(dat), c("time", "status")), drop = FALSE],
    time = as.numeric(dat$time),
    status = as.integer(dat$status),
    dataset_name = "MGUS2"
  )
}

validate_real_survival_data <- function(real_data) {
  required <- c("X", "time", "status")
  missing <- setdiff(required, names(real_data))
  if (length(missing) > 0L) {
    stop(
      "real_data is missing required components: ",
      paste(missing, collapse = ", ")
    )
  }

  X <- as_numeric_matrix(real_data$X)
  time <- as.numeric(real_data$time)
  status <- as.integer(real_data$status)

  if (nrow(X) != length(time) || nrow(X) != length(status)) {
    stop("X, time, and status must contain the same number of observations.")
  }
  if (nrow(X) < 4L) stop("At least four observations are required.")
  if (ncol(X) < 1L) stop("At least one predictor is required.")
  if (anyNA(X) || any(!is.finite(X))) {
    stop("Predictors must not contain missing or non-finite values.")
  }
  if (anyNA(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("time must contain finite, non-negative values without missingness.")
  }
  if (anyNA(status) || !all(status %in% c(0L, 1L))) {
    stop("status must be coded as 0/1 without missingness.")
  }
  if (length(unique(status)) < 2L) {
    stop("status must contain both censored observations and events.")
  }

  list(
    X = as.data.frame(X, check.names = FALSE),
    time = time,
    status = status,
    dataset_name = real_data$dataset_name %||% "real_survival"
  )
}

real_result_list_to_df <- function(results, split_id, dataset_name) {
  out <- result_list_to_df(
    results = results,
    sim_id = split_id,
    outcome = "survival",
    scenario = NA_integer_,
    scenario_name = paste0("Real data: ", dataset_name)
  )
  out$Split <- as.integer(split_id)
  out$Dataset <- as.character(dataset_name)
  out[, c(
    "Split", "Dataset", "Sim", "Outcome", "Scenario", "Scenario_name",
    setdiff(names(out), c(
      "Split", "Dataset", "Sim", "Outcome", "Scenario", "Scenario_name"
    ))
  ), drop = FALSE]
}

run_one_real_survival_split <- function(
    sim_id,
    real_data,
    config,
    dataset_name = real_data$dataset_name %||% "real_survival"
) {
  real_data <- validate_real_survival_data(real_data)
  outcome_type <- "survival"
  metric_name <- "C-index"
  seeds <- make_seed_set(config$base_seed, sim_id)

  split <- make_train_test_split(
    n = nrow(real_data$X),
    train_fraction = config$train_fraction,
    seed = seeds$split,
    strata = NULL
  )

  X_train <- real_data$X[split$train, , drop = FALSE]
  X_test <- real_data$X[split$test, , drop = FALSE]
  time_train <- real_data$time[split$train]
  time_test <- real_data$time[split$test]
  status_train <- real_data$status[split$train]
  status_test <- real_data$status[split$test]
  y_train <- NULL
  y_test <- NULL

  foldid <- make_foldid(
    n = nrow(X_train),
    nfolds = config$nfolds,
    seed = seeds$cv,
    strata = NULL
  )

  pool_error <- NULL
  pool_start <- proc.time()[["elapsed"]]
  pool <- tryCatch(
    generate_bootstrap_rule_pool(
      X = X_train,
      outcome_type = outcome_type,
      y = y_train,
      time = time_train,
      status = status_train,
      B = config$B,
      Mmax = config$Mmax,
      degree = config$degree,
      min_support_prop = config$min_support_prop,
      min_span_prop = config$min_span_prop,
      improvement_tol = config$improvement_tol,
      seed = seeds$brp
    ),
    error = function(e) {
      pool_error <<- conditionMessage(e)
      NULL
    }
  )
  pool_elapsed <- proc.time()[["elapsed"]] - pool_start

  design <- NULL
  if (!is.null(pool)) {
    design <- tryCatch(
      build_brp_design(
        pool, X_train, X_test,
        linear_trim = config$linear_trim
      ),
      error = function(e) {
        pool_error <<- paste("Design construction:", conditionMessage(e))
        NULL
      }
    )
  }

  results <- list()

  results[["Linear"]] <- safe_model(
    "Linear", metric_name,
    fit_linear_baseline(
      X_train = X_train,
      outcome_type = outcome_type,
      y_train = y_train,
      time_train = time_train,
      status_train = status_train,
      X_test = X_test,
      y_test = y_test,
      time_test = time_test,
      status_test = status_test
    )
  )

  results[["GBM"]] <- safe_model(
    "GBM", metric_name,
    fit_gbm_baseline(
      X_train = X_train,
      outcome_type = outcome_type,
      y_train = y_train,
      time_train = time_train,
      status_train = status_train,
      X_test = X_test,
      y_test = y_test,
      time_test = time_test,
      status_test = status_test,
      seed = seeds$gbm
    )
  )

  results[["Random forest"]] <- safe_model(
    "Random forest", metric_name,
    fit_rf_baseline(
      X_train = X_train,
      outcome_type = outcome_type,
      y_train = y_train,
      time_train = time_train,
      status_train = status_train,
      X_test = X_test,
      y_test = y_test,
      time_test = time_test,
      status_test = status_test,
      seed = seeds$rf
    )
  )

  results[["GABLE"]] <- safe_model(
    "GABLE", metric_name,
    fit_gable_baseline(
      X_train = X_train,
      outcome_type = outcome_type,
      y_train = y_train,
      time_train = time_train,
      status_train = status_train,
      X_test = X_test,
      y_test = y_test,
      time_test = time_test,
      status_test = status_test,
      Mmax = config$Mmax,
      degree = config$degree,
      min_support_prop = config$min_support_prop,
      min_span_prop = config$min_span_prop,
      improvement_tol = config$improvement_tol,
      step_k = config$step_k
    )
  )

  results[["RuleFit"]] <- safe_model(
    "RuleFit", metric_name,
    fit_rulefit_default(
      X_train = X_train,
      outcome_type = outcome_type,
      y_train = y_train,
      time_train = time_train,
      status_train = status_train,
      X_test = X_test,
      y_test = y_test,
      time_test = time_test,
      status_test = status_test,
      seed = seeds$rulefit
    )
  )

  if (is.null(design)) {
    msg <- pool_error %||% "BRP-GABLE rule pool/design failed."
    results[["BRP-GABLE rule-only"]] <- failed_model_result(
      "BRP-GABLE rule-only", metric_name, msg, pool_elapsed
    )
    results[["BRP-GABLE full"]] <- failed_model_result(
      "BRP-GABLE full", metric_name, msg, pool_elapsed
    )
  } else {
    results[["BRP-GABLE rule-only"]] <- safe_model(
      "BRP-GABLE rule-only", metric_name,
      fit_brp_lasso(
        design = design,
        outcome_type = outcome_type,
        y_train = y_train,
        time_train = time_train,
        status_train = status_train,
        y_test = y_test,
        time_test = time_test,
        status_test = status_test,
        include_linear = FALSE,
        foldid = foldid,
        lambda_choice = config$lambda_choice
      )
    )
    results[["BRP-GABLE rule-only"]]$elapsed_time <-
      results[["BRP-GABLE rule-only"]]$elapsed_time + pool_elapsed

    results[["BRP-GABLE full"]] <- safe_model(
      "BRP-GABLE full", metric_name,
      fit_brp_lasso(
        design = design,
        outcome_type = outcome_type,
        y_train = y_train,
        time_train = time_train,
        status_train = status_train,
        y_test = y_test,
        time_test = time_test,
        status_test = status_test,
        include_linear = TRUE,
        foldid = foldid,
        lambda_choice = config$lambda_choice
      )
    )
    results[["BRP-GABLE full"]]$elapsed_time <-
      results[["BRP-GABLE full"]]$elapsed_time + pool_elapsed
  }

  real_result_list_to_df(results, sim_id, dataset_name)
}

summarize_real_survival_results <- function(results) {
  groups <- split(results, results$method, drop = TRUE)
  rows <- lapply(groups, function(z) {
    ok <- z$fit_ok & is.finite(z$metric_value)
    values <- z$metric_value[ok]
    data.frame(
      Outcome = "survival",
      Dataset = z$Dataset[1L],
      method = z$method[1L],
      metric_name = z$metric_name[1L],
      mean = if (length(values) > 0L) mean(values) else NA_real_,
      sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
      median = if (length(values) > 0L) stats::median(values) else NA_real_,
      n_success = length(values),
      n_total = nrow(z),
      n_failed = nrow(z) - sum(ok),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$method), , drop = FALSE]
}

run_real_survival <- function(
    real_data = prepare_mgus2_data(),
    n_splits = 100L,
    config = default_config(),
    n_workers = 1L,
    output_dir = file.path("results", "real_data", "experiment1"),
    dataset_name = real_data$dataset_name %||% "real_survival"
) {
  assert_packages(required_packages("survival"))
  real_data <- validate_real_survival_data(real_data)

  n_splits <- as.integer(n_splits)
  n_workers <- as.integer(n_workers)
  if (length(n_splits) != 1L || is.na(n_splits) || n_splits < 1L) {
    stop("n_splits must be a positive integer.")
  }
  if (length(n_workers) != 1L || is.na(n_workers) || n_workers < 1L) {
    stop("n_workers must be a positive integer.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (n_workers <= 1L) {
    rows <- lapply(seq_len(n_splits), function(k) {
      cat(sprintf(
        "[real survival: %s] split %d/%d\n",
        dataset_name, k, n_splits
      ))
      run_one_real_survival_split(
        sim_id = k,
        real_data = real_data,
        config = config,
        dataset_name = dataset_name
      )
    })
  } else {
    assert_packages(c("parallel", "doParallel", "foreach", "doRNG"))
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(doRNG))

    n_workers <- min(n_workers, n_splits)
    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(config$base_seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    rows <- foreach(
      k = seq_len(n_splits),
      .packages = required_packages("survival"),
      .export = ls(envir = .GlobalEnv)
    ) %dopar% {
      run_one_real_survival_split(
        sim_id = k,
        real_data = real_data,
        config = config,
        dataset_name = dataset_name
      )
    }
  }

  results <- do.call(rbind, rows)
  results <- results[order(results$Split, results$method), , drop = FALSE]
  summary_df <- summarize_real_survival_results(results)

  dataset_tag <- gsub("[^A-Za-z0-9_-]+", "_", tolower(dataset_name))
  tag <- paste0(
    dataset_tag,
    "_n", nrow(real_data$X),
    "_splits", n_splits,
    "_B", config$B,
    "_M", config$Mmax,
    "_deg", config$degree,
    "_support", gsub("\\.", "_", config$min_support_prop),
    "_span", gsub("\\.", "_", config$min_span_prop)
  )

  write.csv(
    results,
    file.path(output_dir, paste0("real_survival_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  saveRDS(
    list(
      results = results,
      summary = summary_df,
      config = config,
      dataset_name = dataset_name
    ),
    file.path(output_dir, paste0("real_survival_", tag, ".rds"))
  )
  write.csv(
    summary_df,
    file.path(output_dir, paste0("summary_real_survival_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  invisible(list(
    results = results,
    summary = summary_df,
    config = config,
    dataset_name = dataset_name
  ))
}
