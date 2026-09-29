
default_config <- function() {
  list(
    n = 1000L,
    p = 20L,
    train_fraction = 0.5,
    B = 50L,
    Mmax = 10L,
    degree = 2L,
    min_support_prop = 0.05,
    min_span_prop = 0.20,
    improvement_tol = 1e-8,
    nfolds = 10L,
    lambda_choice = "lambda.1se",
    linear_trim = 0.025,
    step_k = 6,
    x_dist = "normal",
    c0 = 0.36,
    noise_sd = 1,
    baseline_hazard = 0.1,
    censor_upper = 20,
    base_seed = 123L
  )
}

required_packages <- function(outcome_type) {
  unique(c(
    "MASS", "glmnet", "pre", "pROC", "survival",
    "gbm", "randomForest",
    if (outcome_type == "survival") "randomForestSRC" else NULL
  ))
}

run_one_simulation <- function(sim_id, outcome_type, scenario, config) {
  seeds <- make_seed_set(config$base_seed, sim_id)

  dat <- switch(
    outcome_type,
    binary = generate_binary_data(
      n = config$n,
      p = config$p,
      scenario = scenario,
      x_dist = config$x_dist,
      c0 = config$c0,
      seed = seeds$data
    ),
    continuous = generate_continuous_data(
      n = config$n,
      p = config$p,
      scenario = scenario,
      x_dist = config$x_dist,
      c0 = config$c0,
      noise_sd = config$noise_sd,
      seed = seeds$data
    ),
    survival = generate_survival_data(
      n = config$n,
      p = config$p,
      scenario = scenario,
      x_dist = config$x_dist,
      c0 = config$c0,
      baseline_hazard = config$baseline_hazard,
      censor_upper = config$censor_upper,
      seed = seeds$data
    )
  )

  split <- make_train_test_split(
    n = nrow(dat$X),
    train_fraction = config$train_fraction,
    seed = seeds$split,
    strata = NULL
  )

  X_train <- dat$X[split$train, , drop = FALSE]
  X_test <- dat$X[split$test, , drop = FALSE]

  y_train <- if (!is.null(dat$y)) dat$y[split$train] else NULL
  y_test <- if (!is.null(dat$y)) dat$y[split$test] else NULL
  time_train <- if (!is.null(dat$time)) dat$time[split$train] else NULL
  time_test <- if (!is.null(dat$time)) dat$time[split$test] else NULL
  status_train <- if (!is.null(dat$status)) dat$status[split$train] else NULL
  status_test <- if (!is.null(dat$status)) dat$status[split$test] else NULL

  foldid <- make_foldid(
    n = nrow(X_train),
    nfolds = config$nfolds,
    seed = seeds$cv,
    strata = NULL
  )

  metric_name <- switch(
    outcome_type,
    binary = "AUC",
    continuous = "MSE",
    survival = "C-index"
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
      X_train, outcome_type,
      y_train, time_train, status_train,
      X_test, y_test, time_test, status_test
    )
  )

  results[["GBM"]] <- safe_model(
    "GBM", metric_name,
    fit_gbm_baseline(
      X_train, outcome_type,
      y_train, time_train, status_train,
      X_test, y_test, time_test, status_test,
      seed = seeds$gbm
    )
  )

  results[["Random forest"]] <- safe_model(
    "Random forest", metric_name,
    fit_rf_baseline(
      X_train, outcome_type,
      y_train, time_train, status_train,
      X_test, y_test, time_test, status_test,
      seed = seeds$rf
    )
  )

  results[["GABLE"]] <- safe_model(
    "GABLE", metric_name,
    fit_gable_baseline(
      X_train, outcome_type,
      y_train, time_train, status_train,
      X_test, y_test, time_test, status_test,
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
      X_train, outcome_type,
      y_train, time_train, status_train,
      X_test, y_test, time_test, status_test,
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
        design, outcome_type,
        y_train, time_train, status_train,
        y_test, time_test, status_test,
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
        design, outcome_type,
        y_train, time_train, status_train,
        y_test, time_test, status_test,
        include_linear = TRUE,
        foldid = foldid,
        lambda_choice = config$lambda_choice
      )
    )
    results[["BRP-GABLE full"]]$elapsed_time <-
      results[["BRP-GABLE full"]]$elapsed_time + pool_elapsed
  }

  result_list_to_df(
    results,
    sim_id,
    outcome_type,
    scenario = scenario,
    scenario_name = dat$scenario_name
  )
}

run_simulation <- function(
    outcome_type = c("binary", "continuous", "survival"),
    scenarios = 1:4,
    n_sim_each = 100L,
    config = default_config(),
    n_workers = 1L,
    output_dir = "results"
) {
  outcome_type <- match.arg(outcome_type)
  assert_packages(required_packages(outcome_type))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  jobs <- expand.grid(
    Scenario = scenarios,
    Sim = seq_len(n_sim_each),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (n_workers <= 1L) {
    rows <- lapply(seq_len(nrow(jobs)), function(k) {
      cat(
        sprintf(
          "[%s] Scenario=%s, simulation %d/%d\n",
          outcome_type, jobs$Scenario[k], jobs$Sim[k], n_sim_each
        )
      )
      run_one_simulation(
        sim_id = jobs$Sim[k],
        outcome_type = outcome_type,
        scenario = jobs$Scenario[k],
        config = config
      )
    })
  } else {
    assert_packages(c("parallel", "doParallel", "foreach", "doRNG"))
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(doRNG))

    cl <- parallel::makeCluster(n_workers)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(config$base_seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    rows <- foreach(
      k = seq_len(nrow(jobs)),
      .packages = required_packages(outcome_type),
      .export = ls(envir = .GlobalEnv)
    ) %dopar% {
      run_one_simulation(
        sim_id = jobs$Sim[k],
        outcome_type = outcome_type,
        scenario = jobs$Scenario[k],
        config = config
      )
    }
  }

  results <- do.call(rbind, rows)
  results <- results[
    order(results$Scenario, results$Sim, results$method),
  ]

  tag <- paste0(
    outcome_type,
    "_n", config$n,
    "_B", config$B,
    "_M", config$Mmax,
    "_deg", config$degree,
    "_support", gsub("\\.", "_", config$min_support_prop),
    "_span", gsub("\\.", "_", config$min_span_prop)
  )

  write.csv(
    results,
    file.path(output_dir, paste0("simulation_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  saveRDS(
    list(results = results, config = config),
    file.path(output_dir, paste0("simulation_", tag, ".rds"))
  )

  summary_df <- aggregate(
    metric_value ~ Outcome + Scenario + Scenario_name + method,
    data = results[results$fit_ok & is.finite(results$metric_value), ],
    FUN = function(x) c(
      mean = mean(x),
      sd = stats::sd(x),
      median = stats::median(x),
      n = length(x)
    )
  )
  write.csv(
    summary_df,
    file.path(output_dir, paste0("summary_", tag, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  invisible(list(results = results, summary = summary_df, config = config))
}
