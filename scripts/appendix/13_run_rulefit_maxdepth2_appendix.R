script_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) return(sub("^--file=", "", file_arg[[1L]]))
  frame_paths <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(frame_paths)) return(tail(frame_paths, 1L)[[1L]])
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) return(active_path)
  }
  stop("Unable to determine the script location. Run this saved script with Source or Rscript.")
})
setwd(file.path(dirname(script_path), "..", ".."))
rm(script_path)

source(file.path("R", "00_utils.R"))
source(file.path("R", "02_models.R"))
source(file.path("R", "03_data_generators.R"))
source(file.path("R", "04_simulation_engine.R"))

appendix_method <- "RuleFit (max depth 2)"
appendix_output_dir <- file.path("results", "simulation", "rulefit_maxdepth2")

assert_matching_pre_installation <- function() {
  assert_packages("pre")
  if (as.character(packageVersion("pre")) != "1.0.7") {
    stop(
      "This sensitivity analysis must use pre 1.0.7, matching the main ",
      "simulation; installed version is ", as.character(packageVersion("pre")), "."
    )
  }

  f <- formals(pre::pre)
  observed <- list(
    ntrees = as.integer(f[["ntrees"]]),
    maxdepth = as.integer(f[["maxdepth"]]),
    tree.unbiased = as.logical(f[["tree.unbiased"]]),
    removecomplements = as.logical(f[["removecomplements"]]),
    removeduplicates = as.logical(f[["removeduplicates"]]),
    type = as.character(f[["type"]]),
    sampfrac = as.numeric(f[["sampfrac"]]),
    learnrate = as.numeric(f[["learnrate"]]),
    nfolds = as.integer(f[["nfolds"]])
  )
  expected <- list(
    ntrees = 500L,
    maxdepth = 3L,
    tree.unbiased = TRUE,
    removecomplements = TRUE,
    removeduplicates = TRUE,
    type = "both",
    sampfrac = 0.5,
    learnrate = 0.01,
    nfolds = 10L
  )
  if (!isTRUE(all.equal(observed, expected, tolerance = 0))) {
    stop(
      "Installed pre defaults differ from those used in the main simulation.\n",
      "Observed: ",
      paste(names(observed), unlist(observed), sep = "=", collapse = ", "),
      "\nExpected: ",
      paste(names(expected), unlist(expected), sep = "=", collapse = ", ")
    )
  }
  observed
}

formal_specs <- list(
  binary = list(
    base_seed = 123L,
    formal_csv = file.path(
      file.path("results", "simulation", "n1000", "binary"),
      "simulation_binary_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  ),
  continuous = list(
    base_seed = 223L,
    formal_csv = file.path(
      file.path("results", "simulation", "n1000", "continuous"),
      "simulation_continuous_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  ),
  survival = list(
    base_seed = 323L,
    formal_csv = file.path(
      file.path("results", "simulation", "n1000", "survival"),
      "simulation_survival_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  )
)

make_formal_config <- function(base_seed) {
  config <- default_config()
  config$n <- 1000L
  config$p <- 20L
  config$train_fraction <- 0.5
  config$B <- 50L
  config$Mmax <- 10L
  config$degree <- 2L
  config$min_support_prop <- 0.05
  config$min_span_prop <- 0.20
  config$x_dist <- "normal"
  config$c0 <- 0.36
  config$noise_sd <- 1
  config$baseline_hazard <- 0.1
  config$censor_upper <- 20
  config$base_seed <- as.integer(base_seed)
  config
}

generate_formal_dataset <- function(outcome_type, scenario, config, seed) {
  switch(
    outcome_type,
    binary = generate_binary_data(
      n = config$n, p = config$p, scenario = scenario,
      x_dist = config$x_dist, c0 = config$c0, seed = seed
    ),
    continuous = generate_continuous_data(
      n = config$n, p = config$p, scenario = scenario,
      x_dist = config$x_dist, c0 = config$c0,
      noise_sd = config$noise_sd, seed = seed
    ),
    survival = generate_survival_data(
      n = config$n, p = config$p, scenario = scenario,
      x_dist = config$x_dist, c0 = config$c0,
      baseline_hazard = config$baseline_hazard,
      censor_upper = config$censor_upper, seed = seed
    )
  )
}

base_result_key <- function(x) {
  paste(x$Outcome, x$Scenario, x$Sim, sep = "|")
}

run_one_depth2_rulefit <- function(sim_id, outcome_type, scenario, config) {
  seeds <- make_seed_set(config$base_seed, sim_id)
  dat <- generate_formal_dataset(
    outcome_type, scenario, config, seed = seeds$data
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

  metric_name <- switch(
    outcome_type,
    binary = "AUC",
    continuous = "MSE",
    survival = "C-index"
  )

  rulefit_result <- safe_model(
    appendix_method, metric_name,
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
      seed = seeds$rulefit,
      maxdepth = 2L
    )
  )

  linear_check <- fit_linear_baseline(
    X_train, outcome_type,
    y_train, time_train, status_train,
    X_test, y_test, time_test, status_test
  )$metric_value

  row <- result_list_to_df(
    setNames(list(rulefit_result), appendix_method),
    sim_id = sim_id,
    outcome = outcome_type,
    scenario = scenario,
    scenario_name = dat$scenario_name
  )
  row$linear_check_metric <- as.numeric(linear_check)
  row
}

summarize_depth2_results <- function(results) {
  groups <- split(results, results$Scenario)
  rows <- lapply(groups, function(z) {
    ok <- z$fit_ok & is.finite(z$metric_value)
    metric <- z$metric_value[ok]
    terms <- z$n_rules_total[ok & is.finite(z$n_rules_total)]
    degree1 <- z$prop_deg1[ok & is.finite(z$prop_deg1)]
    data.frame(
      Outcome = z$Outcome[1L],
      Scenario = z$Scenario[1L],
      Scenario_name = z$Scenario_name[1L],
      method = z$method[1L],
      metric_name = z$metric_name[1L],
      metric_mean = mean(metric),
      metric_sd = stats::sd(metric),
      metric_median = stats::median(metric),
      final_term_count_mean = mean(terms),
      final_term_count_sd = stats::sd(terms),
      degree1_proportion_mean = mean(degree1),
      degree1_proportion_sd = stats::sd(degree1),
      n_success = sum(ok),
      n_total = nrow(z),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$Scenario), , drop = FALSE]
}

make_depth_comparison <- function(depth2, formal) {
  default <- formal[formal$method == "RuleFit", , drop = FALSE]
  joined <- merge(
    depth2,
    default,
    by = c("Outcome", "Scenario", "Scenario_name", "Sim"),
    suffixes = c("_depth2", "_depth3"),
    sort = TRUE
  )
  if (nrow(joined) != 400L) {
    stop("Could not form all 400 paired depth-2/depth-3 comparisons.")
  }
  data.frame(
    Outcome = joined$Outcome,
    Scenario = joined$Scenario,
    Scenario_name = joined$Scenario_name,
    Sim = joined$Sim,
    metric_name = joined$metric_name_depth2,
    metric_depth2 = joined$metric_value_depth2,
    metric_depth3 = joined$metric_value_depth3,
    metric_difference_depth2_minus_depth3 =
      joined$metric_value_depth2 - joined$metric_value_depth3,
    final_terms_depth2 = joined$n_rules_total_depth2,
    final_terms_depth3 = joined$n_rules_total_depth3,
    final_terms_difference_depth2_minus_depth3 =
      joined$n_rules_total_depth2 - joined$n_rules_total_depth3,
    degree1_proportion_depth2 = joined$prop_deg1_depth2,
    degree1_proportion_depth3 = joined$prop_deg1_depth3,
    degree1_difference_depth2_minus_depth3 =
      joined$prop_deg1_depth2 - joined$prop_deg1_depth3,
    stringsAsFactors = FALSE
  )
}

run_one_outcome <- function(outcome_type, spec, n_workers, overwrite = FALSE) {
  config <- make_formal_config(spec$base_seed)
  formal <- read.csv(
    spec$formal_csv, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (nrow(formal) != 2800L || sum(formal$method == "RuleFit") != 400L) {
    stop("Unexpected main simulation dimensions in ", spec$formal_csv)
  }

  result_csv <- file.path(
    appendix_output_dir,
    paste0("rulefit_maxdepth2_", outcome_type, "_n1000_sim100.csv")
  )
  summary_csv <- file.path(
    appendix_output_dir,
    paste0("summary_rulefit_maxdepth2_", outcome_type, "_n1000_sim100.csv")
  )
  comparison_csv <- file.path(
    appendix_output_dir,
    paste0("paired_rulefit_depth2_vs_depth3_", outcome_type, "_n1000_sim100.csv")
  )
  result_rds <- sub("[.]csv$", ".rds", result_csv)

  output_files <- c(result_csv, summary_csv, comparison_csv, result_rds)
  if (!overwrite && all(file.exists(output_files))) {
    saved <- read.csv(result_csv, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(saved) == 400L && all(saved$fit_ok)) {
      cat("Using completed saved result for ", outcome_type, ".\n", sep = "")
      return(saved)
    }
  }

  jobs <- expand.grid(
    Scenario = 1:4,
    Sim = 1:100,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  packages <- unique(c(
    "pre",
    if (outcome_type == "binary") "pROC" else NULL,
    if (outcome_type == "survival") "survival" else NULL
  ))
  assert_packages(packages)

  cat("\n=== RuleFit maxdepth=2: ", outcome_type, " (400 fits) ===\n", sep = "")
  if (n_workers <= 1L) {
    rows <- lapply(seq_len(nrow(jobs)), function(k) {
      if (k %% 25L == 0L) {
        cat("Completed ", k, "/", nrow(jobs), "\n", sep = "")
      }
      run_one_depth2_rulefit(
        jobs$Sim[k], outcome_type, jobs$Scenario[k], config
      )
    })
  } else {
    assert_packages(c("parallel", "doParallel", "foreach", "doRNG"))
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(doRNG))
    n_workers <- min(as.integer(n_workers), nrow(jobs))
    cl <- parallel::makeCluster(n_workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(config$base_seed + 9300000L)
    rows <- foreach(
      k = seq_len(nrow(jobs)),
      .packages = packages,
      .export = ls(envir = .GlobalEnv)
    ) %dopar% {
      run_one_depth2_rulefit(
        jobs$Sim[k], outcome_type, jobs$Scenario[k], config
      )
    }
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  }

  results <- do.call(rbind, rows)
  results <- results[order(results$Scenario, results$Sim), , drop = FALSE]
  if (nrow(results) != 400L || any(!results$fit_ok) ||
      any(!is.finite(results$metric_value))) {
    stop("Incomplete or failed maxdepth=2 fits for ", outcome_type)
  }

  stored_linear <- formal[formal$method == "Linear", , drop = FALSE]
  linear_match <- match(base_result_key(results), base_result_key(stored_linear))
  if (anyNA(linear_match)) {
    stop("Unable to match stored Linear rows for ", outcome_type)
  }
  max_linear_difference <- max(abs(
    results$linear_check_metric - stored_linear$metric_value[linear_match]
  ))
  if (!is.finite(max_linear_difference) || max_linear_difference > 1e-12) {
    stop(
      "Data/split reproduction failed for ", outcome_type,
      "; maximum Linear metric difference = ", max_linear_difference
    )
  }
  results$linear_check_metric <- NULL

  summary_df <- summarize_depth2_results(results)
  comparison_df <- make_depth_comparison(results, formal)
  write.csv(results, result_csv, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(summary_df, summary_csv, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(
    comparison_df, comparison_csv, row.names = FALSE, fileEncoding = "UTF-8"
  )
  saveRDS(
    list(
      results = results,
      summary = summary_df,
      paired_default_comparison = comparison_df,
      config = config,
      rulefit_override = list(maxdepth = 2L),
      main_formal_csv = spec$formal_csv,
      max_linear_reproduction_difference = max_linear_difference
    ),
    result_rds
  )
  cat(
    "Saved ", outcome_type,
    "; maximum Linear reproduction difference = ",
    format(max_linear_difference, scientific = TRUE), "\n", sep = ""
  )
  results
}

args <- commandArgs(trailingOnly = TRUE)
smoke_test <- "--smoke" %in% args
overwrite <- "--overwrite" %in% args
worker_arg <- grep("^--workers=", args, value = TRUE)
n_workers <- if (length(worker_arg)) {
  as.integer(sub("^--workers=", "", worker_arg[1L]))
} else {
  20L
}
if (!is.finite(n_workers) || n_workers < 1L) stop("Invalid --workers value.")

outcome_arg <- grep("^--outcomes=", args, value = TRUE)
outcomes_to_run <- if (length(outcome_arg)) {
  strsplit(sub("^--outcomes=", "", outcome_arg[1L]), ",", fixed = TRUE)[[1L]]
} else {
  names(formal_specs)
}
if (!all(outcomes_to_run %in% names(formal_specs))) {
  stop("--outcomes must contain binary, continuous, and/or survival.")
}

pre_defaults <- assert_matching_pre_installation()
dir.create(appendix_output_dir, recursive = TRUE, showWarnings = FALSE)

if (smoke_test) {
  smoke_rows <- lapply(outcomes_to_run, function(outcome_type) {
    spec <- formal_specs[[outcome_type]]
    run_one_depth2_rulefit(
      sim_id = 1L,
      outcome_type = outcome_type,
      scenario = 1L,
      config = make_formal_config(spec$base_seed)
    )
  })
  smoke_results <- do.call(rbind, smoke_rows)
  if (any(!smoke_results$fit_ok) ||
      any(!is.finite(smoke_results$metric_value)) ||
      any(!is.finite(smoke_results$linear_check_metric))) {
    stop("RuleFit maxdepth=2 smoke test failed.")
  }
  print(smoke_results[, c(
    "Outcome", "Scenario", "Sim", "method", "metric_name", "metric_value",
    "n_rules_total", "prop_deg1"
  )])
  cat("RuleFit maxdepth=2 smoke test passed.\n")
  quit(save = "no", status = 0L)
}

cat("pre version: ", as.character(packageVersion("pre")), "\n", sep = "")
cat(
  "Defaults retained except maxdepth: ",
  paste(names(pre_defaults), unlist(pre_defaults), sep = "=", collapse = ", "),
  "\nIntentional override: maxdepth=2\n", sep = ""
)
cat("Workers: ", n_workers, "\n", sep = "")

all_results <- lapply(outcomes_to_run, function(outcome_type) {
  run_one_outcome(
    outcome_type = outcome_type,
    spec = formal_specs[[outcome_type]],
    n_workers = n_workers,
    overwrite = overwrite
  )
})
names(all_results) <- outcomes_to_run

combined <- do.call(rbind, all_results)
combined_path <- file.path(
  appendix_output_dir,
  "rulefit_maxdepth2_all_outcomes_n1000_sim100.csv"
)
write.csv(combined, combined_path, row.names = FALSE, fileEncoding = "UTF-8")

writeLines(
  c(
    paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("pre_version=", as.character(packageVersion("pre"))),
    "main_rulefit_maxdepth=3",
    "sensitivity_rulefit_maxdepth=2",
    "all_other_pre_parameters=package_defaults",
    "same_data_split_and_rulefit_seeds_as_main_simulation=TRUE",
    paste0("outcomes=", paste(outcomes_to_run, collapse = ",")),
    "scenarios=1,2,3,4",
    "n_sim_each=100",
    "n=1000",
    "train_fraction=0.5",
    paste0("workers=", n_workers),
    capture.output(sessionInfo())
  ),
  file.path(appendix_output_dir, "manifest.txt")
)

cat("\nCompleted RuleFit maxdepth=2 appendix sensitivity analysis.\n")
cat("Combined result: ", normalizePath(combined_path), "\n", sep = "")
