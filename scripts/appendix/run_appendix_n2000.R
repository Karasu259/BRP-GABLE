
rm(list = ls())

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
source(file.path("R", "01_rule_generation.R"))
source(file.path("R", "02_models.R"))
source(file.path("R", "03_data_generators.R"))
source(file.path("R", "04_simulation_engine.R"))

appendix_specs <- list(
  continuous = list(
    base_seed = 223L,
    output_dir = file.path("results", "simulation", "n2000", "continuous")
  ),
  binary = list(
    base_seed = 123L,
    output_dir = file.path("results", "simulation", "n2000", "binary")
  ),
  survival = list(
    base_seed = 323L,
    output_dir = file.path("results", "simulation", "n2000", "survival")
  )
)

make_appendix_config <- function(base_seed) {
  config <- default_config()
  config$n <- 2000L
  config$p <- 20L
  config$train_fraction <- 0.5
  config$B <- 50L
  config$Mmax <- 10L
  config$degree <- 2L
  config$min_support_prop <- 0.05
  config$min_span_prop <- 0.20
  config$improvement_tol <- 1e-8
  config$nfolds <- 10L
  config$lambda_choice <- "lambda.1se"
  config$linear_trim <- 0.025
  config$step_k <- 6
  config$x_dist <- "normal"
  config$c0 <- 0.36
  config$noise_sd <- 1
  config$baseline_hazard <- 0.1
  config$censor_upper <- 20
  config$base_seed <- as.integer(base_seed)
  config
}

validate_appendix_result <- function(result, outcome_type, n_sim_each) {
  dat <- result$results
  expected_rows <- 4L * n_sim_each * 7L
  expected_per_method <- 4L * n_sim_each

  if (result$config$n != 2000L) {
    stop("Appendix result does not use n=2000.")
  }
  if (nrow(dat) != expected_rows) {
    stop(
      "Unexpected row count for ", outcome_type,
      ": observed ", nrow(dat), ", expected ", expected_rows
    )
  }
  if (any(!dat$fit_ok) || any(!is.finite(dat$metric_value))) {
    failed <- dat[!dat$fit_ok | !is.finite(dat$metric_value), ]
    stop(
      "At least one ", outcome_type, " fit failed: ",
      paste(unique(failed$error_message), collapse = " | ")
    )
  }
  if (anyDuplicated(dat[c("Outcome", "Scenario", "Sim", "method")])) {
    stop("Duplicate result keys for ", outcome_type, ".")
  }
  method_counts <- table(dat$method)
  if (length(method_counts) != 7L || any(method_counts != expected_per_method)) {
    stop("Incomplete method counts for ", outcome_type, ".")
  }
  scenario_counts <- table(dat$Scenario)
  if (length(scenario_counts) != 4L ||
      any(scenario_counts != n_sim_each * 7L)) {
    stop("Incomplete scenario counts for ", outcome_type, ".")
  }
  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
worker_arg <- grep("^--workers=", args, value = TRUE)
n_workers <- if (length(worker_arg)) {
  as.integer(sub("^--workers=", "", worker_arg[1L]))
} else {
  20L
}
if (!is.finite(n_workers) || n_workers < 1L) stop("Invalid --workers value.")

outcome_arg <- grep("^--outcomes=", args, value = TRUE)
outcomes_to_run <- if (length(outcome_arg)) {
  trimws(strsplit(sub("^--outcomes=", "", outcome_arg[1L]), ",")[[1L]])
} else {
  names(appendix_specs)
}
invalid_outcomes <- setdiff(outcomes_to_run, names(appendix_specs))
if (!length(outcomes_to_run) || length(invalid_outcomes)) {
  stop(
    "Invalid --outcomes value. Available outcomes are: ",
    paste(names(appendix_specs), collapse = ", ")
  )
}

dir.create(file.path("results", "simulation", "n2000"), recursive = TRUE, showWarnings = FALSE)
run_started <- Sys.time()
run_summary <- list()

partial_summary_path <- file.path(
  file.path("results", "simulation", "n2000"), "run_summary_partial.csv"
)
if (file.exists(partial_summary_path)) {
  existing_summary <- read.csv(
    partial_summary_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (outcome_type in intersect(
    names(appendix_specs), existing_summary$outcome
  )) {
    run_summary[[outcome_type]] <- existing_summary[
      existing_summary$outcome == outcome_type,
      ,
      drop = FALSE
    ][1L, , drop = FALSE]
  }
}

for (outcome_type in outcomes_to_run) {
  spec <- appendix_specs[[outcome_type]]
  config <- make_appendix_config(spec$base_seed)
  outcome_started <- Sys.time()

  cat(
    "\n=== Starting appendix ", outcome_type,
    " simulations: n=2000, 4 scenarios x 100 ===\n", sep = ""
  )
  result <- run_simulation(
    outcome_type = outcome_type,
    scenarios = 1:4,
    n_sim_each = 100L,
    config = config,
    n_workers = n_workers,
    output_dir = spec$output_dir
  )
  validate_appendix_result(result, outcome_type, 100L)

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), outcome_started, units = "mins")
  )
  run_summary[[outcome_type]] <- data.frame(
    outcome = outcome_type,
    n = config$n,
    scenarios = "1:4",
    simulations_per_scenario = 100L,
    rows = nrow(result$results),
    failed = sum(!result$results$fit_ok),
    elapsed_minutes = elapsed_minutes,
    output_dir = spec$output_dir,
    stringsAsFactors = FALSE
  )
  write.csv(
    do.call(rbind, run_summary),
    file.path("results", "simulation", "n2000", "run_summary_partial.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  cat(
    "Completed and validated appendix ", outcome_type,
    " in ", round(elapsed_minutes, 1), " minutes.\n", sep = ""
  )
}

run_summary_df <- do.call(rbind, run_summary)
run_summary_df <- run_summary_df[
  match(names(appendix_specs), run_summary_df$outcome, nomatch = 0L),
  ,
  drop = FALSE
]

if (all(names(appendix_specs) %in% run_summary_df$outcome)) {
  write.csv(
    run_summary_df,
    file.path("results", "simulation", "n2000", "run_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  writeLines(
    c(
      paste0("started=", format(run_started, "%Y-%m-%d %H:%M:%S %Z")),
      paste0("completed=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "sample_size=2000",
      "predictors=20",
      "train_fraction=0.5",
      "B=50",
      "Mmax=10",
      "degree=2",
      "min_support_prop=0.05",
      "min_span_prop=0.20",
      "n_sim_each=100",
      paste0("workers=", n_workers),
      capture.output(sessionInfo())
    ),
    file.path("results", "simulation", "n2000", "manifest.txt")
  )
}

print(run_summary_df)
if (all(names(appendix_specs) %in% run_summary_df$outcome)) {
  cat("\nAll appendix n=2000 simulations completed successfully.\n")
} else {
  cat("\nSelected appendix n=2000 simulations completed successfully.\n")
}
