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
source(file.path("R", "05_real_data_engine.R"))
source(file.path("R", "10_real_data_explanation_engine.R"))

real_data <- prepare_mgus2_data()

config <- default_config()
config$B <- 50L
config$Mmax <- 10L
config$degree <- 2L
config$min_support_prop <- 0.05
config$min_span_prop <- 0.20
config$base_seed <- 423L

full_result <- fit_full_brp_gable_survival(
  real_data = real_data,
  config = config,
  output_dir = file.path("results", "real_data", "experiment2"),
  dataset_name = "MGUS2",
  full_fit_id = 0L,
  n_ale_bins = 20L,
  n_ale_variables = 3L,
  n_km_rules = 3L
)

cat("Saved full-data result:", attr(full_result, "output_rds"), "\n")
cat(
  "Selected terms:", nrow(full_result$selected_terms),
  "(rules:", nrow(full_result$selected_rules),
  ", linear:", sum(full_result$selected_terms$term_type == "linear"), ")\n"
)
