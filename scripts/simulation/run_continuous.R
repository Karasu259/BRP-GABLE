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

config <- default_config()
config$n <- 1000L
config$p <- 20L
config$B <- 50L
config$Mmax <- 10L
config$degree <- 2L
config$min_support_prop <- 0.05
config$min_span_prop <- 0.20
config$x_dist <- "normal"
config$c0 <- 0.36
config$noise_sd <- 1
config$base_seed <- 223L

result <- run_simulation(
  outcome_type = "continuous",
  scenarios = 1:4,
  n_sim_each = 100L,
  config = config,
  n_workers = 20L,
  output_dir = file.path("results", "simulation", "n1000", "continuous")
)
