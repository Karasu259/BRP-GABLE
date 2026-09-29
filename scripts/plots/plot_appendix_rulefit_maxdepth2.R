
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source(file.path("scripts", "plots", "plot_formal_simulation_results.R"))

depth2_result_files <- c(
  binary = file.path(
    file.path("results", "simulation", "rulefit_maxdepth2"),
    "rulefit_maxdepth2_binary_n1000_sim100.csv"
  ),
  continuous = file.path(
    file.path("results", "simulation", "rulefit_maxdepth2"),
    "rulefit_maxdepth2_continuous_n1000_sim100.csv"
  ),
  survival = file.path(
    file.path("results", "simulation", "rulefit_maxdepth2"),
    "rulefit_maxdepth2_survival_n1000_sim100.csv"
  )
)

figure_output_dir <- file.path(
  "figures", "simulation", "rulefit_maxdepth2"
)

prediction_methods <- c(
  "RuleFit",
  "RuleFit (max depth 2)",
  "BRP-GABLE rule-only",
  "BRP-GABLE full"
)

complexity_methods <- prediction_methods

method_labels <- c(
  "RuleFit" = "RuleFit\ndefault",
  "RuleFit (max depth 2)" = "RuleFit\ndepth 2",
  "BRP-GABLE rule-only" = "BRP-GABLE\nrules",
  "BRP-GABLE full" = "BRP-GABLE\nfull"
)

method_colors <- c(
  "RuleFit\ndefault" = "#D55E00",
  "RuleFit\ndepth 2" = "#EFA57A",
  "BRP-GABLE\nrules" = "#56B4E9",
  "BRP-GABLE\nfull" = "#0072B2"
)

outcome_titles <- c(
  binary = paste(
    "RuleFit maximum-depth sensitivity for binary outcome",
    "across Scenarios 1 to 4"
  ),
  continuous = paste(
    "RuleFit maximum-depth sensitivity for continuous outcome",
    "across Scenarios 1 to 4"
  ),
  survival = paste(
    "RuleFit maximum-depth sensitivity for time-to-event outcome",
    "across Scenarios 1 to 4"
  )
)

read_appendix_depth_results <- function(outcome) {
  formal <- read_formal_results(formal_result_files[[outcome]], outcome)
  formal <- formal[
    formal$method %in% c(
      "RuleFit", "BRP-GABLE rule-only", "BRP-GABLE full"
    ),
    , drop = FALSE
  ]

  depth2 <- read_formal_results(depth2_result_files[[outcome]], outcome)
  combined <- rbind(formal, depth2)

  observed <- aggregate(
    Sim ~ Scenario + method,
    data = combined,
    FUN = function(x) length(unique(x))
  )
  expected <- expand.grid(
    Scenario = 1:4,
    method = prediction_methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  checked <- merge(
    expected, observed,
    by = c("Scenario", "method"),
    all.x = TRUE,
    sort = TRUE
  )
  if (anyNA(checked$Sim) ||
      any(checked$Sim != expected_simulations_per_scenario)) {
    stop(
      "Incomplete appendix plotting data for ", outcome, ": ",
      paste(
        paste0("S", checked$Scenario, "/", checked$method, "=", checked$Sim),
        collapse = "; "
      )
    )
  }

  combined
}

save_appendix_outcome_figure <- function(plot, outcome) {
  dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("appendix_rulefit_maxdepth2_", outcome)
  pdf_path <- file.path(figure_output_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_output_dir, paste0(stem, ".png"))

  ggsave(
    filename = pdf_path,
    plot = plot,
    width = figure_style$export_width,
    height = figure_style$export_height,
    units = "in",
    bg = "white"
  )
  ggsave(
    filename = png_path,
    plot = plot,
    device = grDevices::png,
    width = figure_style$export_width,
    height = figure_style$export_height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  invisible(c(pdf = pdf_path, png = png_path))
}

plot_rulefit_maxdepth2_appendix <- function() {
  missing_files <- c(formal_result_files, depth2_result_files)
  missing_files <- missing_files[!file.exists(missing_files)]
  if (length(missing_files)) {
    stop("Missing result file(s): ", paste(missing_files, collapse = ", "))
  }

  outputs <- lapply(names(formal_result_files), function(outcome) {
    message("Plotting RuleFit maximum-depth sensitivity for ", outcome)
    dat <- read_appendix_depth_results(outcome)
    fig <- make_outcome_figure(dat, outcome)
    save_appendix_outcome_figure(fig, outcome)
  })
  names(outputs) <- names(formal_result_files)
  invisible(outputs)
}

if (sys.nframe() == 0L) {
  appendix_rulefit_depth_outputs <- plot_rulefit_maxdepth2_appendix()
}
