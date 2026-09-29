
source(file.path("scripts", "plots", "plot_formal_simulation_results.R"))

formal_result_files <- c(
  binary = file.path(
    file.path("results", "simulation", "n500"), "binary",
    "simulation_binary_n500_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  continuous = file.path(
    file.path("results", "simulation", "n500"), "continuous",
    "simulation_continuous_n500_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  survival = file.path(
    file.path("results", "simulation", "n500"), "survival",
    "simulation_survival_n500_B50_M10_deg2_support0_05_span0_2.csv"
  )
)

figure_output_dir <- file.path("figures", "simulation", "n500")

figure_style$scenario_text_size <- 3.8

outcome_titles <- c(
  binary = "Appendix results for binary outcome (n = 500)",
  continuous = "Appendix results for continuous outcome (n = 500)",
  survival = "Appendix results for time-to-event outcome (n = 500)"
)

appendix_n500_figure_outputs <- plot_available_formal_results()
