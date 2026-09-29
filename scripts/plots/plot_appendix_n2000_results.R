
source(file.path("scripts", "plots", "plot_formal_simulation_results.R"))

formal_result_files <- c(
  binary = file.path(
    file.path("results", "simulation", "n2000"), "binary",
    "simulation_binary_n2000_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  continuous = file.path(
    file.path("results", "simulation", "n2000"), "continuous",
    "simulation_continuous_n2000_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  survival = file.path(
    file.path("results", "simulation", "n2000"), "survival",
    "simulation_survival_n2000_B50_M10_deg2_support0_05_span0_2.csv"
  )
)

figure_output_dir <- file.path("figures", "simulation", "n2000")

figure_style$scenario_text_size <- 3.8

outcome_titles <- c(
  binary = "Appendix results for binary outcome (n = 2000)",
  continuous = "Appendix results for continuous outcome (n = 2000)",
  survival = "Appendix results for time-to-event outcome (n = 2000)"
)

appendix_n2000_figure_outputs <- plot_available_formal_results()
