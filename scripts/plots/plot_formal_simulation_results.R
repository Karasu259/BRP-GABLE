
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

formal_result_files <- c(
  binary = file.path(
    file.path("results", "simulation", "n1000", "binary"),
    "simulation_binary_n1000_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  continuous = file.path(
    file.path("results", "simulation", "n1000", "continuous"),
    "simulation_continuous_n1000_B50_M10_deg2_support0_05_span0_2.csv"
  ),
  survival = file.path(
    file.path("results", "simulation", "n1000", "survival"),
    "simulation_survival_n1000_B50_M10_deg2_support0_05_span0_2.csv"
  )
)

figure_output_dir <- file.path("figures", "simulation", "main")
expected_simulations_per_scenario <- 100L

figure_style <- list(
  export_width = 10.6,
  export_height = 9.2,
  figure_title_size = 18,
  panel_title_size = 14,
  axis_y_text_size = 11,
  axis_x_text_size = 10.5,
  scenario_text_size = 4.2,
  base_text_size = 12,
  x_text_angle = 45,
  scenario_strip_width = 0.95,
  prediction_panel_width = 2.85,
  term_panel_width = 2.15,
  degree1_panel_width = 2.30,
  scenario_border_linewidth = 0.85,
  panel_border_linewidth = 0.65
)

prediction_methods <- c(
  "Linear",
  "GBM",
  "Random forest",
  "GABLE",
  "RuleFit",
  "BRP-GABLE rule-only",
  "BRP-GABLE full"
)

complexity_methods <- c(
  "GABLE",
  "RuleFit",
  "BRP-GABLE rule-only",
  "BRP-GABLE full"
)

method_labels <- c(
  "Linear" = "Linear",
  "GBM" = "GBM",
  "Random forest" = "RF",
  "GABLE" = "GABLE",
  "RuleFit" = "RuleFit",
  "BRP-GABLE rule-only" = "BRP-GABLE\nrules",
  "BRP-GABLE full" = "BRP-GABLE\nfull"
)

method_colors <- c(
  "Linear" = "#8C8C8C",
  "GBM" = "#E69F00",
  "RF" = "#CC79A7",
  "GABLE" = "#009E73",
  "RuleFit" = "#D55E00",
  "BRP-GABLE\nrules" = "#56B4E9",
  "BRP-GABLE\nfull" = "#0072B2"
)

scenario_labels <- c(
  "1" = "Scenario 1\n(purely\nrule-based)",
  "2" = "Scenario 2\n(purely\nlinear)",
  "3" = "Scenario 3\n(mixed\nlinear and\nrule-based)",
  "4" = "Scenario 4\n(smooth\nnonlinear\nmisspecified)"
)

outcome_titles <- c(
  binary = "Results for binary outcome across Scenarios 1 to 4",
  continuous = "Results for continuous outcome across Scenarios 1 to 4",
  survival = "Results for time-to-event outcome across Scenarios 1 to 4"
)

prediction_titles <- c(
  binary = "AUC",
  continuous = "MSE",
  survival = "C-index"
)

required_result_columns <- c(
  "Sim", "Outcome", "Scenario", "Scenario_name", "method",
  "fit_ok", "metric_value", "n_rules_total", "n_rules_deg1",
  "prop_deg1"
)

read_formal_results <- function(path, outcome) {
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_result_columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "Result file is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat <- dat[dat$Outcome == outcome, , drop = FALSE]
  dat$fit_ok <- as.logical(dat$fit_ok)
  dat <- dat[dat$fit_ok, , drop = FALSE]

  needs_prop <- !is.finite(dat$prop_deg1) &
    is.finite(dat$n_rules_deg1) &
    is.finite(dat$n_rules_total) &
    dat$n_rules_total > 0
  dat$prop_deg1[needs_prop] <-
    dat$n_rules_deg1[needs_prop] / dat$n_rules_total[needs_prop]

  observed <- aggregate(
    Sim ~ Scenario + method,
    data = dat,
    FUN = function(x) length(unique(x))
  )
  incomplete <- observed$Sim != expected_simulations_per_scenario
  if (any(incomplete)) {
    warning(
      "Some Scenario-method combinations do not contain exactly ",
      expected_simulations_per_scenario,
      " simulations: ",
      paste(
        paste0(
          "S", observed$Scenario[incomplete], "/",
          observed$method[incomplete], "=", observed$Sim[incomplete]
        ),
        collapse = "; "
      )
    )
  }

  dat
}

rounded_upper_limit <- function(x, step = 10) {
  x <- x[is.finite(x)]
  if (!length(x)) return(step)
  max(step, ceiling(max(x) / step) * step)
}

prediction_limits <- function(dat, outcome) {
  if (outcome %in% c("binary", "survival")) return(c(0.5, 1.0))
  vals <- dat$metric_value[
    dat$method %in% prediction_methods & is.finite(dat$metric_value)
  ]
  c(0, rounded_upper_limit(vals, step = 1))
}

make_scenario_strip <- function(scenario) {
  label <- unname(scenario_labels[as.character(scenario)])
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = label,
      fontface = "bold",
      size = figure_style$scenario_text_size,
      lineheight = 0.95
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      panel.background = element_rect(
        fill = "#F2F2F2",
        color = "black",
        linewidth = figure_style$scenario_border_linewidth
      ),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_box_panel <- function(
    dat,
    scenario,
    value_column,
    methods,
    title = NULL,
    limits = NULL,
    percent_axis = FALSE,
    show_x = FALSE
) {
  z <- dat[
    dat$Scenario == scenario &
      dat$method %in% methods &
      is.finite(dat[[value_column]]),
    , drop = FALSE
  ]

  z$method_label <- factor(
    unname(method_labels[z$method]),
    levels = unname(method_labels[methods])
  )

  p <- ggplot(
    z,
    aes(x = method_label, y = .data[[value_column]], fill = method_label)
  ) +
    geom_boxplot(
      width = 0.70,
      linewidth = 0.48,
      outlier.size = 0.85,
      outlier.alpha = 0.75,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = figure_style$base_text_size) +
    theme(
      panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        linewidth = figure_style$panel_border_linewidth
      ),
      axis.text.y = element_text(
        size = figure_style$axis_y_text_size,
        color = "black"
      ),
      axis.title = element_blank(),
      plot.title = element_text(
        size = figure_style$panel_title_size,
        face = "bold", hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.margin = margin(2, 3, 2, 3)
    )

  if (show_x) {
    p <- p + theme(
      axis.text.x = element_text(
        size = figure_style$axis_x_text_size,
        color = "black",
        angle = figure_style$x_text_angle,
        hjust = 1, vjust = 1
      ),
      axis.ticks.x = element_line(linewidth = 0.3)
    )
  } else {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  if (!is.null(limits)) {
    p <- p + coord_cartesian(ylim = limits)
  }
  if (percent_axis) {
    p <- p + scale_y_continuous(
      labels = label_percent(accuracy = 1),
      breaks = seq(0, 1, by = 0.25),
      limits = c(0, 1),
      expand = expansion(mult = c(0.01, 0.03))
    )
  }

  p
}

make_outcome_figure <- function(dat, outcome) {
  pred_limits <- prediction_limits(dat, outcome)
  term_values <- dat$n_rules_total[
    dat$method %in% complexity_methods & is.finite(dat$n_rules_total)
  ]
  term_limits <- c(0, rounded_upper_limit(term_values, step = 10))

  rows <- lapply(1:4, function(scenario) {
    show_x <- scenario == 4L
    top_row <- scenario == 1L

    strip <- make_scenario_strip(scenario)
    prediction <- make_box_panel(
      dat = dat,
      scenario = scenario,
      value_column = "metric_value",
      methods = prediction_methods,
      title = if (top_row) unname(prediction_titles[outcome]) else NULL,
      limits = pred_limits,
      show_x = show_x
    )
    term_count <- make_box_panel(
      dat = dat,
      scenario = scenario,
      value_column = "n_rules_total",
      methods = complexity_methods,
      title = if (top_row) "Final term count" else NULL,
      limits = term_limits,
      show_x = show_x
    )
    degree1_prop <- make_box_panel(
      dat = dat,
      scenario = scenario,
      value_column = "prop_deg1",
      methods = complexity_methods,
      title = if (top_row) "Proportion of degree-1 terms" else NULL,
      percent_axis = TRUE,
      show_x = show_x
    )

    strip + prediction + term_count + degree1_prop +
      plot_layout(widths = c(
        figure_style$scenario_strip_width,
        figure_style$prediction_panel_width,
        figure_style$term_panel_width,
        figure_style$degree1_panel_width
      ))
  })

  wrap_plots(rows, ncol = 1, heights = rep(1, 4)) +
    plot_annotation(
      title = unname(outcome_titles[outcome]),
      theme = theme(
        plot.title = element_text(
          size = figure_style$figure_title_size,
          face = "bold", hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}

save_outcome_figure <- function(plot, outcome, output_dir = figure_output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, paste0("simulation_", outcome, ".pdf"))
  png_path <- file.path(output_dir, paste0("simulation_", outcome, ".png"))

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

plot_available_formal_results <- function() {
  outputs <- list()
  for (outcome in names(formal_result_files)) {
    path <- formal_result_files[[outcome]]
    if (!file.exists(path)) {
      message("Skipping ", outcome, ": result file is not available yet.")
      next
    }

    message("Plotting ", outcome, " results from ", path)
    dat <- read_formal_results(path, outcome)
    fig <- make_outcome_figure(dat, outcome)
    outputs[[outcome]] <- save_outcome_figure(fig, outcome)
  }

  if (!length(outputs)) {
    warning("No formal result files were available for plotting.")
  }
  invisible(outputs)
}

if (sys.nframe() == 0L) {
  formal_figure_outputs <- plot_available_formal_results()
}
