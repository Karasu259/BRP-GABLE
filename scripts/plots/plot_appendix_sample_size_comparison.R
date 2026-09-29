
source(file.path("scripts", "plots", "plot_formal_simulation_results.R"))

sample_size_result_files <- list(
  binary = c(
    `500` = file.path(
      file.path("results", "simulation", "n500"), "binary",
      "simulation_binary_n500_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `1000` = file.path(
      file.path("results", "simulation", "n1000", "binary"),
      "simulation_binary_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `2000` = file.path(
      file.path("results", "simulation", "n2000"), "binary",
      "simulation_binary_n2000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  ),
  continuous = c(
    `500` = file.path(
      file.path("results", "simulation", "n500"), "continuous",
      "simulation_continuous_n500_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `1000` = file.path(
      file.path("results", "simulation", "n1000", "continuous"),
      "simulation_continuous_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `2000` = file.path(
      file.path("results", "simulation", "n2000"), "continuous",
      "simulation_continuous_n2000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  ),
  survival = c(
    `500` = file.path(
      file.path("results", "simulation", "n500"), "survival",
      "simulation_survival_n500_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `1000` = file.path(
      file.path("results", "simulation", "n1000", "survival"),
      "simulation_survival_n1000_B50_M10_deg2_support0_05_span0_2.csv"
    ),
    `2000` = file.path(
      file.path("results", "simulation", "n2000"), "survival",
      "simulation_survival_n2000_B50_M10_deg2_support0_05_span0_2.csv"
    )
  )
)

sample_size_titles <- c(
  binary = "Sample-size comparison for binary outcome",
  continuous = "Sample-size comparison for continuous outcome",
  survival = "Sample-size comparison for time-to-event outcome"
)

sample_size_figure_style <- figure_style
sample_size_figure_style$export_width <- 11.8
sample_size_figure_style$export_height <- 9.4
sample_size_figure_style$axis_x_text_size <- 11
sample_size_figure_style$scenario_text_size <- 3.8

method_shapes <- c(
  "Linear" = 16,
  "GBM" = 17,
  "RF" = 15,
  "GABLE" = 18,
  "RuleFit" = 0,
  "BRP-GABLE\nrules" = 1,
  "BRP-GABLE\nfull" = 2
)

read_sample_size_results <- function(outcome) {
  paths <- sample_size_result_files[[outcome]]
  missing_files <- paths[!file.exists(paths)]
  if (length(missing_files)) {
    stop("Missing result file(s): ", paste(missing_files, collapse = ", "))
  }

  pieces <- lapply(names(paths), function(n_value) {
    dat <- read_formal_results(paths[[n_value]], outcome)
    dat$sample_size <- as.numeric(n_value)
    dat
  })
  do.call(rbind, pieces)
}

summarize_mean_ci <- function(dat, value_column, methods) {
  z <- dat[
    dat$method %in% methods & is.finite(dat[[value_column]]),
    c("Scenario", "method", "sample_size", value_column),
    drop = FALSE
  ]
  names(z)[4] <- "value"

  keys <- interaction(z$Scenario, z$method, z$sample_size, drop = TRUE)
  groups <- split(z, keys)
  summaries <- lapply(groups, function(g) {
    n_obs <- nrow(g)
    estimate <- mean(g$value)
    standard_error <- if (n_obs > 1L) stats::sd(g$value) / sqrt(n_obs) else NA_real_
    data.frame(
      Scenario = g$Scenario[1],
      method = g$method[1],
      sample_size = g$sample_size[1],
      n = n_obs,
      estimate = estimate,
      lower = estimate - 1.96 * standard_error,
      upper = estimate + 1.96 * standard_error,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out$sample_size_label <- factor(
    out$sample_size,
    levels = c(500, 1000, 2000),
    labels = c("500", "1,000", "2,000")
  )
  out$method_label <- factor(
    unname(method_labels[out$method]),
    levels = unname(method_labels[prediction_methods])
  )
  out
}

adaptive_prediction_limits <- function(summary_data, scenario, outcome) {
  z <- summary_data[
    summary_data$Scenario == scenario,
    c("estimate", "lower", "upper"),
    drop = FALSE
  ]
  values <- unlist(z, use.names = FALSE)
  values <- values[is.finite(values)]
  if (!length(values)) return(NULL)

  observed_range <- range(values)
  observed_span <- diff(observed_range)

  if (outcome %in% c("binary", "survival")) {
    padding <- max(0.08 * observed_span, 0.005)
    limits <- c(
      max(0, observed_range[1] - padding),
      min(1, observed_range[2] + padding)
    )
    step <- if (diff(limits) <= 0.14) 0.02 else 0.05
    limits <- c(
      max(0, floor(limits[1] / step) * step),
      min(1, ceiling(limits[2] / step) * step)
    )
  } else {
    padding <- max(0.08 * observed_span, 0.05)
    limits <- c(
      max(0, observed_range[1] - padding),
      observed_range[2] + padding
    )
  }

  if (!all(is.finite(limits)) || diff(limits) <= 0) return(NULL)
  limits
}

make_mean_ci_panel <- function(
    summary_data,
    scenario,
    title = NULL,
    limits = NULL,
    percent_axis = FALSE,
    show_x_title = FALSE
) {
  z <- summary_data[summary_data$Scenario == scenario, , drop = FALSE]

  p <- ggplot(
    z,
    aes(
      x = sample_size_label,
      y = estimate,
      color = method_label,
      shape = method_label,
      group = method_label
    )
  ) +
    geom_line(linewidth = 0.62, alpha = 0.9, show.legend = TRUE) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.12,
      linewidth = 0.48,
      show.legend = FALSE
    ) +
    geom_point(size = 2.25, stroke = 0.75, show.legend = TRUE) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = method_shapes, drop = FALSE) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    labs(
      title = title,
      x = if (show_x_title) "Sample size" else NULL,
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    theme_bw(base_size = sample_size_figure_style$base_text_size) +
    theme(
      panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        linewidth = sample_size_figure_style$panel_border_linewidth
      ),
      axis.text.y = element_text(
        size = sample_size_figure_style$axis_y_text_size,
        color = "black"
      ),
      axis.text.x = element_text(
        size = sample_size_figure_style$axis_x_text_size,
        color = "black"
      ),
      axis.title.x = element_text(
        size = sample_size_figure_style$axis_x_text_size,
        color = "black",
        margin = margin(t = 3)
      ),
      plot.title = element_text(
        size = sample_size_figure_style$panel_title_size,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.margin = margin(2, 3, 2, 3),
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.key.width = grid::unit(1.15, "lines"),
      legend.spacing.x = grid::unit(0.15, "lines")
    )

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

make_sample_size_figure <- function(
    dat,
    outcome,
    adaptive_prediction_axis = FALSE
) {
  prediction_summary <- summarize_mean_ci(dat, "metric_value", prediction_methods)
  term_summary <- summarize_mean_ci(dat, "n_rules_total", complexity_methods)
  degree1_summary <- summarize_mean_ci(dat, "prop_deg1", complexity_methods)

  fixed_prediction_limits <- prediction_limits(dat, outcome)
  term_limits <- c(
    0,
    rounded_upper_limit(
      c(term_summary$lower, term_summary$upper),
      step = 10
    )
  )

  rows <- lapply(1:4, function(scenario) {
    top_row <- scenario == 1L
    bottom_row <- scenario == 4L
    pred_limits <- if (isTRUE(adaptive_prediction_axis)) {
      adaptive_prediction_limits(prediction_summary, scenario, outcome)
    } else {
      fixed_prediction_limits
    }

    strip <- make_scenario_strip(scenario)
    prediction <- make_mean_ci_panel(
      prediction_summary,
      scenario,
      title = if (top_row) unname(prediction_titles[outcome]) else NULL,
      limits = pred_limits,
      show_x_title = bottom_row
    )
    term_count <- make_mean_ci_panel(
      term_summary,
      scenario,
      title = if (top_row) "Final term count" else NULL,
      limits = term_limits,
      show_x_title = bottom_row
    )
    degree1_prop <- make_mean_ci_panel(
      degree1_summary,
      scenario,
      title = if (top_row) "Proportion of degree-1 terms" else NULL,
      percent_axis = TRUE,
      show_x_title = bottom_row
    )

    strip + prediction + term_count + degree1_prop +
      plot_layout(widths = c(
        sample_size_figure_style$scenario_strip_width,
        sample_size_figure_style$prediction_panel_width,
        sample_size_figure_style$term_panel_width,
        sample_size_figure_style$degree1_panel_width
      ))
  })

  wrap_plots(rows, ncol = 1, heights = rep(1, 4), guides = "collect") +
    plot_annotation(
      title = unname(sample_size_titles[outcome]),
      theme = theme(
        plot.title = element_text(
          size = sample_size_figure_style$figure_title_size,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom"
      )
    ) &
    theme(legend.position = "bottom")
}

save_sample_size_figure <- function(plot, outcome, file_suffix = "") {
  output_dir <- file.path("figures", "simulation", "sample_size")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(
    output_dir,
    paste0("sample_size_comparison_", outcome, file_suffix, ".pdf")
  )
  png_path <- file.path(
    output_dir,
    paste0("sample_size_comparison_", outcome, file_suffix, ".png")
  )

  ggsave(
    pdf_path,
    plot = plot,
    width = sample_size_figure_style$export_width,
    height = sample_size_figure_style$export_height,
    units = "in",
    bg = "white"
  )
  ggsave(
    png_path,
    plot = plot,
    device = grDevices::png,
    width = sample_size_figure_style$export_width,
    height = sample_size_figure_style$export_height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  invisible(c(pdf = pdf_path, png = png_path))
}

plot_sample_size_comparison <- function(
    outcomes = c("binary", "continuous", "survival"),
    adaptive_prediction_axis = FALSE,
    file_suffix = ""
) {
  outputs <- list()
  for (outcome in outcomes) {
    message("Plotting sample-size comparison for ", outcome)
    dat <- read_sample_size_results(outcome)
    fig <- make_sample_size_figure(
      dat,
      outcome,
      adaptive_prediction_axis = adaptive_prediction_axis
    )
    outputs[[outcome]] <- save_sample_size_figure(
      fig,
      outcome,
      file_suffix = file_suffix
    )
  }
  invisible(outputs)
}

if (sys.nframe() == 0L) {
  sample_size_figure_outputs <- plot_sample_size_comparison(
    adaptive_prediction_axis = TRUE,
    file_suffix = "_adaptive"
  )
}
