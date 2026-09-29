
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

result_file <- file.path(
  file.path("results", "real_data", "experiment1"),
  "real_survival_mgus2_n1338_splits100_B50_M10_deg2_support0_05_span0_2.csv"
)

figure_output_dir <- file.path("figures", "real_data", "experiment1")
expected_splits <- 100L

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

figure_style <- list(
  export_width = 10.6,
  export_height = 3.65,
  figure_title_size = 17,
  panel_title_size = 14,
  axis_y_text_size = 11,
  axis_x_text_size = 10.5,
  base_text_size = 12,
  x_text_angle = 45,
  dataset_strip_width = 0.95,
  prediction_panel_width = 2.85,
  term_panel_width = 2.15,
  degree1_panel_width = 2.30,
  dataset_border_linewidth = 0.85,
  panel_border_linewidth = 0.65
)

required_columns <- c(
  "Split", "Dataset", "method", "fit_ok", "metric_value",
  "n_rules_total", "n_rules_deg1", "prop_deg1"
)

read_real_results <- function(path) {
  if (!file.exists(path)) stop("Result file not found: ", path)

  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(dat))
  if (length(missing_columns) > 0L) {
    stop(
      "Result file is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat$fit_ok <- as.logical(dat$fit_ok)
  expected_methods <- sort(prediction_methods)
  observed_methods <- sort(unique(dat$method))

  if (!identical(observed_methods, expected_methods)) {
    stop(
      "Unexpected method set. Observed: ",
      paste(observed_methods, collapse = ", ")
    )
  }
  if (any(!dat$fit_ok)) {
    stop("At least one fitted model failed; inspect error_message before plotting.")
  }
  if (anyDuplicated(dat[c("Split", "method")])) {
    stop("Duplicate Split-method rows were found.")
  }

  counts <- table(dat$method)
  if (any(counts != expected_splits)) {
    stop(
      "Each method must contain exactly ", expected_splits,
      " successful splits. Observed: ",
      paste(names(counts), counts, sep = "=", collapse = "; ")
    )
  }
  if (!identical(sort(unique(dat$Split)), seq_len(expected_splits))) {
    stop("Split identifiers are not exactly 1 through ", expected_splits, ".")
  }

  needs_prop <- !is.finite(dat$prop_deg1) &
    is.finite(dat$n_rules_deg1) &
    is.finite(dat$n_rules_total) &
    dat$n_rules_total > 0
  dat$prop_deg1[needs_prop] <-
    dat$n_rules_deg1[needs_prop] / dat$n_rules_total[needs_prop]

  dat
}

make_dataset_strip <- function() {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "MGUS2\n(100\nrepeated\n7:3 splits)",
      fontface = "bold", size = 4.2, lineheight = 0.95
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      panel.background = element_rect(
        fill = "#F2F2F2", color = "black",
        linewidth = figure_style$dataset_border_linewidth
      ),
      plot.margin = margin(2, 2, 2, 2)
    )
}

make_box_panel <- function(
    dat,
    value_column,
    methods,
    title,
    limits = NULL,
    percent_axis = FALSE
) {
  z <- dat[
    dat$method %in% methods & is.finite(dat[[value_column]]),
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
        color = "black", linewidth = figure_style$panel_border_linewidth
      ),
      axis.text.y = element_text(
        size = figure_style$axis_y_text_size, color = "black"
      ),
      axis.text.x = element_text(
        size = figure_style$axis_x_text_size,
        color = "black",
        angle = figure_style$x_text_angle,
        hjust = 1,
        vjust = 1
      ),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.title = element_blank(),
      plot.title = element_text(
        size = figure_style$panel_title_size,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.margin = margin(2, 3, 2, 3)
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

make_real_figure <- function(dat) {
  term_values <- dat$n_rules_total[
    dat$method %in% complexity_methods & is.finite(dat$n_rules_total)
  ]
  term_upper <- max(10, ceiling(max(term_values) / 10) * 10)

  strip <- make_dataset_strip()
  prediction <- make_box_panel(
    dat,
    value_column = "metric_value",
    methods = prediction_methods,
    title = "C-index",
    limits = c(0.50, 1.00)
  )
  term_count <- make_box_panel(
    dat,
    value_column = "n_rules_total",
    methods = complexity_methods,
    title = "Final term count",
    limits = c(0, term_upper)
  )
  degree1_prop <- make_box_panel(
    dat,
    value_column = "prop_deg1",
    methods = complexity_methods,
    title = "Proportion of degree-1 terms",
    percent_axis = TRUE
  )

  strip + prediction + term_count + degree1_prop +
    plot_layout(widths = c(
      figure_style$dataset_strip_width,
      figure_style$prediction_panel_width,
      figure_style$term_panel_width,
      figure_style$degree1_panel_width
    )) +
    plot_annotation(
      title = "Results for MGUS2 across 100 repeated 7:3 train-test splits",
      theme = theme(
        plot.title = element_text(
          size = figure_style$figure_title_size,
          face = "bold",
          hjust = 0.5,
          margin = margin(b = 5)
        ),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}

save_real_figure <- function(plot, output_dir = figure_output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, "real_data_experiment1_mgus2.pdf")
  png_path <- file.path(output_dir, "real_data_experiment1_mgus2.png")

  ggsave(
    pdf_path,
    plot = plot,
    width = figure_style$export_width,
    height = figure_style$export_height,
    units = "in",
    bg = "white"
  )
  ggsave(
    png_path,
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

real_data <- read_real_results(result_file)
real_figure <- make_real_figure(real_data)
real_figure_outputs <- save_real_figure(real_figure)
print(real_figure_outputs)
