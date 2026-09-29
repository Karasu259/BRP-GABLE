
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

result_files <- list.files(
  file.path("results", "real_data", "experiment2"),
  pattern = paste0(
    "^brp_gable_mgus2_full_n[0-9]+_B50_M10_deg2_",
    "support0_05_span0_2[.]rds$"
  ),
  full.names = TRUE
)
figure_output_dir <- file.path("figures", "real_data", "experiment2")

if (!length(result_files)) {
  stop(
    "No full-data Experiment 2 result was found in ",
    "results/real_data/experiment2.",
    "\nRun run_real_survival_full.R first."
  )
}
result_file <- result_files[which.max(file.info(result_files)$mtime)]
bundle <- readRDS(result_file)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.65),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

save_figure <- function(plot, name, width, height) {
  pdf_path <- file.path(figure_output_dir, paste0(name, ".pdf"))
  png_path <- file.path(figure_output_dir, paste0(name, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in", bg = "white")
  ggsave(
    png_path, plot,
    width = width, height = height, units = "in", dpi = 300, bg = "white"
  )
  invisible(c(pdf = pdf_path, png = png_path))
}

wrap_rule_text <- function(x, width = 34L) {
  vapply(
    x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

make_support_effect_plot <- function(bundle) {
  z <- bundle$selected_rules
  if (!nrow(z)) return(NULL)
  z <- z[order(-z$importance, z$rule_id), , drop = FALSE]
  z$degree_label <- factor(
    paste0("Degree ", z$degree),
    levels = paste0("Degree ", sort(unique(z$degree)))
  )
  importance_levels <- c("0-10", "10-20", "20-30", "Above 30")
  z$importance_band <- cut(
    z$importance_relative,
    breaks = c(-Inf, 10, 20, 30, Inf),
    labels = importance_levels,
    right = FALSE
  )
  z$rule_label <- paste0(
    z$rule_id, ": ",
    gsub(" AND ", "\n& ", z$term, fixed = TRUE)
  )

  label_offset_x <- c(-0.18, -0.18, -0.18, 0.08, 0.08, -0.18, 0.08)
  label_multiplier_y <- c(0.985, 0.985, 1.008, 0.970, 1.020, 1.012, 1.012)
  z$label_x <- pmin(
    pmax(z$support + rep(label_offset_x, length.out = nrow(z)), 0.02),
    0.80
  )
  z$label_y <- z$hazard_ratio * rep(
    label_multiplier_y, length.out = nrow(z)
  )
  importance_legend_data <- data.frame(
    support = NA_real_,
    hazard_ratio = NA_real_,
    importance_band = factor(importance_levels, levels = importance_levels)
  )

  colors <- c("Degree 1" = "#009E73", "Degree 2" = "#0072B2")
  ggplot(
    z,
    aes(
      x = support,
      y = hazard_ratio,
      size = importance_band,
      color = degree_label
    )
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
    geom_segment(
      data = z,
      aes(
        x = support, y = hazard_ratio,
        xend = label_x, yend = label_y
      ),
      inherit.aes = FALSE,
      color = "#777777", linewidth = 0.35,
      show.legend = FALSE
    ) +
    geom_point(alpha = 0.85) +
    geom_point(
      data = importance_legend_data,
      aes(x = support, y = hazard_ratio, size = importance_band),
      inherit.aes = FALSE,
      color = "black", na.rm = TRUE, show.legend = TRUE
    ) +
    geom_label(
      aes(x = label_x, y = label_y, label = rule_label),
      color = "black",
      fill = "white", size = 3.2, lineheight = 0.95,
      linewidth = 0.18, hjust = 0,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = label_percent(accuracy = 1),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      limits = c(0.83, 1.03),
      labels = label_number(accuracy = 0.01)
    ) +
    scale_color_manual(values = colors, drop = FALSE, name = "Rule degree") +
    scale_size_manual(
      values = c("0-10" = 3.5, "10-20" = 5.5,
                 "20-30" = 7.5, "Above 30" = 9.5),
      limits = importance_levels,
      drop = FALSE,
      name = "Relative importance"
    ) +
    labs(
      title = "Rule support-effect map of selected BRP-GABLE rules",
      x = "Rule support in the full MGUS2 data",
      y = "Model-based hazard ratio"
    ) +
    base_theme +
    guides(
      color = guide_legend(
        override.aes = list(color = "#0072B2", size = 4)
      )
    )
}

make_generation_frequency_plot <- function(bundle) {
  z <- bundle$term_table[
    bundle$term_table$term_type == "rule",
    , drop = FALSE
  ]
  ord <- order(-z$count_bootstrap, -z$count, z$rule_key)
  z <- z[ord, , drop = FALSE]
  z$rank <- seq_len(nrow(z))
  z$selection_label <- ifelse(z$selected, "Selected by lasso", "Not selected")
  selected_index <- which(z$selected)
  z$label_y <- z$count_bootstrap
  label_level <- ave(
    z$rank[selected_index],
    z$count_bootstrap[selected_index],
    FUN = seq_along
  )
  z$label_y[selected_index] <- pmin(
    z$count_bootstrap[selected_index] + 0.8 + 1.2 * (label_level - 1L),
    bundle$config$B - 0.3
  )

  ggplot(z, aes(x = rank, y = count_bootstrap)) +
    geom_point(
      data = z[!z$selected, , drop = FALSE],
      color = "#8C8C8C", alpha = 0.55, size = 1.5
    ) +
    geom_point(
      data = z[z$selected, , drop = FALSE],
      color = "#0072B2", size = 3.2
    ) +
    geom_segment(
      data = z[z$selected, , drop = FALSE],
      aes(xend = rank, yend = label_y - 0.2),
      color = "#0072B2", linewidth = 0.35
    ) +
    geom_text(
      data = z[z$selected, , drop = FALSE],
      aes(y = label_y, label = rule_id),
      color = "#0072B2", fontface = "bold",
      size = 3.2, check_overlap = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, bundle$config$B),
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    labs(
      title = "Bootstrap generation frequency before lasso",
      x = "Candidate-rule rank",
      y = paste0("Number of bootstrap fits (out of ", bundle$config$B, ")")
    ) +
    base_theme
}

rule_importance_order <- function(bundle) {
  z <- bundle$selected_rules
  if (!nrow(z)) return(character(0))
  z$rule_id[order(-z$importance, z$rule_id)]
}

make_activation_heatmap <- function(bundle, rule_order) {
  A <- bundle$selected_rule_activation
  if (!ncol(A)) return(NULL)
  A <- A[, rule_order, drop = FALSE]

  z <- data.frame(
    patient_index = rep(seq_len(nrow(A)), times = ncol(A)),
    rule_id = rep(colnames(A), each = nrow(A)),
    activated = factor(
      as.vector(A), levels = c(0, 1),
      labels = c("Not activated", "Activated")
    ),
    stringsAsFactors = FALSE
  )
  z$rule_id <- factor(z$rule_id, levels = rule_order)

  ggplot(z, aes(x = rule_id, y = patient_index, fill = activated)) +
    geom_raster() +
    geom_vline(
      xintercept = seq(1.5, ncol(A) - 0.5, by = 1),
      color = "white", linewidth = 0.35
    ) +
    scale_fill_manual(
      values = c("Not activated" = "#E6E6E6", "Activated" = "#0072B2"),
      name = NULL
    ) +
    scale_y_continuous(
      breaks = pretty_breaks(n = 4),
      expand = c(0, 0)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(
      title = "Patient-rule activation heatmap",
      x = "Selected rule",
      y = "Patient index"
    ) +
    base_theme +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.key = element_rect(color = "#777777", linewidth = 0.4)
    )
}

make_similarity_heatmap <- function(bundle, rule_order) {
  sim <- bundle$jaccard_similarity
  if (!length(sim)) return(NULL)
  sim <- sim[rule_order, rule_order, drop = FALSE]
  z <- as.data.frame(as.table(sim), stringsAsFactors = FALSE)
  names(z) <- c("rule_row", "rule_col", "jaccard")
  z$rule_row <- factor(z$rule_row, levels = rev(rule_order))
  z$rule_col <- factor(z$rule_col, levels = rule_order)

  p <- ggplot(z, aes(x = rule_col, y = rule_row, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.45) +
    scale_fill_gradient(
      low = "#FFF7BC", high = "#D7301F",
      limits = c(0, 1), name = "Jaccard\nsimilarity"
    ) +
    coord_equal() +
    labs(
      title = "Rule activation-set similarity",
      x = "Selected rule",
      y = "Selected rule"
    ) +
    base_theme +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  if (length(rule_order) <= 15L) {
    p <- p + geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3)
  }
  p
}

variable_axis_labels <- c(
  age = "Age",
  sex = "Sex",
  dxyr = "Year of diagnosis",
  hgb = "Hemoglobin",
  creat = "Serum creatinine",
  mspike = "M-spike"
)

variable_label <- function(v) {
  if (v %in% names(variable_axis_labels)) {
    unname(variable_axis_labels[v])
  } else {
    v
  }
}

make_variable_effect_plots <- function(bundle) {
  plots <- list()
  if (nrow(bundle$ale_curves)) {
    for (v in unique(bundle$ale_curves$variable)) {
      z <- bundle$ale_curves[bundle$ale_curves$variable == v, , drop = FALSE]
      display_limits <- as.numeric(stats::quantile(
        bundle$data$X[[v]], c(0.01, 0.99), na.rm = TRUE, names = FALSE
      ))
      plots[[length(plots) + 1L]] <- ggplot(z, aes(x = x, y = ale)) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
        geom_step(direction = "hv", color = "#0072B2", linewidth = 0.9) +
        geom_point(color = "#0072B2", size = 1.6) +
        coord_cartesian(xlim = display_limits) +
        labs(
          title = variable_label(v),
          x = variable_label(v),
          y = "Centered log relative hazard"
        ) +
        base_theme
    }
  }

  if (nrow(bundle$binary_effects)) {
    for (v in unique(bundle$binary_effects$variable)) {
      z <- bundle$binary_effects[
        bundle$binary_effects$variable == v,
        , drop = FALSE
      ]
      level_labels <- as.character(z$level)
      if (v == "sex" && identical(as.numeric(z$level), c(0, 1))) {
        level_labels <- c("Female", "Male")
      }
      z$level_label <- factor(level_labels, levels = level_labels)
      plots[[length(plots) + 1L]] <- ggplot(
        z, aes(x = level_label, y = centered_effect, group = 1)
      ) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
        geom_line(color = "#0072B2", linewidth = 0.8) +
        geom_point(color = "#0072B2", size = 3) +
        labs(
          title = variable_label(v),
          x = variable_label(v),
          y = "Centered log relative hazard"
        ) +
        base_theme
    }
  }

  if (!length(plots)) return(NULL)
  wrap_plots(plots, ncol = if (length(plots) <= 3L) length(plots) else 2L) +
    plot_annotation(
      title = "Model-based variable-effect summaries",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
      )
    )
}

make_km_panel <- function(bundle, rule_id) {
  curves <- bundle$km_curves[bundle$km_curves$rule_id == rule_id, , drop = FALSE]
  lookup <- bundle$representative_rules[
    bundle$representative_rules$rule_id == rule_id,
    , drop = FALSE
  ]
  colors <- c("Not satisfied" = "#8C8C8C", "Satisfied" = "#0072B2")

  km_plot <- ggplot(
    curves,
    aes(x = time, y = survival, color = group, group = group)
  ) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = colors, name = NULL) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = rule_id,
      subtitle = wrap_rule_text(lookup$term, width = 30L),
      x = "Follow-up time (months)",
      y = "Overall survival"
    ) +
    base_theme +
    theme(legend.position = "bottom")

  km_plot
}

make_km_figure <- function(bundle) {
  if (!nrow(bundle$representative_rules)) return(NULL)
  panels <- lapply(
    bundle$representative_rules$rule_id,
    function(rid) make_km_panel(bundle, rid)
  )
  wrap_plots(panels, ncol = length(panels), guides = "collect") +
    plot_annotation(
      title = "Kaplan-Meier curves for representative selected rules",
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
      )
    ) &
    theme(legend.position = "bottom")
}

outputs <- list()

support_effect <- make_support_effect_plot(bundle)
if (!is.null(support_effect)) {
  outputs$support_effect <- save_figure(
    support_effect, "rule_support_effect_map", 8.5, 6.2
  )
}

generation_frequency <- make_generation_frequency_plot(bundle)
outputs$generation_frequency <- save_figure(
  generation_frequency, "bootstrap_rule_generation_frequency", 8.5, 5.4
)

rule_order <- rule_importance_order(bundle)
activation_heatmap <- make_activation_heatmap(bundle, rule_order)
similarity_heatmap <- make_similarity_heatmap(bundle, rule_order)
if (!is.null(activation_heatmap)) {
  outputs$activation_heatmap <- save_figure(
    activation_heatmap, "patient_rule_activation_heatmap", 8.5, 6.5
  )
}
if (!is.null(similarity_heatmap)) {
  outputs$similarity_heatmap <- save_figure(
    similarity_heatmap, "rule_similarity_heatmap", 7.2, 6.4
  )
}
if (!is.null(activation_heatmap) && !is.null(similarity_heatmap)) {
  combined <- activation_heatmap + similarity_heatmap +
    plot_layout(widths = c(1.25, 1))
  outputs$activation_similarity <- save_figure(
    combined, "rule_activation_and_similarity", 12, 6.5
  )
}

variable_effects <- make_variable_effect_plots(bundle)
if (!is.null(variable_effects)) {
  n_panels <- length(unique(bundle$ale_curves$variable)) +
    length(unique(bundle$binary_effects$variable))
  outputs$variable_effects <- save_figure(
    variable_effects,
    "variable_effect_summaries",
    width = if (n_panels <= 3L) 10.5 else 9,
    height = if (n_panels <= 3L) 4.3 else 7.5
  )
}

km_figure <- make_km_figure(bundle)
if (!is.null(km_figure)) {
  outputs$km <- save_figure(
    km_figure, "representative_rule_km_curves", 11, 4.6
  )
}

print(outputs)
