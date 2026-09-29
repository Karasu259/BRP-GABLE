
user_library <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_library)) {
  user_library <- path.expand(user_library)
  if (dir.exists(user_library) && !user_library %in% .libPaths()) {
    .libPaths(c(user_library, .libPaths()))
  }
}

source(file.path("R", "00_utils.R"))
source(file.path("R", "01_rule_generation.R"))
source(file.path("R", "02_models.R"))
source(file.path("R", "05_real_data_engine.R"))
source(file.path("R", "10_real_data_explanation_engine.R"))

assert_packages(c(
  "survival", "glmnet", "pre", "gbm", "randomForestSRC",
  "ggplot2", "patchwork", "scales"
))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

result_dir <- file.path("results", "real_data", "experiment2")
figure_dir <- file.path("figures", "real_data", "appendix_ale")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

brp_files <- list.files(
  result_dir,
  pattern = paste0(
    "^brp_gable_mgus2_full_n[0-9]+_B50_M10_deg2_",
    "support0_05_span0_2[.]rds$"
  ),
  full.names = TRUE
)
if (!length(brp_files)) {
  stop("No saved full-data BRP-GABLE Experiment 2 fit was found.")
}
brp_file <- brp_files[which.max(file.info(brp_files)$mtime)]
brp_bundle <- readRDS(brp_file)

X <- as.data.frame(brp_bundle$data$X, check.names = FALSE)
time <- brp_bundle$data$time
status <- brp_bundle$data$status

continuous_variables <- names(X)[
  vapply(X, function(z) length(unique(z)) > 2L, logical(1))
]
variable_order <- c("age", "dxyr", "hgb", "creat", "mspike")
missing_variables <- setdiff(variable_order, continuous_variables)
if (length(missing_variables)) {
  stop(
    "Expected continuous MGUS2 predictors were not found: ",
    paste(missing_variables, collapse = ", ")
  )
}

seed_set <- make_seed_set(brp_bundle$config$base_seed, 0L)
rulefit_cache <- file.path(
  result_dir,
  paste0(
    "rulefit_full_mgus2_default_n", nrow(X),
    "_seed", seed_set$rulefit, ".rds"
  )
)

if (file.exists(rulefit_cache)) {
  rulefit_saved <- readRDS(rulefit_cache)
  rulefit_fit <- rulefit_saved$fit
  message("Using cached full-data RuleFit fit: ", rulefit_cache)
} else {
  message("Fitting package-default RuleFit to the full MGUS2 data...")
  rulefit_result <- fit_rulefit_default(
    X_train = X,
    outcome_type = "survival",
    time_train = time,
    status_train = status,
    X_test = X,
    time_test = time,
    status_test = status,
    seed = seed_set$rulefit
  )
  rulefit_fit <- rulefit_result$fit
  saveRDS(
    list(
      fit = rulefit_fit,
      seed = seed_set$rulefit,
      n = nrow(X),
      predictors = names(X),
      pre_version = as.character(utils::packageVersion("pre"))
    ),
    rulefit_cache
  )
}

gbm_cache <- file.path(
  result_dir,
  paste0(
    "gbm_full_mgus2_default_n", nrow(X),
    "_seed", seed_set$gbm, ".rds"
  )
)

if (file.exists(gbm_cache)) {
  gbm_saved <- readRDS(gbm_cache)
  gbm_fit <- gbm_saved$fit
  message("Using cached full-data GBM fit: ", gbm_cache)
} else {
  message("Fitting package-default GBM to the full MGUS2 data...")
  gbm_result <- fit_gbm_baseline(
    X_train = X,
    outcome_type = "survival",
    time_train = time,
    status_train = status,
    X_test = X,
    time_test = time,
    status_test = status,
    seed = seed_set$gbm
  )
  gbm_fit <- gbm_result$fit
  saveRDS(
    list(
      fit = gbm_fit,
      seed = seed_set$gbm,
      n = nrow(X),
      predictors = names(X),
      gbm_version = as.character(utils::packageVersion("gbm"))
    ),
    gbm_cache
  )
}

rf_cache <- file.path(
  result_dir,
  paste0(
    "random_forest_full_mgus2_default_n", nrow(X),
    "_seed", seed_set$rf, ".rds"
  )
)

if (file.exists(rf_cache)) {
  rf_saved <- readRDS(rf_cache)
  rf_fit <- rf_saved$fit
  message("Using cached full-data random survival forest fit: ", rf_cache)
} else {
  message("Fitting package-default random survival forest to the full MGUS2 data...")
  rf_result <- fit_rf_baseline(
    X_train = X,
    outcome_type = "survival",
    time_train = time,
    status_train = status,
    X_test = X,
    time_test = time,
    status_test = status,
    seed = seed_set$rf
  )
  rf_fit <- rf_result$fit
  saveRDS(
    list(
      fit = rf_fit,
      seed = seed_set$rf,
      n = nrow(X),
      predictors = names(X),
      randomForestSRC_version = as.character(
        utils::packageVersion("randomForestSRC")
      )
    ),
    rf_cache
  )
}

predict_brp_raw <- function(X_new) {
  predict_brp_from_raw_X(brp_bundle, X_new, type = "link")
}

predict_rulefit_raw <- function(X_new) {
  as.numeric(stats::predict(
    rulefit_fit,
    newdata = as.data.frame(X_new, check.names = FALSE),
    type = "link",
    penalty.par.val = "lambda.1se"
  ))
}

predict_gbm_raw <- function(X_new) {
  as.numeric(stats::predict(
    gbm_fit,
    newdata = as.data.frame(X_new, check.names = FALSE),
    n.trees = gbm_fit$n.trees,
    type = "link"
  ))
}

predict_rf_raw <- function(X_new) {
  prediction <- stats::predict(
    rf_fit,
    newdata = as.data.frame(X_new, check.names = FALSE)
  )
  as.numeric(prediction$predicted)
}

raw_prediction_functions <- list(
  "BRP-GABLE" = predict_brp_raw,
  "RuleFit" = predict_rulefit_raw,
  "GBM" = predict_gbm_raw,
  "RF" = predict_rf_raw
)

prediction_sd <- vapply(
  raw_prediction_functions,
  function(predict_fun) stats::sd(predict_fun(X)),
  numeric(1)
)
if (any(!is.finite(prediction_sd) | prediction_sd <= 0)) {
  stop(
    "Non-positive or non-finite fitted-score standard deviation for: ",
    paste(names(prediction_sd)[!is.finite(prediction_sd) | prediction_sd <= 0],
          collapse = ", ")
  )
}

standardized_prediction_functions <- Map(
  function(predict_fun, scale_value) {
    force(predict_fun)
    force(scale_value)
    function(X_new) predict_fun(X_new) / scale_value
  },
  raw_prediction_functions,
  prediction_sd
)

make_quantile_breaks <- function(x, n_bins = 20L) {
  breaks <- unique(as.numeric(stats::quantile(
    x,
    probs = seq(0, 1, length.out = as.integer(n_bins) + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )))
  if (length(breaks) < 3L) {
    stop("Too few unique quantile boundaries for ALE calculation.")
  }
  breaks
}

compute_ale_with_breaks <- function(X, variable, breaks, predict_link, method) {
  x <- as.numeric(X[[variable]])
  n_intervals <- length(breaks) - 1L
  interval <- findInterval(
    x,
    breaks,
    all.inside = TRUE,
    rightmost.closed = TRUE
  )
  interval <- pmin(pmax(interval, 1L), n_intervals)

  local_effect <- numeric(n_intervals)
  n_in_interval <- integer(n_intervals)
  for (ell in seq_len(n_intervals)) {
    ids <- which(interval == ell)
    n_in_interval[ell] <- length(ids)
    if (!length(ids)) next

    X_low <- X[ids, , drop = FALSE]
    X_high <- X[ids, , drop = FALSE]
    X_low[[variable]] <- breaks[ell]
    X_high[[variable]] <- breaks[ell + 1L]
    local_effect[ell] <- mean(
      predict_link(X_high) - predict_link(X_low)
    )
  }

  accumulated <- cumsum(local_effect)
  center <- stats::weighted.mean(accumulated, n_in_interval)
  data.frame(
    variable = variable,
    method = method,
    x = breaks,
    ale = c(-center, accumulated - center),
    interval_lower = c(NA_real_, breaks[-length(breaks)]),
    interval_upper = c(NA_real_, breaks[-1L]),
    n_in_interval = c(NA_integer_, n_in_interval),
    stringsAsFactors = FALSE
  )
}

ale_rows <- list()
row_id <- 1L
for (variable in variable_order) {
  shared_breaks <- make_quantile_breaks(X[[variable]], n_bins = 20L)
  for (method in names(standardized_prediction_functions)) {
    ale_rows[[row_id]] <- compute_ale_with_breaks(
      X,
      variable,
      shared_breaks,
      standardized_prediction_functions[[method]],
      method
    )
    row_id <- row_id + 1L
  }
}
ale_comparison <- do.call(rbind, ale_rows)
ale_comparison$method <- factor(
  ale_comparison$method,
  levels = c("BRP-GABLE", "RuleFit", "GBM", "RF")
)

data_path <- file.path(
  result_dir,
  "ale_comparison_brp_gable_rulefit_gbm_rf_mgus2_full.csv"
)
utils::write.csv(
  ale_comparison,
  data_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

scale_path <- file.path(
  result_dir,
  "ale_prediction_score_scales_mgus2_full.csv"
)
utils::write.csv(
  data.frame(
    method = names(prediction_sd),
    fitted_score_sd = as.numeric(prediction_sd),
    stringsAsFactors = FALSE
  ),
  scale_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

variable_labels <- c(
  age = "Age",
  dxyr = "Year of diagnosis",
  hgb = "Hemoglobin",
  creat = "Serum creatinine",
  mspike = "Mspike"
)
method_colors <- c(
  "BRP-GABLE" = "#0072B2",
  "RuleFit" = "#D55E00",
  "GBM" = "#E69F00",
  "RF" = "#CC79A7"
)
method_linetypes <- c(
  "BRP-GABLE" = "solid",
  "RuleFit" = "dashed",
  "GBM" = "dotdash",
  "RF" = "twodash"
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "#E3E3E3", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.65),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 11),
    plot.background = element_rect(fill = "white", color = NA)
  )

panels <- lapply(variable_order, function(variable) {
  z <- ale_comparison[
    ale_comparison$variable == variable,
    , drop = FALSE
  ]
  display_limits <- as.numeric(stats::quantile(
    X[[variable]],
    probs = c(0.01, 0.99),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))

  ggplot(
    z,
    aes(x = x, y = ale, color = method, linetype = method, group = method)
  ) +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.45) +
    geom_step(direction = "hv", linewidth = 0.95) +
    geom_point(size = 1.35, alpha = 0.85) +
    coord_cartesian(xlim = display_limits) +
    scale_color_manual(values = method_colors, name = NULL) +
    scale_linetype_manual(values = method_linetypes, name = NULL) +
    labs(
      title = unname(variable_labels[variable]),
      x = unname(variable_labels[variable]),
      y = "Standardized accumulated local effect"
    ) +
    base_theme
})

comparison_figure <- wrap_plots(
  panels,
  design = paste(
    "AABBCC",
    "#DDEE#",
    sep = "\n"
  ),
  guides = "collect"
) +
  plot_annotation(
    title = "Comparison of accumulated local effects",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17)
    )
  ) &
  theme(legend.position = "bottom")

pdf_path <- file.path(
  figure_dir,
  "ale_comparison_brp_gable_rulefit_gbm_rf.pdf"
)
png_path <- file.path(
  figure_dir,
  "ale_comparison_brp_gable_rulefit_gbm_rf.png"
)
ggsave(
  pdf_path,
  comparison_figure,
  width = 11,
  height = 7.8,
  units = "in",
  bg = "white"
)
grDevices::png(
  filename = png_path,
  width = 11,
  height = 7.8,
  units = "in",
  res = 300,
  bg = "white"
)
print(comparison_figure)
grDevices::dev.off()

cat("Saved ALE comparison data:", data_path, "\n")
cat("Saved fitted-score scales:", scale_path, "\n")
cat("Saved ALE comparison PDF:", pdf_path, "\n")
cat("Saved ALE comparison PNG:", png_path, "\n")
