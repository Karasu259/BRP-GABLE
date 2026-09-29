
scenario_name <- function(scenario) {
  labels <- c(
    "Purely rule-based",
    "Purely linear",
    "Mixed linear and rule-based",
    "Smooth nonlinear misspecified"
  )
  scenario <- as.integer(scenario)
  if (length(scenario) != 1L || is.na(scenario) ||
      scenario < 1L || scenario > 4L) {
    stop("scenario must be one of 1, 2, 3, or 4.")
  }
  labels[scenario]
}

generate_predictors <- function(
    n,
    p = 20L,
    x_dist = c("normal", "uniform", "bernoulli"),
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  x_dist <- match.arg(x_dist)

  X <- switch(
    x_dist,
    normal = matrix(stats::rnorm(n * p), nrow = n, ncol = p),
    uniform = matrix(stats::runif(n * p, -1, 1), nrow = n, ncol = p),
    bernoulli = matrix(stats::rbinom(n * p, 1, 0.5), nrow = n, ncol = p)
  )
  colnames(X) <- paste0("X", seq_len(p))
  as.data.frame(X, check.names = FALSE)
}

true_function <- function(X, scenario, c0 = 0.36) {
  X <- as.data.frame(X, check.names = FALSE)
  scenario <- as.integer(scenario)

  if (scenario %in% 1:3 && ncol(X) < 9L) {
    stop("Scenarios 1-3 require at least 9 predictors.")
  }
  if (scenario == 4L && ncol(X) < 4L) {
    stop("Scenario 4 requires at least 4 predictors.")
  }

  linear_component <- if (ncol(X) >= 6L) {
    rowSums(as.matrix(X[, 1:6, drop = FALSE]))
  } else {
    rep(NA_real_, nrow(X))
  }

  rule_component <- if (ncol(X) >= 9L) {
    4 * rowSums(as.matrix(X[, 4:9, drop = FALSE]) > 0) +
      4 * as.numeric(X$X4 > 0 & X$X5 > 0) +
      4 * as.numeric(X$X6 > 0 & X$X7 > 0) -
      14
  } else {
    rep(NA_real_, nrow(X))
  }

  eta <- switch(
    as.character(scenario),
    "1" = c0 * rule_component,
    "2" = linear_component,
    "3" = 0.5 * linear_component + 0.5 * c0 * rule_component,
    "4" = 4 * (
      (sin(X$X1) + 1) / 2 +
        log1p(abs(X$X2)) -
        sqrt(abs(X$X3 * X$X4))
    ),
    stop("scenario must be one of 1, 2, 3, or 4.")
  )

  list(
    eta = as.numeric(eta),
    linear_component = as.numeric(linear_component),
    rule_component = as.numeric(rule_component),
    c0 = as.numeric(c0),
    scenario = scenario,
    scenario_name = scenario_name(scenario)
  )
}

generate_binary_data <- function(
    n,
    p = 20L,
    scenario,
    x_dist = "normal",
    c0 = 0.36,
    seed = NULL
) {
  X <- generate_predictors(n, p, x_dist, seed)
  sig <- true_function(X, scenario, c0)
  probability <- stats::plogis(sig$eta)
  y <- stats::rbinom(n, size = 1L, prob = probability)

  list(
    X = X,
    y = y,
    eta = sig$eta,
    probability = probability,
    scenario = sig$scenario,
    scenario_name = sig$scenario_name,
    c0 = sig$c0
  )
}

generate_continuous_data <- function(
    n,
    p = 20L,
    scenario,
    x_dist = "normal",
    c0 = 0.36,
    noise_sd = 1,
    seed = NULL
) {
  X <- generate_predictors(n, p, x_dist, seed)
  sig <- true_function(X, scenario, c0)
  error <- stats::rnorm(n, mean = 0, sd = noise_sd)
  y <- sig$eta + error

  list(
    X = X,
    y = y,
    eta = sig$eta,
    error = error,
    scenario = sig$scenario,
    scenario_name = sig$scenario_name,
    c0 = sig$c0
  )
}

generate_survival_data <- function(
    n,
    p = 20L,
    scenario,
    x_dist = "normal",
    c0 = 0.36,
    baseline_hazard = 0.1,
    censor_upper = 20,
    seed = NULL
) {
  X <- generate_predictors(n, p, x_dist, seed)
  sig <- true_function(X, scenario, c0)

  individual_hazard <- baseline_hazard * exp(sig$eta)
  true_time <- stats::rexp(n, rate = individual_hazard)
  censor_time <- stats::runif(n, min = 0, max = censor_upper)

  time <- pmin(true_time, censor_time)
  status <- as.integer(true_time <= censor_time)

  list(
    X = X,
    time = time,
    status = status,
    true_time = true_time,
    censor_time = censor_time,
    eta = sig$eta,
    individual_hazard = individual_hazard,
    scenario = sig$scenario,
    scenario_name = sig$scenario_name,
    c0 = sig$c0
  )
}
