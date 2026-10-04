# Dynamic Statistical Arbitrage with GARCH-EVT Tail-Risk Control
# =============================================================================
# Research question:
# Does a time-varying ES/NQ hedge plus an EVT risk overlay improve out-of-sample
# downside performance relative to static and unfiltered trading rules?
#
# Yahoo's ES=F and NQ=F series are convenient continuous daily proxies. For a
# final study, set local_data_file to a roll-adjusted CSV with Date, ES, NQ.

required_packages <- c("quantmod", "tseries", "urca", "rugarch", "xts", "zoo")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install these packages first: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(quantmod)
  library(tseries)
  library(urca)
  library(rugarch)
  library(xts)
  library(zoo)
})
set.seed(4205)

# Configuration ----------------------------------------------------------------
config <- list(
  symbols = c(ES = "ES=F", NQ = "NQ=F"),
  start_date = as.Date("2018-01-01"),
  end_date = Sys.Date(),
  local_data_file = NULL,
  training_fraction = 0.60,
  minimum_training_observations = 500L,
  kalman_initial_window = 60L,
  entry_z = 2.0,
  exit_z = 0.50,
  stop_z = 4.0,
  transaction_cost_bps = 1.0,
  annualization = 252,
  evt_window = 500L,
  evt_threshold_probability = 0.90,
  evt_target_probability = 0.99,
  evt_min_exceedances = 25L,
  evt_refit_every = 5L,
  risk_history = 252L,
  risk_cutoff_probability = 0.90,
  stressed_position_scale = 0.25
)
stopifnot(
  config$entry_z > config$exit_z,
  config$evt_target_probability > config$evt_threshold_probability
)

# Data and small helpers --------------------------------------------------------
safe_test <- function(expr) {
  tryCatch(expr, error = function(e) {
    structure(list(error = conditionMessage(e)), class = "test_error")
  })
}

load_prices <- function(cfg) {
  if (!is.null(cfg$local_data_file)) {
    raw <- read.csv(cfg$local_data_file, stringsAsFactors = FALSE)
    if (!all(c("Date", "ES", "NQ") %in% names(raw))) {
      stop("Local data must contain Date, ES, and NQ columns.")
    }
    raw$Date <- as.Date(raw$Date)
    prices <- xts(raw[, c("ES", "NQ")], order.by = raw$Date)
  } else {
    market_env <- new.env(parent = emptyenv())
    downloaded_names <- suppressWarnings(getSymbols(
      unname(cfg$symbols), src = "yahoo", from = cfg$start_date,
      to = cfg$end_date, env = market_env, auto.assign = TRUE,
      warnings = FALSE
    ))
    stored_names <- ls(envir = market_env)
    if (length(stored_names) != 2L) {
      stop(
        "Yahoo did not return both requested series. Returned: ",
        paste(downloaded_names, collapse = ", ")
      )
    }
    es_name <- stored_names[grepl("^ES", stored_names)][1L]
    nq_name <- stored_names[grepl("^NQ", stored_names)][1L]
    es_object <- get(es_name, envir = market_env)
    nq_object <- get(nq_name, envir = market_env)
    prices <- merge(Ad(es_object), Ad(nq_object), join = "inner")
    names(prices) <- c("ES", "NQ")
  }
  prices <- na.omit(prices)
  prices <- prices[
    is.finite(prices$ES) & is.finite(prices$NQ) &
      prices$ES > 0 & prices$NQ > 0
  ]
  if (NROW(prices) < 750L) {
    stop("At least 750 aligned observations are required; found ", NROW(prices))
  }
  prices
}

descriptive_statistics <- function(x) {
  x <- x[is.finite(x)]
  centered <- x - mean(x)
  s <- sd(x)
  data.frame(
    n = length(x), mean = mean(x), sd = s,
    skewness = mean(centered^3) / s^3,
    excess_kurtosis = mean(centered^4) / s^4 - 3,
    minimum = min(x), maximum = max(x), row.names = NULL
  )
}

tls_hedge <- function(nq, es) {
  m <- cbind(NQ = as.numeric(nq), ES = as.numeric(es))
  fit <- prcomp(m, center = TRUE, scale. = FALSE)
  v <- fit$rotation[, ncol(fit$rotation)]
  if (abs(v[["NQ"]]) < sqrt(.Machine$double.eps)) {
    stop("TLS normalization failed: the NQ loading is effectively zero.")
  }
  beta <- -v[["ES"]] / v[["NQ"]]
  list(
    alpha = mean(m[, "NQ"] - beta * m[, "ES"]),
    beta = beta,
    fit = fit
  )
}

# Kalman model: NQ_t = alpha_t + beta_t ES_t + error_t -------------------------
kalman_filter <- function(log_variances, x, y, initial_theta,
                          initial_covariance, start_index, keep_path = TRUE) {
  variances <- exp(log_variances)
  measurement_variance <- variances[1L]
  state_variance <- diag(variances[2:3], 2L)
  identity_matrix <- diag(2L)
  n <- length(y)
  theta <- matrix(initial_theta, 2L)
  covariance <- initial_covariance
  log_likelihood <- 0

  innovation <- innovation_variance <- rep(NA_real_, n)
  alpha_prior <- beta_prior <- alpha_filtered <- beta_filtered <-
    rep(NA_real_, n)

  for (t in seq.int(start_index, n)) {
    theta_prior <- theta
    covariance_prior <- covariance + state_variance
    h <- matrix(c(1, x[t]), 1L)
    error <- as.numeric(y[t] - h %*% theta_prior)
    error_variance <- as.numeric(
      h %*% covariance_prior %*% t(h) + measurement_variance
    )
    if (!is.finite(error_variance) || error_variance <= 0) {
      return(list(negative_log_likelihood = 1e20))
    }
    gain <- covariance_prior %*% t(h) / error_variance
    theta <- theta_prior + gain * error
    kh <- identity_matrix - gain %*% h
    covariance <- kh %*% covariance_prior %*% t(kh) +
      gain %*% matrix(measurement_variance, 1L, 1L) %*% t(gain)
    covariance <- (covariance + t(covariance)) / 2
    log_likelihood <- log_likelihood - 0.5 * (
      log(2 * pi) + log(error_variance) + error^2 / error_variance
    )

    if (keep_path) {
      innovation[t] <- error
      innovation_variance[t] <- error_variance
      alpha_prior[t] <- theta_prior[1L]
      beta_prior[t] <- theta_prior[2L]
      alpha_filtered[t] <- theta[1L]
      beta_filtered[t] <- theta[2L]
    }
  }
  list(
    negative_log_likelihood = -log_likelihood,
    innovation = innovation,
    innovation_variance = innovation_variance,
    z_score = innovation / sqrt(innovation_variance),
    alpha_prior = alpha_prior,
    beta_prior = beta_prior,
    alpha_filtered = alpha_filtered,
    beta_filtered = beta_filtered
  )
}

fit_kalman <- function(x, y, training_end, initial_window) {
  initial_fit <- lm(y[seq_len(initial_window)] ~ x[seq_len(initial_window)])
  initial_theta <- unname(coef(initial_fit))
  residual_variance <- max(var(residuals(initial_fit)), 1e-8)
  x_variance <- max(var(x[seq_len(initial_window)]), 1e-8)
  initial_covariance <- diag(c(residual_variance, residual_variance / x_variance))
  start_index <- initial_window + 1L

  objective <- function(par) {
    kalman_filter(
      par, x[seq_len(training_end)], y[seq_len(training_end)],
      initial_theta, initial_covariance, start_index, FALSE
    )$negative_log_likelihood
  }
  starts <- list(
    log(c(residual_variance, residual_variance * 1e-4,
          residual_variance * 1e-6)),
    log(c(residual_variance, residual_variance * 1e-3,
          residual_variance * 1e-5)),
    log(c(residual_variance * 0.5, residual_variance * 1e-5,
          residual_variance * 1e-7))
  )
  fits <- lapply(starts, function(start) {
    optim(
      start, objective, method = "L-BFGS-B",
      lower = log(c(1e-10, 1e-12, 1e-12)),
      upper = log(c(10, 1, 1)),
      control = list(maxit = 1000)
    )
  })
  best <- fits[[which.min(vapply(fits, function(z) z$value, numeric(1)))]]
  if (best$convergence != 0L) {
    warning("Kalman optimization warning: ", best$message)
  }
  path <- kalman_filter(
    best$par, x, y, initial_theta, initial_covariance, start_index, TRUE
  )
  list(
    optimizer = best,
    parameters = setNames(exp(best$par), c("R", "Q_alpha", "Q_beta")),
    initial_theta = setNames(initial_theta, c("alpha", "beta")),
    path = path
  )
}

# Conditional volatility: choose sGARCH or GJR-GARCH by formation-sample BIC ----
fit_garch_candidates <- function(training_innovations) {
  models <- c("sGARCH", "gjrGARCH")
  fits <- lapply(models, function(model_name) {
    specification <- ugarchspec(
      variance.model = list(model = model_name, garchOrder = c(1L, 1L)),
      mean.model = list(armaOrder = c(0L, 0L), include.mean = FALSE),
      distribution.model = "std"
    )
    tryCatch(
      ugarchfit(
        specification, training_innovations, solver = "hybrid",
        solver.control = list(trace = 0)
      ),
      error = function(e) NULL
    )
  })
  names(fits) <- models
  bic <- vapply(fits, function(fit) {
    if (is.null(fit)) Inf else as.numeric(infocriteria(fit)[2L])
  }, numeric(1))
  if (all(!is.finite(bic))) stop("Both GARCH candidate fits failed.")
  selected <- names(which.min(bic))
  list(model = selected, fit = fits[[selected]], candidate_bic = bic)
}

filter_garch <- function(innovations, garch_result) {
  coefficients <- coef(garch_result$fit)
  omega <- unname(coefficients[["omega"]])
  alpha <- unname(coefficients[["alpha1"]])
  beta <- unname(coefficients[["beta1"]])
  gamma <- if ("gamma1" %in% names(coefficients)) {
    unname(coefficients[["gamma1"]])
  } else 0

  n <- length(innovations)
  variance <- next_variance <- rep(NA_real_, n)
  valid <- which(is.finite(innovations))
  persistence <- alpha + beta + gamma / 2
  initial_variance <- if (persistence < 0.999) {
    omega / (1 - persistence)
  } else {
    var(innovations[valid])
  }
  initial_variance <- max(initial_variance, 1e-10)
  variance[valid[1L]] <- initial_variance

  if (length(valid) > 1L) {
    for (t in valid[-1L]) {
      previous <- t - 1L
      if (!is.finite(innovations[previous]) || !is.finite(variance[previous])) {
        variance[t] <- initial_variance
      } else {
        shock <- innovations[previous]
        variance[t] <- omega + alpha * shock^2 +
          gamma * (shock < 0) * shock^2 + beta * variance[previous]
      }
      variance[t] <- max(variance[t], 1e-12)
    }
  }
  for (t in valid) {
    shock <- innovations[t]
    next_variance[t] <- omega + alpha * shock^2 +
      gamma * (shock < 0) * shock^2 + beta * variance[t]
    next_variance[t] <- max(next_variance[t], 1e-12)
  }
  list(
    coefficients = coefficients,
    sigma = sqrt(variance),
    sigma_next = sqrt(next_variance),
    standardized_residual = innovations / sqrt(variance)
  )
}

# EVT: rolling two-sided POT/GPD on absolute standardized innovations -----------
fit_gpd <- function(sample, threshold_probability, minimum_exceedances) {
  sample <- sample[is.finite(sample)]
  threshold <- as.numeric(
    quantile(sample, threshold_probability, names = FALSE)
  )
  excess <- sample[sample > threshold] - threshold
  if (length(excess) < minimum_exceedances) return(NULL)

  objective <- function(par) {
    scale <- exp(par[1L])
    shape <- par[2L]
    support <- 1 + shape * excess / scale
    if (any(!is.finite(support)) || any(support <= 0)) return(1e20)
    if (abs(shape) < 1e-6) {
      length(excess) * log(scale) + sum(excess) / scale
    } else {
      length(excess) * log(scale) +
        (1 / shape + 1) * sum(log(support))
    }
  }
  fit <- optim(
    c(log(max(mean(excess), 1e-6)), 0.05), objective,
    method = "L-BFGS-B",
    lower = c(log(1e-8), -0.45),
    upper = c(log(max(sample) * 100 + 1), 0.95),
    control = list(maxit = 500)
  )
  if (fit$convergence != 0L || !is.finite(fit$value)) return(NULL)
  list(
    threshold = threshold,
    scale = exp(fit$par[1L]),
    shape = fit$par[2L],
    exceedance_probability = length(excess) / length(sample),
    exceedances = length(excess)
  )
}

gpd_risk_measures <- function(fit, target_probability) {
  if (is.null(fit)) return(c(VaR = NA_real_, ES = NA_real_))
  tail_probability <- 1 - target_probability
  if (tail_probability >= fit$exceedance_probability) {
    return(c(VaR = NA_real_, ES = NA_real_))
  }
  ratio <- fit$exceedance_probability / tail_probability
  value_at_risk <- if (abs(fit$shape) < 1e-6) {
    fit$threshold + fit$scale * log(ratio)
  } else {
    fit$threshold + fit$scale / fit$shape * (ratio^fit$shape - 1)
  }
  expected_shortfall <- if (fit$shape < 1) {
    (value_at_risk + fit$scale - fit$shape * fit$threshold) /
      (1 - fit$shape)
  } else Inf
  c(VaR = value_at_risk, ES = expected_shortfall)
}

rolling_evt <- function(standardized_residuals, sigma_next, cfg) {
  n <- length(standardized_residuals)
  standardized_var <- standardized_es <- shape <- threshold <-
    rep(NA_real_, n)
  exceedances <- rep(NA_integer_, n)
  current_fit <- NULL
  current_risk <- c(VaR = NA_real_, ES = NA_real_)
  last_refit <- -Inf

  for (t in seq_len(n)) {
    available <- which(seq_len(n) <= t & is.finite(standardized_residuals))
    available <- tail(available, cfg$evt_window)
    enough <- length(available) >= max(100L, cfg$evt_min_exceedances * 4L)
    if (enough && (is.null(current_fit) ||
                   t - last_refit >= cfg$evt_refit_every)) {
      new_fit <- fit_gpd(
        abs(standardized_residuals[available]),
        cfg$evt_threshold_probability,
        cfg$evt_min_exceedances
      )
      if (!is.null(new_fit)) {
        current_fit <- new_fit
        current_risk <- gpd_risk_measures(
          current_fit, cfg$evt_target_probability
        )
        last_refit <- t
      }
    }
    if (!is.null(current_fit)) {
      standardized_var[t] <- current_risk[["VaR"]]
      standardized_es[t] <- current_risk[["ES"]]
      shape[t] <- current_fit$shape
      threshold[t] <- current_fit$threshold
      exceedances[t] <- current_fit$exceedances
    }
  }
  list(
    standardized_var = standardized_var,
    standardized_es = standardized_es,
    conditional_var = standardized_var * sigma_next,
    conditional_es = standardized_es * sigma_next,
    shape = shape,
    threshold = threshold,
    exceedances = exceedances
  )
}

make_risk_multiplier <- function(conditional_es, training_end, cfg) {
  n <- length(conditional_es)
  multiplier <- rep(1, n)
  for (t in seq.int(training_end, n)) {
    if (!is.finite(conditional_es[t]) || conditional_es[t] <= 0) next
    history_start <- max(1L, t - cfg$risk_history)
    history <- conditional_es[history_start:(t - 1L)]
    history <- history[is.finite(history) & history > 0]
    if (length(history) < 30L) next
    normal_risk <- median(history)
    stress_cutoff <- as.numeric(
      quantile(history, cfg$risk_cutoff_probability, names = FALSE)
    )
    multiplier[t] <- min(1, normal_risk / conditional_es[t])
    if (conditional_es[t] > stress_cutoff) {
      multiplier[t] <- min(multiplier[t], cfg$stressed_position_scale)
    }
  }
  multiplier
}

# Trading and evaluation --------------------------------------------------------
make_positions <- function(z_score, start_index, entry_z, exit_z, stop_z) {
  n <- length(z_score)
  position <- rep(0, n)
  for (t in seq.int(start_index, n)) {
    previous <- if (t == 1L) 0 else position[t - 1L]
    z <- z_score[t]
    if (!is.finite(z)) {
      position[t] <- previous
    } else if (previous == 0) {
      if (z <= -entry_z) position[t] <- 1
      if (z >= entry_z) position[t] <- -1
    } else if (previous == 1) {
      position[t] <- if (z >= -exit_z || z <= -stop_z) 0 else 1
    } else {
      position[t] <- if (z <= exit_z || z >= stop_z) 0 else -1
    }
  }
  position
}

backtest_pair <- function(position, beta, position_scale, es_return, nq_return,
                          cost_bps, dates, test_start) {
  n <- length(position)
  beta[!is.finite(beta)] <- 0
  position_scale[!is.finite(position_scale)] <- 1
  denominator <- 1 + abs(beta)
  target_nq <- position * position_scale / denominator
  target_es <- -position * position_scale * beta / denominator
  applied_nq <- c(0, head(target_nq, -1L))
  applied_es <- c(0, head(target_es, -1L))
  target_turnover <- c(
    abs(target_nq[1L]) + abs(target_es[1L]),
    abs(diff(target_nq)) + abs(diff(target_es))
  )
  applied_turnover <- c(0, head(target_turnover, -1L))
  gross_return <- applied_nq * nq_return + applied_es * es_return
  gross_return[!is.finite(gross_return)] <- 0
  trading_cost <- applied_turnover * cost_bps / 10000
  net_return <- gross_return - trading_cost
  net_return[seq_len(test_start - 1L)] <- 0
  data.frame(
    Date = dates, position = position, position_scale = position_scale,
    beta = beta, weight_ES = target_es, weight_NQ = target_nq,
    turnover = applied_turnover, gross_return = gross_return,
    trading_cost = trading_cost, net_return = net_return
  )
}

performance_metrics <- function(backtest, annualization, test_start) {
  selection <- seq.int(test_start, nrow(backtest))
  r <- backtest$net_return[selection]
  wealth <- cumprod(1 + r)
  drawdown <- wealth / cummax(c(1, wealth))[-1L] - 1
  losses <- -r
  var_99 <- as.numeric(quantile(losses, 0.99, names = FALSE))
  es_99 <- mean(losses[losses >= var_99])
  annual_return <- prod(1 + r)^(annualization / length(r)) - 1
  annual_volatility <- sd(r) * sqrt(annualization)
  c(
    annual_return = annual_return,
    annual_volatility = annual_volatility,
    sharpe = if (annual_volatility > 0) {
      mean(r) / sd(r) * sqrt(annualization)
    } else NA_real_,
    max_drawdown = min(drawdown),
    VaR_99 = var_99,
    ES_99 = es_99,
    total_return = tail(wealth, 1L) - 1,
    average_turnover = mean(backtest$turnover[selection]),
    active_days = sum(abs(backtest$position[selection]) > 0)
  )
}

# Run the study -----------------------------------------------------------------
prices_xts <- load_prices(config)
log_prices_xts <- log(prices_xts)
analysis_data <- data.frame(
  Date = as.Date(index(log_prices_xts)),
  ES = as.numeric(prices_xts$ES),
  NQ = as.numeric(prices_xts$NQ),
  log_ES = as.numeric(log_prices_xts$ES),
  log_NQ = as.numeric(log_prices_xts$NQ)
)
analysis_data$return_ES <- c(NA_real_, diff(analysis_data$log_ES))
analysis_data$return_NQ <- c(NA_real_, diff(analysis_data$log_NQ))

n_observations <- nrow(analysis_data)
training_end <- max(
  config$minimum_training_observations,
  floor(config$training_fraction * n_observations)
)
training_end <- min(training_end, n_observations - config$annualization)
if (training_end <= config$kalman_initial_window + 100L) {
  stop("Formation sample is too short for the requested models.")
}
test_start <- training_end + 1L
formation <- seq_len(training_end)

return_statistics <- rbind(
  ES = descriptive_statistics(analysis_data$return_ES),
  NQ = descriptive_statistics(analysis_data$return_NQ)
)
eda_tests <- list(
  ES_jarque_bera = safe_test(jarque.bera.test(na.omit(analysis_data$return_ES))),
  NQ_jarque_bera = safe_test(jarque.bera.test(na.omit(analysis_data$return_NQ))),
  ES_squared_ljung_box = safe_test(Box.test(
    na.omit(analysis_data$return_ES)^2, lag = 12L, type = "Ljung-Box"
  )),
  NQ_squared_ljung_box = safe_test(Box.test(
    na.omit(analysis_data$return_NQ)^2, lag = 12L, type = "Ljung-Box"
  )),
  ES_price_adf = safe_test(adf.test(analysis_data$log_ES[formation])),
  NQ_price_adf = safe_test(adf.test(analysis_data$log_NQ[formation])),
  ES_return_adf = safe_test(adf.test(na.omit(analysis_data$return_ES[formation]))),
  NQ_return_adf = safe_test(adf.test(na.omit(analysis_data$return_NQ[formation])))
)

formation_matrix <- cbind(
  NQ = analysis_data$log_NQ[formation],
  ES = analysis_data$log_ES[formation]
)
tls_fit <- tls_hedge(formation_matrix[, "NQ"], formation_matrix[, "ES"])
ols_fit <- lm(log_NQ ~ log_ES, data = analysis_data[formation, ])
static_spread <- analysis_data$log_NQ - tls_fit$alpha -
  tls_fit$beta * analysis_data$log_ES
static_z <- (static_spread - mean(static_spread[formation])) /
  sd(static_spread[formation])

cointegration_tests <- list(
  phillips_ouliaris = safe_test(
    ca.po(formation_matrix, demean = "constant", type = "Pu")
  ),
  johansen_trace = safe_test(ca.jo(
    formation_matrix, type = "trace", ecdet = "const", K = 2L,
    spec = "transitory"
  )),
  static_spread_adf = safe_test(adf.test(static_spread[formation]))
)

kalman_fit <- fit_kalman(
  analysis_data$log_ES, analysis_data$log_NQ,
  training_end, config$kalman_initial_window
)
kalman_path <- kalman_fit$path
training_innovations <- kalman_path$innovation[formation]
training_innovations <- training_innovations[is.finite(training_innovations)]
garch_fit <- fit_garch_candidates(training_innovations)
garch_path <- filter_garch(kalman_path$innovation, garch_fit)
evt_path <- rolling_evt(
  garch_path$standardized_residual, garch_path$sigma_next, config
)
risk_multiplier <- make_risk_multiplier(
  evt_path$conditional_es, training_end, config
)

dynamic_spread_tests <- list(
  formation_adf = safe_test(adf.test(na.omit(
    kalman_path$innovation[formation]
  ))),
  test_adf = safe_test(adf.test(na.omit(
    kalman_path$innovation[test_start:n_observations]
  ))),
  standardized_residual_ljung_box = safe_test(Box.test(
    na.omit(garch_path$standardized_residual[formation]),
    lag = 12L, type = "Ljung-Box"
  )),
  standardized_squared_ljung_box = safe_test(Box.test(
    na.omit(garch_path$standardized_residual[formation])^2,
    lag = 12L, type = "Ljung-Box"
  ))
)

static_position <- make_positions(
  static_z, training_end, config$entry_z, config$exit_z, config$stop_z
)
kalman_position <- make_positions(
  kalman_path$z_score, training_end,
  config$entry_z, config$exit_z, config$stop_z
)
static_beta <- rep(tls_fit$beta, n_observations)
dynamic_beta <- kalman_path$beta_filtered
dynamic_beta[!is.finite(dynamic_beta)] <- tls_fit$beta

static_backtest <- backtest_pair(
  static_position, static_beta, rep(1, n_observations),
  analysis_data$return_ES, analysis_data$return_NQ,
  config$transaction_cost_bps, analysis_data$Date, test_start
)
kalman_backtest <- backtest_pair(
  kalman_position, dynamic_beta, rep(1, n_observations),
  analysis_data$return_ES, analysis_data$return_NQ,
  config$transaction_cost_bps, analysis_data$Date, test_start
)
evt_backtest <- backtest_pair(
  kalman_position, dynamic_beta, risk_multiplier,
  analysis_data$return_ES, analysis_data$return_NQ,
  config$transaction_cost_bps, analysis_data$Date, test_start
)
performance <- rbind(
  Static_TLS = performance_metrics(
    static_backtest, config$annualization, test_start
  ),
  Kalman = performance_metrics(
    kalman_backtest, config$annualization, test_start
  ),
  Kalman_EVT = performance_metrics(
    evt_backtest, config$annualization, test_start
  )
)

# Console results and four core plots ------------------------------------------
cat("\nDynamic Statistical Arbitrage with EVT Tail-Risk Control\n")
cat(
  "Formation:", as.character(analysis_data$Date[1L]), "to",
  as.character(analysis_data$Date[training_end]), "(", training_end, "rows )\n"
)
cat(
  "Test:", as.character(analysis_data$Date[test_start]), "to",
  as.character(tail(analysis_data$Date, 1L)), "(",
  n_observations - training_end, "rows )\n\n"
)
cat("Return stylized facts:\n")
print(round(return_statistics, 6))
cat("\nStatic TLS beta (NQ on ES):", round(tls_fit$beta, 6), "\n")
cat("Static OLS beta:", round(unname(coef(ols_fit)[2L]), 6), "\n")
cat("Kalman variance parameters:\n")
print(signif(kalman_fit$parameters, 5))
cat("Selected volatility model:", garch_fit$model, "(Student-t innovations)\n")
cat("Candidate BIC:\n")
print(garch_fit$candidate_bic)
cat("\nOut-of-sample performance after transaction costs:\n")
print(round(performance, 5))

test_dates <- analysis_data$Date[test_start:n_observations]
static_wealth <- cumprod(
  1 + static_backtest$net_return[test_start:n_observations]
)
kalman_wealth <- cumprod(
  1 + kalman_backtest$net_return[test_start:n_observations]
)
evt_wealth <- cumprod(
  1 + evt_backtest$net_return[test_start:n_observations]
)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(2L, 2L), mar = c(4, 4, 3, 1))
matplot(
  analysis_data$Date,
  cbind(analysis_data$log_ES, analysis_data$log_NQ),
  type = "l", lty = 1, col = c("steelblue", "darkorange"),
  xlab = "Date", ylab = "Log price", main = "ES and NQ log prices"
)
legend(
  "topleft", c("ES", "NQ"), col = c("steelblue", "darkorange"),
  lty = 1, bty = "n"
)
abline(v = analysis_data$Date[training_end], lty = 2, col = "grey40")
plot(
  analysis_data$Date, dynamic_beta, type = "l", col = "navy",
  xlab = "Date", ylab = expression(beta[t]),
  main = "Filtered Kalman hedge ratio"
)
abline(h = tls_fit$beta, col = "firebrick", lty = 2)
abline(v = analysis_data$Date[training_end], lty = 2, col = "grey40")
legend(
  "topleft", c("Kalman", "Static TLS"),
  col = c("navy", "firebrick"), lty = c(1, 2), bty = "n"
)
plot(
  analysis_data$Date, kalman_path$z_score, type = "l", col = "grey25",
  xlab = "Date", ylab = "Innovation z-score",
  main = "Past-only trading signal"
)
abline(
  h = c(-config$entry_z, 0, config$entry_z),
  col = c("firebrick", "grey50", "firebrick"), lty = c(2, 1, 2)
)
abline(v = analysis_data$Date[training_end], lty = 2, col = "grey40")
matplot(
  test_dates, cbind(static_wealth, kalman_wealth, evt_wealth),
  type = "l", lty = 1, col = c("grey45", "steelblue", "darkgreen"),
  xlab = "Date", ylab = "Growth of $1",
  main = "Out-of-sample net wealth"
)
legend(
  "topleft", c("Static TLS", "Kalman", "Kalman + EVT"),
  col = c("grey45", "steelblue", "darkgreen"), lty = 1, bty = "n"
)

# All numerical outputs remain available after source("Project.R").
project_results <- list(
  config = config,
  data = analysis_data,
  split = list(training_end = training_end, test_start = test_start),
  return_statistics = return_statistics,
  eda_tests = eda_tests,
  tls_fit = tls_fit,
  ols_fit = ols_fit,
  cointegration_tests = cointegration_tests,
  kalman_fit = kalman_fit,
  dynamic_spread_tests = dynamic_spread_tests,
  garch_fit = garch_fit,
  garch_path = garch_path,
  evt_path = evt_path,
  risk_multiplier = risk_multiplier,
  backtests = list(
    Static_TLS = static_backtest,
    Kalman = kalman_backtest,
    Kalman_EVT = evt_backtest
  ),
  performance = performance
)
invisible(project_results)
