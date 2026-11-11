# Synthetic data and reporting helpers for the separate simulation example.
# These functions are not part of the installed package API.

simulate_cohort <- function(seed = 20261110L, n_per_group = 4000L,
                            a2 = 0, b2 = 1.4,
                            design = c("z2", "balanced"), gz = 0.5,
                            known_noncase_fraction = 0.1,
                            validation_fraction = 0.1) {
  design <- match.arg(design)
  set.seed(seed)
  n <- 2L * n_per_group
  group <- rep(c("1", "2"), each = n_per_group)
  latent_risk <- rbeta(n, rep(c(2, 4), each = n_per_group), 8)
  latent_logit <- qlogis(latent_risk)
  true_outcome <- rbinom(n, 1, latent_risk)
  risk <- plogis(rep(c(0, a2), each = n_per_group) +
                   rep(c(1, b2), each = n_per_group) * latent_logit +
                   rnorm(n, sd = 0.1))
  Z1 <- latent_logit + rnorm(n, sd = sqrt(0.1))
  Z2 <- Z1 + rnorm(n, sd = 1)
  Z3 <- Z2 + rnorm(n, sd = sqrt(5))
  sensitivity <- c(0.8, 0.9)
  specificity <- c(0.9, 0.8)
  phenotype <- integer(n)
  for (g in 1:2) {
    cases <- which(group == as.character(g) & true_outcome == 1L)
    noncases <- which(group == as.character(g) & true_outcome == 0L)
    phenotype[cases] <- rbinom(length(cases), 1, sensitivity[g])
    phenotype[noncases] <- rbinom(length(noncases), 1, 1 - specificity[g])
  }

  candidate <- rep(1L, n)
  validated <- integer(n)
  true_validation_probability <- numeric(n)
  risk_stratum <- rep(NA_character_, n)
  known_count <- as.integer(round(known_noncase_fraction * n_per_group))
  validation_budget <- validation_fraction * n_per_group
  diagnostics <- vector("list", 2L)
  for (g in 1:2) {
    subgroup <- which(group == as.character(g))
    noncases <- subgroup[true_outcome[subgroup] == 0L]
    # Fixed known-noncase counts are selected separately in each subgroup.
    known <- noncases[sample.int(length(noncases), known_count)]
    candidate[known] <- 0L
    ii <- subgroup[candidate[subgroup] == 1L]
    breaks <- quantile(risk[ii], c(0, 1 / 3, 2 / 3, 1), names = FALSE)
    risk_stratum[ii] <- as.character(cut(
      risk[ii], breaks, include.lowest = TRUE,
      labels = c("low", "middle", "high")
    ))
    if (design == "z2") {
      standardized_z2 <- as.numeric(scale(Z2[ii]))
      intercept <- uniroot(function(x) {
        sum(plogis(x + gz * standardized_z2)) - validation_budget
      }, interval = c(-40, 40), tol = 1e-11)$root
      probabilities <- plogis(intercept + gz * standardized_z2)
    } else {
      # Near-equal expected allocation across nonempty risk x phenotype cells.
      # Bounds are imposed while solving for the requested subgroup budget.
      cell <- interaction(risk_stratum[ii], phenotype[ii], drop = TRUE)
      inverse_cell_size <- 1 / ave(rep(1, length(ii)), cell, FUN = length)
      scale_log <- uniroot(function(x) {
        sum(pmin(0.95, pmax(0.03, exp(x) * inverse_cell_size))) -
          validation_budget
      }, interval = c(-45, 45), tol = 1e-11)$root
      probabilities <- pmin(0.95, pmax(0.03, exp(scale_log) * inverse_cell_size))
    }
    true_validation_probability[ii] <- probabilities
    validated[ii] <- rbinom(length(ii), 1, probabilities)
    diagnostics[[g]] <- data.frame(
      group = as.character(g), n = length(subgroup),
      known_noncases = length(known), candidates = length(ii),
      expected_validation = sum(probabilities),
      observed_validation = sum(validated[ii])
    )
  }
  outcome <- rep(NA_integer_, n)
  outcome[candidate == 0L] <- 0L
  outcome[validated == 1L] <- true_outcome[validated == 1L]
  observed <- data.frame(
    group, risk, outcome, candidate, validated, Z1, Z2, Z3, phenotype,
    risk_stratum = factor(risk_stratum, levels = c("low", "middle", "high"))
  )
  list(observed = observed, true_outcome = true_outcome,
       true_validation_probability = true_validation_probability,
       sampling = do.call(rbind, diagnostics), design = design)
}

fit_simulation_methods <- function(simulation, density = "llogit", tau = 0.2) {
  methods <- c("complete", "contaminated", "AIPW-C", "AIPW-P", "full")
  validation_formula <- if (simulation$design == "balanced") {
    ~ risk_stratum * phenotype
  } else {
    ~ Z1 + Z2 + Z3
  }
  fits <- setNames(vector("list", length(methods)), methods)
  rows <- lapply(methods, function(method) {
    warnings <- character()
    error_message <- ""
    analysis_data <- simulation$observed
    # Only the explicitly requested full-data benchmark receives hidden labels.
    if (method == "full") analysis_data$outcome <- simulation$true_outcome
    result <- tryCatch(withCallingHandlers(
      aTPRAIPW::fit_atpr(
        data = analysis_data, risk = "risk", group = "group",
        outcome = "outcome", candidate = "candidate", validated = "validated",
        validation_formula = validation_formula,
        outcome_formula = ~ Z1 + Z2 + Z3, phenotype = "phenotype",
        method = method, reference = "1", tau = tau, density = density,
        prob_floor = 0.005, se = method %in% c("AIPW-C", "AIPW-P"),
        variance = "stratified"
      ), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        error_message <<- conditionMessage(e)
        list(estimate = NA_real_, std.error = NA_real_)
      })
    fits[method] <<- list(result)
    warnings <- c(warnings, result$diagnostics_details$warnings)
    estimate <- result$estimate
    standard_error <- result$std.error
    data.frame(method = method, estimate = estimate, if_se = standard_error,
               point_success = is.finite(estimate),
               if_success = is.finite(standard_error), error = error_message,
               warnings = paste(unique(warnings), collapse = " | "))
  })
  estimates <- do.call(rbind, rows)
  benchmark <- estimates$estimate[estimates$method == "full"]
  estimates$benchmark <- benchmark
  estimates$paired_difference <- estimates$estimate - benchmark
  list(estimates = estimates, fits = fits)
}

summarize_replicates <- function(replicates) {
  keys <- interaction(replicates$a2, replicates$b2, replicates$method, drop = TRUE)
  summaries <- lapply(split(replicates, keys), function(rows) {
    point_ok <- rows$point_success
    paired_ok <- point_ok & is.finite(rows$benchmark)
    if_ok <- point_ok & rows$if_success
    point <- rows$estimate[point_ok]
    empirical_sd <- sd(point)
    mean_if_se <- mean(rows$if_se[if_ok])
    if_empirical_sd <- sd(rows$estimate[if_ok])
    data.frame(
      a2 = rows$a2[1L], b2 = rows$b2[1L], method = rows$method[1L],
      attempted = nrow(rows), point_successes = sum(point_ok),
      paired_successes = sum(paired_ok), if_successes = sum(if_ok),
      failures = sum(!point_ok),
      mean = mean(point),
      bias = mean(rows$paired_difference[paired_ok]),
      empirical_sd = empirical_sd, mean_if_se = mean_if_se,
      if_subset_empirical_sd = if_empirical_sd,
      ratio = mean_if_se / if_empirical_sd,
      q025 = quantile(point, 0.025, names = FALSE),
      q975 = quantile(point, 0.975, names = FALSE),
      benchmark_mean = mean(rows$benchmark[paired_ok])
    )
  })
  result <- do.call(rbind, summaries)
  method_order <- c("complete", "contaminated", "AIPW-C", "AIPW-P", "full")
  result <- result[order(result$a2, result$b2, match(result$method, method_order)), ]
  rownames(result) <- NULL
  result
}

plot_simulation_summary <- function(summary, file) {
  methods <- c("complete", "contaminated", "AIPW-C", "AIPW-P", "full")
  labels <- c("CC", "CL", "AIPW-C", "AIPW-P", "Full")
  colors <- c("#377DFF", "#F3A02E", "#9467BD", "#E41A1C", "black")
  a_values <- sort(unique(summary$a2))
  b_values <- sort(unique(summary$b2))
  grDevices::pdf(file, width = max(6, 3.5 * length(b_values)),
                 height = max(4.5, 3 * length(a_values)))
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(a_values), length(b_values)),
                mar = c(4, 4.5, 2, 0.8), mgp = c(2.6, 0.7, 0))
  for (a in a_values) for (b in b_values) {
    rows <- summary[summary$a2 == a & summary$b2 == b, ]
    rows <- rows[match(methods, rows$method), ]
    limits <- range(c(rows$q025, rows$q975, rows$mean), na.rm = TRUE)
    padding <- max(0.02, diff(limits) * 0.08)
    graphics::plot(1:5, rows$mean, type = "n", xaxt = "n", xlab = "Estimator",
                   ylab = "Difference in aTPR", xlim = c(0.6, 5.4),
                   ylim = limits + c(-padding, padding))
    graphics::axis(1, at = 1:5, labels = labels)
    benchmark <- rows$mean[rows$method == "full"]
    graphics::abline(h = benchmark, lty = 2, col = "grey50")
    graphics::segments(1:5, rows$q025, 1:5, rows$q975, col = colors)
    graphics::segments((1:5) - 0.08, rows$q025, (1:5) + 0.08, rows$q025, col = colors)
    graphics::segments((1:5) - 0.08, rows$q975, (1:5) + 0.08, rows$q975, col = colors)
    graphics::points(1:5, rows$mean, pch = 19, col = colors)
    if (length(a_values) * length(b_values) > 1L) {
      graphics::mtext(sprintf("a2 = %g, b2 = %g", a, b), side = 3, line = 0.4)
    }
  }
  invisible(file)
}
