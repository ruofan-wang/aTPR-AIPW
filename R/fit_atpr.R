#' Estimate a risk-standardized true positive rate difference
#'
#' Fit covariate-based or phenotype-assisted augmented inverse probability
#' weighted (AIPW) estimators using a partially labeled cohort. Complete-case,
#' contaminated-label, and full-data point estimators are also available.
#'
#' @param data A data frame containing one row per cohort member.
#' @param risk,group,outcome,candidate,validated Column names, supplied as strings.
#'   `risk` is the original prediction score in `(0, 1)`; `group` contains exactly
#'   two subgroup labels; `outcome` is the binary gold-standard outcome;
#'   `candidate` is 1 for the candidate pool and 0 for known noncases; and
#'   `validated` is 1 for candidates with observed gold-standard outcomes.
#'   Outcomes of unvalidated candidates may be `NA` and are not used by AIPW or
#'   complete-case estimation. Known noncases must have outcome 0 or `NA`, and
#'   validation indicator 0.
#' @param validation_formula One-sided formula for the logistic validation model,
#'   fitted separately in each subgroup using all candidates.
#' @param outcome_formula One-sided formula for the candidate outcome model,
#'   fitted separately in each subgroup using only validated candidates.
#'   Formula predictors must be observed for all candidates. Offsets are not
#'   supported. For balanced sampling, precompute risk strata separately within
#'   each subgroup's candidate pool and use, for example, `~ risk_bin * phenotype`.
#' @param phenotype Name of the binary algorithmic-phenotype column. Required for
#'   `method = "AIPW-P"` or `"contaminated"`. AIPW-P adds its main effect to
#'   `outcome_formula`; it is not added to the validation formula automatically.
#' @param method Estimation method. `"AIPW-C"` uses the supplied outcome formula;
#'   `"AIPW-P"` additionally uses the phenotype; `"complete"` uses known noncases
#'   and validated candidates at every stage; `"contaminated"` uses the phenotype
#'   as the outcome throughout the cohort; and `"full"` requires observed true
#'   outcomes for every cohort member.
#' @param reference Label of the reference subgroup. Defaults to the first
#'   subgroup appearing in `data`; specifying it explicitly is recommended.
#' @param tau Threshold applied to the original score: a positive prediction is
#'   defined by `risk > tau`.
#' @param density Density-ratio specification. `"llogit"` uses a linear term in
#'   the logit of the calibrated risk; `"qlogit"` additionally uses its square.
#' @param prob_floor Lower bound on fitted validation probabilities used in the
#'   AIPW residual correction. The default 0.005 limits inverse weights to 200.
#'   Set to 0 to disable this bound.
#' @param risk_floor Numerical bound applied to the original score before its
#'   logit is used in AIPW and complete-case calibration. Full-data and
#'   contaminated-label calibration use the unbounded logit, matching the
#'   original simulation implementation. Set to 0 to disable bounding.
#' @param se Calculate the working influence-function standard error for AIPW.
#'   Standard errors are not implemented for the three comparator methods, for
#'   which `std.error` is always `NA`.
#' @param variance `"stratified"` centers influence contributions within each
#'   subgroup, treating subgroup sample sizes as fixed; `"iid"` uses uncentered
#'   contributions under independent identically distributed sampling.
#'
#' @details
#' Risk calibration is logistic and linear in the logit of the original score.
#' Density estimation models reference-group membership using the chosen basis
#' in the logit of the calibrated risk. Its odds, multiplied by the comparison
#' to reference sample-size ratio, estimate the reference density divided by the
#' comparison density. Each subgroup's aTPR averages `I(risk > tau)` weighted by
#' calibrated risk and density ratio; the reported contrast is comparison minus
#' reference.
#'
#' The package extends `fairRisk`: ordinary calibration uses
#' [fairRisk::calibrate()], density ratios use [fairRisk::estDensityRatio()]
#' with `transform = TRUE`, and aTPRs use [fairRisk::get_aTPR()]. AIPW
#' pseudo-outcomes and augmented calibration are implemented here because
#' pseudo-outcomes can lie outside the binary-response range. For inference,
#' the equivalent density model is also fitted in a raw polynomial basis to
#' obtain coefficients not returned by the upstream density function.
#'
#' The stacked working influence function includes fitted validation, outcome,
#' calibration, and density-ratio parameters. Generated strata, known-noncase
#' selection, and other supplied design features are treated as fixed. It does
#' not add a separate correction for fixed-count noncase sampling or adaptive
#' validation designs. Fitted probability bounding is included in the local
#' derivatives away from the bound.
#'
#' @return An object of class `atpr_fit` with `estimate`, `std.error`, subgroup
#'   `atpr`, per-row `calibrated_risk`, `density_ratio`, and
#'   `validation_probability`, subgroup `diagnostics`, and `diagnostics_details`.
#'   Per-row outputs follow the original data order; excluded complete cases
#'   have `NA`. The `engine` component contains the parameter vector, score
#'   function, sensitivity matrix, and gradient used for AIPW inference and is
#'   intended for numerical audits. Methods include [print()], [summary()],
#'   [stats::coef()], [stats::vcov()], and [stats::confint()].
#' @export
fit_atpr <- function(data, risk, group, outcome, candidate, validated,
                     validation_formula = ~ Z1 + Z2 + Z3,
                     outcome_formula = ~ Z1 + Z2 + Z3,
                     phenotype = NULL,
                     method = c("AIPW-C", "AIPW-P", "complete", "contaminated", "full"),
                     reference = NULL, tau = 0.2,
                     density = c("llogit", "qlogit"),
                     prob_floor = 0.005, risk_floor = 1e-4, se = TRUE,
                     variance = c("stratified", "iid")) {
  call <- match.call()
  method <- match.arg(method)
  density <- match.arg(density)
  variance <- match.arg(variance)
  score <- data[[risk]]
  labels <- as.character(data[[group]])
  y <- data[[outcome]]
  q <- as.numeric(data[[candidate]])
  v <- as.numeric(data[[validated]])
  if (is.null(reference)) reference <- labels[1L]
  levels <- c(as.character(reference), setdiff(unique(labels), as.character(reference)))
  groups <- match(labels, levels)
  n <- nrow(data)
  # Hidden candidate outcomes do not enter AIPW or complete-case estimation.
  observed <- q == 0 | v == 1
  y_work <- numeric(n)
  if (method %in% c("AIPW-C", "AIPW-P", "complete")) {
    y_work[q == 1 & v == 1] <- as.numeric(y[q == 1 & v == 1])
  }
  phen <- NULL
  if (method %in% c("AIPW-P", "contaminated")) {
    phen <- data[[phenotype]]
  }
  degree <- if (density == "llogit") 1L else 2L
  rows <- if (method == "complete") which(observed) else seq_len(n)
  xrisk <- if (method %in% c("AIPW-C", "AIPW-P", "complete")) {
    stats::qlogis(pmin(pmax(score, risk_floor), 1 - risk_floor))
  } else stats::qlogis(score)
  gmat <- cbind("(Intercept)" = 1, logit_risk = xrisk)
  engine <- NULL
  standard_error <- NA_real_
  fitted_probability <- rep(NA_real_, n)
  models <- list()

  if (method %in% c("AIPW-C", "AIPW-P")) {
    if (method == "AIPW-P") {
      # Appending a main effect is safe when it is already present: terms()
      # removes duplicates. Merely appearing in an interaction is insufficient.
      rhs <- call("+", outcome_formula[[2L]], as.name(phenotype))
      outcome_formula <- stats::as.formula(call("~", rhs), env = environment(outcome_formula))
    }
    nuisance <- vector("list", 2L)
    theta <- vector("list", 2L)
    eta <- numeric(n)
    fitted_probability[q == 0] <- 1
    for (s in seq_len(2L)) {
      all_rows <- which(groups == s)
      candidates <- which(groups == s & q == 1)
      validated_rows <- which(groups == s & q == 1 & v == 1)
      a <- .atpr_model_matrix(validation_formula, data[candidates, , drop = FALSE])
      w <- .atpr_model_matrix(outcome_formula, data[candidates, , drop = FALSE])
      gf <- .atpr_logistic(a, v[candidates])
      validated_positions <- match(validated_rows, candidates)
      mf <- .atpr_logistic(w[validated_positions, , drop = FALSE], y_work[validated_rows])
      fitted_probability[candidates] <- stats::plogis(drop(a %*% gf$coefficients))
      outcome_probability <- stats::plogis(drop(w %*% mf$coefficients))
      bounded <- pmax(fitted_probability[validated_rows], prob_floor)
      pseudo <- numeric(n)
      pseudo[candidates] <- outcome_probability
      pseudo[validated_rows] <- outcome_probability[validated_positions] +
        (y_work[validated_rows] - outcome_probability[validated_positions]) / bounded
      theta[[s]] <- .atpr_calibrate(gmat[all_rows, , drop = FALSE], pseudo[all_rows])
      eta[all_rows] <- drop(gmat[all_rows, , drop = FALSE] %*% theta[[s]])
      afull <- matrix(0, n, ncol(a), dimnames = list(NULL, colnames(a)))
      wfull <- matrix(0, n, ncol(w), dimnames = list(NULL, colnames(w)))
      afull[candidates, ] <- a
      wfull[candidates, ] <- w
      nuisance[[s]] <- list(A = afull, W = wfull, gamma = gf$coefficients,
                            beta = mf$coefficients, validation = gf, outcome = mf)
    }
    dens <- .atpr_density(eta, groups, degree, score, tau, need_fit = se)
    if (se) {
      engine <- .atpr_influence(nuisance, theta, dens$fit$coefficients,
                                gmat, groups, q, v, y_work, score > tau,
                                degree, prob_floor)
      standard_error <- if (variance == "stratified") engine$se_stratified else engine$se_iid
    }
    models <- list(validation = lapply(nuisance, `[[`, "validation"),
                   outcome = lapply(nuisance, `[[`, "outcome"),
                   calibration = theta, density = dens$fit)
  } else {
    if (method == "full") {
      y_work <- as.numeric(y)
    } else if (method == "contaminated") {
      y_work <- as.numeric(phen)
    }
    eta <- rep(NA_real_, n)
    calibration_models <- vector("list", 2L)
    for (s in seq_len(2L)) {
      subgroup_rows <- rows[groups[rows] == s]
      fit <- .atpr_ordinary_calibrate(xrisk[subgroup_rows], y_work[subgroup_rows])
      eta[subgroup_rows] <- fit$eta
      calibration_models[[s]] <- fit
    }
    dens <- .atpr_density(eta[rows], groups[rows], degree, score[rows], tau)
    models <- list(calibration = calibration_models, density = dens$fit)
  }
  calibrated <- density_ratio <- rep(NA_real_, n)
  calibrated[rows] <- dens$calibrated
  density_ratio[rows] <- dens$ratio
  names(dens$atpr) <- levels
  diagnostics <- do.call(rbind, lapply(seq_len(2L), function(s) {
    cand <- groups == s & q == 1
    val <- cand & v == 1
    pp <- fitted_probability[cand]
    weights <- 1 / pmax(fitted_probability[val], prob_floor)
    data.frame(group = levels[s], n = sum(groups == s),
               n_known_noncase = sum(groups == s & q == 0),
               n_candidates = sum(cand), n_validated = sum(val),
               min_validation_probability = min(pp),
               max_validation_probability = max(pp),
               n_candidates_below_floor = sum(pp < prob_floor),
               n_truncated = sum(fitted_probability[val] < prob_floor),
               max_inverse_weight = max(weights))
  }))
  warnings <- unlist(lapply(models$validation, `[[`, "warnings"), use.names = FALSE)
  warnings <- unique(c(warnings, unlist(lapply(models$outcome, `[[`, "warnings"),
                                           use.names = FALSE), dens$fit$warnings))
  if (!method %in% c("AIPW-C", "AIPW-P")) {
    warnings <- unique(c(warnings, unlist(lapply(models$calibration, `[[`, "warnings"),
                                              use.names = FALSE)))
  }
  result <- list(call = call, method = method, density = density, tau = tau,
                 reference = levels[1L], comparison = levels[2L],
                 estimate = dens$estimate, std.error = standard_error,
                 atpr = dens$atpr, calibrated_risk = calibrated,
                 density_ratio = density_ratio,
                 validation_probability = fitted_probability,
                 diagnostics = diagnostics,
                 diagnostics_details = list(
                   warnings = warnings,
                   score_residual = engine$score_residual,
                   scaled_score_residual = engine$scaled_score_residual,
                   rcond = engine$rcond,
                   max_density_ratio = max(density_ratio, na.rm = TRUE)),
                 prob_floor = prob_floor, risk_floor = risk_floor, variance = variance,
                 models = models, engine = engine)
  class(result) <- "atpr_fit"
  result
}
