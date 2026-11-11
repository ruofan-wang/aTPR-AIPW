# Adapters to the exported fairRisk calibration and standardization functions.

.atpr_upstream <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(messages))
}

.atpr_ordinary_calibrate <- function(logit_risk, outcome) {
  dat <- data.frame(y = outcome, lp.gX = logit_risk)
  result <- .atpr_upstream(fairRisk::calibrate(
    train = dat, test = dat, method = "logit", args = list(1L)))
  list(eta = stats::qlogis(result$value$rs.gX),
       fitted.values = result$value$rs.gX, warnings = result$warnings,
       source = "fairRisk::calibrate")
}

.atpr_density <- function(eta, groups, degree, risk, tau, need_fit = FALSE) {
  calibrated <- stats::plogis(eta)
  dat <- data.frame(s = as.character(groups), rs.gX = calibrated,
                    gX = risk, row_id = seq_along(groups))
  result <- .atpr_upstream(fairRisk::estDensityRatio(
    train = dat, test = dat, method = "logit", args = list(degree),
    refgp = "1", transform = TRUE))
  dr <- result$value[order(result$value$row_id), , drop = FALSE]
  # Inject column symbols into fairRisk's data-masking interface.
  rates <- .atpr_upstream(suppressMessages(rlang::inject(
    fairRisk::get_aTPR(dr, orig_risk = !!rlang::sym("gX"),
                      cal_risk = !!rlang::sym("rs.gX"),
                      dens_ratio = !!rlang::sym("w_s"),
                      groupvar = !!rlang::sym("s"), taus = tau))))
  atpr <- rates$value$aTPR[match(c("1", "2"), rates$value$s)]
  fit <- list(source = "fairRisk::estDensityRatio")
  if (need_fit) {
    # fairRisk returns weights but not coefficients. Refit the equivalent raw
    # polynomial basis solely to assemble the stacked influence function.
    fit <- .atpr_logistic(.atpr_density_basis(eta, degree),
                          as.numeric(groups == 1L))
    fit$source <- "fairRisk weights; equivalent coefficient fit for inference"
  }
  fit$warnings <- unique(c(result$warnings, rates$warnings, fit$warnings))
  list(fit = fit, calibrated = calibrated, ratio = dr$w_s, atpr = atpr,
       estimate = unname(atpr[2L] - atpr[1L]))
}
