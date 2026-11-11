# Internal model helpers.

.atpr_model_matrix <- function(formula, data) {
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass,
                              drop.unused.levels = TRUE)
  stats::model.matrix(formula, frame)
}

.atpr_logistic <- function(x, y) {
  messages <- character()
  fit <- withCallingHandlers(
    stats::glm.fit(x, y, family = stats::binomial()),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  fit$warnings <- unique(messages)
  fit
}

.atpr_calibrate <- function(x, y) {
  n <- length(y)
  objective <- function(theta) {
    eta <- drop(x %*% theta)
    mean(pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta)
  }
  polish <- function(theta) {
    for (iteration in seq_len(100L)) {
      p <- stats::plogis(drop(x %*% theta))
      score <- drop(crossprod(x, y - p)) / n
      if (max(abs(score)) < 1e-10) return(theta)
      information <- crossprod(x, x * (p * (1 - p))) / n
      step <- drop(solve(information, score))
      value <- objective(theta)
      fraction <- 1
      for (halving in seq_len(35L)) {
        candidate <- theta + fraction * step
        candidate_value <- objective(candidate)
        if (candidate_value <= value + 1e-12) break
        fraction <- fraction / 2
      }
      theta <- candidate
    }
    theta
  }
  unname(polish(c(0, 1)))
}

.atpr_density_basis <- function(eta, degree) {
  ans <- if (degree == 1L) cbind(1, eta) else cbind(1, eta, eta^2)
  colnames(ans) <- c("(Intercept)", "logit_risk", "logit_risk_squared")[
    seq_len(degree + 1L)]
  ans
}

