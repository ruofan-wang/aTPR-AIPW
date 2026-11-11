#' Methods for fitted adjusted-TPR contrasts
#'
#' Print, summarize, or extract a point estimate and its working
#' influence-function uncertainty.
#' @param x,object An `atpr_fit` object returned by [fit_atpr()].
#' @param digits Number of digits printed.
#' @param level Confidence level for a working normal-approximation interval.
#' @param parm Unused; the object contains the single contrast delta.
#' @param ... Additional arguments; currently unused.
#' @return `coef()` returns the named contrast; `vcov()` returns its estimated
#'   variance; `confint()` returns a two-column matrix; and `summary()` returns a
#'   list with estimates, diagnostics, and a working normal interval.
#'   Print methods invisibly return their input.
#' @name atpr_fit_methods
NULL

#' @rdname atpr_fit_methods
#' @export
print.atpr_fit <- function(x, digits = 4, ...) {
  cat("Adjusted TPR difference: ", x$comparison, " - ", x$reference, "\n", sep = "")
  cat("Method: ", x$method, "; density: ", x$density,
      "; threshold: ", format(x$tau), "\n", sep = "")
  cat("Estimate: ", format(x$estimate, digits = digits), sep = "")
  cat("; IF standard error: ", format(x$std.error, digits = digits), sep = "")
  cat("\n")
  invisible(x)
}

#' @rdname atpr_fit_methods
#' @export
coef.atpr_fit <- function(object, ...) c(delta = object$estimate)

#' @rdname atpr_fit_methods
#' @export
vcov.atpr_fit <- function(object, ...) {
  matrix(object$std.error^2, nrow = 1L, dimnames = list("delta", "delta"))
}

#' @rdname atpr_fit_methods
#' @export
confint.atpr_fit <- function(object, parm = "delta", level = 0.95, ...) {
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  matrix(object$estimate + stats::qnorm(probs) * object$std.error,
         nrow = 1L, dimnames = list("delta", paste0(format(100 * probs, trim = TRUE), " %")))
}

#' @rdname atpr_fit_methods
#' @export
summary.atpr_fit <- function(object, level = 0.95, ...) {
  interval <- unname(stats::confint(object, level = level)[1L, ])
  ans <- list(method = object$method, density = object$density, tau = object$tau,
              reference = object$reference, comparison = object$comparison,
              coefficients = matrix(c(object$estimate, object$std.error, interval),
                                    nrow = 1L,
                                    dimnames = list("delta", c("Estimate", "IF SE", "Lower", "Upper"))),
              level = level, atpr = object$atpr, diagnostics = object$diagnostics,
              diagnostics_details = object$diagnostics_details)
  class(ans) <- "summary.atpr_fit"
  ans
}

#' @rdname atpr_fit_methods
#' @export
print.summary.atpr_fit <- function(x, digits = 4, ...) {
  cat("Adjusted TPR difference: ", x$comparison, " - ", x$reference, "\n", sep = "")
  cat("Method: ", x$method, "; density: ", x$density, "; threshold: ",
      format(x$tau), "\n", sep = "")
  print(x$coefficients, digits = digits)
  cat("Lower/Upper: ", format(100 * x$level), "% working normal interval.\n", sep = "")
  cat("\nSubgroup aTPR estimates:\n")
  print(x$atpr, digits = digits)
  cat("\nValidation diagnostics:\n")
  print(x$diagnostics, row.names = FALSE, digits = digits)
  invisible(x)
}
