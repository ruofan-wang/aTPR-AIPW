central_jacobian <- function(fun, parameters) {
  vapply(seq_along(parameters), function(j) {
    step <- 1e-5 * max(1, abs(parameters[j]))
    upper <- lower <- parameters
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    (fun(upper) - fun(lower)) / (2 * step)
  }, numeric(length(fun(parameters))))
}

test_that("stacked IF derivatives agree with independent finite differences", {
  data <- read_fixture()
  # A deliberately active probability floor also checks its derivative.
  for (floor in c(0.005, 0.15)) {
    fit <- fixture_fit(data, method = "AIPW-P", density = "qlogit",
                       prob_floor = floor)
    engine <- fit$engine
    numerical_j <- -central_jacobian(
      function(z) colMeans(engine$components(z)$L), engine$z
    )
    numerical_gradient <- as.numeric(central_jacobian(
      function(z) engine$components(z, with_scores = FALSE)$delta, engine$z
    ))
    expect_lt(max(abs(engine$J - numerical_j)), 2e-6)
    expect_lt(max(abs(engine$gradient - numerical_gradient)), 2e-6)
  }
})
