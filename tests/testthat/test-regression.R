test_that("the packaged estimators reproduce the verified research implementation", {
  reference <- read.csv(test_path("fixtures", "reference-results.csv"))
  for (design in c("z2", "balanced")) {
    data <- read_fixture(design)
    expected <- reference[reference$design == design, ]
    for (j in seq_len(nrow(expected))) {
      row <- expected[j, ]
      info <- paste(row$design, row$density, row$method)
      fit <- fixture_fit(data, design, row$method, row$density)
      expect_equal(unname(fit$estimate), row$estimate, tolerance = 1e-6, info = info)
      if (!is.na(row$std.error)) {
        expect_equal(unname(fit$std.error), row$std.error,
                     tolerance = 1e-6, info = info)
        iid <- fixture_fit(data, design, row$method, row$density, variance = "iid")
        expect_equal(unname(iid$std.error), row$se_iid,
                     tolerance = 1e-6, info = paste(info, "iid"))
      }
    }
  }
})

test_that("unvalidated outcomes cannot change AIPW fits", {
  data <- read_fixture()
  unobserved <- data$candidate == 1 & data$validated == 0
  altered <- data
  altered$outcome[unobserved] <- rep(c(0, 1), length.out = sum(unobserved))
  altered$true_outcome <- 1 - altered$true_outcome
  for (method in c("AIPW-C", "AIPW-P")) {
    fit <- fixture_fit(data, method = method)
    hidden <- fixture_fit(altered, method = method)
    expect_identical(hidden$estimate, fit$estimate)
    expect_identical(hidden$std.error, fit$std.error)
    expect_identical(hidden$validation_probability, fit$validation_probability)
  }
})

test_that("density weights use logit calibrated risk and the correct prior odds", {
  data <- read_fixture()
  # Unequal group sizes reveal an omitted or inverted n_comparison/n_reference.
  data <- data[data$group == "1" | seq_len(nrow(data)) %% 2L == 0L, ]
  for (density in c("llogit", "qlogit")) {
    fit <- fixture_fit(data, method = "full", density = density)
    x <- qlogis(fit$calibrated_risk)
    reference <- as.numeric(data$group == "1")
    model <- if (density == "llogit") {
      glm(reference ~ x, family = binomial())
    } else {
      glm(reference ~ x + I(x^2), family = binomial())
    }
    weights <- exp(predict(model, type = "link")) *
      sum(reference == 0) / sum(reference == 1)
    expect_equal(fit$density_ratio[reference == 1], rep(1, sum(reference == 1)))
    expect_equal(unname(fit$density_ratio[reference == 0]),
                 unname(weights[reference == 0]), tolerance = 1e-6)
  }
})

test_that("complete-case fitting uses only known noncases and validated candidates", {
  data <- read_fixture()
  observed <- data$candidate == 0 | data$validated == 1
  fit <- fixture_fit(data, method = "complete", density = "qlogit")
  retained <- fixture_fit(data[observed, ], method = "full", density = "qlogit")
  expect_equal(fit$estimate, retained$estimate, tolerance = 1e-7)
})

test_that("column names and subgroup labels do not encode the estimator", {
  data <- read_fixture()
  original <- fixture_fit(data, method = "AIPW-P", density = "qlogit")
  renamed <- data
  columns <- c(risk = "score", group = "stratum", outcome = "gold",
               candidate = "pool", validated = "reviewed", phenotype = "proxy",
               Z1 = "age_proxy", Z2 = "clinical_proxy", Z3 = "noisy_proxy")
  names(renamed)[match(names(columns), names(renamed))] <- unname(columns)
  renamed$stratum <- factor(ifelse(renamed$stratum == "1", "reference", "comparison"),
                             levels = c("comparison", "reference"))
  fit <- fit_atpr(
    renamed, risk = "score", group = "stratum", outcome = "gold",
    candidate = "pool", validated = "reviewed", phenotype = "proxy",
    validation_formula = ~ age_proxy + clinical_proxy + noisy_proxy,
    outcome_formula = ~ age_proxy + clinical_proxy + noisy_proxy,
    method = "AIPW-P", reference = "reference", density = "qlogit"
  )
  expect_equal(fit$estimate, original$estimate, tolerance = 1e-10)
  expect_equal(fit$std.error, original$std.error, tolerance = 1e-10)
})
