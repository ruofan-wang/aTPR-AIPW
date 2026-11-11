read_fixture <- function(design = "z2") {
  readRDS(testthat::test_path("fixtures", paste0(design, ".rds")))
}

fixture_fit <- function(data, design = "z2", method = "AIPW-C",
                        density = "llogit", ...) {
  if (method == "full") data$outcome <- data$true_outcome
  validation <- if (design == "balanced") {
    ~ risk_stratum * phenotype
  } else {
    ~ Z1 + Z2 + Z3
  }
  fit_atpr(
    data, risk = "risk", group = "group", outcome = "outcome",
    candidate = "candidate", validated = "validated",
    validation_formula = validation, outcome_formula = ~ Z1 + Z2 + Z3,
    phenotype = "phenotype", method = method, reference = "1",
    density = density, se = method %in% c("AIPW-C", "AIPW-P"), ...
  )
}
