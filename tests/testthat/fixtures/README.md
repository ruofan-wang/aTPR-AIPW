# Frozen synthetic regression fixtures

These two cohorts contain simulated data only. They do not contain any EHR
records or other patient data. Each has 4,000 individuals per subgroup, exactly
400 known noncases per subgroup, and a validation design with an expected 400
validated candidates per subgroup. Unvalidated outcomes are missing in the
`outcome` column. `true_outcome` is retained solely for the simulated full-data
benchmark and leakage checks.

The scenario is `a2 = 0`, `b2 = 1.4`, replicate 1. The research cohort seed is
60001; the validation seed is 202610081. Auxiliary error variances are 0.1, 1,
and 5. `z2.rds` uses standardized Z2-dependent sampling with slope 0.5;
`balanced.rds` uses risk-tertile by phenotype allocation.

`reference-results.csv` was generated from the verified research code on
2026-11-10, before substituting the package implementation:

- `logit_density/z2_linear_definitions.R`: cohort and validation draws;
- `logit_density/definitions.R`: logit-risk density and point estimates;
- `logit_density/estimated_if.R`: stacked influence-function standard errors;
- `simulation_validation/calibration_solver_recovery.R`: calibration solver.

Point estimates were independently compared with the IF implementation to
within 1e-7 when generating these fixtures. The CSV stores both standard errors
centered within fixed-size subgroups and the uncentered IID counterpart. Tests
compare the package against these values under both density specifications and
both sampling designs. The fixtures are test data, not a simulation result table
for the manuscript.
