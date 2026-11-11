# Source mapping

This package refactors the current, verified research implementation as an
extension of `fairRisk`. The original research scripts and archived simulation
results are retained separately.

## Upstream dependency

The original scripts sourced `fairRisk` 1.0.1.2 and called its calibration,
density-ratio, and aTPR functions. The package now uses exported functions from
`fairRisk` 1.1.1.1, which supports logit-transformed calibrated risk directly
through `transform = TRUE` in density-ratio estimation.

| Calculation | Implementation |
| --- | --- |
| Complete-case, contaminated-label, and full-data calibration | `fairRisk::calibrate()` |
| Linear and quadratic logit density ratios for all five methods | `fairRisk::estDensityRatio(transform = TRUE)` |
| Subgroup aTPRs and their comparison | `fairRisk::get_aTPR()` followed by comparison minus reference |
| Fitted validation and outcome models, AIPW pseudo-outcomes, augmented calibration | `aTPRAIPW` |
| Stacked influence-function standard errors | `aTPRAIPW` |

The upstream density function returns weights rather than model coefficients.
For AIPW inference, the same logistic model is refitted in a raw polynomial
basis to obtain the coefficients used by the stacked equations. Point estimates
and reported density ratios use the upstream outputs.

The upstream simulation repository,
[aTPR-estimation](https://github.com/sarahhegarty/aTPR-estimation), explicitly
uses the separate [fairRisk R package](https://github.com/sarahhegarty/fairRisk).
This package extends that same dependency. The tested installation is
`fairRisk` 1.1.1.1 at revision
`e8e0f8b6bc58f29fec68acb070064940ac88d13c`, pinned in `DESCRIPTION`.
The upstream README and the pinned revision were verified on GitHub on
2026-11-10.

## Research source mapping

| Research source | Role in the package |
| --- | --- |
| `logit_density/definitions.R` | Logistic subgroup membership in logit calibrated risk; aTPR standardization. |
| `logit_density/estimated_if.R` | Stacked influence function including fitted validation probabilities, outcome regressions, calibration, and density parameters. |
| `simulation_validation/calibration_solver_recovery.R` | Newton updates and step halving for the calibration estimating equations. |
| `simulation_validation/balanced_logistic.R` | Candidate-only validation fitting, validated-only outcome fitting, known-noncase contribution, and balanced allocation. |
| `estimated_validation/definitions.R` | Separate subgroup validation models using clinical covariates. |
| `logit_density/z2_linear_definitions.R` | Primary Z2-dependent design, slope 0.5, and current point-estimator audits. |
| `logit_density/z2_quadratic_definitions.R` | Quadratic logit density sensitivity. |
| `simulation_validation/balanced_logistic_var0p1.R` | Auxiliary noise variances 0.1, 1, and 5. |
| `sim_9settings_fairrisk.R` / `dataGen.R` | Outcome/risk generation and subgroup-specific known-noncase sampling. |

## What changed in the refactoring

- Replaced chains of `source()` and global settings with explicit function
  arguments, per-subgroup design matrices, and an R package namespace.
- Retained upstream functionality through a declared `fairRisk` dependency and
  direct calls to its exported functions.
- Made column names and regression predictors configurable.
- Kept the estimated-validation-probability calculation and the logit risk
  density models used in the current manuscript.
- Separated synthetic cohort generation and Monte Carlo loops from estimation.
- Added fitted-object methods, installed help, fit diagnostics, and numerical
  regression tests.

The example generates fresh synthetic random draws with explicit seeds. It
matches the current model and sampling setup but does not promise to reproduce
the archived 500-replicate manuscript tables bit for bit. The frozen regression
tests, in contrast, use selected original synthetic draws and numerical
reference outputs from the original implementation.

## Reference fixtures

The files in `tests/testthat/fixtures/` contain two synthetic cohorts of 8,000
observations, one per validation design, at `(a2, b2) = (0, 1.4)`. They contain
no clinical records. Expected point estimates cover the five methods under
linear and quadratic logit density estimation. Expected standard errors cover
both AIPW methods under each design/specification, including the subgroup-
centered and uncentered versions.

Fixtures should not be regenerated from the new package to make a failing
comparison pass. A change in methodology should instead be identified and
reviewed explicitly.

## Attribution

The original research code used `fairRisk` version 1.0.1.2, authored by Sarah
Hegarty and released under MIT. This package depends on version 1.1.1.1 or later,
uses its risk-standardization functions, and adds the partially observed-outcome
augmentation and inference from this research project. See `inst/NOTICE`.
