# aTPR-AIPW

aTPR-AIPW is an R package for comparing adjusted true positive rates (aTPRs)
between two subgroups when gold-standard outcomes are partially observed.
It extends [fairRisk](https://github.com/sarahhegarty/fairRisk) with
covariate-based and phenotype-assisted augmented inverse probability weighted
(AIPW) estimators and influence-function standard errors.

## Method Overview

Both AIPW estimators follow the same workflow:

1. Fit a logistic validation model within each subgroup using all candidates.
2. Fit an outcome model using validated candidates and construct AIPW
   pseudo-outcomes. Known noncases contribute outcome zero.
3. Estimate subgroup-specific calibration functions, linear in the logit of
   the original risk score.
4. Estimate the density ratio using subgroup membership and calibrated risk.
5. Standardize the comparison subgroup to the reference subgroup's
   calibrated-risk distribution and calculate the aTPR difference.

| Method | Outcome model |
|---|---|
| `"AIPW-C"` | Uses the covariates specified in `outcome_formula`. |
| `"AIPW-P"` | Adds the algorithmic phenotype to the same outcome model. |

The estimated difference is **comparison subgroup minus reference subgroup**.
Validation probabilities are estimated from the observed data. Both methods
use `fairRisk` for density-ratio estimation and aTPR calculation; the AIPW
augmentation and influence-function calculation are implemented here.

## Installation

Install the package from GitHub:

```r
# If needed: install.packages("remotes")
remotes::install_github("ruofan-wang/aTPR-AIPW", dependencies = TRUE)
```

The package dependencies, including `fairRisk`, are installed automatically.
The repository is named `aTPR-AIPW`; the installed R package is named
`aTPRAIPW`. R version 4.1.0 or later is required.

## Example Data

The examples below use a simulated cohort with two subgroups, three clinical
auxiliary variables, and a binary algorithmic phenotype. Download or clone
the repository, then run them from its root directory after installation.

```r
library(aTPRAIPW)
source("examples/simulation_helpers.R")

simulation <- simulate_cohort(
  seed = 20261110,
  n_per_group = 4000,
  design = "z2"
)

dat <- simulation$observed
```

Each subgroup contains 4,000 observations. Exactly 400 known noncases are
sampled within each subgroup, and the expected validation sample size is
400 per subgroup. Validation depends on standardized Z2.

The analysis data contain one row per observation:

| Column | Description |
|---|---|
| `risk` | Original risk score, strictly between 0 and 1. |
| `group` | Subgroup label: `"1"` or `"2"` in this example. |
| `candidate` | 0 for known noncases; 1 for candidates. |
| `validated` | 1 for validated candidates; 0 otherwise, including known noncases. |
| `outcome` | Gold-standard outcome for validated candidates, 0 for known noncases, and `NA` for unvalidated candidates. |
| `Z1`, `Z2`, `Z3` | Clinical auxiliary variables. |
| `phenotype` | Binary algorithmic phenotype. |

Column names are supplied as strings. Other observed covariates can be used
through the model formulas. The simulated outcomes of unvalidated candidates
are stored separately and are not used by either AIPW estimator.

## Covariate-Based AIPW

Use `method = "AIPW-C"` to fit the outcome model with clinical covariates:

```r
fit_c <- fit_atpr(
  data = dat,
  risk = "risk",
  group = "group",
  outcome = "outcome",
  candidate = "candidate",
  validated = "validated",
  validation_formula = ~ Z1 + Z2 + Z3,
  outcome_formula = ~ Z1 + Z2 + Z3,
  method = "AIPW-C",
  reference = "1",
  tau = 0.2,
  density = "llogit"
)

fit_c$estimate
fit_c$std.error
confint(fit_c)
```

The positive prediction rule is `risk > tau`. Setting `reference = "1"`
standardizes subgroup 2 to subgroup 1 and reports aTPR(2) minus aTPR(1).

## Phenotype-Assisted AIPW

Use `method = "AIPW-P"` and specify the phenotype column to include it as an
additional main effect in the outcome model:

```r
fit_p <- fit_atpr(
  data = dat,
  risk = "risk",
  group = "group",
  outcome = "outcome",
  candidate = "candidate",
  validated = "validated",
  validation_formula = ~ Z1 + Z2 + Z3,
  outcome_formula = ~ Z1 + Z2 + Z3,
  phenotype = "phenotype",
  method = "AIPW-P",
  reference = "1",
  tau = 0.2,
  density = "llogit"
)

fit_p$estimate
fit_p$std.error
confint(fit_p)
```

The phenotype is added to the outcome model automatically. To include it in
the validation model, specify it separately in `validation_formula`.

## Density-Ratio Models

Both specifications model subgroup membership using the logit of calibrated
risk:

| `density` | Specification |
|---|---|
| `"llogit"` | Linear term in the logit of calibrated risk. |
| `"qlogit"` | Linear and quadratic terms in the logit of calibrated risk. |

Use `density = "qlogit"` in either call above to fit the quadratic logit
specification. Risk calibration remains linear in the logit of the original
risk score.

## Arguments

| Argument | Description |
|---|---|
| `data` | Data frame with one row per cohort member. |
| `risk`, `group`, `outcome`, `candidate`, `validated` | Names of the corresponding columns. |
| `validation_formula` | One-sided logistic model formula, fitted to candidates separately by subgroup. |
| `outcome_formula` | One-sided logistic model formula, fitted to validated candidates separately by subgroup. |
| `phenotype` | Binary phenotype column name, required for AIPW-P and contaminated-label estimation. |
| `method` | `"AIPW-C"`, `"AIPW-P"`, `"complete"`, `"contaminated"`, or `"full"`. |
| `reference` | Reference subgroup label. Defaults to the first subgroup appearing in the data. |
| `tau` | Threshold applied to the original risk score; default is 0.2. |
| `density` | `"llogit"` or `"qlogit"`; default is `"llogit"`. |
| `prob_floor` | Lower bound for fitted validation probabilities used in weighting; default is 0.005. |
| `risk_floor` | Bound applied before taking the original score's logit in AIPW and complete-case calibration; default is 0.0001. |
| `se` | Calculate AIPW influence-function standard errors; default is `TRUE`. |
| `variance` | `"stratified"` for fixed subgroup sizes, or `"iid"`; default is `"stratified"`. |

The complete-case comparator uses known noncases and validated candidates.
The contaminated-label comparator uses the phenotype as its outcome throughout
the cohort. The full-data comparator requires gold-standard outcomes for all
observations. Standard errors are available for the two AIPW methods.

## Returned Values

`fit_atpr()` returns an object of class `atpr_fit` containing:

| Value | Description |
|---|---|
| `estimate` | Estimated aTPR difference: comparison minus reference. |
| `std.error` | Influence-function standard error for an AIPW estimator. |
| `atpr` | Adjusted true positive rate in each subgroup. |
| `calibrated_risk` | Calibrated risk for each observation. |
| `density_ratio` | Density weight for each observation; 1 in the reference subgroup. |
| `validation_probability` | Fitted validation probabilities for AIPW estimation, before bounding; 1 for known noncases. |
| `diagnostics` | Subgroup sample sizes, validation counts, and probability truncation summaries. |

Use `summary()` to inspect a fit and `confint()` for a normal-approximation
confidence interval. Standard errors account for fitted validation, outcome,
calibration, and density-ratio models. Supplied sampling strata and
known-noncase allocations are treated as fixed.

## Simulation Example

The separate [`examples/`](examples/README.md) directory compares all five
methods. From the repository root:

```sh
Rscript examples/quickstart.R
Rscript examples/run_simulation.R --reps=10 --design=z2 --density=llogit
```

The Monte Carlo example uses 10 replicates of one scenario by default.
The [example documentation](examples/README.md) explains how to use balanced
sampling, quadratic logit density estimation, or the nine-scenario grid. These
examples generate fresh synthetic data and do not reproduce the manuscript's
archived random draws.

## Documentation

The complete function reference is available in R with `?fit_atpr`.
See [PROVENANCE.md](PROVENANCE.md) for the upstream dependency and research
source mapping, and [DEVELOPMENT.md](DEVELOPMENT.md) for package checks.

The risk-standardization framework follows Hegarty et al. (2025),
*Assessing Algorithm Fairness Requires Adjustment for Risk Distribution
Differences: Re-Considering the Equal Opportunity Criterion*,
doi: [10.1101/2025.01.31.25321489](https://doi.org/10.1101/2025.01.31.25321489).
Upstream attribution is provided in [`inst/NOTICE`](inst/NOTICE).

Maintainer: Ruofan Wang <ruofanwa@upenn.edu>. License: MIT.
