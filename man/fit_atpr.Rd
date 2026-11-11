\name{fit_atpr}
\alias{fit_atpr}
\title{Estimate a risk-standardized true positive rate difference}
\usage{
fit_atpr(
  data,
  risk,
  group,
  outcome,
  candidate,
  validated,
  validation_formula = ~Z1 + Z2 + Z3,
  outcome_formula = ~Z1 + Z2 + Z3,
  phenotype = NULL,
  method = c("AIPW-C", "AIPW-P", "complete", "contaminated", "full"),
  reference = NULL,
  tau = 0.2,
  density = c("llogit", "qlogit"),
  prob_floor = 0.005,
  risk_floor = 1e-04,
  se = TRUE,
  variance = c("stratified", "iid")
)
}
\arguments{
\item{data}{A data frame containing one row per cohort member.}

\item{risk, group, outcome, candidate, validated}{Column names, supplied as strings.
\code{risk} is the original prediction score in \verb{(0, 1)}; \code{group} contains exactly
two subgroup labels; \code{outcome} is the binary gold-standard outcome;
\code{candidate} is 1 for the candidate pool and 0 for known noncases; and
\code{validated} is 1 for candidates with observed gold-standard outcomes.
Outcomes of unvalidated candidates may be \code{NA} and are not used by AIPW or
complete-case estimation. Known noncases must have outcome 0 or \code{NA}, and
validation indicator 0.}

\item{validation_formula}{One-sided formula for the logistic validation model,
fitted separately in each subgroup using all candidates.}

\item{outcome_formula}{One-sided formula for the candidate outcome model,
fitted separately in each subgroup using only validated candidates.
Formula predictors must be observed for all candidates. Offsets are not
supported. For balanced sampling, precompute risk strata separately within
each subgroup's candidate pool and use, for example, \code{~ risk_bin * phenotype}.}

\item{phenotype}{Name of the binary algorithmic-phenotype column. Required for
\code{method = "AIPW-P"} or \code{"contaminated"}. AIPW-P adds its main effect to
\code{outcome_formula}; it is not added to the validation formula automatically.}

\item{method}{Estimation method. \code{"AIPW-C"} uses the supplied outcome formula;
\code{"AIPW-P"} additionally uses the phenotype; \code{"complete"} uses known noncases
and validated candidates at every stage; \code{"contaminated"} uses the phenotype
as the outcome throughout the cohort; and \code{"full"} requires observed true
outcomes for every cohort member.}

\item{reference}{Label of the reference subgroup. Defaults to the first
subgroup appearing in \code{data}; specifying it explicitly is recommended.}

\item{tau}{Threshold applied to the original score: a positive prediction is
defined by \code{risk > tau}.}

\item{density}{Density-ratio specification. \code{"llogit"} uses a linear term in
the logit of the calibrated risk; \code{"qlogit"} additionally uses its square.}

\item{prob_floor}{Lower bound on fitted validation probabilities used in the
AIPW residual correction. The default 0.005 limits inverse weights to 200.
Set to 0 to disable this bound.}

\item{risk_floor}{Numerical bound applied to the original score before its
logit is used in AIPW and complete-case calibration. Full-data and
contaminated-label calibration use the unbounded logit, matching the
original simulation implementation. Set to 0 to disable bounding.}

\item{se}{Calculate the working influence-function standard error for AIPW.
Standard errors are not implemented for the three comparator methods, for
which \code{std.error} is always \code{NA}.}

\item{variance}{\code{"stratified"} centers influence contributions within each
subgroup, treating subgroup sample sizes as fixed; \code{"iid"} uses uncentered
contributions under independent identically distributed sampling.}
}
\value{
An object of class \code{atpr_fit} with \code{estimate}, \code{std.error}, subgroup
\code{atpr}, per-row \code{calibrated_risk}, \code{density_ratio}, and
\code{validation_probability}, subgroup \code{diagnostics}, and \code{diagnostics_details}.
Per-row outputs follow the original data order; excluded complete cases
have \code{NA}. The \code{engine} component contains the parameter vector, score
function, sensitivity matrix, and gradient used for AIPW inference and is
intended for numerical audits. Methods include \code{\link[=print]{print()}}, \code{\link[=summary]{summary()}},
\code{\link[stats:coef]{stats::coef()}}, \code{\link[stats:vcov]{stats::vcov()}}, and \code{\link[stats:confint]{stats::confint()}}.
}
\description{
Fit covariate-based or phenotype-assisted augmented inverse probability
weighted (AIPW) estimators using a partially labeled cohort. Complete-case,
contaminated-label, and full-data point estimators are also available.
}
\details{
Risk calibration is logistic and linear in the logit of the original score.
Density estimation models reference-group membership using the chosen basis
in the logit of the calibrated risk. Its odds, multiplied by the comparison
to reference sample-size ratio, estimate the reference density divided by the
comparison density. Each subgroup's aTPR averages \code{I(risk > tau)} weighted by
calibrated risk and density ratio; the reported contrast is comparison minus
reference.

The package extends \code{fairRisk}: ordinary calibration uses
\code{\link[fairRisk:calibrate]{fairRisk::calibrate()}}, density ratios use \code{\link[fairRisk:estDensityRatio]{fairRisk::estDensityRatio()}}
with \code{transform = TRUE}, and aTPRs use \code{\link[fairRisk:get_aTPR]{fairRisk::get_aTPR()}}. AIPW
pseudo-outcomes and augmented calibration are implemented here because
pseudo-outcomes can lie outside the binary-response range. For inference,
the equivalent density model is also fitted in a raw polynomial basis to
obtain coefficients not returned by the upstream density function.

The stacked working influence function includes fitted validation, outcome,
calibration, and density-ratio parameters. Generated strata, known-noncase
selection, and other supplied design features are treated as fixed. It does
not add a separate correction for fixed-count noncase sampling or adaptive
validation designs. Fitted probability bounding is included in the local
derivatives away from the bound.
}
