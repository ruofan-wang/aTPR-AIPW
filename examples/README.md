# Small simulation example

These scripts generate synthetic observations and call the installed `aTPRAIPW`
package. They are separate from the package estimator and do not read any EHR
files or previous simulation outputs. Run the commands below from the repository
root.

## Install and run one cohort

First install the package and its dependencies as described in the main README.
The estimator calls `fairRisk` calibration and risk-standardization functions
directly.

```sh
R CMD INSTALL .
Rscript examples/quickstart.R
```

`quickstart.R` compares complete-case (CC), contaminated-label (CL), AIPW-C,
AIPW-P, and full-data analyses for one generated cohort. Gold-standard outcomes
for unvalidated candidates are `NA` in the analysis data. Their simulated true
outcomes are kept separately and are supplied only to the full-data benchmark.
Both AIPW estimators fit validation probabilities from the observed candidates;
they do not use the probabilities that generated the validation sample.

## Run a small Monte Carlo example

```sh
Rscript examples/run_simulation.R
```

The default is **10 attempted replicates of one scenario**, with `a2 = 0`,
`b2 = 1.4`, Z2-dependent validation, and a linear logit density-ratio model. This
is a runnable example, not a precision study; 10 replicates are insufficient for
stable empirical SDs or 2.5th/97.5th percentiles. No 500-replicate run starts
unless requested explicitly.

The data-generating setup matches the revised simulation design:

- 4,000 observations in each of two subgroups, with latent risks drawn from
  Beta(2, 8) and Beta(4, 8), respectively.
- Risk scores satisfy `logit(R) = a_s + b_s * logit(pi*) + epsilon`, with
  `epsilon` SD 0.1, `(a1, b1) = (0, 1)`, and configurable `(a2, b2)`.
- Clinical auxiliaries have independent incremental noise variances 0.1, 1, and
  5 for Z1, Z2, and Z3, respectively.
- Phenotype sensitivity/specificity are (0.8, 0.9) in subgroup 1 and (0.9, 0.8)
  in subgroup 2.
- Exactly 400 true noncases are sampled without replacement **within each
  subgroup**. Validation is Bernoulli sampling with **expected size 400 within
  each subgroup**; its realized size is random.
- Z2-dependent sampling uses within-candidate-pool standardized Z2 and coefficient
  `gz = 0.5`. Balanced sampling uses near-equal expected allocation across risk
  tertile-by-phenotype cells, subject to probability bounds [0.03, 0.95].
- Validation models are fitted separately by subgroup: `~ Z1 + Z2 + Z3` for
  Z2-dependent sampling or `~ risk_stratum * phenotype` for balanced sampling.
  Outcome models use Z1, Z2, Z3, with phenotype added for AIPW-P.
- Calibration is linear in `logit(R)`. Density-ratio estimation uses
  `logit(calibrated risk)` for `llogit`, and adds its square for `qlogit`.
- The decision threshold is 0.2, subgroup 1 is the reference, and fitted
  validation probabilities have a lower bound of 0.005 for weighting.

The example uses fresh, explicit seeds. It illustrates the same statistical
setup but does **not** reproduce the archived manuscript's exact random draws
or printed numerical results. The same seed sequence is used across density
specifications so those comparisons use the same cohorts and validation samples.

## Change the design or density model

```sh
Rscript examples/run_simulation.R --reps=10 --design=balanced --density=qlogit
Rscript examples/run_simulation.R --reps=10 --a2=-0.75 --b2=1.8
Rscript examples/run_simulation.R --help
```

`--n` changes the subgroup size and scales the known-noncase count and expected
validation budget to 10% of that size. Model warnings and errors are recorded
in the replicate output. No failed replicate is replaced with a new sample.

For a longer study, these four **optional** commands run the nine combinations
of `a2 = -0.75, 0, 0.75` and `b2 = 1, 1.4, 1.8`, with 500 attempts per scenario:

```sh
Rscript examples/run_simulation.R --grid --reps=500 --design=z2 --density=llogit
Rscript examples/run_simulation.R --grid --reps=500 --design=z2 --density=qlogit
Rscript examples/run_simulation.R --grid --reps=500 --design=balanced --density=llogit
Rscript examples/run_simulation.R --grid --reps=500 --design=balanced --density=qlogit
```

## Outputs and interpretation

Each configuration writes to `examples/results/<design>_<density>/`, or a
directory selected with `--outdir=...`:

- `replicates.csv`: every attempted method fit, including errors and warnings.
- `summary.csv`: Monte Carlo means, paired benchmark differences (Bias),
  empirical SDs, mean IF SEs, empirical percentiles, and success counts.
- `failures.csv`: failed point estimates or AIPW standard errors; no failed
  replicate is replaced.
- `sampling.csv`: known-noncase counts and expected/realized validation counts
  for each subgroup and replicate.
- `estimates.pdf`: means and 2.5th--97.5th percentile intervals, with the full-data
  benchmark mean shown as a dashed line. These are sampling-distribution
  intervals, not confidence intervals for individual estimates.
- `configuration.R` and `session-info.txt`: analysis settings and R environment.

Empirical SD uses all successful point estimates. Mean IF SE uses successful
point-and-SE pairs. Ratio divides Mean IF SE by `if_subset_empirical_sd`, computed
on that same subset; when every IF fit succeeds this equals the main empirical
SD. Bias uses paired successful estimator and full-data fits. The explicit count
columns identify these subsets if a fit fails. IF SEs use the package's working
stacked influence function with contributions centered within each subgroup;
empirical risk cutpoints and known-noncase selection are treated as fixed.

Generated outputs are ignored by Git. They may be regenerated at any time.
