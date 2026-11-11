#!/usr/bin/env Rscript
# Run from the repository root. Only a small example runs by default.
library(aTPRAIPW)
source("examples/simulation_helpers.R")

arguments <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% arguments) {
  cat(paste(
    "Usage: Rscript examples/run_simulation.R [options]",
    "  --reps=10          Number of attempted replicates per scenario",
    "  --design=z2        z2 or balanced validation",
    "  --density=llogit   llogit or qlogit density model",
    "  --seed=20261110    Base simulation seed",
    "  --n=4000          Cohort size per subgroup",
    "  --a2=0 --b2=1.4   Single-scenario parameters",
    "  --grid             Run all nine manuscript parameter combinations",
    "  --outdir=PATH      Output folder (default: examples/results/DESIGN_DENSITY)",
    sep = "\n"
  ), "\n")
  quit(status = 0L)
}
get_option <- function(name, default) {
  values <- arguments[startsWith(arguments, paste0("--", name, "="))]
  if (!length(values)) return(default)
  sub(paste0("^--", name, "="), "", values)
}
integer_option <- function(name, default) {
  as.integer(get_option(name, as.character(default)))
}
replicates <- integer_option("reps", 10L)
n_per_group <- integer_option("n", 4000L)
base_seed <- integer_option("seed", 20261110L)
design <- match.arg(get_option("design", "z2"), c("z2", "balanced"))
density <- match.arg(get_option("density", "llogit"), c("llogit", "qlogit"))
outdir <- get_option("outdir", file.path("examples", "results", paste(design, density, sep = "_")))
if ("--grid" %in% arguments) {
  scenarios <- expand.grid(b2 = c(1, 1.4, 1.8), a2 = c(-0.75, 0, 0.75))
} else {
  scenarios <- data.frame(
    b2 = as.numeric(get_option("b2", "1.4")),
    a2 = as.numeric(get_option("a2", "0"))
  )
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
raw_rows <- sampling_rows <- list()
counter <- 0L
for (scenario in seq_len(nrow(scenarios))) {
  a2 <- scenarios$a2[scenario]
  b2 <- scenarios$b2[scenario]
  for (replicate in seq_len(replicates)) {
    # Shared seeds permit paired comparisons when rerunning another density model.
    seed <- base_seed + replicate - 1L
    simulation <- simulate_cohort(seed = seed, n_per_group = n_per_group,
                                  a2 = a2, b2 = b2, design = design)
    results <- fit_simulation_methods(simulation, density = density)
    counter <- counter + 1L
    raw_rows[[counter]] <- cbind(
      data.frame(a2 = a2, b2 = b2, replicate = replicate, seed = seed,
                 design = design, density = density), results$estimates
    )
    sampling_rows[[counter]] <- cbind(
      data.frame(a2 = a2, b2 = b2, replicate = replicate), simulation$sampling
    )
  }
  cat(sprintf("Completed (%g, %g): %d attempted replicates.\n", a2, b2, replicates))
}
raw <- do.call(rbind, raw_rows)
summary <- summarize_replicates(raw)
failures <- raw[!raw$point_success |
                 (raw$method %in% c("AIPW-C", "AIPW-P") & !raw$if_success), ]
write.csv(raw, file.path(outdir, "replicates.csv"), row.names = FALSE)
write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)
write.csv(failures, file.path(outdir, "failures.csv"), row.names = FALSE)
write.csv(do.call(rbind, sampling_rows), file.path(outdir, "sampling.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(outdir, "session-info.txt"))
configuration <- list(
  design = design, gz = if (design == "z2") 0.5 else NA_real_, density = density,
  tau = 0.2, replicates = replicates, base_seed = base_seed,
  n_per_group = n_per_group, known_noncase_fraction = 0.1,
  expected_validation_fraction = 0.1, e1_variance = 0.1,
  probability_floor = 0.005, scenarios = scenarios,
  seed_policy = "Fresh illustrative seeds; not the archived manuscript random draws."
)
dput(configuration, file = file.path(outdir, "configuration.R"))
plot_simulation_summary(summary, file.path(outdir, "estimates.pdf"))
print(summary, row.names = FALSE)
cat("Outputs:", normalizePath(outdir), "\n")
