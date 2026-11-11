# Run from the repository root after installing aTPRAIPW:
# Rscript examples/quickstart.R
library(aTPRAIPW)
source("examples/simulation_helpers.R")

simulation <- simulate_cohort(seed = 20261110L, design = "z2")
print(simulation$sampling, row.names = FALSE)

results <- fit_simulation_methods(simulation, density = "llogit")
print(results$estimates, row.names = FALSE)

