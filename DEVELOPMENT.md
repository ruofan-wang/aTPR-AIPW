# Development

Run from the repository root. The checked-in `man/` files and `NAMESPACE` allow
installation without roxygen2.

Install the local package and its dependencies before building or testing.
The `Remotes` field pins the tested `fairRisk` revision for dependency
installation and GitHub Actions.

```r
# If needed: install.packages("remotes")
remotes::install_local(".", dependencies = TRUE)
install.packages(c("testthat", "numDeriv", "roxygen2"))
roxygen2::roxygenise()
testthat::test_local()
```

From the parent directory:

```sh
R CMD build aTPR-AIPW
R CMD check --no-manual aTPRAIPW_0.1.0.tar.gz
```

The regression fixtures contain synthetic data only. Tests compare estimates
and influence-function standard errors with the archived research
implementation; they also check missing-label handling and IF derivatives.
GitHub Actions checks the package on Linux, macOS, and Windows.

The `examples/` directory is a companion simulation example and is excluded from
the installed package. Running it requires a clone or download of this
repository and an installed copy of aTPRAIPW. Generated example outputs are
ignored by Git.

The maintainer is Ruofan Wang <ruofanwa@upenn.edu>. No manuscript data or private
clinical data are needed to run the package tests or synthetic examples.
