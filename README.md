# BMA-BOIN-ET

Endpoint-specific Bayesian partition model averaging for efficacy-toxicity
basket dose finding. The R package is named **BMABOINET** (R package names
cannot contain hyphens).

## Install

```r
install.packages("remotes")
remotes::install_github("haohaostats/BMA-BOIN-ET")
library(BMABOINET)
```

## Reproduce the manuscript example

```r
design <- setup_design()
fit <- update_trial(design, example_manuscript_data())
select_od(fit)
```

`example_manuscript_data()` contains the deterministic marginal counts in
manuscript Table 7. The final BMA recommendations reproduce Table 8; these are
constructed illustrative data, not observations from a clinical trial.

## Prospective workflow

```r
design <- setup_design(
  basket_names = paste0("Basket ", 1:4), dose_levels = 1:5,
  target_tox = 0.30, efficacy_null = 0.30, efficacy_alt = 0.50,
  cohort_size = 3, max_sample = 24
)
get_boundaries(design)
fit <- update_trial(design, example_interim_data(), current_dose = rep(3, 4))
recommend_dose(fit)
borrowing_diagnostics(fit)
plot_borrowing(fit, endpoint = "toxicity", type = "coclustering", dose = 3)
select_od(fit)
```

Provide cumulative `basket`, `dose`, `n`, `tox`, and `eff` counts. At subsequent
reviews, pass `previous_analysis = fit` and the actual `current_dose` when it
differs from the previous recommendation. This preserves irreversible
elimination and closure. Archive the design before enrollment. Maximum sample
sizes are validated; stop enrollment when a basket reaches its planned total.

## Locked method

| Component | Implementation |
|:--|:--|
| Sharing | Exact partition averaging separately for each endpoint and dose |
| Prior | Beta(1,1); uniform number of groups, then uniform partitions within each group count |
| Permanent toxicity elimination | Local tail probability greater than 0.95, after at least 3 local patients; removes that dose and higher doses |
| Permanent futility elimination | Both local and BMA futility probabilities greater than 0.50, after at least 3 local patients |
| Allocation | Local or BMA toxicity action takes priority; BOIN-ET interval decisions and BMA efficacy ranking |
| Lowest-dose toxicity action | Closes the basket and deactivates all doses |
| Final selection | Active dose, at least 6 local patients, local and BMA toxicity tails at most 0.70, BMA futility at most 0.50 |
| Final ranking | Local-sample-weighted isotonic toxicity, utility 40 - 40 pT + 60 pE; ties use lower toxicity then lower dose |

For four baskets, all 15 partitions are enumerated. The basket set remains
fixed after closure. Exact enumeration scales with the Bell number and can
become expensive for larger basket counts. Posterior tails use the full beta
mixture, not a beta approximation fitted to its mean.

Diagnostics report partition probabilities and symmetric, dose-specific
co-clustering probabilities. They are not directional donor weights or BESS.

## Migration from 1.0.1

This is a major method change, not a renamed version of the old algorithm.
Use `library(BMABOINET)` and recreate saved design and analysis objects.
The public workflow function names remain familiar, but output classes and
diagnostics have changed. `plot_borrowing()` now accepts `"coclustering"` or
`"partitions"`; decision tables contain `bma_mean_tox` and `bma_mean_eff`.

The historical SA-BOIN-ETB package remains reproducible from its original tag:

```r
remotes::install_github("haohaostats/BMA-BOIN-ET@v1.0.1")
library(SAboinETB)
```

## Shiny application

The Shiny application is deployed separately from the public package.
`run_app()` opens the existing hosted URL; inspect the software name and version
in its audit panel before use. Installing this package does not update the
hosted application. The private application source is not included here.

## Validation

The test suite checks exact partition calculations, endpoint-specific learning,
ordered dose actions, irreversible closure, cumulative data validation, and
reproduction of the manuscript example. The package intentionally ships no function help pages or vignettes; all
usage documentation is maintained in this README. R CMD check therefore
reports missing documentation for exported functions.

## Development checks

Run `testthat::test_local()` from the package source directory. The release
was also compared with the locked manuscript simulation engine on 100 seeded
random states (posterior means, eliminations, next doses, and final selections).
The private Shiny server was exercised for design setup, example analysis,
both diagnostic displays, and the audit panel using BMABOINET 2.0.0.

Code contains no comments; this README is the sole maintained documentation.
The private Shiny deployment manifest must be regenerated before redeployment.

License: GPL-3.
