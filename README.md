# BMA-BOIN-ET

Endpoint-specific Bayesian partition model averaging for efficacy-toxicity
basket dose finding.

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

## Shiny application

The Shiny application is deployed separately from the public package.
`run_app()` opens the existing hosted URL; inspect the software name and version
in its audit panel before use. Installing this package does not update the
hosted application. The private application source is not included here.
