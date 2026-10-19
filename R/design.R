setup_design <- function(
    basket_names = paste0("Basket ", seq_len(4)),
    dose_levels = seq_len(5),
    target_tox = 0.30,
    efficacy_null = 0.30,
    efficacy_alt = 0.50,
    cohort_size = 3L,
    max_sample = rep(24L, length(basket_names))) {
  if (anyNA(basket_names) || anyNA(dose_levels) ||
      any(!nzchar(as.character(basket_names))) ||
      any(!nzchar(as.character(dose_levels)))) stop("Labels must be non-missing and non-empty.")
  if (!is.numeric(cohort_size) || length(cohort_size) != 1L ||
      !is.finite(cohort_size) || cohort_size < 1 || cohort_size != floor(cohort_size)) {
    stop("cohort_size must be a positive integer.")
  }
  if (!is.numeric(max_sample) || !length(max_sample) ||
      any(!is.finite(max_sample)) || any(max_sample != floor(max_sample))) {
    stop("max_sample must contain finite integer sample sizes.")
  }
  for (value in list(target_tox, efficacy_null, efficacy_alt)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("Clinical probabilities must be finite numeric scalars.")
    }
  }
  K <- length(basket_names)
  J <- length(dose_levels)
  tox_low <- 0.60 * target_tox
  tox_high <- 1.40 * target_tox
  prior_tox <- prior_eff <- c(1, 1)
  toxicity_elimination_cutoff <- 0.95
  efficacy_elimination_cutoff <- 0.50
  toxicity_selection_cutoff <- 0.70
  min_eliminate_t <- 3L
  min_eliminate_e <- 3L
  min_select <- 6L
  utility <- c(100, 40, 60, 0)

  if (K < 2L || J < 2L) stop("At least two baskets and two ordered doses are required.")
  if (anyDuplicated(basket_names) || anyDuplicated(dose_levels)) {
    stop("Basket names and dose levels must be unique.")
  }
  probs <- c(target_tox, tox_low, tox_high, efficacy_null, efficacy_alt,
             toxicity_elimination_cutoff, efficacy_elimination_cutoff,
             toxicity_selection_cutoff)
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("Probability inputs must lie strictly between zero and one.")
  }
  if (!(tox_low < target_tox && target_tox < tox_high)) {
    stop("Require tox_low < target_tox < tox_high.")
  }
  if (efficacy_null >= efficacy_alt) {
    stop("efficacy_null must be below efficacy_alt.")
  }
  if (length(max_sample) == 1L) max_sample <- rep(max_sample, K)
  if (length(max_sample) != K || any(max_sample < cohort_size)) {
    stop("max_sample must provide one valid value per basket.")
  }

  partitions <- .set_partitions(K)
  group_counts <- vapply(partitions, max, integer(1))
  stirling_counts <- tabulate(group_counts, nbins = K)
  partition_labels <- vapply(partitions, paste, collapse = "", character(1))
  design <- list(
    K = K, J = J,
    basket_names = as.character(basket_names), dose_levels = dose_levels,
    target_tox = target_tox, tox_low = tox_low, tox_high = tox_high,
    efficacy_null = efficacy_null, efficacy_alt = efficacy_alt,
    cohort_size = as.integer(cohort_size), max_sample = as.integer(max_sample),
    prior_tox = prior_tox, prior_eff = prior_eff,
    partitions = partitions, partition_labels = partition_labels,
    stirling_counts = stirling_counts,
    partition_prior = "uniform number of groups, uniform partitions within group count",
    toxicity_elimination_cutoff = toxicity_elimination_cutoff,
    efficacy_elimination_cutoff = efficacy_elimination_cutoff,
    toxicity_selection_cutoff = toxicity_selection_cutoff,
    min_eliminate_t = min_eliminate_t, min_eliminate_e = min_eliminate_e,
    min_select = min_select, utility = utility,
    implementation_profile = "BMA-BOIN-ET-2.0",
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = "2.0.0"
  )
  class(design) <- "BMABOINET_design"
  design
}

print.BMABOINET_design <- function(x, ...) {
  cat("BMA-BOIN-ET Bayesian partition-model-averaged design\n")
  cat("  Baskets:", x$K, " | Doses:", x$J,
      " | Partitions per endpoint-dose:", length(x$partitions), "\n")
  cat("  Target toxicity:", format(x$target_tox, digits = 3),
      " | Efficacy null:", format(x$efficacy_null, digits = 3), "\n")
  cat("  Partition prior:", x$partition_prior, "\n")
  invisible(x)
}

get_boundaries <- function(design) {
  stopifnot(inherits(design, "BMABOINET_design"))
  lambda_e <- log((1 - design$tox_low) / (1 - design$target_tox)) /
    log(design$target_tox * (1 - design$tox_low) /
          (design$tox_low * (1 - design$target_tox)))
  lambda_d <- log((1 - design$target_tox) / (1 - design$tox_high)) /
    log(design$tox_high * (1 - design$target_tox) /
          (design$target_tox * (1 - design$tox_high)))
  lambda_eff <- log((1 - design$efficacy_alt) / (1 - design$efficacy_null)) /
    log(design$efficacy_null * (1 - design$efficacy_alt) /
          (design$efficacy_alt * (1 - design$efficacy_null)))
  c(escalate = lambda_e, deescalate = lambda_d, efficacy = lambda_eff)
}
