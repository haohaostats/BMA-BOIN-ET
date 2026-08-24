#' Set up a similarity-adaptive basket dose-optimization design
#'
#' @param basket_names Character vector of basket labels.
#' @param dose_levels Numeric or character vector of ordered dose labels.
#' @param target_tox Target toxicity probability.
#' @param efficacy_null Clinically uninteresting efficacy probability.
#' @param efficacy_alt Clinically desirable efficacy probability.
#' @param cohort_size Planned cohort size.
#' @param max_sample Maximum sample size in each basket.
#' @param toxicity_similarity Endpoint-specific prior basket-similarity
#'   probability or symmetric matrix for toxicity. A scalar of `0.5` is
#'   uninformative.
#' @param efficacy_similarity Endpoint-specific prior basket-similarity
#'   probability or symmetric matrix for efficacy. A scalar of `0.5` is
#'   uninformative.
#' @details The public interface deliberately accepts only clinically interpretable
#'   trial inputs. BOIN reference probabilities, beta priors, borrowing limits,
#'   stopping cutoffs, reliability discounts, and utility scores are supplied by
#'   the versioned standard implementation profile and stored in the returned
#'   object for review and audit.
#' @return An object of class `SAboinETB_design`.
#' @export
setup_design <- function(
    basket_names = paste0("Basket ", seq_len(4)),
    dose_levels = seq_len(5),
    target_tox = 0.30,
    efficacy_null = 0.30,
    efficacy_alt = 0.50,
    cohort_size = 3L,
    max_sample = rep(24L, length(basket_names)),
    toxicity_similarity = 0.50,
    efficacy_similarity = 0.50) {
  K <- length(basket_names)
  J <- length(dose_levels)
  ## Versioned standard implementation profile. These are intentionally not
  ## user-facing arguments: changing them defines a different calibrated design.
  tox_low <- 0.60 * target_tox
  tox_high <- 1.40 * target_tox
  prior_tox <- prior_eff <- c(1, 1)
  bess_cap_t <- 12
  bess_cap_e <- 18
  bess_ratio_t <- 1.0
  bess_ratio_e <- 1.5
  escalation_gate_cutoff <- 0.80
  toxicity_elimination_cutoff <- 0.95
  efficacy_elimination_cutoff <- 0.95
  local_toxicity_selection_cutoff <- 0.95
  working_toxicity_selection_cutoff <- 0.95
  working_efficacy_selection_cutoff <- 0.85
  min_eliminate <- 3L
  min_select <- 6L
  utility <- c(100, 40, 60, 0)
  if (K < 2L || J < 2L) stop("At least two baskets and two ordered doses are required.")
  if (anyDuplicated(basket_names) || anyDuplicated(dose_levels)) {
    stop("Basket names and dose levels must be unique.")
  }
  probs <- c(target_tox, tox_low, tox_high, efficacy_null, efficacy_alt,
             escalation_gate_cutoff, toxicity_elimination_cutoff,
             efficacy_elimination_cutoff,
             local_toxicity_selection_cutoff,
             working_toxicity_selection_cutoff,
             working_efficacy_selection_cutoff)
  if (any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("Probability inputs must lie strictly between zero and one.")
  }
  if (!(tox_low < target_tox && target_tox < tox_high)) {
    stop("Require tox_low < target_tox < tox_high.")
  }
  if (efficacy_null >= efficacy_alt) stop("efficacy_null must be below efficacy_alt.")
  if (length(max_sample) == 1L) max_sample <- rep(max_sample, K)
  if (length(max_sample) != K || any(max_sample < cohort_size)) {
    stop("max_sample must provide one valid value per basket.")
  }
  make_similarity <- function(x, label) {
    if (length(x) == 1L) {
      if (!is.finite(x) || x <= 0 || x >= 1) stop(label, " must be in (0, 1).")
      out <- matrix(x, K, K)
      diag(out) <- 1
      return(out)
    }
    out <- as.matrix(x)
    if (!identical(dim(out), c(K, K)) || any(!is.finite(out)) ||
        any(out < 0 | out > 1) || max(abs(out - t(out))) > 1e-10) {
      stop(label, " must be a symmetric K by K matrix with entries in [0, 1].")
    }
    diag(out) <- 1
    out
  }
  design <- list(
    K = K, J = J, basket_names = as.character(basket_names),
    dose_levels = dose_levels, target_tox = target_tox, tox_low = tox_low,
    tox_high = tox_high, efficacy_null = efficacy_null,
    efficacy_alt = efficacy_alt, cohort_size = as.integer(cohort_size),
    max_sample = as.integer(max_sample), prior_tox = prior_tox,
    prior_eff = prior_eff,
    prior_similarity_t = make_similarity(toxicity_similarity, "toxicity_similarity"),
    prior_similarity_e = make_similarity(efficacy_similarity, "efficacy_similarity"),
    bess_cap_t = bess_cap_t, bess_cap_e = bess_cap_e,
    bess_ratio_t = bess_ratio_t, bess_ratio_e = bess_ratio_e,
    escalation_gate_cutoff = escalation_gate_cutoff,
    toxicity_elimination_cutoff = toxicity_elimination_cutoff,
    efficacy_elimination_cutoff = efficacy_elimination_cutoff,
    local_toxicity_selection_cutoff = local_toxicity_selection_cutoff,
    working_toxicity_selection_cutoff = working_toxicity_selection_cutoff,
    working_efficacy_selection_cutoff = working_efficacy_selection_cutoff,
    # Compatibility aliases retained for early development design records.
    safety_cutoff = toxicity_elimination_cutoff,
    activity_cutoff = working_efficacy_selection_cutoff,
    min_eliminate = as.integer(min_eliminate), min_select = as.integer(min_select),
    min_escalate = 3L,
    utility = utility, curve_saturation = 3, min_borrow = 3L,
    conflict_scale_t = 2.0, conflict_scale_e = 3.0,
    prior_odds_cap = stats::qlogis(0.95), prior_conflict_scale = 2.0,
    exploration_toxicity_cutoff = 0.24,
    exploration_toxicity_tail_cutoff = 1.00,
    exploration_efficacy_cutoff = 0.55,
    implementation_profile = "SA-BOIN-ETB-1.0",
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = "1.0.0"
  )
  class(design) <- "SAboinETB_design"
  design
}

#' Print an SA-BOIN-ETB object
#'
#' Compact print methods for locked designs and completed interim analyses.
#'
#' @param x An `SAboinETB_design` or `SAboinETB_analysis` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
#' @rdname print.SAboinETB
#' @export
print.SAboinETB_design <- function(x, ...) {
  cat("SA-BOIN-ETB trial design\n")
  cat("  Baskets:", x$K, " | Doses:", x$J, " | Cohort size:", x$cohort_size, "\n")
  cat("  Target toxicity:", format(x$target_tox, digits = 3),
      " | Efficacy null:", format(x$efficacy_null, digits = 3), "\n")
  cat("  Incoming BESS caps (toxicity/efficacy):",
      x$bess_cap_t, "/", x$bess_cap_e, "\n")
  cat("  Frontier gate / elimination cutoffs:",
      x$escalation_gate_cutoff, "/", x$toxicity_elimination_cutoff, "\n")
  invisible(x)
}

#' Obtain the locked interval boundaries
#' @param design A design returned by [setup_design()].
#' @return Named numeric vector containing escalation, de-escalation, and
#'   efficacy boundaries.
#' @export
get_boundaries <- function(design) {
  stopifnot(inherits(design, "SAboinETB_design"))
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
