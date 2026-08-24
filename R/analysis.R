#' Validate cumulative interim data
#'
#' @param data Data frame with columns `basket`, `dose`, `n`, `tox`, and `eff`.
#' @param design A design returned by [setup_design()].
#' @return A validated and ordered data frame.
#' @export
validate_interim_data <- function(data, design) {
  stopifnot(inherits(design, "SAboinETB_design"))
  required <- c("basket", "dose", "n", "tox", "eff")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("Interim data must contain basket, dose, n, tox, and eff columns.")
  }
  x <- data[required]
  if (anyDuplicated(x[c("basket", "dose")])) {
    stop("Each basket-dose combination must appear at most once.")
  }
  if (any(!x$basket %in% design$basket_names)) stop("Unknown basket label detected.")
  if (any(!x$dose %in% design$dose_levels)) stop("Unknown dose label detected.")
  for (nm in c("n", "tox", "eff")) {
    if (any(!is.finite(x[[nm]])) || any(x[[nm]] < 0) ||
        any(abs(x[[nm]] - round(x[[nm]])) > 1e-8)) {
      stop(nm, " must contain non-negative integer counts.")
    }
    x[[nm]] <- as.integer(x[[nm]])
  }
  if (any(x$tox > x$n) || any(x$eff > x$n)) {
    stop("Toxicity and efficacy events cannot exceed the treated sample size.")
  }
  totals <- tapply(x$n, x$basket, sum)
  if (any(totals > design$max_sample[match(names(totals), design$basket_names)])) {
    stop("At least one basket exceeds its prespecified maximum sample size.")
  }
  x <- x[order(match(x$basket, design$basket_names),
               match(x$dose, design$dose_levels)), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.as_matrices <- function(data, design) {
  dn <- list(design$basket_names, as.character(design$dose_levels))
  n <- tox <- eff <- matrix(0, design$K, design$J, dimnames = dn)
  for (i in seq_len(nrow(data))) {
    k <- match(data$basket[i], design$basket_names)
    j <- match(data$dose[i], design$dose_levels)
    n[k, j] <- data$n[i]
    tox[k, j] <- data$tox[i]
    eff[k, j] <- data$eff[i]
  }
  list(n = n, tox = tox, eff = eff)
}

.expected_utility <- function(p_t, p_e, utility) {
  u01 <- utility[1]; u00 <- utility[2]
  u11 <- utility[3]; u10 <- utility[4]
  u00 + (u10 - u00) * p_t + (u01 - u00) * p_e
}

.weighted_pava <- function(x, w) {
  stopifnot(length(x) == length(w), all(w > 0))
  level <- as.numeric(x); weight <- as.numeric(w)
  start <- end <- seq_along(x)
  m <- length(level); i <- 1L
  while (i < m) {
    if (level[i] <= level[i + 1L] + 1e-14) {
      i <- i + 1L
    } else {
      new_w <- weight[i] + weight[i + 1L]
      level[i] <- (weight[i] * level[i] +
        weight[i + 1L] * level[i + 1L]) / new_w
      weight[i] <- new_w
      end[i] <- end[i + 1L]
      if (i + 1L < m) {
        level[(i + 1L):(m - 1L)] <- level[(i + 2L):m]
        weight[(i + 1L):(m - 1L)] <- weight[(i + 2L):m]
        start[(i + 1L):(m - 1L)] <- start[(i + 2L):m]
        end[(i + 1L):(m - 1L)] <- end[(i + 2L):m]
      }
      m <- m - 1L
      level <- level[seq_len(m)]; weight <- weight[seq_len(m)]
      start <- start[seq_len(m)]; end <- end[seq_len(m)]
      if (i > 1L) i <- i - 1L
    }
  }
  out <- numeric(length(x))
  for (g in seq_len(m)) out[start[g]:end[g]] <- level[g]
  out
}

.working_rate <- function(post, k, j) {
  if (post$n[k, j] <= 0) return(0)
  post$y[k, j] / post$n[k, j]
}

.current_indices <- function(design, mat, current_dose, previous_analysis) {
  if (is.null(current_dose) && !is.null(previous_analysis)) {
    previous <- previous_analysis$decisions$recommended_dose
    names(previous) <- previous_analysis$decisions$basket
    current_dose <- previous[design$basket_names]
    missing <- is.na(current_dose)
    if (any(missing)) {
      fallback <- previous_analysis$decisions$current_dose
      names(fallback) <- previous_analysis$decisions$basket
      current_dose[missing] <- fallback[design$basket_names][missing]
    }
  }
  if (is.null(current_dose)) {
    return(vapply(seq_len(design$K), function(k) {
      tried <- which(mat$n[k, ] > 0)
      if (length(tried)) max(tried) else 1L
    }, integer(1)))
  }
  if (is.null(names(current_dose))) {
    if (length(current_dose) != design$K) {
      stop("current_dose must have one value per basket.")
    }
    names(current_dose) <- design$basket_names
  }
  current_dose <- current_dose[design$basket_names]
  out <- match(current_dose, design$dose_levels)
  if (anyNA(out)) stop("current_dose contains missing or unknown dose labels.")
  as.integer(out)
}

.initial_active <- function(design, previous_analysis) {
  if (is.null(previous_analysis)) {
    return(matrix(TRUE, design$K, design$J,
                  dimnames = list(design$basket_names,
                                  as.character(design$dose_levels))))
  }
  if (!inherits(previous_analysis, "SAboinETB_analysis")) {
    stop("previous_analysis must be returned by update_trial().")
  }
  previous_design <- previous_analysis$design
  if (!identical(previous_design$basket_names, design$basket_names) ||
      !identical(previous_design$dose_levels, design$dose_levels)) {
    stop("previous_analysis does not use the locked design.")
  }
  previous_analysis$active
}

.apply_elimination <- function(active, mat, bt, be, design) {
  open <- rowSums(active) > 0
  for (k in seq_len(design$K)) {
    if (!open[k]) next
    for (j in seq_len(design$J)) {
      if (!active[k, j]) next
      local_t_bad <- mat$n[k, j] >= design$min_eliminate &&
        .beta_prob_gt(mat$tox[k, j], mat$n[k, j], design$target_tox,
                      design$prior_tox) > design$toxicity_elimination_cutoff
      if (local_t_bad) {
        active[k, j:design$J] <- FALSE
        next
      }
      local_e_bad <- mat$n[k, j] >= design$min_eliminate &&
        .beta_prob_lt(mat$eff[k, j], mat$n[k, j], design$efficacy_null,
                      design$prior_eff) > design$efficacy_elimination_cutoff
      working_e_bad <- be$n[k, j] >= design$min_eliminate &&
        .beta_prob_lt(be$y[k, j], be$n[k, j], design$efficacy_null,
                      design$prior_eff) > design$efficacy_elimination_cutoff
      if (local_e_bad && working_e_bad) active[k, j] <- FALSE
    }
    open[k] <- any(active[k, ])
  }
  list(active = active, open = open)
}

.local_escalation_ok <- function(mat, k, j, design) {
  mat$n[k, j] >= design$min_escalate &&
    .beta_prob_gt(mat$tox[k, j], mat$n[k, j], design$target_tox,
                  design$prior_tox) <= design$escalation_gate_cutoff
}

.choose_next_dose <- function(k, current, active, mat, bt, be, design) {
  j <- current[k]
  p_t <- .working_rate(bt, k, j)
  bounds <- get_boundaries(design)
  tried_doses <- which(mat$n[k, ] > 0)
  frontier <- if (length(tried_doses)) max(tried_doses) else j

  if (p_t > bounds["deescalate"]) {
    active_lower <- which(active[k, ] & seq_len(design$J) < j)
    if (length(active_lower)) {
      return(list(
        index = as.integer(max(active_lower)), frontier = frontier,
        rule = paste(
          "Working toxicity crossed the de-escalation boundary;",
          "selected the highest active lower dose"
        )
      ))
    }
    if (isTRUE(active[k, j])) {
      return(list(
        index = as.integer(j), frontier = frontier,
        rule = paste(
          "Working toxicity crossed the de-escalation boundary, but the",
          "current dose is the lowest active dose"
        )
      ))
    }
    return(list(
      index = NA_integer_, frontier = frontier,
      rule = paste(
        "Working toxicity crossed the de-escalation boundary and no",
        "active lower dose remains; close the basket"
      )
    ))
  }

  frontier_gate <- .local_escalation_ok(mat, k, frontier, design)
  if (frontier < design$J && active[k, frontier + 1L] && frontier_gate) {
    frontier_t <- .working_rate(bt, k, frontier)
    frontier_e <- .working_rate(be, k, frontier)
    frontier_t_tail <- .beta_prob_gt(
      bt$y[k, frontier], bt$n[k, frontier], design$target_tox,
      design$prior_tox
    )
    if (frontier_t <= design$exploration_toxicity_cutoff &&
        frontier_t_tail <= design$exploration_toxicity_tail_cutoff &&
        frontier_e <= design$exploration_efficacy_cutoff) {
      return(list(index = as.integer(frontier + 1L), frontier = frontier,
                  rule = "Safety-gated extension of the exploration frontier"))
    }
  }

  tried <- which(active[k, ] & mat$n[k, ] > 0)
  if (!length(tried)) {
    return(list(index = NA_integer_, frontier = frontier,
                rule = "No active previously tried dose; close the basket"))
  }
  locally_safe <- vapply(tried, function(x) {
    .beta_prob_gt(mat$tox[k, x], mat$n[k, x], design$target_tox,
                  design$prior_tox) <=
      design$local_toxicity_selection_cutoff
  }, logical(1))
  candidates <- tried[locally_safe]
  if (!length(candidates)) {
    return(list(
      index = NA_integer_, frontier = frontier,
      rule = paste(
        "No locally toxicity-admissible active previously tried dose;",
        "close the basket"
      )
    ))
  }
  pt <- vapply(candidates, function(x) {
    .beta_mean(bt$y[k, x], bt$n[k, x], design$prior_tox)
  }, numeric(1))
  pe <- vapply(candidates, function(x) {
    .beta_mean(be$y[k, x], be$n[k, x], design$prior_eff)
  }, numeric(1))
  utility <- .expected_utility(pt, pe, design$utility)
  best <- candidates[utility == max(utility)]
  if (length(best) > 1L) {
    best_pt <- pt[candidates %in% best]
    best <- best[best_pt == min(best_pt)]
  }
  list(index = as.integer(min(best)), frontier = frontier,
       rule = "Utility selection among locally toxicity-admissible tried doses")
}

#' Run a cumulative interim analysis
#'
#' @param design A design returned by [setup_design()].
#' @param data Cumulative basket-dose counts.
#' @param current_dose Optional named vector giving the dose currently assigned
#'   in each basket. If omitted, the preceding recommendation is used when
#'   `previous_analysis` is supplied; otherwise the highest observed dose is
#'   used.
#' @param previous_analysis Optional preceding analysis. Supplying it preserves
#'   irreversible dose eliminations and the prospective dose trajectory.
#' @details Recommendations are restricted to the post-elimination active set.
#'   Immediate de-escalation selects the highest active dose below the current
#'   dose. If no active lower dose exists, the basket stays only when the
#'   current dose remains active; otherwise it closes. When frontier extension
#'   is unavailable, utility is maximized only over previously tried, active,
#'   locally toxicity-admissible doses. An empty candidate set closes the
#'   basket.
#' @return An object of class `SAboinETB_analysis`.
#' @export
update_trial <- function(design, data, current_dose = NULL,
                         previous_analysis = NULL) {
  data <- validate_interim_data(data, design)
  mat <- .as_matrices(data, design)
  bt <- .borrow_endpoint(mat$tox, mat$n, "toxicity", design)
  be <- .borrow_endpoint(mat$eff, mat$n, "efficacy", design)
  current <- .current_indices(design, mat, current_dose, previous_analysis)
  active <- .initial_active(design, previous_analysis)
  elimination <- .apply_elimination(active, mat, bt, be, design)
  active <- elimination$active; open <- elimination$open

  decisions <- vector("list", design$K)
  for (k in seq_len(design$K)) {
    j <- current[k]
    tried <- which(mat$n[k, ] > 0)
    frontier <- if (length(tried)) max(tried) else j
    if (!open[k]) {
      target <- NA_integer_
      rule <- "All doses are inactive; close the basket"
    } else {
      choice <- .choose_next_dose(k, current, active, mat, bt, be, design)
      target <- choice$index; frontier <- choice$frontier; rule <- choice$rule
    }
    action <- if (is.na(target)) "Close basket" else if (target > j) {
      "Escalate"
    } else if (target < j) "De-escalate" else "Stay"
    local_tox_tail <- .beta_prob_gt(
      mat$tox[k, j], mat$n[k, j], design$target_tox, design$prior_tox)
    working_tox_tail <- .beta_prob_gt(
      bt$y[k, j], bt$n[k, j], design$target_tox, design$prior_tox)
    local_futility <- .beta_prob_lt(
      mat$eff[k, j], mat$n[k, j], design$efficacy_null, design$prior_eff)
    working_futility <- .beta_prob_lt(
      be$y[k, j], be$n[k, j], design$efficacy_null, design$prior_eff)
    decisions[[k]] <- data.frame(
      basket = design$basket_names[k],
      current_dose = design$dose_levels[j],
      frontier_dose = design$dose_levels[frontier],
      local_n = mat$n[k, j], local_tox = mat$tox[k, j],
      local_eff = mat$eff[k, j],
      posterior_tox = .beta_mean(bt$y[k, j], bt$n[k, j], design$prior_tox),
      posterior_eff = .beta_mean(be$y[k, j], be$n[k, j], design$prior_eff),
      working_rate_tox = .working_rate(bt, k, j),
      working_rate_eff = .working_rate(be, k, j),
      toxicity_bess = bt$bess[k, j], efficacy_bess = be$bess[k, j],
      local_tox_tail = local_tox_tail,
      working_tox_tail = working_tox_tail,
      local_futility = local_futility,
      working_futility = working_futility,
      frontier_gate = .local_escalation_ok(mat, k, frontier, design),
      action = action,
      recommended_dose = if (is.na(target)) NA else design$dose_levels[target],
      reason = rule, stringsAsFactors = FALSE
    )
  }
  analysis <- list(
    design = design, data = data, matrices = mat,
    toxicity = bt, efficacy = be, active = active, open = open,
    current_index = current, decisions = do.call(rbind, decisions),
    boundaries = get_boundaries(design),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  class(analysis) <- "SAboinETB_analysis"
  analysis
}

#' @rdname print.SAboinETB
#' @export
print.SAboinETB_analysis <- function(x, ...) {
  cat("SA-BOIN-ETB interim analysis\n")
  cat("  Timestamp:", x$timestamp, "\n\n")
  print(x$decisions[c("basket", "current_dose", "frontier_dose", "local_n",
                      "action", "recommended_dose", "toxicity_bess",
                      "efficacy_bess")], row.names = FALSE)
  invisible(x)
}

#' Obtain basket-specific interim recommendations
#' @param analysis An analysis returned by [update_trial()].
#' @return Audit-ready recommendation table.
#' @export
recommend_dose <- function(analysis) {
  stopifnot(inherits(analysis, "SAboinETB_analysis"))
  analysis$decisions
}

#' Select final optimal doses
#' @param analysis An analysis returned by [update_trial()].
#' @return One-row-per-basket final selection table.
#' @export
select_od <- function(analysis) {
  stopifnot(inherits(analysis, "SAboinETB_analysis"))
  d <- analysis$design; mat <- analysis$matrices
  rows <- lapply(seq_len(d$K), function(k) {
    admissible <- integer()
    for (j in seq_len(d$J)) {
      if (!analysis$active[k, j] || mat$n[k, j] < d$min_select) next
      local_t <- .beta_prob_gt(mat$tox[k, j], mat$n[k, j],
                               d$target_tox, d$prior_tox)
      working_t <- .beta_prob_gt(analysis$toxicity$y[k, j],
        analysis$toxicity$n[k, j], d$target_tox, d$prior_tox)
      local_e <- .beta_prob_lt(mat$eff[k, j], mat$n[k, j],
                               d$efficacy_null, d$prior_eff)
      working_e <- .beta_prob_lt(analysis$efficacy$y[k, j],
        analysis$efficacy$n[k, j], d$efficacy_null, d$prior_eff)
      safe <- local_t <= d$local_toxicity_selection_cutoff &&
        working_t <= d$working_toxicity_selection_cutoff
      active <- local_e <= d$efficacy_elimination_cutoff &&
        working_e <= d$working_efficacy_selection_cutoff
      if (safe && active) admissible <- c(admissible, j)
    }
    if (!length(admissible)) {
      return(data.frame(
        basket = d$basket_names[k], selected_dose = NA,
        posterior_tox = NA_real_, posterior_eff = NA_real_,
        expected_utility = NA_real_, status = "No admissible active dose"))
    }
    tox_means <- vapply(seq_len(d$J), function(j) {
      .beta_mean(analysis$toxicity$y[k, j], analysis$toxicity$n[k, j],
                 d$prior_tox)
    }, numeric(1))
    tox_iso <- .weighted_pava(
      tox_means, sum(d$prior_tox) + analysis$toxicity$n[k, ])
    eff_means <- vapply(seq_len(d$J), function(j) {
      .beta_mean(analysis$efficacy$y[k, j], analysis$efficacy$n[k, j],
                 d$prior_eff)
    }, numeric(1))
    utility <- .expected_utility(tox_iso[admissible],
                                 eff_means[admissible], d$utility)
    best <- admissible[utility == max(utility)]
    if (length(best) > 1L) best <- best[tox_iso[best] == min(tox_iso[best])]
    best <- min(best)
    data.frame(
      basket = d$basket_names[k], selected_dose = d$dose_levels[best],
      posterior_tox = tox_iso[best], posterior_eff = eff_means[best],
      expected_utility = .expected_utility(tox_iso[best], eff_means[best],
                                           d$utility),
      status = "Selected")
  })
  do.call(rbind, rows)
}

#' Extract borrowing diagnostics
#' @param analysis An analysis returned by [update_trial()].
#' @return List of endpoint-specific similarity, BESS, and donor weight arrays.
#' @export
borrowing_diagnostics <- function(analysis) {
  stopifnot(inherits(analysis, "SAboinETB_analysis"))
  list(
    toxicity_similarity = analysis$toxicity$similarity,
    efficacy_similarity = analysis$efficacy$similarity,
    toxicity_bess = analysis$toxicity$bess,
    efficacy_bess = analysis$efficacy$bess,
    toxicity_weights = analysis$toxicity$weights,
    efficacy_weights = analysis$efficacy$weights)
}
