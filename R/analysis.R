validate_interim_data <- function(data, design) {
  stopifnot(inherits(design, "BMABOINET_design"))
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

.local_rate <- function(events, n, k, j) {
  if (n[k, j] <= 0) return(NA_real_)
  events[k, j] / n[k, j]
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
  if (!inherits(previous_analysis, "BMABOINET_analysis")) {
    stop("previous_analysis must be returned by update_trial().")
  }
  previous_design <- previous_analysis$design
  locked <- setdiff(names(design), "created")
  if (!identical(previous_design[locked], design[locked])) {
    stop("previous_analysis does not use the locked design.")
  }
  if (!identical(previous_design$basket_names, design$basket_names) ||
      !identical(previous_design$dose_levels, design$dose_levels)) {
    stop("previous_analysis does not use the locked design.")
  }
  previous_analysis$active
}

.apply_elimination <- function(active, mat, toxicity, efficacy, design) {
  for (k in seq_len(design$K)) {
    for (j in seq_len(design$J)) {
      if (!active[k, j]) next
      local_t_bad <- mat$n[k, j] >= design$min_eliminate_t &&
        .beta_prob_gt(mat$tox[k, j], mat$n[k, j], design$target_tox,
                      design$prior_tox) > design$toxicity_elimination_cutoff
      if (local_t_bad) {
        active[k, j:design$J] <- FALSE
        next
      }
      local_e_bad <- mat$n[k, j] >= design$min_eliminate_e &&
        .beta_prob_lt(mat$eff[k, j], mat$n[k, j], design$efficacy_null,
                      design$prior_eff) > design$efficacy_elimination_cutoff
      bma_e_bad <- .mixture_prob_lt(
        efficacy, k, j, design$efficacy_null
      ) > design$efficacy_elimination_cutoff
      if (local_e_bad && bma_e_bad) active[k, j] <- FALSE
    }
  }
  list(active = active, open = rowSums(active) > 0)
}

.best_by_efficacy <- function(candidates, k, toxicity, efficacy) {
  if (!length(candidates)) return(NA_integer_)
  pe <- efficacy$mean[k, candidates]
  best <- candidates[pe == max(pe)]
  if (length(best) > 1L) {
    pt <- toxicity$mean[k, best]
    best <- best[pt == min(pt)]
  }
  as.integer(min(best))
}

.deescalate <- function(k, j, active) {
  lower <- which(active[k, ] & seq_len(ncol(active)) < j)
  if (length(lower)) {
    return(list(index = as.integer(max(lower)),
                rule = "Toxicity de-escalation to the highest active lower dose"))
  }
  if (isTRUE(active[k, j])) {
    return(list(index = as.integer(j),
                rule = "Toxicity signal at the lowest active dose; stay pending elimination"))
  }
  list(index = NA_integer_,
       rule = "Toxicity action required and no active lower dose remains; close the basket")
}

.choose_next_dose <- function(k, current, active, mat, toxicity, efficacy,
                              design) {
  j <- current[k]
  bounds <- get_boundaries(design)
  local_t <- .local_rate(mat$tox, mat$n, k, j)
  bma_t <- toxicity$mean[k, j]
  bma_e <- efficacy$mean[k, j]

  if ((!is.na(local_t) && local_t >= bounds["deescalate"]) ||
      bma_t >= bounds["deescalate"]) {
    return(.deescalate(k, j, active))
  }

  if (!active[k, j]) {
    higher <- which(active[k, ] & seq_len(design$J) > j)
    if (length(higher)) {
      return(list(index = as.integer(min(higher)),
                  rule = "Current dose was eliminated for futility; move to the next active higher dose"))
    }
    lower_tried <- which(active[k, ] & seq_len(design$J) < j & mat$n[k, ] > 0)
    target <- .best_by_efficacy(lower_tried, k, toxicity, efficacy)
    if (!is.na(target)) {
      return(list(index = target,
                  rule = "Current dose was eliminated for futility; select the best active tried lower dose"))
    }
    return(list(index = NA_integer_,
                rule = "Current dose was eliminated and no assignable active dose remains; close the basket"))
  }

  if (bma_e > bounds["efficacy"]) {
    return(list(index = as.integer(j),
                rule = "BOIN-ET: efficacy above its boundary and toxicity below de-escalation"))
  }

  if (bma_t <= bounds["escalate"]) {
    if (j < design$J && active[k, j + 1L]) {
      return(list(index = as.integer(j + 1L),
                  rule = "BOIN-ET: low efficacy and toxicity; escalate"))
    }
    return(list(index = as.integer(j),
                rule = "BOIN-ET escalation unavailable; stay"))
  }

  if (j < design$J && active[k, j + 1L] && mat$n[k, j + 1L] == 0) {
    return(list(index = as.integer(j + 1L),
                rule = "BOIN-ET middle cell: explore the untried adjacent higher dose"))
  }
  neighbours <- intersect(seq.int(max(1L, j - 1L), min(design$J, j + 1L)),
                          which(active[k, ] & mat$n[k, ] > 0))
  target <- .best_by_efficacy(neighbours, k, toxicity, efficacy)
  if (is.na(target)) {
    return(list(index = NA_integer_,
                rule = "BOIN-ET middle cell has no active observed neighbour; close the basket"))
  }
  list(index = target,
       rule = "BOIN-ET middle cell: maximize BMA efficacy over active observed neighbours")
}

update_trial <- function(design, data, current_dose = NULL,
                         previous_analysis = NULL) {
  data <- validate_interim_data(data, design)
  mat <- .as_matrices(data, design)
  toxicity <- .partition_bma_endpoint(mat$tox, mat$n, design$prior_tox, design)
  efficacy <- .partition_bma_endpoint(mat$eff, mat$n, design$prior_eff, design)
  current <- .current_indices(design, mat, current_dose, previous_analysis)
  active <- .initial_active(design, previous_analysis)
  elimination <- .apply_elimination(active, mat, toxicity, efficacy, design)
  active <- elimination$active
  open <- elimination$open
  if (!is.null(previous_analysis)) {
    old <- previous_analysis$matrices
    if (any(mat$n < old$n) || any(mat$tox < old$tox) || any(mat$eff < old$eff) ||
        any(mat$tox - old$tox > mat$n - old$n) ||
        any(mat$eff - old$eff > mat$n - old$n)) {
      stop("Cumulative counts must extend the previous analysis without changing past outcomes.")
    }
    open <- open & previous_analysis$open
  }

  decisions <- vector("list", design$K)
  for (k in seq_len(design$K)) {
    j <- current[k]

    lowest_toxicity <- j == 1L && mat$n[k, j] > 0 &&
      (mat$tox[k, j] / mat$n[k, j] >= get_boundaries(design)["deescalate"] ||
       toxicity$mean[k, j] >= get_boundaries(design)["deescalate"])
    if (lowest_toxicity) {
      active[k, ] <- FALSE
      open[k] <- FALSE
    }
    if (!open[k]) {
      target <- NA_integer_
      rule <- if (lowest_toxicity) "Toxicity action at the lowest dose; close the basket" else
        "Basket is closed; no further enrollment"
    } else if (sum(mat$n[k, ]) == 0) {
      target <- 1L
      rule <- "Start the basket at the lowest dose"
    } else {
      choice <- .choose_next_dose(k, current, active, mat, toxicity,
                                  efficacy, design)
      target <- choice$index
      rule <- choice$rule
    }
    if (is.na(target)) open[k] <- FALSE
    action <- if (is.na(target)) "Close basket" else if (target > j) {
      "Escalate"
    } else if (target < j) "De-escalate" else "Stay"
    decisions[[k]] <- data.frame(
      basket = design$basket_names[k], current_dose = design$dose_levels[j],
      local_n = mat$n[k, j], local_tox = mat$tox[k, j],
      local_eff = mat$eff[k, j],
      bma_mean_tox = toxicity$mean[k, j],
      bma_mean_eff = efficacy$mean[k, j],
      local_tox_tail = .beta_prob_gt(
        mat$tox[k, j], mat$n[k, j], design$target_tox, design$prior_tox),
      bma_tox_tail = .mixture_prob_gt(toxicity, k, j, design$target_tox),
      local_futility = .beta_prob_lt(
        mat$eff[k, j], mat$n[k, j], design$efficacy_null, design$prior_eff),
      bma_futility = .mixture_prob_lt(efficacy, k, j, design$efficacy_null),
      action = action,
      recommended_dose = if (is.na(target)) NA else design$dose_levels[target],
      reason = rule, stringsAsFactors = FALSE
    )
  }
  analysis <- list(
    design = design, data = data, matrices = mat,
    toxicity = toxicity, efficacy = efficacy,
    active = active, open = open, current_index = current,
    decisions = do.call(rbind, decisions), boundaries = get_boundaries(design),
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  class(analysis) <- "BMABOINET_analysis"
  analysis
}

print.BMABOINET_analysis <- function(x, ...) {
  cat("BMA-BOIN-ET partition-model-averaged interim analysis\n")
  cat("  Timestamp:", x$timestamp, "\n\n")
  print(x$decisions[c("basket", "current_dose", "local_n", "action",
                      "recommended_dose", "bma_mean_tox", "bma_mean_eff")],
        row.names = FALSE)
  invisible(x)
}

recommend_dose <- function(analysis) {
  stopifnot(inherits(analysis, "BMABOINET_analysis"))
  analysis$decisions
}

select_od <- function(analysis) {
  stopifnot(inherits(analysis, "BMABOINET_analysis"))
  d <- analysis$design
  mat <- analysis$matrices
  rows <- lapply(seq_len(d$K), function(k) {
    observed <- which(mat$n[k, ] > 0)
    tox_iso <- rep(NA_real_, d$J)
    if (length(observed)) {
      tox_iso[observed] <- .weighted_pava(
        analysis$toxicity$mean[k, observed], mat$n[k, observed]
      )
    }
    admissible <- integer()
    for (j in seq_len(d$J)) {
      if (!analysis$active[k, j] || mat$n[k, j] < d$min_select) next
      local_t <- .beta_prob_gt(mat$tox[k, j], mat$n[k, j],
                               d$target_tox, d$prior_tox)
      bma_t <- .mixture_prob_gt(analysis$toxicity, k, j, d$target_tox)
      bma_e <- .mixture_prob_lt(analysis$efficacy, k, j, d$efficacy_null)
      if (max(local_t, bma_t) <= d$toxicity_selection_cutoff &&
          bma_e <= d$efficacy_elimination_cutoff) {
        admissible <- c(admissible, j)
      }
    }
    if (!length(admissible)) {
      return(data.frame(
        basket = d$basket_names[k], selected_dose = NA,
        posterior_tox = NA_real_, posterior_eff = NA_real_,
        expected_utility = NA_real_, status = "No admissible active dose"
      ))
    }
    utility <- .expected_utility(
      tox_iso[admissible], analysis$efficacy$mean[k, admissible], d$utility
    )
    best <- admissible[utility == max(utility)]
    if (length(best) > 1L) best <- best[tox_iso[best] == min(tox_iso[best])]
    best <- min(best)
    data.frame(
      basket = d$basket_names[k], selected_dose = d$dose_levels[best],
      posterior_tox = tox_iso[best],
      posterior_eff = analysis$efficacy$mean[k, best],
      expected_utility = .expected_utility(
        tox_iso[best], analysis$efficacy$mean[k, best], d$utility
      ),
      status = "Selected"
    )
  })
  do.call(rbind, rows)
}

borrowing_diagnostics <- function(analysis) {
  stopifnot(inherits(analysis, "BMABOINET_analysis"))
  list(
    toxicity_partition_weights = analysis$toxicity$partition_weights,
    efficacy_partition_weights = analysis$efficacy$partition_weights,
    toxicity_coclustering = analysis$toxicity$coclustering,
    efficacy_coclustering = analysis$efficacy$coclustering
  )
}
