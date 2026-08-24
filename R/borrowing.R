.beta_mean <- function(y, n, prior) (prior[1] + y) / (sum(prior) + n)

.beta_prob_gt <- function(y, n, threshold, prior) {
  stats::pbeta(threshold, prior[1] + y, prior[2] + n - y, lower.tail = FALSE)
}

.beta_prob_lt <- function(y, n, threshold, prior) {
  stats::pbeta(threshold, prior[1] + y, prior[2] + n - y)
}

.log_bf_equal <- function(y1, n1, y2, n2, prior) {
  a <- prior[1]; b <- prior[2]
  lbeta(a + y1 + y2, b + n1 + n2 - y1 - y2) + lbeta(a, b) -
    lbeta(a + y1, b + n1 - y1) - lbeta(a + y2, b + n2 - y2)
}

.curve_similarity <- function(y, n, xi, design, prior) {
  K <- nrow(n); J <- ncol(n); out <- diag(1, K)
  for (k in seq_len(K - 1L)) for (kp in (k + 1L):K) {
    prior_odds <- log(xi[k, kp] / max(1e-8, 1 - xi[k, kp]))
    prior_odds <- max(-design$prior_odds_cap,
                      min(design$prior_odds_cap, prior_odds))
    observed <- which(n[k, ] > 0 & n[kp, ] > 0)
    if (!length(observed)) {
      score <- 0
    } else {
      evidence <- 0; common_information <- 0
      for (j in observed) {
        m <- min(n[k, j], n[kp, j])
        reliability <- min(1, m / 6)
        lbf <- .log_bf_equal(y[k, j], n[k, j], y[kp, j], n[kp, j], prior)
        evidence <- evidence + reliability * max(-5, min(5, lbf))
        common_information <- common_information + m
      }
      curve_reliability <- min(1, common_information / design$curve_saturation)
      weighted_evidence <- curve_reliability * evidence
      prior_discount <- if (prior_odds * weighted_evidence < 0) {
        exp(-abs(weighted_evidence) / design$prior_conflict_scale)
      } else 1
      score <- stats::plogis(prior_discount * prior_odds + weighted_evidence) *
        curve_reliability
    }
    out[k, kp] <- out[kp, k] <- score
  }
  out
}

.local_compatibility <- function(y1, n1, y2, n2, scale, prior) {
  if (n1 <= 0 || n2 <= 0) return(0)
  m1 <- .beta_mean(y1, n1, prior); m2 <- .beta_mean(y2, n2, prior)
  v1 <- m1 * (1 - m1) / (n1 + sum(prior) + 1)
  v2 <- m2 * (1 - m2) / (n2 + sum(prior) + 1)
  z <- abs(m1 - m2) / sqrt(v1 + v2 + 1e-10)
  exp(-0.5 * (z / scale)^2)
}

.borrow_endpoint <- function(y, n, endpoint, design) {
  is_tox <- identical(endpoint, "toxicity")
  prior <- if (is_tox) design$prior_tox else design$prior_eff
  xi <- if (is_tox) design$prior_similarity_t else design$prior_similarity_e
  cap_abs <- if (is_tox) design$bess_cap_t else design$bess_cap_e
  cap_ratio <- if (is_tox) design$bess_ratio_t else design$bess_ratio_e
  conflict_scale <- if (is_tox) design$conflict_scale_t else design$conflict_scale_e
  similarity <- .curve_similarity(y, n, xi, design, prior)
  K <- nrow(n); J <- ncol(n); yb <- y; nb <- n
  weights <- array(0, c(K, K, J),
                   dimnames = list(design$basket_names, design$basket_names,
                                   as.character(design$dose_levels)))
  for (k in seq_len(K)) for (j in seq_len(J)) {
    if (n[k, j] < design$min_borrow) next
    raw <- numeric(K)
    for (kp in setdiff(seq_len(K), k)) {
      if (n[kp, j] <= 0) next
      raw[kp] <- similarity[k, kp] *
        .local_compatibility(y[k, j], n[k, j], y[kp, j], n[kp, j],
                             conflict_scale, prior)
    }
    raw_ess <- sum(raw * n[, j])
    cap <- min(cap_abs, cap_ratio * n[k, j])
    scaling <- if (raw_ess > 0) min(1, cap / raw_ess) else 0
    weights[k, , j] <- raw * scaling
    yb[k, j] <- y[k, j] + sum(weights[k, , j] * y[, j])
    nb[k, j] <- n[k, j] + sum(weights[k, , j] * n[, j])
  }
  dimnames(similarity) <- list(design$basket_names, design$basket_names)
  list(y = yb, n = nb, weights = weights, similarity = similarity,
       bess = nb - n)
}

