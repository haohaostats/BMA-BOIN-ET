.beta_mean <- function(y, n, prior) (prior[1] + y) / (sum(prior) + n)

.beta_prob_gt <- function(y, n, threshold, prior) {
  stats::pbeta(threshold, prior[1] + y, prior[2] + n - y,
               lower.tail = FALSE)
}

.beta_prob_lt <- function(y, n, threshold, prior) {
  stats::pbeta(threshold, prior[1] + y, prior[2] + n - y)
}

.set_partitions <- function(K) {
  if (K < 1L) return(list(integer()))
  out <- list()
  labels <- integer(K)
  labels[1] <- 1L
  visit <- function(i, largest) {
    if (i > K) {
      out[[length(out) + 1L]] <<- labels
      return(invisible(NULL))
    }
    for (g in seq_len(largest + 1L)) {
      labels[i] <<- g
      visit(i + 1L, max(largest, g))
    }
    invisible(NULL)
  }
  if (K == 1L) out[[1L]] <- labels else visit(2L, 1L)
  out
}

.log_sum_exp <- function(x) {
  anchor <- max(x)
  anchor + log(sum(exp(x - anchor)))
}

.partition_bma_endpoint <- function(y, n, prior, design) {
  K <- nrow(n)
  J <- ncol(n)
  partitions <- design$partitions
  group_counts <- vapply(partitions, max, integer(1))
  log_prior <- -log(K) - log(design$stirling_counts[group_counts])
  means <- matrix(NA_real_, K, J, dimnames = dimnames(n))
  mixtures <- vector("list", J)
  partition_weights <- matrix(
    NA_real_, length(partitions), J,
    dimnames = list(design$partition_labels, as.character(design$dose_levels))
  )
  coclustering <- array(
    0, c(K, K, J),
    dimnames = list(design$basket_names, design$basket_names,
                    as.character(design$dose_levels))
  )

  for (j in seq_len(J)) {
    log_marginal <- numeric(length(partitions))
    component_alpha <- matrix(NA_real_, length(partitions), K)
    component_beta <- matrix(NA_real_, length(partitions), K)
    for (r in seq_along(partitions)) {
      membership <- partitions[[r]]
      for (g in seq_len(max(membership))) {
        members <- which(membership == g)
        total_y <- sum(y[members, j])
        total_n <- sum(n[members, j])
        alpha <- prior[1] + total_y
        beta <- prior[2] + total_n - total_y
        log_marginal[r] <- log_marginal[r] +
          lbeta(alpha, beta) - lbeta(prior[1], prior[2])
        component_alpha[r, members] <- alpha
        component_beta[r, members] <- beta
      }
    }
    log_posterior <- log_prior + log_marginal
    weights <- exp(log_posterior - .log_sum_exp(log_posterior))
    partition_weights[, j] <- weights
    mixtures[[j]] <- lapply(seq_len(K), function(k) {
      list(weight = weights, alpha = component_alpha[, k],
           beta = component_beta[, k])
    })
    means[, j] <- vapply(mixtures[[j]], function(mix) {
      sum(mix$weight * mix$alpha / (mix$alpha + mix$beta))
    }, numeric(1))
    for (k in seq_len(K)) for (kp in seq_len(K)) {
      together <- vapply(partitions, function(partition) {
        partition[k] == partition[kp]
      }, logical(1))
      coclustering[k, kp, j] <- sum(weights[together])
    }
  }

  list(mean = means, mixtures = mixtures,
       partition_weights = partition_weights,
       coclustering = coclustering)
}

.mixture_prob_gt <- function(endpoint, k, j, threshold) {
  mix <- endpoint$mixtures[[j]][[k]]
  sum(mix$weight * stats::pbeta(threshold, mix$alpha, mix$beta,
                                lower.tail = FALSE))
}

.mixture_prob_lt <- function(endpoint, k, j, threshold) {
  mix <- endpoint$mixtures[[j]][[k]]
  sum(mix$weight * stats::pbeta(threshold, mix$alpha, mix$beta))
}

.mixture_quantile <- function(endpoint, k, j, probability,
                              tolerance = 1e-10) {
  mix <- endpoint$mixtures[[j]][[k]]
  cdf <- function(x) sum(mix$weight * stats::pbeta(x, mix$alpha, mix$beta))
  stats::uniroot(function(x) cdf(x) - probability, c(0, 1),
                 tol = tolerance)$root
}
