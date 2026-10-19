test_that("the fixed four-basket design enumerates all partitions", {
  design <- setup_design()
  expect_equal(length(design$partitions), 15L)
  expect_equal(design$stirling_counts, c(1L, 7L, 6L, 1L))
  expect_equal(design$implementation_profile, "BMA-BOIN-ET-2.0")
  expect_false(any(c("bess_cap_t", "conflict_scale_t", "min_borrow",
                     "prior_conflict_scale") %in% names(design)))
})

test_that("the partition prior is uniform over the number of groups", {
  design <- setup_design()
  zero <- matrix(0, design$K, design$J)
  endpoint <- BMABOINET:::.partition_bma_endpoint(
    zero, zero, design$prior_tox, design
  )
  groups <- vapply(design$partitions, max, integer(1))
  mass <- vapply(seq_len(design$K), function(b) {
    sum(endpoint$partition_weights[groups == b, 1])
  }, numeric(1))
  expect_equal(mass, rep(1 / design$K, design$K), tolerance = 1e-12)
})

test_that("two-basket posterior partition weights match direct calculation", {
  design <- setup_design(basket_names = c("A", "B"), dose_levels = 1:2)
  n <- matrix(c(6, 0, 6, 0), 2, 2, byrow = TRUE)
  y <- matrix(c(3, 0, 3, 0), 2, 2, byrow = TRUE)
  endpoint <- BMABOINET:::.partition_bma_endpoint(y, n, c(1, 1), design)
  pooled <- exp(lbeta(7, 7) - lbeta(1, 1))
  separate <- exp(2 * (lbeta(4, 4) - lbeta(1, 1)))
  expected <- c(pooled, separate) / (pooled + separate)
  expect_equal(unname(endpoint$partition_weights[, 1]), expected, tolerance = 1e-12)
  expect_equal(endpoint$coclustering[1, 2, 1], expected[1], tolerance = 1e-12)
})

test_that("endpoint-specific data yield endpoint-specific co-clustering", {
  design <- setup_design(max_sample = rep(80, 4))
  n <- matrix(0, 4, 5); n[, 2] <- 40
  tox <- eff <- matrix(0, 4, 5)
  tox[, 2] <- c(4, 4, 28, 28)
  eff[, 2] <- c(8, 32, 8, 32)
  bt <- BMABOINET:::.partition_bma_endpoint(tox, n, c(1, 1), design)
  be <- BMABOINET:::.partition_bma_endpoint(eff, n, c(1, 1), design)
  expect_gt(bt$coclustering[1, 2, 2], bt$coclustering[1, 3, 2])
  expect_gt(be$coclustering[1, 3, 2], be$coclustering[1, 2, 2])
})

test_that("synthetic interim data produce bounded auditable recommendations", {
  design <- setup_design()
  analysis <- update_trial(design, example_interim_data())
  recommendations <- recommend_dose(analysis)
  expect_s3_class(analysis, "BMABOINET_analysis")
  expect_equal(nrow(recommendations), design$K)
  expect_true(all(recommendations$action %in%
                    c("Close basket", "De-escalate", "Escalate", "Stay")))
  expect_true(all(analysis$toxicity$mean >= 0 & analysis$toxicity$mean <= 1))
  expect_true(all(analysis$efficacy$mean >= 0 & analysis$efficacy$mean <= 1))
  expect_true(all(nzchar(recommendations$reason)))
})

test_that("all eliminations are completed before dose action", {
  design <- setup_design(
    basket_names = c("A", "B"), dose_levels = 1:3,
    max_sample = c(30, 30)
  )
  x <- expand.grid(basket = c("A", "B"), dose = 1:3,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  x <- x[order(x$basket, x$dose), ]
  x$n <- rep(c(3, 9, 3), 2)
  x$tox <- rep(c(0, 0, 2), 2)
  x$eff <- rep(c(1, 0, 1), 2)
  analysis <- update_trial(design, x, current_dose = c(A = 3, B = 3))
  expect_false(analysis$active[1, 2])
  expect_equal(analysis$decisions$recommended_dose[1], 1)
  expect_equal(analysis$decisions$action[1], "De-escalate")
})

test_that("de-escalation uses inclusive upper boundary", {
  design <- setup_design(basket_names = c("A", "B"), dose_levels = 1:3)
  active <- matrix(TRUE, 2, 3)
  mat <- list(n = matrix(c(3, 3, 0, 3, 0, 0), 2, 3, byrow = TRUE),
              tox = matrix(0, 2, 3), eff = matrix(0, 2, 3))
  toxicity <- efficacy <- list(mean = matrix(0, 2, 3))
  toxicity$mean[1, 2] <- get_boundaries(design)["deescalate"]
  choice <- BMABOINET:::.choose_next_dose(
    1, c(2L, 1L), active, mat, toxicity, efficacy, design
  )
  expect_equal(choice$index, 1L)
})

test_that("middle-cell neighbours include the current dose", {
  design <- setup_design(basket_names = c("A", "B"), dose_levels = 1:3)
  active <- matrix(TRUE, 2, 3)
  mat <- list(n = matrix(c(3, 3, 3, 3, 0, 0), 2, 3, byrow = TRUE),
              tox = matrix(0, 2, 3), eff = matrix(0, 2, 3))
  toxicity <- list(mean = matrix(0.30, 2, 3))
  efficacy <- list(mean = matrix(0.20, 2, 3))
  efficacy$mean[1, ] <- c(0.20, 0.29, 0.21)
  choice <- BMABOINET:::.choose_next_dose(
    1, c(2L, 1L), active, mat, toxicity, efficacy, design
  )
  expect_equal(choice$index, 2L)
})

test_that("utility equals the additive BOIN12 marginal utility", {
  p_t <- c(0, 0.25, 1)
  p_e <- c(0, 0.50, 1)
  expect_equal(
    BMABOINET:::.expected_utility(p_t, p_e, c(100, 40, 60, 0)),
    40 - 40 * p_t + 60 * p_e
  )
})

test_that("irreversible elimination is retained across analyses", {
  design <- setup_design()
  x <- example_interim_data()
  row <- x$basket == "Basket 1" & x$dose == 3
  x$n[row] <- 3L; x$tox[row] <- 3L; x$eff[row] <- 0L
  first <- update_trial(design, x, current_dose = rep(3, 4))
  expect_false(any(first$active[1, 3:5]))
  second <- update_trial(design, x, previous_analysis = first)
  expect_false(any(second$active[1, 3:5]))
})

test_that("invalid cumulative counts are rejected", {
  design <- setup_design()
  bad <- example_interim_data()
  bad$tox[1] <- bad$n[1] + 1L
  expect_error(update_trial(design, bad), "cannot exceed")
})

test_that("the hosted web app link remains available", {
  app_url <- run_app(launch.browser = FALSE)
  expect_true(grepl("share.connect.posit.cloud", app_url, fixed = TRUE))
})

test_that("manuscript Table 7 reproduces Table 8 BMA recommendations", {
  d <- setup_design()
  fit <- update_trial(d, example_manuscript_data())
  expect_equal(select_od(fit)$selected_dose, c(3, 4, 2, 2))
  expect_equal(d$toxicity_selection_cutoff, 0.70)
  expect_equal(d$efficacy_elimination_cutoff, 0.50)
  expect_equal(unname(fit$toxicity$mean[, 2]),
               c(0.1455586, 0.1455586, 0.2867173, 0.2867173), tolerance = 1e-6)
})

test_that("lowest-dose toxicity closes a basket irreversibly", {
  d <- setup_design(basket_names = c("A", "B"), dose_levels = 1:3)
  x <- data.frame(basket = c("A", "B"), dose = 1, n = 3, tox = 2, eff = 3)
  first <- update_trial(d, x, current_dose = c(1, 1))
  expect_false(any(first$open))
  expect_false(any(first$active))
  expect_true(all(first$decisions$action == "Close basket"))
  x$n <- 12
  x$eff <- 12
  later <- update_trial(d, x, previous_analysis = first)
  expect_false(any(later$open))
  expect_true(all(is.na(select_od(later)$selected_dose)))
})

test_that("an empty trial starts at dose one", {
  d <- setup_design()
  x <- example_interim_data(); x$n <- x$tox <- x$eff <- 0
  fit <- update_trial(d, x)
  expect_equal(fit$decisions$recommended_dose, rep(1, 4))
  expect_true(all(fit$open))
})

test_that("previous counts and design cannot silently change", {
  d <- setup_design()
  x <- example_interim_data()
  first <- update_trial(d, x)
  x$n[1] <- x$n[1] - 1
  expect_error(update_trial(d, x, previous_analysis = first), "Cumulative")
  other <- setup_design(target_tox = 0.25)
  expect_error(update_trial(other, example_interim_data(), previous_analysis = first), "locked design")
})

test_that("partition mixtures are normalized and co-clustering is symmetric", {
  fit <- update_trial(setup_design(), example_manuscript_data())
  for (endpoint in list(fit$toxicity, fit$efficacy)) {
    expect_equal(unname(colSums(endpoint$partition_weights)), rep(1, 5))
    for (j in 1:5) {
      q <- endpoint$coclustering[, , j]
      expect_equal(q, t(q))
      expect_equal(unname(diag(q)), rep(1, 4))
    }
  }
})
