test_that("synthetic interim data produce bounded, auditable recommendations", {
  design <- setup_design()
  analysis <- update_trial(design, example_interim_data())
  recommendations <- recommend_dose(analysis)

  expect_s3_class(analysis, "SAboinETB_analysis")
  expect_equal(nrow(recommendations), design$K)
  expect_true(all(recommendations$action %in%
                    c("Close basket", "De-escalate", "Escalate", "Stay")))
  expect_lte(max(analysis$toxicity$bess), design$bess_cap_t + 1e-8)
  expect_lte(max(analysis$efficacy$bess), design$bess_cap_e + 1e-8)
  expect_true(all(nzchar(recommendations$reason)))
})

test_that("invalid cumulative counts are rejected", {
  design <- setup_design()
  bad <- example_interim_data()
  bad$tox[1] <- bad$n[1] + 1L
  expect_error(update_trial(design, bad), "cannot exceed")
})

test_that("local evidence is retained when baskets conflict", {
  design <- setup_design()
  x <- example_interim_data()
  analysis <- update_trial(design, x)
  expect_equal(analysis$matrices$n[1, 1], 3)
  expect_gte(analysis$toxicity$n[1, 1], analysis$matrices$n[1, 1])
  expect_gte(analysis$efficacy$n[1, 1], analysis$matrices$n[1, 1])
})

test_that("the hosted web app is linked without bundling Shiny source", {
  app_url <- run_app(launch.browser = FALSE)
  expect_equal(
    app_url,
    paste0(
      "https://01a03253-f738-56f5-30c7-6829683c4a56",
      ".share.connect.posit.cloud/"
    )
  )
  bundled_app <- system.file(
    "shiny", "SAboinETB", "app.R", package = "SAboinETB"
  )
  expect_false(file.exists(bundled_app))
})

test_that("the locked profile matches the standard SA-BOIN-ETB algorithm", {
  design <- setup_design()
  expect_equal(design$curve_saturation, 3)
  expect_equal(design$conflict_scale_t, 2)
  expect_equal(design$conflict_scale_e, 3)
  expect_equal(design$escalation_gate_cutoff, 0.80)
  expect_equal(design$toxicity_elimination_cutoff, 0.95)
  expect_equal(design$efficacy_elimination_cutoff, 0.95)
  expect_equal(design$local_toxicity_selection_cutoff, 0.95)
  expect_equal(design$working_efficacy_selection_cutoff, 0.85)
  expect_equal(design$exploration_toxicity_cutoff, 0.24)
  expect_equal(design$exploration_efficacy_cutoff, 0.55)
})

test_that("frontier extension is gated at the historical frontier", {
  design <- setup_design()
  x <- example_interim_data()
  x$n[x$dose == 3] <- 3L
  x$tox[x$dose == 3] <- 0L
  x$eff[x$dose == 3] <- 0L
  analysis <- update_trial(design, x, current_dose = rep(2, 4))

  expect_equal(analysis$decisions$current_dose, rep(2, 4))
  expect_equal(analysis$decisions$frontier_dose, rep(3, 4))
  expect_true(all(analysis$decisions$frontier_gate))
  expect_equal(analysis$decisions$recommended_dose, rep(4, 4))
  expect_true(all(grepl("exploration frontier", analysis$decisions$reason)))
})

test_that("reference recommendations match the formal simulation engine", {
  design <- setup_design()
  analysis <- update_trial(design, example_interim_data(),
                           current_dose = rep(3, 4))

  expect_equal(analysis$decisions$recommended_dose, c(4, 3, 2, 3))
  expect_equal(analysis$decisions$action,
               c("Escalate", "Stay", "De-escalate", "Stay"))
  expect_equal(select_od(analysis)$selected_dose, c(3, 2, 2, 2))
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

test_that("de-escalation skips inactive intervening doses", {
  design <- setup_design(
    basket_names = c("B1", "B2"), dose_levels = 1:4,
    max_sample = c(24, 24)
  )
  mat <- list(
    n = rbind(c(3, 3, 3, 3), c(3, 0, 0, 0)),
    tox = matrix(0, 2, 4), eff = matrix(0, 2, 4)
  )
  active <- matrix(TRUE, 2, 4)
  active[1, 3:4] <- FALSE
  bt <- list(y = matrix(0, 2, 4), n = matrix(1, 2, 4))
  bt$y[1, 4] <- 1
  be <- list(y = matrix(0, 2, 4), n = matrix(1, 2, 4))

  choice <- SAboinETB:::.choose_next_dose(
    1, current = c(4L, 1L), active, mat, bt, be, design
  )

  expect_equal(choice$index, 2L)
  expect_true(active[1, choice$index])
  expect_match(choice$rule, "highest active lower dose")
})

test_that("de-escalation stays only when the current dose is lowest active", {
  design <- setup_design(
    basket_names = c("B1", "B2"), dose_levels = 1:3,
    max_sample = c(24, 24)
  )
  mat <- list(
    n = rbind(c(3, 3, 0), c(3, 0, 0)),
    tox = matrix(0, 2, 3), eff = matrix(0, 2, 3)
  )
  active <- matrix(TRUE, 2, 3)
  active[1, 1] <- FALSE
  bt <- list(y = matrix(0, 2, 3), n = matrix(1, 2, 3))
  bt$y[1, 2] <- 1
  be <- list(y = matrix(0, 2, 3), n = matrix(1, 2, 3))

  choice <- SAboinETB:::.choose_next_dose(
    1, current = c(2L, 1L), active, mat, bt, be, design
  )

  expect_equal(choice$index, 2L)
  expect_match(choice$rule, "lowest active dose")

  active[1, 2] <- FALSE
  closed <- SAboinETB:::.choose_next_dose(
    1, current = c(2L, 1L), active, mat, bt, be, design
  )
  expect_true(is.na(closed$index))
  expect_match(closed$rule, "close the basket")
})

test_that("an empty local-toxicity candidate set closes the basket", {
  design <- setup_design(
    basket_names = c("B1", "B2"), dose_levels = 1:3,
    max_sample = c(24, 24)
  )
  mat <- list(
    n = rbind(c(0, 2, 0), c(0, 0, 0)),
    tox = rbind(c(0, 2, 0), c(0, 0, 0)),
    eff = matrix(0, 2, 3)
  )
  active <- matrix(TRUE, 2, 3)
  bt <- list(y = matrix(0, 2, 3), n = matrix(1, 2, 3))
  bt$y[1, 2] <- 1
  bt$n[1, 2] <- 4
  be <- list(y = matrix(0, 2, 3), n = matrix(1, 2, 3))

  expect_gt(
    SAboinETB:::.beta_prob_gt(2, 2, design$target_tox, design$prior_tox),
    design$local_toxicity_selection_cutoff
  )
  choice <- SAboinETB:::.choose_next_dose(
    1, current = c(2L, 1L), active, mat, bt, be, design
  )

  expect_true(is.na(choice$index))
  expect_match(choice$rule, "No locally toxicity-admissible")
})
