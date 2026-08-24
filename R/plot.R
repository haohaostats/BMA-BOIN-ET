#' Plot endpoint-specific borrowing diagnostics
#' @param analysis An analysis returned by [update_trial()].
#' @param endpoint Either `"toxicity"` or `"efficacy"`.
#' @param type Either `"similarity"` or `"bess"`.
#' @param dose Dose label used when `type = "bess"`.
#' @param ... Additional graphical parameters.
#' @return The plotted matrix, invisibly.
#' @export
plot_borrowing <- function(analysis, endpoint = c("toxicity", "efficacy"),
                           type = c("similarity", "bess"), dose = NULL, ...) {
  stopifnot(inherits(analysis, "SAboinETB_analysis"))
  endpoint <- match.arg(endpoint); type <- match.arg(type)
  object <- if (endpoint == "toxicity") analysis$toxicity else analysis$efficacy
  if (type == "similarity") {
    z <- object$similarity
    main <- paste(tools::toTitleCase(endpoint), "similarity")
  } else {
    if (is.null(dose)) dose <- analysis$design$dose_levels[1]
    j <- match(dose, analysis$design$dose_levels)
    if (is.na(j)) stop("Unknown dose label.")
    z <- object$weights[, , j] *
      matrix(analysis$matrices$n[, j], analysis$design$K, analysis$design$K,
             byrow = TRUE)
    main <- paste(tools::toTitleCase(endpoint), "donor BESS at dose", dose)
  }
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old))
  graphics::par(mar = c(6, 7, 3, 2))
  graphics::image(seq_len(ncol(z)), seq_len(nrow(z)), t(z[nrow(z):1, ]),
                  col = grDevices::gray(seq(1, 0.12, length.out = 20)),
                  axes = FALSE, xlab = "Donor basket", ylab = "",
                  main = main, ...)
  graphics::axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2)
  graphics::axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 2)
  z_max <- max(z, na.rm = TRUE)
  for (i in seq_len(nrow(z))) for (j in seq_len(ncol(z))) {
    value <- z[i, j]
    if (!is.finite(value)) next
    shade <- if (z_max > 0) value / z_max else 0
    label <- if (type == "similarity") sprintf("%.2f", value) else sprintf("%.1f", value)
    graphics::text(j, nrow(z) - i + 1, labels = label,
                   col = if (shade > 0.58) "white" else "black",
                   cex = 1.12, font = 2)
  }
  graphics::mtext("Recipient basket", side = 2, line = 5.2)
  graphics::box()
  invisible(z)
}

#' Synthetic cumulative interim data
#'
#' This data set is supplied solely to demonstrate the prospective software
#' workflow and can be replaced directly with trial-specific cumulative data.
#' @export
example_interim_data <- function() {
  data.frame(
    basket = rep(paste0("Basket ", 1:4), each = 3),
    dose = rep(1:3, 4),
    n = c(3, 6, 6, 3, 6, 3, 3, 6, 3, 3, 6, 3),
    tox = c(0, 1, 1, 0, 1, 0, 0, 2, 2, 0, 1, 1),
    eff = c(0, 2, 4, 1, 3, 2, 0, 1, 0, 1, 2, 1),
    stringsAsFactors = FALSE
  )
}
