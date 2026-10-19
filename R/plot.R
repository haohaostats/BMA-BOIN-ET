plot_borrowing <- function(analysis, endpoint = c("toxicity", "efficacy"),
                           type = c("coclustering", "partitions"),
                           dose = NULL, ...) {
  stopifnot(inherits(analysis, "BMABOINET_analysis"))
  endpoint <- match.arg(endpoint); type <- match.arg(type)
  object <- if (endpoint == "toxicity") analysis$toxicity else analysis$efficacy
  if (is.null(dose)) dose <- analysis$design$dose_levels[1]
  j <- match(dose, analysis$design$dose_levels)
  if (is.na(j)) stop("Unknown dose label.")
  if (type == "partitions") {
    weights <- object$partition_weights[, j]
    graphics::barplot(
      weights, names.arg = rownames(object$partition_weights), las = 2,
      ylab = "Posterior probability", xlab = "Partition label",
      main = paste(tools::toTitleCase(endpoint),
                   "partition posterior at dose", dose), ...
    )
    return(invisible(weights))
  }
  z <- object$coclustering[, , j]
  main <- paste(tools::toTitleCase(endpoint),
                "posterior co-clustering at dose", dose)
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old))
  graphics::par(mar = c(6, 7, 3, 2))
  graphics::image(seq_len(ncol(z)), seq_len(nrow(z)), t(z[nrow(z):1, ]),
                  col = grDevices::gray(seq(1, 0.12, length.out = 20)),
                  axes = FALSE, xlab = "Basket", ylab = "",
                  main = main, ...)
  graphics::axis(1, at = seq_len(ncol(z)), labels = colnames(z), las = 2)
  graphics::axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 2)
  z_max <- max(z, na.rm = TRUE)
  for (i in seq_len(nrow(z))) for (j in seq_len(ncol(z))) {
    value <- z[i, j]
    if (!is.finite(value)) next
    shade <- if (z_max > 0) value / z_max else 0
    label <- sprintf("%.2f", value)
    graphics::text(j, nrow(z) - i + 1, labels = label,
                   col = if (shade > 0.58) "white" else "black",
                   cex = 1.12, font = 2)
  }
  graphics::mtext("Basket", side = 2, line = 5.2)
  graphics::box()
  invisible(z)
}

example_interim_data <- function() {
  example_manuscript_data()
}

example_manuscript_data <- function() {
  data.frame(
    basket = rep(paste0("Basket ", 1:4), each = 5),
    dose = rep(1:5, 4),
    n = rep(c(3, 6, 6, 6, 3), 4),
    tox = c(rep(c(0, 0, 1, 1, 1), 2), rep(c(0, 2, 3, 4, 2), 2)),
    eff = rep(c(1, 3, 4, 3, 1, 0, 2, 2, 4, 2), 2),
    stringsAsFactors = FALSE
  )
}
