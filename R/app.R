run_app <- function(launch.browser = interactive(), ...) {
  app_url <- paste0(
    "https://01a03253-f738-56f5-30c7-6829683c4a56",
    ".share.connect.posit.cloud/"
  )

  if (!is.logical(launch.browser) || length(launch.browser) != 1L ||
      is.na(launch.browser)) {
    stop("`launch.browser` must be TRUE or FALSE.", call. = FALSE)
  }

  if (isTRUE(launch.browser)) {
    utils::browseURL(app_url, ...)
  }

  invisible(app_url)
}
