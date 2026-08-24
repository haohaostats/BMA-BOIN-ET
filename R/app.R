#' Open the hosted SA-BOIN-ETB web application
#'
#' Opens the hosted trial implementation console in a web browser. The hosted
#' interface is deployed separately and its source code is not bundled with
#' the public R package.
#'
#' @param launch.browser Logical; whether to open the application in the system
#'   web browser. Set to `FALSE` to return the URL without opening it.
#' @param ... Additional arguments passed to [utils::browseURL()] when
#'   `launch.browser = TRUE`.
#' @return Invisibly returns the hosted application URL.
#' @export
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
