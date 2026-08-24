\name{run_app}
\alias{run_app}
\title{Open the hosted SA-BOIN-ETB web application}
\usage{
run_app(launch.browser = interactive(), ...)
}
\arguments{
\item{launch.browser}{Logical; whether to open the application in the system
web browser. Set to \code{FALSE} to return the URL without opening it.}

\item{...}{Additional arguments passed to \code{utils::browseURL()} when
\code{launch.browser = TRUE}.}
}
\value{
Invisibly returns the hosted application URL.
}
\description{
Opens the hosted trial implementation console in a web browser. The hosted
interface is deployed separately and its source code is not bundled with
the public R package.
}
