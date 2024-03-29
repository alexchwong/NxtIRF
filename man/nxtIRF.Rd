% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dash.R
\name{nxtIRF}
\alias{nxtIRF}
\title{Launches the NxtIRF Graphics User Interface using Shiny Dashboard}
\usage{
nxtIRF(mode = c("dialog", "browser"), width = 1600, height = 900)
}
\arguments{
\item{mode}{(default \code{"dialog"}) \code{"dialog"} displays NxtIRF in a dialog box
with specified width and height. \code{"browser"} opens NxtIRF in a browser-
like resizable window.}

\item{width, height}{If \code{mode} is set to \code{"dialog"}, the specified width
and height of the NxtIRF app.}
}
\value{
An interactive shinydashboard NxtIRF app runs.
}
\description{
This function launches the NxtIRF interactive app using Shiny Dashboard
}
\examples{
if(interactive()) {

# Launches interactive ShinyDashboard NxtIRF app as fixed-size dialog box
    nxtIRF(mode = "dialog", width = 1600, height = 900) 

# Launches interactive ShinyDashboard NxtIRF app as browser window
    nxtIRF(mode = "browser")

}

}
