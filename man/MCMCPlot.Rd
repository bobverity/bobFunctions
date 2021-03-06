% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{MCMCPlot}
\alias{MCMCPlot}
\title{MCMCPlot}
\usage{
MCMCPlot(x, y = NULL, col = grey(0.7), pch = 4, cex = 0.4,
  xmin = NULL, xmax = NULL, ymin = min(x, na.rm = T), ymax = max(x,
  na.rm = T), xlab = NULL, main = "", ylab = main)
}
\description{
Produces simple plot useful for visualising MCMC chains. Enter x vector to plot single chain, or x and y vectors to visualise correlation between chains.
}
