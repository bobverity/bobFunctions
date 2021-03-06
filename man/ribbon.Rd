% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{ribbon}
\alias{ribbon}
\title{ribbon}
\usage{
ribbon(y1, y2, x = 1:length(y1), upperLower = TRUE, density = NA,
  border = NA, col = "#FF000020")
}
\arguments{
\item{y1}{lower limit of ribbon, or alternatively middle value if upperLower=FALSE}

\item{y2}{upper limit of ribbon, or alternatively thickness if upperLower=FALSE}

\item{x}{x-values}

\item{upperLower}{whether y-values represent upper and lower values of the ribbon}

\item{density}{density of shading (leave NA for no shading)}

\item{border}{colour of border (leave NA for no border)}

\item{col}{colour of ribbon, or colour of shading lines if density>0}
}
\description{
Adds ribbon to plot. Set upper and lower limits of ribbon, or middle values and thickness.
}
