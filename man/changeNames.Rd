% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{changeNames}
\alias{changeNames}
\title{changeNames}
\usage{
changeNames(df, oldNames = names(df), newNames)
}
\description{
Change names of selected variables in a data frame. Note that for large data frames the raw code of this function should be copied into script rather than using this function, as there will be a memory cost in copying the data frame within the function.
}
