% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cell.R
\name{bh_defineCell}
\alias{bh_defineCell}
\title{Define a new cell}
\usage{
bh_defineCell(
  name = NULL,
  cytoplasm = NULL,
  nucleus = NULL,
  organelle = NULL,
  markers = NULL
)
}
\arguments{
\item{cytoplasm}{genericShape.}

\item{nucleus}{genericShape.}

\item{organelle}{genericSape. (not implemented yet)}

\item{markers}{list, makers.}
}
\value{
An object of class cellPrototype
}
\description{
helper function to create a new cell prototype.
}
