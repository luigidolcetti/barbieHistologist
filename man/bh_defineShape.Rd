% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genericShape.R
\name{bh_defineShape}
\alias{bh_defineShape}
\title{Define a new shape}
\usage{
bh_defineShape(
  majorAxis = NULL,
  minorAxis = NULL,
  roundness = NULL,
  nArms = NULL,
  fixedArms = NULL,
  orientation = NULL,
  armExt = NULL,
  armElbow = NULL,
  armSwing = NULL,
  armTrig = NULL
)
}
\arguments{
\item{majorAxis}{numeric, mean and sd pair.}

\item{minorAxis}{numeric, mean and sd pair.}

\item{nArms}{numeric, number of arms to be produced. This value is passed to
rpois to obtain the actual number of arms.}

\item{armExt}{numeric, mean and sd pair. Proportion of the major axis.
E.g. if major axis is 10 AU and arm extent is 1.5 the total lenght of
the arm will be approx 15 AU with the branching starting after 10 AU.}

\item{armSwing}{numeric, arm swinging extent in degree.}

\item{armTrig}{numeric, mean and sd pair. Determine the behaviour of branches.
A positive value determine the branches to get longer towards the edge,
a negative value makes branches to get smaller.}

\item{armElow}{numeric, number of branching points.}
}
\value{
An object of class genericShape
}
\description{
helper function to create a new shape prototype.
}
