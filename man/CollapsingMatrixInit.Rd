% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collapsingmatrix.R
\name{CollapsingMatrixInit}
\alias{CollapsingMatrixInit}
\title{initializes data frame with discrete genes as rows and samples as columns.
each element is an integer representing how many qualifying variants fall
in that gene for the particular sample in question.}
\usage{
CollapsingMatrixInit(geneset, samples, qualvar = NULL)
}
\arguments{
\item{geneset}{set of genes to build collapsing matrix rows from.}

\item{samples}{set of samples to build collapsing matrix columns from.}

\item{qualvar}{qualvar data frame to code collapsing matrix on.}
}
\description{
initializes data frame with discrete genes as rows and samples as columns.
each element is an integer representing how many qualifying variants fall
in that gene for the particular sample in question.
}

