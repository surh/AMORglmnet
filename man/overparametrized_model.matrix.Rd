% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/overparametrized_model.matrix.r
\name{overparametrized_model.matrix}
\alias{overparametrized_model.matrix}
\title{Obtain overparametrized design matrix}
\usage{
overparametrized_model.matrix(formula, data, remove.constant = TRUE,
  intercept = TRUE)
}
\arguments{
\item{formula}{Righ hand side of formula for which a design matrix is desired.}

\item{data}{\code{data.frame} containing the variables used in the formula}

\item{remove.constant}{Logical indicating whether columns that never vary should be removed.
Most commonly these are due to specific combinations of factors that never happen on the
dataset.}

\item{intercept}{Logical indicating whether or not to add an intercept}
}
\description{
This functions takes a formula and data frame and creates a design
matrix that has a column for every level of every factor.
}
\examples{
data(Rhizo.map)

X <- overparametrized_model.matrix(~ fraction * accession,
                                   data = Rhizo.map)
image(X)
}
\author{
Sur from Dangl lab
}

