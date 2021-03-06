% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spec.l1.solution.R
\name{spec.l1.solution}
\alias{spec.l1.solution}
\title{Get the solution that is a (partial) minimal (or/and) maximal l1 norm solution for some lambda, when the lasso penalty is included in the penalty matrix.}
\usage{
spec.l1.solution(
  object,
  lambda.seq,
  min.indices = c(),
  max.indices = c(),
  tol = 1e-10
)
}
\arguments{
\item{object}{ESPgenlasso object.}

\item{lambda.seq}{lambda seq.}

\item{min.indices}{beta coef indices that would expected to have the minimum l1 norm.}

\item{max.indices}{beta coef indices that would expected to have the maximum l1 norm.}

\item{tol}{tolerance.}
}
\value{
solutions for each lambda seq.
}
\description{
This function gets the solution that is a (partial) minimal (or/and) maximal l1 norm solution for some lambda, when the lasso penalty is included in the penalty matrix.
}
