% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.u.R
\name{get.u}
\alias{get.u}
\title{Get the dual solution of the genlasso problem at a fixed lambda.}
\usage{
get.u(object, lambda)
}
\arguments{
\item{object}{ESPgenlasso object.}

\item{lambda}{a specific lambda.}
}
\value{
the dual solution at lambda.
}
\description{
This function gets the dual solution of the genlasso problem at a fixed lambda.
}
\examples{
y <- matrix(c(1000), nrow = 1)
X <- matrix(c(1,1,0),nrow = 1)
D <- matrix(c(1,0,0,0,1,1,0,0,-1), nrow = 3)
object <- ESPgenlasso(y,X,D,genlasso.option = FALSE)
get.u(object, 700)
}
