% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TLC.R
\name{TLC}
\alias{TLC}
\title{Solve Truncated L1 Constrained Problem}
\usage{
TLC(
  Y,
  X,
  K,
  tau = 1e-05,
  tlc_weight = rep(1, ncol(X)),
  maxit_tlc = 100,
  tol_tlc = 1e-05
)
}
\arguments{
\item{Y}{Response vector of length n.}

\item{X}{Design matrix with n rows and p columns.}

\item{K}{Constraint parameter in TLC.}

\item{tau}{Parameter tau in TLC.}

\item{tlc_weight}{Length p vector of weights corresponding to p variables,
each element is either 1 (with constraint) or 0 (without constraint).
Default are all 1's.}

\item{maxit_tlc}{Maximum number of DC iteration.}

\item{tol_tlc}{Convergence tolerance for TLC.}
}
\value{
beta: a vector of length p containing estimated coefficients.
}
\description{
Using Difference Convex (DC) method to solve minimization problem with Truncated L1 Constraint.
}
