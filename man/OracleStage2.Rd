% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OracleStage2.R
\name{OracleStage2}
\alias{OracleStage2}
\title{Fit Second Stage of Oracle-2SLS}
\usage{
OracleStage2(
  gamma.hat.stage1,
  cor.Y2Z2,
  Estimated.Sigma,
  n1,
  n2,
  p,
  ind.stage2,
  posdef = 1e-04,
  Est.Sigma1Square,
  Est.Sigma2Square
)
}
\arguments{
\item{gamma.hat.stage1}{A vector of length p, containing estimated effects from IVs to exposure
obtained in the first stage model.}

\item{cor.Y2Z2}{A vector of length p, containing sample correlations (summary data) between
outcome Y and p IVs calculated with the second sample.}

\item{Estimated.Sigma}{Estimated correlation matrix with reference panel.}

\item{n1}{Sample size in first stage.}

\item{n2}{Sample size in second stage.}

\item{p}{Number of instruments.}

\item{ind.stage2}{Indices of truly invalid IVs (as this is the oracle model).}

\item{posdef}{A small positive constant to be added on diagonal of covariance matrix
to ensure it is positive definite.}

\item{Est.Sigma1Square}{A number, estimated variance of random error in the first stage.}

\item{Est.Sigma2Square}{A number, estimated variance of random error in the second stage.}
}
\value{
This function returns a list with two elements,
(1) betaalpha.hat.stage2: a vector of length (p+1), the first element is the estimated beta
(causal effect from exposure to outcome), the rest p elements are estimated alpha's
(direct effects from IVs to outcome);
(2) Asympt.Var.BetaHat: the uncorrected estimated variance for estimated beta.
}
\description{
This function is for fitting the second stage model of oracle-2SLS method, assuming the
truly set of invalid IVs is known.
}
