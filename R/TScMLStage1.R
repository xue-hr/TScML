#' Fit First Stage Model for 2ScML
#'
#' This function is used for fitting first stage model of 2ScML.
#'
#' @param cor.D1Z1 A vector of length p, containing sample correlations  
#' (calculated with the first sample) between exposure D and p instruments.
#' @param Cap.Sigma.stage1 Sample correlation matrix between p instruments, 
#' could be calculated either with the first sample individual-level instruments data,
#' or with the reference panel.
#' @param n1 First stage sample size.
#' @param p Number of instruments.
#' @param ind.stage1 Indices of instruments relevant to exposure D.
#' @param posdef A small positive constant to be added on diagonal of covariance matrix 
#' to ensure it is positive definite.
#'
#' @return  A vector of length p, containing estimated effects from IVs to exposure D.
#' @export
#'
#' @examples
TScMLStage1 <- function(cor.D1Z1,
                        Cap.Sigma.stage1,
                        n1,p,
                        ind.stage1,
                        posdef = 1e-7)
{
  ### Stage 1
  Cap.Sigma.stage1 = Cap.Sigma.stage1 + diag(posdef,nrow(Cap.Sigma.stage1))
  Eigen.Decompose.stage1 = eigen(Cap.Sigma.stage1)
  Cap.Sigma.sqrt.stage1 = 
    Eigen.Decompose.stage1$vectors%*%
    diag(sqrt(Eigen.Decompose.stage1$values))%*%
    t(Eigen.Decompose.stage1$vectors)
  Response.stage1 = solve(Cap.Sigma.sqrt.stage1, tol = 0)%*%cor.D1Z1
  Predictor.stage1 = Cap.Sigma.sqrt.stage1
  
  gamma.hat.stage1 = SubsetLM(Y = Response.stage1,
                              X = Predictor.stage1,
                              ind = ind.stage1)
  
  
  return(gamma.hat.stage1)
}
