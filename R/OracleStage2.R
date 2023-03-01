#' Fit Second Stage of Oracle-2SLS
#' 
#' This function is for fitting the second stage model of oracle-2SLS method, assuming the 
#' truly set of invalid IVs is known.
#'
#' @param gamma.hat.stage1 A vector of length p, containing estimated effects from IVs to exposure
#' obtained in the first stage model.
#' @param cor.Y2Z2 A vector of length p, containing sample correlations (summary data) between 
#' outcome Y and p IVs calculated with the second sample.
#' @param Estimated.Sigma Estimated correlation matrix with reference panel.
#' @param n1 Sample size in first stage.
#' @param n2 Sample size in second stage.
#' @param p Number of instruments.
#' @param ind.stage2 Indices of truly invalid IVs (as this is the oracle model).
#' @param posdef A small positive constant to be added on diagonal of covariance matrix 
#' to ensure it is positive definite.
#' @param Est.Sigma1Square A number, estimated variance of random error in the first stage. 
#' @param Est.Sigma2Square A number, estimated variance of random error in the second stage.
#'
#' @return This function returns a list with two elements, 
#' (1) betaalpha.hat.stage2: a vector of length (p+1), the first element is the estimated beta 
#' (causal effect from exposure to outcome), the rest p elements are estimated alpha's 
#' (direct effects from IVs to outcome); 
#' (2) Asympt.Var.BetaHat: the uncorrected estimated variance for estimated beta.
#' @export
#'
#' @examples
OracleStage2 <- function(gamma.hat.stage1,cor.Y2Z2,
                         Estimated.Sigma,
                         n1,n2,p,
                         ind.stage2,
                         posdef = 0.0001,
                         Est.Sigma1Square,
                         Est.Sigma2Square)
{

  vec.D2hatZ2Y2 = 
    c(sum(cor.Y2Z2*gamma.hat.stage1),
      cor.Y2Z2)
  
  vec.D2hatZ2 = as.numeric(Estimated.Sigma%*%gamma.hat.stage1)
  
  num.D2hatD2hat = as.numeric(gamma.hat.stage1%*%Estimated.Sigma%*%gamma.hat.stage1)
  mat.D2hatZ2 = 
    rbind(c(num.D2hatD2hat,vec.D2hatZ2),
          cbind(vec.D2hatZ2,Estimated.Sigma))
  mat.D2hatZ2 = mat.D2hatZ2 + diag(p+1)*posdef
  
  Eigen.Decompose.stage2 = eigen(mat.D2hatZ2)
  Cap.Sigma.sqrt.stage2 = 
    Eigen.Decompose.stage2$vectors%*%
    diag(sqrt(Eigen.Decompose.stage2$values))%*%
    t(Eigen.Decompose.stage2$vectors)
  Response.stage2 = 
    as.numeric(solve(Cap.Sigma.sqrt.stage2, tol = 0)%*%vec.D2hatZ2Y2)
  Predictor.stage2 = Cap.Sigma.sqrt.stage2
  
  
  betaalpha.hat.stage2 = SubsetLM(Y = Response.stage2,
                                  X = Predictor.stage2,
                                  ind = ind.stage2)
  
    
  ### Calculate Variance
  stage1.nonzero.ind = which(gamma.hat.stage1!=0)
  stage2.nonzero.ind = which(betaalpha.hat.stage2!=0)[-1]-1
  
  Cap.Sigma.Var = 
    rbind(c(gamma.hat.stage1%*%Estimated.Sigma%*%gamma.hat.stage1,
            gamma.hat.stage1%*%Estimated.Sigma[,stage2.nonzero.ind]),
          cbind(c(gamma.hat.stage1%*%Estimated.Sigma[,stage2.nonzero.ind]),
                Estimated.Sigma[stage2.nonzero.ind,stage2.nonzero.ind]))
  Inv.Cap.Sigma.Var = solve(Cap.Sigma.Var, tol = 0)
  if(length(stage1.nonzero.ind) == 1)
  {
    Cap.Psi.Part1 = cbind(Estimated.Sigma[stage1.nonzero.ind,]%*%gamma.hat.stage1,
                          t(Estimated.Sigma[stage1.nonzero.ind,stage2.nonzero.ind]))
  } else {
    Cap.Psi.Part1 = cbind(Estimated.Sigma[stage1.nonzero.ind,]%*%gamma.hat.stage1,
                          Estimated.Sigma[stage1.nonzero.ind,stage2.nonzero.ind])
    
  }
  
  Cap.Psi.Part2 = solve(Estimated.Sigma[stage1.nonzero.ind,stage1.nonzero.ind], tol = 0)
  Cap.Psi.Var = t(Cap.Psi.Part1)%*%Cap.Psi.Part2%*%Cap.Psi.Part1
  
  Asympt.Var.BetaHat = 
    Inv.Cap.Sigma.Var[1,1]*Est.Sigma2Square + 
    n2/n1*betaalpha.hat.stage2[1]^2*Est.Sigma1Square*
    (Inv.Cap.Sigma.Var%*%Cap.Psi.Var%*%Inv.Cap.Sigma.Var)[1,1]
  
  return(list(betaalpha.hat.stage2 = betaalpha.hat.stage2,
              Asympt.Var.BetaHat = Asympt.Var.BetaHat/n2
  ))
}
