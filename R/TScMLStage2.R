#' Fit Second Stage of 2ScML
#'
#' This function fitting the second stage of 2ScML is the main function of our proposed method.
#' It applies Truncated L1 Constraint in the second stage regression to select invalid instruments.
#'
#' @param gamma.hat.stage1 A vector of length p, containing estimated effects from IVs to exposure
#' obtained in the first stage model.
#' @param cor.Y2Z2 A vector of length p, containing sample correlations (summary data) between 
#' outcome Y and p IVs calculated with the second sample.
#' @param Estimated.Sigma Estimated correlation matrix with reference panel.
#' @param n1 Sample size in first stage.
#' @param n2 Sample size in second stage.
#' @param p Number of instruments.
#' @param K.vec.stage2 Set of candidate K's being used in the constraint.
#' @param posdef A small positive constant to be added on diagonal of covariance matrix 
#' to ensure it is positive definite.
#' @param Est.Sigma1Square A number, estimated variance of random error in the first stage. 
#' @param Est.Sigma2Square A number, estimated variance of random error in the second stage. 
#'
#' @return This function returns a list with three elements, 
#' (1) betaalpha.hat.stage2: a vector of length (p+1), the first element is the estimated beta 
#' (causal effect from exposure to outcome), the rest p elements are estimated alpha's 
#' (direct effects from IVs to outcome); 
#' (2) Asympt.Var.BetaHat: the uncorrected estimated variance for estimated beta.
#' (3) BIC.vec.stage2: the vector of BIC, having the same length as K.vec.stage2.
#' @export
#'
#' @examples
TScMLStage2 <- function(gamma.hat.stage1,cor.Y2Z2,
                        Estimated.Sigma,
                        n1,n2,p,
                        K.vec.stage2 = 0:(p/2),
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
  
  Est.BetaAlpha.mat = NULL
  BIC.vec.stage2 = NULL

  for(K.stage2 in K.vec.stage2)
  {
    Est.BetaAlpha = 
      TLC(Y = Response.stage2,
          X = Predictor.stage2,
          K = K.stage2,
          tlc_weight = c(0,rep(1,p)))
    Est.BetaAlpha.mat = 
      rbind(Est.BetaAlpha.mat,Est.BetaAlpha)
    BIC = as.numeric(1 - 2*sum(Est.BetaAlpha*vec.D2hatZ2Y2) + 
                       Est.BetaAlpha%*%mat.D2hatZ2%*%Est.BetaAlpha)
    BIC.vec.stage2 = c(BIC.vec.stage2,BIC)
  }
  BIC.vec.stage2 = n2*log(BIC.vec.stage2) + log(n2)*K.vec.stage2
  K.hat.stage2 = which.min(BIC.vec.stage2)
  betaalpha.hat.stage2 = Est.BetaAlpha.mat[K.hat.stage2,]
  
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
              Asympt.Var.BetaHat = Asympt.Var.BetaHat/n2,
              BIC.vec.stage2 = BIC.vec.stage2))
}