#' Two Sample Estimate and Standard Error of Causal Effect
#'
#' @param Y Response vector of length n.
#' @param Z IV matrix, n rows and p columns, p columns correspond to p IVs.
#' @param stage1_ind Indices of IVs selected as relevant in stage 1.
#' @param gamma_hat Estimated gamma in the first stage.
#' @param Theta0_hat Asymptotic covariance matrix of gamma_hat as in Assumption 5.
#' @param sigma2_hat Estimated sigma2 square as in Assumption 5.
#' @param n2 Sample size of the second sample used for estimating gamma.
#' @param stage2_ind Indices of IVs selected as invalid in stage 2.
#'
#' @return A vector of two elements, 
#' the first one is estimated beta,
#' the second one is its standard error.
#' @export
#'
#' @examples
TwoSampleEstSE <- function(Y,Z,stage1_ind,gamma_hat,
                           Theta0_hat,sigma2_hat,n2,
                           stage2_ind)
{
  n = length(Y)
  
  Dhat = Z[,stage1_ind]%*%gamma_hat
  lm_stage2 = summary(lm(Y~cbind(Dhat,Z[,stage2_ind])))
  beta_est = lm_stage2$coef[2,1]
  
  X = cbind(Dhat,Z[,stage2_ind])
  inv_Cap_Sigma = solve(t(X)%*%X/n)
  
  MM = t(X)%*%Z[,stage1_ind]/n
  Cap_Psi = MM %*% Theta0_hat %*% t(MM)
  
  v1hat = sum((Y - X%*%(lm_stage2$coefficients[-1,1]))^2)/n
  vhat = 
    v1hat*inv_Cap_Sigma + 
    n/n2*beta_est^2*sigma2_hat*(inv_Cap_Sigma%*%Cap_Psi%*%inv_Cap_Sigma)
  return(c(beta_est,sqrt(vhat[1,1]/n)))
}