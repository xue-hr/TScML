#' 2ScML with Two Sample
#'
#' Main function perform 2ScML with two sample, 
#' using BIC for model selection in the second stage.
#'
#' @param Y Response vector of length n.
#' @param Z IV matrix, n rows and p columns, p columns correspond to p IVs.
#' @param stage1_ind Indices of relevant IVs in the first stage.
#' @param gamma_hat Estimated gamma in the first stage.
#' @param Theta0_hat Asymptotic covariance matrix of gamma_hat as in Assumption 5.
#' @param sigma2_hat Estimated sigma2 square as in Assumption 5.
#' @param n2 Sample size of the second sample used for estimating gamma.
#' @param K2_vec Vector of candidate K2's for stage 2.
#' @param tau2 Parameter tau in TLC for stage 2.
#'
#' @return A list contains following elements. 
#' K2: BIC selected K2 in stage 2;
#' stage2_ind: indices of IVs selected as invalid in stage 2;
#' beta_est: estimate of causal effect of exposure to outcome;
#' beta_se: standard error of beta_est.
#' @export
#'
#' @examples
TwoSample_2ScML <- function(Y,Z,stage1_ind,gamma_hat,
                            Theta0_hat,sigma2_hat,n2,
                            K2_vec,tau2 = 1e-5)
{
  
  n = length(Y)
  p = ncol(Z)
  
  Y = scale(Y,scale = F)
  Z = scale(Z,scale = F)
  
  Dhat = Z[,stage1_ind]%*%gamma_hat
  
  ### Stage2
  BIC_vec = NULL
  for(K in K2_vec)
  {
    Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K,
                     tlc_weight = c(0,rep(1,p)),
                     tau = tau2)
    nonzero_ind = which(Stage2_Fit!=0)
    lm_stage2 = lm(Y ~ 0 + cbind(Dhat,Z)[,nonzero_ind])
    BIC_vec = c(BIC_vec,
                n*log(sum((lm_stage2$residuals)^2)) + log(n)*length(nonzero_ind))
  }
  K2 = K2_vec[which.min(BIC_vec)]
  Stage2_Fit = TLC(Y = Y,X = cbind(Dhat,Z),K = K2,
                   tlc_weight = c(0,rep(1,p)),
                   tau = tau2)
  ind2 = which(Stage2_Fit!=0)[-1]-1
  
  ### result
  TLP_Wald_test = 
    TwoSampleEstSE(Y,Z,stage1_ind,gamma_hat,
                   Theta0_hat,sigma2_hat,n2,
                   ind2)
  
  return(list(K2 = K2, 
              stage2_ind = ind2,
              beta_est = TLP_Wald_test[1],
              beta_se = TLP_Wald_test[2]))
  
}