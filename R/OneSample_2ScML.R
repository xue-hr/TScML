#' 2ScML with One Sample
#'
#' Main function perform 2ScML with one sample, 
#' using BIC for model selection in both stages.
#'
#' @param Y Response vector of length n.
#' @param D Exposure vector of length n.
#' @param Z IV matrix, n rows and p columns, p columns correspond to p IVs.
#' @param K1_vec Vector of candidate K1's for stage 1.
#' @param K2_vec Vector of candidate K2's for stage 2.
#' @param tau1 Parameter tau in TLC for stage 1.
#' @param tau2 Parameter tau in TLC for stage 2.
#'
#' @return A list contains following elements. 
#' K1: BIC selected K1 in stage 1;
#' K2: BIC selected K2 in stage 2;
#' stage1_ind: indices of IVs selected as relevant in stage 1;
#' stage2_ind: indices of IVs selected as invalid in stage 2;
#' beta_est: estimate of causal effect of exposure to outcome;
#' beta_se: standard error of beta_est.
#' @export
#'
#' @examples
OneSample_2ScML <- function(Y,D,Z,
                            K1_vec,K2_vec,
                            tau1 = 1e-5,
                            tau2 = 1e-5)
{
  n = length(Y)
  p = ncol(Z)
  
  Y = scale(Y,scale = F)
  D = scale(D,scale = F)
  Z = scale(Z,scale = F)
  
  ### Stage1
  BIC_vec = NULL
  for(K in K1_vec)
  {
    Stage1_Fit = TLC(Y = D,X = Z,K = K,
                     tau = tau1)
    nonzero_ind = which(Stage1_Fit!=0)
    lm_stage1 = lm(D ~ 0 + Z[,nonzero_ind])
    BIC_vec = c(BIC_vec,
                n*log(sum((lm_stage1$residuals)^2)) + log(n)*length(nonzero_ind))
  }
  K1 = K1_vec[which.min(BIC_vec)]
  Stage1_Fit = TLC(Y = D,X = Z,K = K1,
                   tau = tau1)
  ind1 = which(Stage1_Fit!=0)
  lm_stage1 = lm(D~0+Z[,ind1])
  Dhat = predict(lm_stage1)
  
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
    OneSampleEstSE(Y,D,Z,ind1,ind2)
  return(list(K1 = K1, K2 = K2, 
              stage1_ind = ind1,
              stage2_ind = ind2,
              beta_est = TLP_Wald_test[1],
              beta_se = TLP_Wald_test[2]))
  
}