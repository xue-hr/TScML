#' One Sample Estimate and Standard Error of Causal Effect
#' 
#' For one sample case, given selected relevant IVs in stage 1 and 
#' selected invalid IVs in stage 2,
#' get estimate of beta, which is the causal effect from exposure to outcome, 
#' and its standard error.
#'
#' @param Y Response vector of length n.
#' @param D Exposure vector of length n.
#' @param Z IV matrix, n rows and p columns, p columns correspond to p IVs.
#' @param stage1_ind Indices of IVs selected as relevant in stage 1.
#' @param stage2_ind Indices of IVs selected as invalid in stage 2.
#'
#' @return A vector of two elements, 
#' the first one is estimated beta,
#' the second one is its standard error.
#' @export
#'
#' @examples
OneSampleEstSE <- function(Y,D,Z,stage1_ind,stage2_ind)
{
  n = length(Y)
  lm_stage1 = lm(D~Z[,stage1_ind])
  Dhat = predict(lm_stage1)
  
  lm_stage2 = summary(lm(Y~cbind(Dhat,Z[,stage2_ind])))
  
  X = cbind(Dhat,Z[,stage2_ind])
  inv_Cap_Sigma = solve((t(X)%*%X)/n)
  ZA = Z[,stage1_ind]
  PZA = ZA%*%solve(t(ZA)%*%ZA)%*%t(ZA)
  Cap_Psi = t(X)%*%PZA%*%X/n
  
  Xstar = cbind(D,Z[,stage2_ind])
  v1hat = sum((Y - X%*%(lm_stage2$coefficients[-1,1]))^2)/n
  v2hat = sum((Y - Xstar%*%(lm_stage2$coefficients[-1,1]))^2)/n
  vhat = v1hat*inv_Cap_Sigma - 
    (v1hat-v2hat)*((inv_Cap_Sigma%*%Cap_Psi%*%inv_Cap_Sigma))
  return(c(lm_stage2$coefficients[2,1],sqrt(vhat[1,1]/n)))
}


