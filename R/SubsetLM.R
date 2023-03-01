#' Fit Linear Model on A Subset of Design Matrix
#'
#' This function is for fitting a linear model of response Y on a subset of design matrix 
#' X.
#'
#' @param Y Response vector.
#' @param X Full design matrix.
#' @param ind Set of columns of X to be used in the linear model.
#'
#' @return beta: a vector of estimated coefficients, with non-zero elements for variables being 
#' used and 0's for others.
#' @export
#'
#' @examples
SubsetLM <- function(Y,X,ind)
{
  n = length(Y)
  p = ncol(X)
  
  nonzero_ind = ind
  if(length(nonzero_ind) == 0)
  {
    lm1 = lm(Y ~ 0)
  } else{
    lm1 = lm(Y ~ 0 + X[,nonzero_ind])
  }
  beta1 = lm1$coefficients
  names(beta1) = NULL
  beta = rep(0,p)
  beta[nonzero_ind] = beta1
  
  return(beta)
  
}
