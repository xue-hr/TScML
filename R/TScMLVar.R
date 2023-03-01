#' Corrected Variance for Estimated Causal Effect
#' 
#' This function is for estimating corrected variance of estimated causal effect (i.e. beta) 
#' obtained from second stage. As described in our draft, the correction accounts for the 
#' effect of using reference panel.
#'
#' @param Z.ref.original The original reference panel, which is a n.ref (sample size of reference panel)
#' by p (number of instruments) matrix.
#' @param Stage1FittedModel A vector of length p, containing estimated effects from IVs to exposure
#' obtained in the first stage model.
#' @param betaalpha.hat.stage2 A vector of length (p+1), the first element is the estimated beta 
#' (causal effect from exposure to outcome), the rest p elements are estimated alpha's 
#' (direct effects from IVs to outcome).
#' @param Est.Sigma1Square A number, estimated variance of random error in the first stage. 
#' @param Est.Sigma2Square A number, estimated variance of random error in the second stage. 
#' @param n1 Sample size in first stage.
#' @param n2 Sample size in second stage.
#' @param n.ref Sample size of reference panel.
#'
#' @return A number, which is the estimated corrected variance of estimated causal effect.
#' @export
#'
#' @examples
TScMLVar <- function(Z.ref.original,
                     Stage1FittedModel,
                     betaalpha.hat.stage2,
                     Est.Sigma1Square,
                     Est.Sigma2Square,
                     n1,n2,n.ref)
{
  scale.Z.ref = scale(Z.ref.original)
  set.A = which(Stage1FittedModel!=0)
  set.B = which(betaalpha.hat.stage2!=0)[-1]-1
  hat.gamma.A = Stage1FittedModel[set.A]
  U.AA = cov(scale.Z.ref[,set.A])
  if(length(set.B) == 1)
  {
    U.BB = diag(1)
  } else {
    U.BB = cov(scale.Z.ref[,set.B])
  }
  U.BA = cov(scale.Z.ref[,set.B],scale.Z.ref[,set.A])
  U.AB = t(U.BA)
  Est.Theta = solve(cov(scale.Z.ref[,set.A]) + diag(0.00001,length(set.A)))
  Est.beta = betaalpha.hat.stage2[1]
  Est.alpha = (betaalpha.hat.stage2[-1])[set.B]
  
  Est.Sigma = rbind(cbind(hat.gamma.A%*%U.AA%*%hat.gamma.A,hat.gamma.A%*%U.AB),
                    cbind(U.BA%*%hat.gamma.A,U.BB))
  Est.Sigma = Est.Sigma + diag(0.00001,nrow(Est.Sigma))
  inv.Est.Sigma = solve(Est.Sigma)
  
  Matrix.A = (t(t(c(Est.beta,Est.alpha)))%*%t(c(Est.beta,Est.alpha)))
  Matrix.B = 
    n2*Est.Sigma2Square*Est.Sigma + 
    Est.beta^2*Est.Sigma1Square/n1*
    ((n2^2+n2)*
       rbind(cbind(hat.gamma.A%*%U.AA%*%Est.Theta%*%U.AA%*%hat.gamma.A,
                   hat.gamma.A%*%U.AA%*%Est.Theta%*%U.AB),
             cbind(U.BA%*%Est.Theta%*%U.AA%*%hat.gamma.A,U.BA%*%Est.Theta%*%U.AB))+
       n2*sum(diag(U.AA%*%Est.Theta))*
       rbind(cbind(hat.gamma.A%*%U.AA%*%hat.gamma.A,hat.gamma.A%*%U.AB),
             cbind(U.BA%*%hat.gamma.A,U.BB)))
  
  const.p = length(set.B)+1
  constant.d = (n.ref-const.p)*(n.ref-const.p-1)*(n.ref-const.p-3)
  
  Est.Cov.Mat = 
    ((n.ref - const.p - 1)/constant.d*(n2+n2^2)*Matrix.A + 
       n2*(n2+n.ref)/constant.d*sum(diag(Matrix.A%*%Est.Sigma))*inv.Est.Sigma)*n.ref^2/n2^2 +
    n.ref^2/n2^2*(
      1/constant.d*sum(diag(inv.Est.Sigma%*%Matrix.B))*inv.Est.Sigma+
        (n.ref - const.p - 1)/constant.d*(inv.Est.Sigma%*%Matrix.B%*%inv.Est.Sigma)) - 
    n.ref^2/(n.ref - const.p - 1)^2*Matrix.A
  return(Est.Cov.Mat[1,1])
}

