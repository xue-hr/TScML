
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TScML

<!-- badges: start -->

<!-- badges: end -->

R package for Two-Stage Constrained Maximum Likelihood (2ScML) method.

## Installation

Install the package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xue-hr/TScML")
```

## Example of One-Sample Case

Here we show how to apply the method to one-sample case. First we
generate simulated data.

``` r
library(TScML)
library(MASS)
n = 300 # sample size
pz = 100 # number of IVs
Beta_DY = 0 # true causal effect

# Covariance matrix of Z, AR(0.5)
Sigma_Z = matrix(0,pz,pz)
for(i in 1:(pz))
{
  for(j in 1:(pz))
  {
    Sigma_Z[i,j] = 0.5^(abs(i-j))
  }
}

# Covariance matrix of Errors
Sigma_error = matrix(c(1.5,0.75,0.75,1.5),nrow = 2)

# Effect Sizes
Gamma_ZD = c(0,rep(1,7),rep(0,pz-8))*0.5 # IV 2 to 8 are relevant
Pi_ZY = c(rep(0,6),rep(1,2),rep(0,pz-8))*0.5 # IV 7, 8 are invalid

set.seed(1)
# generate IV
Z = mvrnorm(n, mu = rep(0,pz), Sigma = Sigma_Z)

# randome error terms
Random_Error = mvrnorm(n, mu = c(0,0), Sigma = Sigma_error)

# generate exposure D and outcome Y
D = Z%*%Gamma_ZD  + Random_Error[,1]
Y = Beta_DY*D + Z%*%Pi_ZY  + Random_Error[,2]

Z = scale(Z,scale = F)
D = scale(D,scale = F)
Y = scale(Y,scale = F)
```

Now we perform 2ScML with function `OneSample_2ScML()`. As 7 IVs are
relevant, we select ![K\_1](https://latex.codecogs.com/png.latex?K_1
"K_1") from 5 to 10; 2 IVs are invalid, we select
![K\_2](https://latex.codecogs.com/png.latex?K_2 "K_2") from 0 to 5.

``` r
OneSample_2ScML(Y = Y,D = D,Z = Z,
                K1_vec = 5:10,K2_vec = 0:5,
                tau1 = 1e-5,
                tau2 = 1e-5)
#> $K1
#> [1] 7
#> 
#> $K2
#> [1] 2
#> 
#> $stage1_ind
#> [1] 2 3 4 5 6 7 8
#> 
#> $stage2_ind
#> [1] 7 8
#> 
#> $beta_est
#> [1] 0.01877348
#> 
#> $beta_se
#> [1] 0.04184548
```

We can see, in the first stage, BIC chooses the correct ![K\_1
= 7](https://latex.codecogs.com/png.latex?K_1%20%3D%207 "K_1 = 7"), and
2ScML correctly select the 7 relevant IVs,
![2^{nd}](https://latex.codecogs.com/png.latex?2%5E%7Bnd%7D "2^{nd}") to
![8^{th}](https://latex.codecogs.com/png.latex?8%5E%7Bth%7D "8^{th}");
in the second stage, BIC chooses the correct ![K\_2
= 2](https://latex.codecogs.com/png.latex?K_2%20%3D%202 "K_2 = 2"), and
2ScML correctly select the 2 invalid IVs,
![7^{th}](https://latex.codecogs.com/png.latex?7%5E%7Bth%7D "7^{th}")
and ![8^{th}](https://latex.codecogs.com/png.latex?8%5E%7Bth%7D
"8^{th}"). The estimated causal effect ![\\hat{\\beta}
= 0.0188](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbeta%7D%20%3D%200.0188
"\\hat{\\beta} = 0.0188"), with standard error ![se(\\hat{\\beta})
= 0.0418](https://latex.codecogs.com/png.latex?se%28%5Chat%7B%5Cbeta%7D%29%20%3D%200.0418
"se(\\hat{\\beta}) = 0.0418").

## Example of Two-Sample Case

Now we show how to apply the method to two-sample case. First we
generate simulated data.

``` r
n = 20000 # sample size
p = 30 # number of IVs
num_invalid = 0.3*p # 30% IVs are invalid
MAF = 0.3 # MAF for SNPs
theta_XY = -0.1 # causal effect from X to Y
sd_phi_ZU = 0.1 # standard deviation of phi_ZU, correlated direct (pleiotropic) effects


set.seed(1)
# generate gamma_ZX from truncated normal at 0.1
gamma_ZX = NULL
i = 0
while (i<p) {
  g = rnorm(1,0,0.2)
  if(abs(g)>0.1)
  {
    gamma_ZX = c(gamma_ZX,g)
    i = i + 1
  }
}

# generate alpha_ZY from normal
alpha_ZY = rnorm(num_invalid,0.5,0.075)

# generate phi_ZU from normal
phi_ZU = rnorm(num_invalid,0,sd_phi_ZU)

# generate first sample
Z1 = matrix(rbinom(n*p,2,MAF),nrow = n)
U1 = Z1[,1:num_invalid]%*%phi_ZU + rnorm(n,0,1)
X1 = Z1%*%gamma_ZX + U1 + rnorm(n,0,1)
Y1 = theta_XY*X1 + Z1[,1:num_invalid]%*%alpha_ZY + U1 + rnorm(n,0,1)

Z1 = scale(Z1,scale = F)
X1 = scale(X1,scale = F)
Y1 = scale(Y1,scale = F)

# generate second sample
Z2 = matrix(rbinom(n*p,2,MAF),nrow = n)
U2 = Z2[,1:num_invalid]%*%phi_ZU + rnorm(n,0,1)
X2 = Z2%*%gamma_ZX + U2 + rnorm(n,0,1)
Y2 = theta_XY*X2 + Z2[,1:num_invalid]%*%alpha_ZY + U2 + rnorm(n,0,1)

Z2 = scale(Z2,scale = F)
X2 = scale(X2,scale = F)
Y2 = scale(Y2,scale = F)
```

From the data generation we can see, all 30 IVs are relevant, and the
![1^{st}](https://latex.codecogs.com/png.latex?1%5E%7Bst%7D "1^{st}") to
![9^{th}](https://latex.codecogs.com/png.latex?9%5E%7Bth%7D "9^{th}")
are invalid. Now we perform first stage with linear regression.

``` r
FirstStage = lm(X1~0+Z1)

# estimated gamma
gamma_hat = as.numeric(FirstStage$coef) 

# covariance matrix of estimated gamma, as in Assumption 5
Theta0_hat = solve(t(Z1)%*%Z1/n) 

# estimated sigma2 square, as in Assumption 5
sigma2_hat = (summary(FirstStage)$sigma)^2 
```

Now we perform the second stage with function `TwoSample_2ScML()`.

``` r
TwoSample_2ScML(Y = Y2,Z = Z2,
                stage1_ind = 1:30,
                gamma_hat = gamma_hat,
                Theta0_hat = Theta0_hat,
                sigma2_hat = sigma2_hat,
                n2 = n,
                K2_vec = 7:11,tau2 = 1e-5)
#> $K2
#> [1] 9
#> 
#> $stage2_ind
#> [1] 1 2 3 4 5 6 7 8 9
#> 
#> $beta_est
#> [1] -0.09919649
#> 
#> $beta_se
#> [1] 0.01577231
```

We can see, in the second stage, BIC select correct ![K\_2
= 9](https://latex.codecogs.com/png.latex?K_2%20%3D%209 "K_2 = 9"), and
2ScML correctly select the 9 invalid IVs,
![1^{st}](https://latex.codecogs.com/png.latex?1%5E%7Bst%7D "1^{st}") to
![9^{th}](https://latex.codecogs.com/png.latex?9%5E%7Bth%7D "9^{th}").
The estimated causal effect is ![\\hat{\\beta} =
-0.0992](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cbeta%7D%20%3D%20-0.0992
"\\hat{\\beta} = -0.0992"), with standard error ![se(\\hat{\\beta})
= 0.0158](https://latex.codecogs.com/png.latex?se%28%5Chat%7B%5Cbeta%7D%29%20%3D%200.0158
"se(\\hat{\\beta}) = 0.0158").
