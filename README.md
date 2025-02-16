
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TScML

<!-- badges: start -->
<!-- badges: end -->

The goal of TScML is to perform the *Two-stage Constrained Maximum
Likelihood (2ScML)* method. In the following example, we apply 2ScML in
TWAS simulation to illustrate how to use our software.

## Installation

You can install the development version of TScML from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xue-hr/TScML")
```

## Load Packages

We need three R packages:

- **MASS**: This package is available on CRAN. We use `MASS::mvrnorm()`
  function to generate random variables.
- **devtools**: This package is available on CRAN. We use
  `devtools::install_github` to install package from GitHub.
- **lasso2**: This package is not available on CRAN for R 4.2.2, it is
  availale on GitHub at `https://github.com/cran/lasso2`. We use
  `lasso2::l1ce` to solve constrained lasso problem. Note that, for
  recent **R** versions (for example R 4.4.2), the `Sint` data type has
  been removed, which was used in `lasso2` package. Therefore, the
  `lasso2` package cannot be directly installed from GitHub for recent
  versions of R. To solve this issue, you can download the source file
  at
  `https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz`,
  then replace `Sint` with `int` in the two files `src/lasso.h` and
  `src/lasso.c`, and install the package from the source file.

``` r
if(!require("MASS"))
{
  install.packages("MASS")
  library(MASS)
}
#> Loading required package: MASS
if(!require("devtools"))
{
  install.packages("devtools")
  library(devtools)
}
#> Loading required package: devtools
#> Loading required package: usethis
if(!require("lasso2"))
{
  devtools::install_github("cran/lasso2")
  library(lasso2)
}
#> Loading required package: lasso2
#> R Package to solve regression problems while imposing
#>   an L1 constraint on the parameters. Based on S-plus Release 2.1
#> Copyright (C) 1998, 1999
#> Justin Lokhorst   <jlokhors@stats.adelaide.edu.au>
#> Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
#> Bill Venables     <wvenable@stats.adelaide.edu.au>
#> 
#> Copyright (C) 2002
#> Martin Maechler <maechler@stat.math.ethz.ch>
library(TScML)
```

## Simulation Setup

Parameters in our simulation are set as below:

``` r
load("MAFB_SNP_May16.Rdata") # data used in simulation
p = 56 # number of SNPs
n1 = 500 # first sample size
n2 = 50000 # second sample size
n.ref = 10000 # reference panel size
random.error.size = 2 # random error size
beta.true = 0 # true causal effect
gamma.vec = c(0,rep(1,7),rep(0,p-8))*1 # effects from SNPs to exposure
alpha.vec = c(1,rep(0,5),rep(1,3),rep(0,p-9))*1 # direct effects from SNPs to outcome
sigma.12 = 0.5 # correlation between error terms
Sigma.Err = 
  matrix(c(1,sigma.12,sigma.12,1),
         nrow = 2)*random.error.size # covariance matrix of errors
```

## Generate Simulated Data

We first generate individual level $Z$’s:

``` r
set.seed(123)
### generate all Z
Z.Stage1 = SNP_BED[sample(1:408339,n1,replace = T),]
Z.Stage2 = SNP_BED[sample(1:408339,n2,replace = T),]
```

Generate the first sample:

``` r
### generate stage 1 sample
Z = Z.Stage1
Random_Error = mvrnorm(n1, mu = c(0,0), Sigma = Sigma.Err)
D = Z%*%gamma.vec  + Random_Error[,1]
Y = beta.true*D + Z%*%alpha.vec + Random_Error[,2]
Z1 = scale(Z,scale = F)
D1 = scale(D,scale = F)
Y1 = scale(Y,scale = F)
```

Generate the second sample:

``` r
### generate stage 2 sample
Z = Z.Stage2
Random_Error = mvrnorm(n2, mu = c(0,0), Sigma = Sigma.Err)
D = Z%*%gamma.vec  + Random_Error[,1]
Y = beta.true*D + Z%*%alpha.vec + Random_Error[,2]
Z2 = scale(Z,scale = F)
D2 = scale(D,scale = F)
Y2 = scale(Y,scale = F)
```

Generate reference panel:

``` r
### generate reference panel
Z.ref = sample(1:408339,n.ref,replace = T)
Z.ref.original = SNP_BED[Z.ref,]
cor.Z.ref.original = cor(Z.ref.original) + diag(0.00001,p)
```

Calculate summary data:

``` r
### Input
cor.D1Z1.original = as.numeric(cor(D1,Z1))
cor.Y2Z2.original = as.numeric(cor(Y2,Z2))
```

## Apply 2ScML and Oracle Methods

We apply **2ScML** and **Oracle** methods to simulated data. We fit the
model in first stage:

``` r
### Stage1 with individual-level data
Stage1FittedModel = 
  TScMLStage1(cor.D1Z1 = cor.D1Z1.original,
              Cap.Sigma.stage1 = cor(Z1),
              n1 = n1,
              p = p,
              ind.stage1 = 2:8)
```

Apply the **Oracle** model in second stage:

``` r
### Oracle Stage 2 with summary data and reference panel
Est.Sigma1Square = 
  as.numeric(1 - cor.D1Z1.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.D1Z1.original)
Est.Sigma2Square = 
  as.numeric(1 - cor.Y2Z2.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.Y2Z2.original)
if(Est.Sigma1Square<=0)
{
  Est.Sigma1Square = 1
}
if(Est.Sigma2Square<=0)
{
  Est.Sigma2Square = 1
}
OracleStage2.ref =
  OracleStage2(gamma.hat.stage1 = Stage1FittedModel,
               cor.Y2Z2 = cor.Y2Z2.original,
               Estimated.Sigma = cor.Z.ref.original,
               n1 = n1,
               n2 = n2,
               p = p,
               ind.stage2 = c(1,2,8,9,10),
               Est.Sigma1Square = Est.Sigma1Square,
               Est.Sigma2Square = Est.Sigma2Square)
Oracle.Summary.Var = 
  TScMLVar(Z.ref.original = Z.ref.original,
           Stage1FittedModel = Stage1FittedModel,
           betaalpha.hat.stage2 = OracleStage2.ref$betaalpha.hat.stage2,
           Est.Sigma1Square = Est.Sigma1Square,
           Est.Sigma2Square = Est.Sigma2Square,
           n1 = n1,n2 = n2,n.ref = n.ref)
```

Apply 2ScML in the second stage:

``` r
### TScML Stage2 with summary data and reference panel
Est.Sigma1Square = 
  as.numeric(1 - cor.D1Z1.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.D1Z1.original)
Est.Sigma2Square = 
  as.numeric(1 - cor.Y2Z2.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.Y2Z2.original)
if(Est.Sigma1Square<=0)
{
  Est.Sigma1Square = 1
}
if(Est.Sigma2Square<=0)
{
  Est.Sigma2Square = 1
}
start.time = Sys.time()
TScMLStage2.Ref =
  TScMLStage2(gamma.hat.stage1 = Stage1FittedModel,
              cor.Y2Z2 = cor.Y2Z2.original,
              Estimated.Sigma = cor.Z.ref.original,
              n1 = n1,
              n2 = n2,
              p = p,
              K.vec.stage2 = 0:10,
              Est.Sigma1Square = Est.Sigma1Square,
              Est.Sigma2Square = Est.Sigma2Square)
end.time = Sys.time()
TScML.Summary.Var = 
  TScMLVar(Z.ref.original = Z.ref.original,
           Stage1FittedModel = Stage1FittedModel,
           betaalpha.hat.stage2 = TScMLStage2.Ref$betaalpha.hat.stage2,
           Est.Sigma1Square = Est.Sigma1Square,
           Est.Sigma2Square = Est.Sigma2Square,
           n1 = n1,n2 = n2,n.ref = n.ref)
```

We record the run time of our main function `TScMLStage2`. Different
from methods that use individual-level data, sample size does not
influence run time of our main function `TScMLStage2` as we only use
summary data. The run time of `TScMLStage2` mainly depends on two
values, the first one is the number of instruments, i.e. $p$; the second
one is the set of candidate $K$’s. Here we have $p = 56$ and
$K = 0,\cdots,10$, and the run time is about 0.3 second, which should be
adequately efficient:

``` r
run.time = end.time - start.time
run.time
#> Time difference of 0.231663 secs
```

Now we can show the results from **Oracle** and **2ScML**.

``` r
# Oracle Estimate
OracleStage2.ref$betaalpha.hat.stage2[1] 
#> [1] -0.01266214
# Uncorrected Variance of Oracle Estimate
OracleStage2.ref$Asympt.Var.BetaHat 
#> [1] 0.0001017406
# Corrected Variance of Oracle Estimate
Oracle.Summary.Var 
#> [1] 0.000291072
# 2ScML Estimate
TScMLStage2.Ref$betaalpha.hat.stage2[1] 
#> [1] -0.01266214
# Uncorrected Variance of 2ScML Estimate
TScMLStage2.Ref$Asympt.Var.BetaHat 
#> [1] 0.0001017406
# Corrected Variance of 2ScML Estimate
TScML.Summary.Var 
#> [1] 0.000291072
```

We can see, in this simulation, the proposed 2ScML method gives same
result as oracle method.

## R Session Information

Here is the R session information. We performed the example in the
latest R release 4.2.2.

``` r
sessionInfo()
#> R version 4.4.2 (2024-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Asia/Shanghai
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] TScML_0.0.0.9000 lasso2_1.2-22    devtools_2.4.5   usethis_3.1.0   
#> [5] MASS_7.3-61     
#> 
#> loaded via a namespace (and not attached):
#>  [1] miniUI_0.1.1.1    compiler_4.4.2    promises_1.3.2    Rcpp_1.0.14      
#>  [5] later_1.4.1       yaml_2.3.10       fastmap_1.2.0     mime_0.12        
#>  [9] R6_2.5.1          knitr_1.49        htmlwidgets_1.6.4 profvis_0.4.0    
#> [13] shiny_1.10.0      rlang_1.1.5       cachem_1.1.0      httpuv_1.6.15    
#> [17] xfun_0.50         fs_1.6.5          pkgload_1.4.0     memoise_2.0.1    
#> [21] cli_3.6.3         magrittr_2.0.3    digest_0.6.37     rstudioapi_0.17.1
#> [25] xtable_1.8-4      remotes_2.5.0     lifecycle_1.0.4   vctrs_0.6.5      
#> [29] evaluate_1.0.3    glue_1.8.0        urlchecker_1.0.1  sessioninfo_1.2.2
#> [33] pkgbuild_1.4.6    rmarkdown_2.29    purrr_1.0.2       tools_4.4.2      
#> [37] ellipsis_0.3.2    htmltools_0.5.8.1
```

## Reference

Haoran Xue, Xiaotong Shen & Wei Pan (2023) Causal Inference in
Transcriptome-Wide Association Studies with Invalid Instruments and GWAS
Summary Data, Journal of the American Statistical Association, DOI:
10.1080/01621459.2023.2183127

## Contact

Feel free to contact the author at <xuexx268@umn.edu> for any comments!
