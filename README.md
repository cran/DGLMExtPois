
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DGLMExtPois

<!-- badges: start -->

<!-- badges: end -->

DGLMExtPois is a package that contains statistical functions for the
model estimation, dispersion testing and diagnosis of hyper-Poisson and
Conway-Maxwell-Poisson regression models.

## Installation

You can install the released version of DGLMExtPois from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("DGLMExtPois")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(DGLMExtPois)
library(Ecdat)
#> Loading required package: Ecfun
#> 
#> Attaching package: 'Ecfun'
#> The following object is masked from 'package:base':
#> 
#>     sign
#> 
#> Attaching package: 'Ecdat'
#> The following object is masked from 'package:datasets':
#> 
#>     Orange
Bids$size.sq <- Bids$size ^ 2
hP.bids <- glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest +
                    whtknght + bidprem + insthold + size + size.sq +
                    regulatn,
                    formula.gamma = numbids ~ 1,
                    data = Bids)
summary(hP.bids)
#> 
#> Call:
#> glm.hP(formula.mu = numbids ~ leglrest + rearest + finrest + 
#>     whtknght + bidprem + insthold + size + size.sq + regulatn, 
#>     formula.gamma = numbids ~ 1, data = Bids)
#> 
#> Mean model coefficients (with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.042145   0.387027   2.693  0.00709 ** 
#> leglrest     0.240887   0.109638   2.197  0.02801 *  
#> rearest     -0.268646   0.144930  -1.854  0.06379 .  
#> finrest      0.104245   0.163049   0.639  0.52260    
#> whtknght     0.487929   0.110133   4.430 9.41e-06 ***
#> bidprem     -0.709086   0.273832  -2.589  0.00961 ** 
#> insthold    -0.363993   0.304749  -1.194  0.23232    
#> size         0.173023   0.048291   3.583  0.00034 ***
#> size.sq     -0.007371   0.002479  -2.973  0.00295 ** 
#> regulatn    -0.008751   0.118167  -0.074  0.94097    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Dispersion model coefficients (with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -2.6219     0.3489  -7.516 5.67e-14 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> AIC: 362.3
AIC(hP.bids)
#> [1] 362.3072
BIC(hP.bids)
#> [1] 393.5063
coef(hP.bids)
#> $mean_model
#>  (Intercept)     leglrest      rearest      finrest     whtknght 
#>  1.042145307  0.240886951 -0.268645993  0.104245100  0.487928599 
#>      bidprem     insthold         size      size.sq     regulatn 
#> -0.709086033 -0.363993485  0.173023482 -0.007370863 -0.008751012 
#> 
#> $dispersion_model
#> (Intercept) 
#>   -2.621855
confint(hP.bids)
#>                   2.5 %       97.5 %
#> (Intercept)  0.28358635  1.800704268
#> leglrest     0.02600000  0.455773898
#> rearest     -0.55270282  0.015410837
#> finrest     -0.21532507  0.423815275
#> whtknght     0.27207160  0.703785596
#> bidprem     -1.24578669 -0.172385379
#> insthold    -0.96129006  0.233303088
#> size         0.07837414  0.267672828
#> size.sq     -0.01223058 -0.002511151
#> regulatn    -0.24035456  0.222852539
head(fitted(hP.bids))
#>        1        2        3        4        5        6 
#> 2.733621 1.331997 2.196977 1.176840 1.231121 2.088129
head(residuals(hP.bids))
#>          1          2          3          4          5          6 
#> -0.5421896 -1.7307986 -1.0424542 -0.2501029 -0.3175111  0.8264357
```

``` r
CMP.bids <- glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest +
                    whtknght + bidprem + insthold + size + size.sq +
                    regulatn,
                    formula.nu = numbids ~ 1,
                    data = Bids)
summary(CMP.bids)
#> 
#> Call:
#> glm.CMP(formula.mu = numbids ~ leglrest + rearest + finrest + 
#>     whtknght + bidprem + insthold + size + size.sq + regulatn, 
#>     formula.nu = numbids ~ 1, data = Bids)
#> 
#> Mean model coefficients (with log link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.990004   0.435140   2.275 0.022898 *  
#> leglrest     0.267903   0.122808   2.181 0.029148 *  
#> rearest     -0.173273   0.154704  -1.120 0.262703    
#> finrest      0.067916   0.174300   0.390 0.696794    
#> whtknght     0.481172   0.131654   3.655 0.000257 ***
#> bidprem     -0.685007   0.307470  -2.228 0.025889 *  
#> insthold    -0.367923   0.346620  -1.061 0.288481    
#> size         0.179279   0.047604   3.766 0.000166 ***
#> size.sq     -0.007580   0.002483  -3.052 0.002270 ** 
#> regulatn    -0.037561   0.130235  -0.288 0.773031    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Dispersion model coefficients (with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.56208    0.04809   11.69   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> AIC: 382.2
AIC(CMP.bids)
#> [1] 382.1753
BIC(CMP.bids)
#> [1] 413.3744
coef(CMP.bids)
#> $mean_model
#>  (Intercept)     leglrest      rearest      finrest     whtknght 
#>  0.990003872  0.267902606 -0.173272510  0.067916284  0.481171553 
#>      bidprem     insthold         size      size.sq     regulatn 
#> -0.685006796 -0.367923118  0.179279126 -0.007580393 -0.037561467 
#> 
#> $dispersion_model
#> (Intercept) 
#>   0.5620821
confint(CMP.bids)
#>                   2.5 %       97.5 %
#> (Intercept)  0.13714473  1.842863012
#> leglrest     0.02720360  0.508601611
#> rearest     -0.47648689  0.129941873
#> finrest     -0.27370514  0.409537708
#> whtknght     0.22313460  0.739208506
#> bidprem     -1.28763746 -0.082376130
#> insthold    -1.04728529  0.311439049
#> size         0.08597638  0.272581868
#> size.sq     -0.01244781 -0.002712976
#> regulatn    -0.29281835  0.217695412
head(fitted(CMP.bids))
#>        1        2        3        4        5        6 
#> 2.736143 1.296889 2.140255 1.189666 1.205770 2.092416
head(residuals(CMP.bids))
#>          1          2          3          4          5          6 
#> -0.5646967 -1.3669047 -0.9753079 -0.2069397 -0.2233074  0.7839429
```
