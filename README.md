
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BOP2-FE

<!-- badges: start -->
<!-- badges: end -->

The primary purpose of an oncology single-arm trial is to evaluate the
effectiveness of anticancer agents and make a go / no-go decision while
maintaining patient safety. This R-package implements a flexible
Bayesian optimal phase II design with futility and efficacy stopping
boundaries for single-arm clinical trials, named the BOP2-FE design,
proposed by Xu et al (under review). The proposed BOP2-FE design allows
for early stopping of efficacy when the observed antitumor effect is
sufficiently higher than the null hypothesis value in the interim looks
and retains the benefits of the original BOP2 design, such as explicitly
controlling the type I error rate while maximizing power, accommodating
different types of endpoint, flexible number of interim looks, and
stopping boundaries calculated before the start of the trial. The
package handles multiple endpoints including binary, nested, co-primary
and joint efficacy and toxicity.

## Installation

You can install the development version of BOP2FE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("belayb/BOP2FE")
```

## Example

This is a basic example which shows you how to use BOP2FE:

``` r
library(BOP2FE)
## basic example code
```

### Binary endpoint

``` r
test_binary <- BOP2FE_binary(
 H0=0.2, H1= 0.4,
 n = c(10, 5, 5, 5, 5, 5, 5),
 nsim = 1000, t1e = 0.1, method = "power",
 lambda1 = 0, lambda2 = 1, grid1 = 11,
 gamma1 = 0, gamma2 = 1, grid2 = 11,
 eta1 = 0, eta2 = 3, grid3 = 31,
 seed = 123
)
summary(test_binary)
#> $design_pars
#> $design_pars$H0
#> [1] 0.2
#> 
#> $design_pars$H1
#> [1] 0.4
#> 
#> $design_pars$n
#> [1] 10  5  5  5  5  5  5
#> 
#> $design_pars$nsim
#> [1] 1000
#> 
#> $design_pars$t1e
#> [1] 0.1
#> 
#> $design_pars$method
#> [1] "power"
#> 
#> 
#> $opt_pars
#>   lambda gamma eta
#> 1    0.9     1 2.1
#> 
#> $boundary
#>                   IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary   1   2   3   5   7   9  11
#> Efficacy boundary   6   8   9  10  10  11  12
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.868    0.095
#> Early stop for efficacy    0.090    0.882
#> Average sample size       20.195   22.085
#> Reject null                0.099    0.897
#plot(test_binary)
```

### Nested endpoint

``` r
test_nested <- search_optimal_pars_nested(
 H0=c(0.05,0.05, 0.15, 0.75),
 H1= c(0.15,0.15, 0.20, 0.50),
 n = c(10, 5, 5, 5, 5, 5, 5),
 nsim = 1000, t1e = 0.1, method = "power",
 lambda1 = 0, lambda2 = 1, grid1 = 11,
 gamma1 = 0, gamma2 = 1, grid2 = 11,
 eta1 = 0, eta2 = 3, grid3 = 31,
 seed = 123
)
summary(test_nested)
#>  fut_boundary_CR1 fut_boundary_CR2 fut_boundary_CR3 fut_boundary_CR4
#>  Min.   : 5.000   Min.   : 8.00    Min.   :12.00    Min.   :16.00   
#>  1st Qu.: 6.000   1st Qu.:10.00    1st Qu.:14.00    1st Qu.:18.00   
#>  Median : 6.000   Median :11.00    Median :15.00    Median :19.00   
#>  Mean   : 6.472   Mean   :10.45    Mean   :14.54    Mean   :18.59   
#>  3rd Qu.: 7.000   3rd Qu.:11.00    3rd Qu.:15.00    3rd Qu.:20.00   
#>  Max.   :10.000   Max.   :15.00    Max.   :20.00    Max.   :25.00   
#>  fut_boundary_CR5 fut_boundary_CR6 fut_boundary_CR7 fut_boundary_CR/PR1
#>  Min.   :20.00    Min.   :24.00    Min.   :28.00    Min.   : 5.000     
#>  1st Qu.:21.00    1st Qu.:25.00    1st Qu.:29.00    1st Qu.: 6.000     
#>  Median :23.00    Median :27.00    Median :31.00    Median : 7.000     
#>  Mean   :22.72    Mean   :26.86    Mean   :31.02    Mean   : 6.923     
#>  3rd Qu.:24.00    3rd Qu.:28.00    3rd Qu.:33.00    3rd Qu.: 7.750     
#>  Max.   :30.00    Max.   :35.00    Max.   :40.00    Max.   :10.000     
#>  fut_boundary_CR/PR2 fut_boundary_CR/PR3 fut_boundary_CR/PR4
#>  Min.   : 9.00       Min.   :14.00       Min.   :18.00      
#>  1st Qu.:10.00       1st Qu.:15.00       1st Qu.:19.00      
#>  Median :11.00       Median :16.00       Median :20.00      
#>  Mean   :11.18       Mean   :15.58       Mean   :19.92      
#>  3rd Qu.:12.00       3rd Qu.:16.00       3rd Qu.:21.00      
#>  Max.   :15.00       Max.   :20.00       Max.   :25.00      
#>  fut_boundary_CR/PR5 fut_boundary_CR/PR6 fut_boundary_CR/PR7 sup_boundary_CR1
#>  Min.   :22.00       Min.   :26.00       Min.   :30          Min.   : 7.000  
#>  1st Qu.:23.00       1st Qu.:27.00       1st Qu.:31          1st Qu.: 9.000  
#>  Median :24.00       Median :29.00       Median :33          Median :10.000  
#>  Mean   :24.22       Mean   :28.61       Mean   :33          Mean   : 9.562  
#>  3rd Qu.:25.00       3rd Qu.:30.00       3rd Qu.:35          3rd Qu.:10.000  
#>  Max.   :30.00       Max.   :35.00       Max.   :40          Max.   :10.000  
#>  sup_boundary_CR2 sup_boundary_CR3 sup_boundary_CR4 sup_boundary_CR5
#>  Min.   :10.00    Min.   :14.00    Min.   :18.00    Min.   :22.00   
#>  1st Qu.:13.00    1st Qu.:17.00    1st Qu.:21.00    1st Qu.:24.00   
#>  Median :14.00    Median :18.00    Median :22.00    Median :25.00   
#>  Mean   :13.76    Mean   :17.69    Mean   :21.47    Mean   :25.18   
#>  3rd Qu.:15.00    3rd Qu.:19.00    3rd Qu.:22.00    3rd Qu.:26.00   
#>  Max.   :15.00    Max.   :20.00    Max.   :25.00    Max.   :30.00   
#>  sup_boundary_CR6 sup_boundary_CR7 sup_boundary_CR/PR1 sup_boundary_CR/PR2
#>  Min.   :25.0     Min.   :29.00    Min.   : 7.000      Min.   :11.00      
#>  1st Qu.:28.0     1st Qu.:30.00    1st Qu.:10.000      1st Qu.:14.00      
#>  Median :29.0     Median :32.00    Median :10.000      Median :14.00      
#>  Mean   :28.8     Mean   :32.01    Mean   : 9.733      Mean   :14.23      
#>  3rd Qu.:30.0     3rd Qu.:34.00    3rd Qu.:10.000      3rd Qu.:15.00      
#>  Max.   :35.0     Max.   :40.00    Max.   :10.000      Max.   :15.00      
#>  sup_boundary_CR/PR3 sup_boundary_CR/PR4 sup_boundary_CR/PR5
#>  Min.   :15.0        Min.   :19.00       Min.   :23.00      
#>  1st Qu.:18.0        1st Qu.:22.00       1st Qu.:26.00      
#>  Median :19.0        Median :23.00       Median :26.00      
#>  Mean   :18.4        Mean   :22.52       Mean   :26.52      
#>  3rd Qu.:19.0        3rd Qu.:23.00       3rd Qu.:27.00      
#>  Max.   :20.0        Max.   :25.00       Max.   :30.00      
#>  sup_boundary_CR/PR6 sup_boundary_CR/PR7 earlystopfuti_mean_h0
#>  Min.   :27.0        Min.   :31.00       Min.   :0.9390       
#>  1st Qu.:29.0        1st Qu.:32.00       1st Qu.:1.0000       
#>  Median :30.0        Median :34.00       Median :1.0000       
#>  Mean   :30.4        Mean   :33.99       Mean   :0.9989       
#>  3rd Qu.:31.0        3rd Qu.:36.00       3rd Qu.:1.0000       
#>  Max.   :35.0        Max.   :40.00       Max.   :1.0000       
#>  earlystopsupe_mean_h0   ss_mean_h0    rejectnull_mean_h0 earlystopfuti_mean_h1
#>  Min.   :0.000000      Min.   :10.00   Min.   :0.000000   Min.   :0.5800       
#>  1st Qu.:0.000000      1st Qu.:10.02   1st Qu.:0.000000   1st Qu.:0.9650       
#>  Median :0.000000      Median :10.07   Median :0.000000   Median :0.9830       
#>  Mean   :0.001131      Mean   :10.16   Mean   :0.001131   Mean   :0.9662       
#>  3rd Qu.:0.000000      3rd Qu.:10.30   3rd Qu.:0.000000   3rd Qu.:0.9910       
#>  Max.   :0.061000      Max.   :11.03   Max.   :0.061000   Max.   :0.9910       
#>  earlystopsupe_mean_h1   ss_mean_h1    rejectnull_mean_h1     lambda      
#>  Min.   :0.00000       Min.   :10.00   Min.   :0.00000    Min.   :0.1000  
#>  1st Qu.:0.00900       1st Qu.:10.34   1st Qu.:0.00900    1st Qu.:0.2000  
#>  Median :0.01300       Median :11.38   Median :0.01700    Median :0.4000  
#>  Mean   :0.03215       Mean   :11.66   Mean   :0.03261    Mean   :0.4644  
#>  3rd Qu.:0.02100       3rd Qu.:12.67   3rd Qu.:0.02700    3rd Qu.:0.7000  
#>  Max.   :0.42000       Max.   :16.36   Max.   :0.42000    Max.   :1.0000  
#>      gamma             eta       
#>  Min.   :0.0000   Min.   :0.000  
#>  1st Qu.:0.2000   1st Qu.:0.400  
#>  Median :0.5000   Median :1.100  
#>  Mean   :0.4844   Mean   :1.192  
#>  3rd Qu.:0.8000   3rd Qu.:1.800  
#>  Max.   :1.0000   Max.   :3.000
#plot(test_nested)
```

### Co-primary endpoint

``` r
test_coprimary <- BOP2FE_coprimary(
 H0=c(0.05,0.05, 0.15, 0.75),
 H1= c(0.15,0.15, 0.20, 0.50),
 n = c(10, 5, 5, 5, 5, 5, 5),
 nsim = 1000, t1e = 0.1, method = "power",
 lambda1 = 0, lambda2 = 1, grid1 = 11,
 gamma1 = 0, gamma2 = 1, grid2 = 11,
 eta1 = 0, eta2 = 3, grid3 = 31,
 seed = 123
)
summary(test_coprimary)
#> $design_pars
#> $design_pars$H0
#> [1] 0.05 0.05 0.15 0.75
#> 
#> $design_pars$H1
#> [1] 0.15 0.15 0.20 0.50
#> 
#> $design_pars$n
#> [1] 10  5  5  5  5  5  5
#> 
#> $design_pars$nsim
#> [1] 1000
#> 
#> $design_pars$t1e
#> [1] 0.1
#> 
#> $design_pars$method
#> [1] "power"
#> 
#> 
#> $opt_pars
#>   lambda gamma eta
#> 1    0.9   0.4 2.6
#> 
#> $boundary
#>                          IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary (OR)     1   2   2   3   4   5   6
#> Futility boundary (PFS6)   2   3   5   6   8   9  11
#> Efficacy boundary (OR)     5   6   6   7   7   7   7
#> Efficacy boundary (PFS6)   7   8   9  10  11  11  12
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.880    0.088
#> Early stop for efficacy    0.082    0.904
#> Average sample size       17.250   19.225
#> Reject null                0.095    0.910
#plot(test_coprimary)
```

### Joint endpoint

``` r
test_joint <- BOP2FE_jointefftox(
 H0=c(0.15,0.30, 0.15, 0.40),
 H1= c(0.18,0.42, 0.02, 0.38),
 n = c(10, 5, 5, 5, 5, 5, 5),
 nsim = 1000, t1e = 0.1, method = "power",
 lambda1 = 0, lambda2 = 1, grid1 = 11,
 gamma1 = 0, gamma2 = 1, grid2 = 11,
 eta1 = 0, eta2 = 3, grid3 = 31,
 seed = 123
)
summary(test_joint)
#> $design_pars
#> $design_pars$H0
#> [1] 0.15 0.30 0.15 0.40
#> 
#> $design_pars$H1
#> [1] 0.18 0.42 0.02 0.38
#> 
#> $design_pars$n
#> [1] 10  5  5  5  5  5  5
#> 
#> $design_pars$nsim
#> [1] 1000
#> 
#> $design_pars$t1e
#> [1] 0.1
#> 
#> $design_pars$method
#> [1] "power"
#> 
#> 
#> $opt_pars
#>   lambda gamma eta
#> 1    0.6   0.9 2.7
#> 
#> $boundary
#>                         IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary (OR)    3   5   8  10  13  16  18
#> Futility boundary (Tox)   5   6   8   9  10  11  12
#> Efficacy boundary (OR)    9  11  13  15  17  18  19
#> Efficacy boundary (Tox)   0   1   3   4   6   9  11
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.881    0.258
#> Early stop for efficacy    0.067    0.675
#> Average sample size       19.365   27.450
#> Reject null                0.099    0.726
#plot(test_joint)
```
