
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
#> 1    0.9     1   2
#> 
#> $boundary
#>                   IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary   1   2   3   5   7   9  11
#> Efficacy boundary   6   7   9   9  10  11  12
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.873    0.119
#> Early stop for efficacy    0.087    0.850
#> Average sample size       20.445   20.865
#> Reject null                0.096    0.868
#plot(test_binary)
```

``` r
plot(test_binary)
```

<img src="man/figures/README-binary2-1.png" width="100%" />

### Nested endpoint

``` r
test_nested <- BOP2FE_nested(
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
#> 1    0.1   0.9   0
#> 
#> $boundary
#>                           IA1 IA2 IA3 IA4 IA5 IA6 IA7 IA8 IA9 IA10 IA11 IA12
#> Futility boundary (CR)      5   9  12  16  20  24  28   5   9   14   18   22
#> Futility boundary (CR/PR)   5   9  14  18  22  26  30   5   9   14   18   22
#> Efficacy boundary (CR)      7  10  14  18  22  25  29   7  11   15   19   23
#> Efficacy boundary (CR/PR)   7  11  15  19  23  27  31   7  11   15   19   23
#>                           IA13 IA14
#> Futility boundary (CR)      26   30
#> Futility boundary (CR/PR)   26   30
#> Efficacy boundary (CR)      27   31
#> Efficacy boundary (CR/PR)   27   31
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.946    0.604
#> Early stop for efficacy    0.054    0.396
#> Average sample size       10.640   11.725
#> Reject null                0.054    0.396
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
#> 1    0.9   0.3 2.6
#> 
#> $boundary
#>                          IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary (OR)     1   2   3   4   4   5   6
#> Futility boundary (PFS6)   2   3   5   6   8   9  11
#> Efficacy boundary (OR)     5   6   6   7   7   7   7
#> Efficacy boundary (PFS6)   7   8   9  10  11  11  12
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.882    0.097
#> Early stop for efficacy    0.084    0.896
#> Average sample size       16.620   18.940
#> Reject null                0.099    0.902
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
#> 1    0.6   0.9 1.2
#> 
#> $boundary
#>                         IA1 IA2 IA3 IA4 IA5 IA6 IA7
#> Futility boundary (OR)    3   5   8  10  13  16  18
#> Futility boundary (Tox)   5   6   8   9  10  11  12
#> Efficacy boundary (OR)    7  10  12  14  16  17  19
#> Efficacy boundary (Tox)   1   2   4   5   7   9  11
#> 
#> $oc
#>                         Under H0 Under H1
#> Early stop for futility    0.888    0.242
#> Early stop for efficacy    0.093    0.733
#> Average sample size       18.150   22.880
#> Reject null                0.098    0.751
#plot(test_joint)
```
