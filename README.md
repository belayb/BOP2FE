
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BOP2FE

<!-- badges: start -->
<!-- badges: end -->

The primary purpose of an oncology single-arm trial is to evaluate the
effectiveness of anticancer agents and make a go / no-go decision while
maintaining patient safety. This R-package implements a flexible
Bayesian optimal phase II design with futility and efficacy stopping
boundaries for single-arm clinical trials, named the BOP2-FE design,
proposed by Xu et al (). The proposed BOP2-FE design allows for early
stopping of efficacy when the observed antitumor effect is sufficiently
higher than the null hypothesis value in the interim looks and retains
the benefits of the original BOP2 design, such as explicitly controlling
the type I error rate while maximizing power, accommodating different
types of endpoint, flexible number of interim looks, and stopping
boundaries calculated before the start of the trial.

## Installation

You can install the development version of BOP2FE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("belayb/BOP2FE")
```

## Example

This is a basic example which shows you how to use BOPFE:

``` r
library(BOP2FE)
## basic example code
```

### Binary endpoint

``` r
BOP2FE_binary(H0 =0.2, n = c(10,5,5,5,5,5,5), lambda = 0.909, gamma=1, eta=NULL, nsim = 10000, seed = 1234)
```

<img src="man/figures/README-binary-1.png" width="100%" />

### Nested endpoint

``` r
BOP2FE_nested(CR0 = 0.15, CRPR0 = 0.30, n=c(10,5,5,5,5,5,5), lambda = 0.95, gamma=1, eta=NULL,  method = "obrian",seed = 123)
```

<img src="man/figures/README-nested-1.png" width="100%" />

### Co-primary endpoint

``` r
BOP2FE_coprimary(H0 = c(0.05,0.05,0.15,0.75), n=c(10,5,5,5,5,5,5), lambda = 0.95, gamma=1, eta=NULL,  method = NULL,seed = 123)
```

<img src="man/figures/README-coprimary-1.png" width="100%" />

### Joint endpoint

``` r
BOP2FE_jointefftox(H0 = c(0.15, 0.30, 0.15, 0.40), n=c(10,5,5,5,5,5,5), lambda = 0.7, gamma=1, eta=NULL,  method = "O’Brien-Fleming",seed = 123)
```

<img src="man/figures/README-joint-1.png" width="100%" />
