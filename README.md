
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
BOP2FE_binary(H0 =0.2, H1=0.4, n = c(10,5,5,5,5,5,5), lambda = 0.909, gamma=1, method = "OF", nsim = 10000, seed = 1234)
#> $boundary
#> # A tibble: 2 × 8
#>   `Interim analysis`   `1`   `2`   `3`   `4`   `5`   `6`   `7`  
#>   <chr>                <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 Futility boundary    1     2     3     5     7     9     11   
#> 2 Superiority boundary 7     8     9     10    11    11    12   
#> 
#> $Oc
#> # A tibble: 4 × 3
#>   Statistic                          Under_H0 Under_H1
#>   <chr>                                 <dbl>    <dbl>
#> 1 Early stopping for Futility (%)       87.7      10.3
#> 2 Early stopping for Superiority (%)     7.22     85.9
#> 3 Average sample size                   20.5      23.6
#> 4 Null rejection (%)                     8.58     88.5
#> 
#> $plot
```

<img src="man/figures/README-binary-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "BOP2FE"

### Nested endpoint

``` r
BOP2FE_nested(H0 = c(0.15,0.15,0.70), H1 = c(0.25,0.25,0.50), n=c(10,5,5,5,5,5,5), lambda = 0.95, gamma=1, method = "OF", seed = 123)
#> $boundary
#> # A tibble: 4 × 8
#>   `Interim analysis`         `1`   `2`   `3`   `4`   `5`   `6`   `7`  
#>   <chr>                      <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 Futility boundary CR       0     1     3     4     5     7     10   
#> 2 Futility boundary CR/PR    2     3     6     8     10    13    17   
#> 3 Superiority boundary CR    8     8     8     9     10    10    11   
#> 4 Superiority boundary CR/PR 9     11    12    13    15    16    18   
#> 
#> $Oc
#> # A tibble: 4 × 3
#>   Statistic                          Under_H0 Under_H1
#>   <chr>                                 <dbl>    <dbl>
#> 1 Early stopping for Futility (%)       81.8      9.95
#> 2 Early stopping for Superiority (%)     6.06    78.0 
#> 3 Average sample size                   24.3     26.9 
#> 4 Null rejection (%)                     7.15    82.2 
#> 
#> $plot
```

<img src="man/figures/README-nested-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "BOP2FE"

### Co-primary endpoint

``` r
BOP2FE_coprimary(H0 = c(0.05,0.05,0.15,0.75), H1 =c(0.15,0.15,0.20,0.50) , n=c(10,5,5,5,5,5,5), lambda = 0.95, gamma=1, method = "OF", seed = 123)
#> $boundary
#> # A tibble: 4 × 8
#>   `Interim analysis`        `1`   `2`   `3`   `4`   `5`   `6`   `7`  
#>   <chr>                     <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 Futility boundary ORR     0     1     2     3     4     5     7    
#> 2 Futility boundary PFS6    1     2     4     5     7     9     12   
#> 3 Superiority boundary ORR  7     7     7     7     7     8     8    
#> 4 Superiority boundary PFS6 8     9     10    11    11    12    13   
#> 
#> $Oc
#> # A tibble: 4 × 3
#>   Statistic                          Under_H0 Under_H1
#>   <chr>                                 <dbl>    <dbl>
#> 1 Early stopping for Futility (%)       80.8      7.79
#> 2 Early stopping for Superiority (%)     6.23    82.1 
#> 3 Average sample size                   24.1     25.9 
#> 4 Null rejection (%)                     8.39    88.5 
#> 
#> $plot
```

<img src="man/figures/README-coprimary-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "BOP2FE"

### Joint endpoint

``` r
BOP2FE_jointefftox(H0 = c(0.15, 0.30, 0.15, 0.40), H1=c(0.18, 0.42, 0.02, 0.38), n=c(10,5,5,5,5,5,5), lambda = 0.7, gamma=1, method = "OF", seed = 123)
#> $boundary
#> # A tibble: 4 × 8
#>   `Interim analysis`            `1`   `2`   `3`   `4`   `5`   `6`   `7`  
#>   <chr>                         <chr> <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 Futility boundary response    3     5     8     10    13    16    19   
#> 2 Futility boundary toxicity    5     6     7     8     9     10    11   
#> 3 Superiority boundary response 8     10    12    14    16    18    20   
#> 4 Superiority boundary toxicity 0     2     4     5     7     8     10   
#> 
#> $Oc
#> # A tibble: 4 × 3
#>   Statistic                          Under_H0 Under_H1
#>   <chr>                                 <dbl>    <dbl>
#> 1 Early stopping for Futility (%)       90.2      30.0
#> 2 Early stopping for Superiority (%)     7.44     65.8
#> 3 Average sample size                   17.9      22.9
#> 4 Null rejection (%)                     8.28     68.8
#> 
#> $plot
```

<img src="man/figures/README-joint-1.png" width="100%" />

    #> 
    #> attr(,"class")
    #> [1] "BOP2FE"

### Optimizing Design Parameters

The operating characteristics of the BOP2-FE design with futility and
efficacy stopping rules depend on the specification of the probability
cutoffs for futility and efficacy. To maximize statistical power while
controlling the type I error rate at a certain pre-specified level, the
tuning parameters, (lambda, gamma, eta) when using the power function
for the probability cutoffs for efficacy, and (lambda, gamma) when using
the O’Brien-Fleming type function for the probability cutoffs for
efficacy, should be identified through simulation.
*search_optimal_pars_binary()*, *search_optimal_pars_nested()*,
*search_optimal_pars_coprimary()*, and
*search_optimal_pars_jointefftox()* can facilitate one to do a grid
search. An example for a nested outcome is given below.

``` r
Sim_res <- search_optimal_pars_nested(H0 = c(0.15, 0.15, 0.70), H1=c(0.25, 0.25,0.50), n=c(10,5,5,5,5,5,5), nsim=1000, t1e=0.1, method="power",
                              lambda1=0.85, lambda2=1, grid1=4, gamma1=0.9, gamma2=1, grid2=5, eta1=0.9, eta2=1, grid3=5)
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |===                                                                   |   4%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |=================                                                     |  24%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |===============================                                       |  44%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |=============================================                         |  64%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |===========================================================           |  84%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%

head(Sim_res)
#>   fut_boundary_CR_1 fut_boundary_CR_2 fut_boundary_CR_3 fut_boundary_CR_4
#> 1                 0                 1                 3                 4
#> 2                 0                 1                 3                 4
#> 3                 0                 1                 3                 4
#> 4                 0                 1                 3                 4
#> 5                 0                 1                 3                 4
#> 6                 0                 1                 3                 4
#>   fut_boundary_CR_5 fut_boundary_CR_6 fut_boundary_CR_7 fut_boundary_CR.PR_1
#> 1                 5                 7                10                    2
#> 2                 5                 7                10                    2
#> 3                 5                 7                10                    2
#> 4                 5                 7                10                    2
#> 5                 5                 7                10                    2
#> 6                 5                 7                10                    2
#>   fut_boundary_CR.PR_2 fut_boundary_CR.PR_3 fut_boundary_CR.PR_4
#> 1                    4                    6                    8
#> 2                    4                    6                    8
#> 3                    4                    6                    8
#> 4                    4                    6                    8
#> 5                    4                    6                    8
#> 6                    4                    6                    8
#>   fut_boundary_CR.PR_5 fut_boundary_CR.PR_6 fut_boundary_CR.PR_7
#> 1                   10                   13                   17
#> 2                   10                   13                   17
#> 3                   10                   13                   17
#> 4                   10                   13                   17
#> 5                   10                   13                   17
#> 6                   10                   13                   17
#>   Sup_boundary_CR_1 Sup_boundary_CR_2 Sup_boundary_CR_3 Sup_boundary_CR_4
#> 1                 5                 6                 7                 8
#> 2                 5                 6                 7                 8
#> 3                 5                 6                 7                 8
#> 4                 5                 6                 7                 8
#> 5                 5                 6                 7                 8
#> 6                 5                 6                 7                 8
#>   Sup_boundary_CR_5 Sup_boundary_CR_6 Sup_boundary_CR_7 Sup_boundary_CR.PR_1
#> 1                 9                10                11                    7
#> 2                 9                10                11                    7
#> 3                 9                10                11                    7
#> 4                 9                10                11                    7
#> 5                 9                10                11                    7
#> 6                 9                10                11                    7
#>   Sup_boundary_CR.PR_2 Sup_boundary_CR.PR_3 Sup_boundary_CR.PR_4
#> 1                    9                   11                   13
#> 2                    9                   11                   13
#> 3                    9                   11                   13
#> 4                    9                   11                   13
#> 5                    9                   11                   13
#> 6                    9                   11                   13
#>   Sup_boundary_CR.PR_5 Sup_boundary_CR.PR_6 Sup_boundary_CR.PR_7
#> 1                   15                   16                   18
#> 2                   15                   16                   18
#> 3                   15                   16                   18
#> 4                   15                   16                   18
#> 5                   15                   16                   18
#> 6                   15                   16                   18
#>   earlystopfuti_mean_h0 earlystopsupe_mean_h0 ss_mean_h0 rejectnull_mean_h0
#> 1                 0.811                 0.063     23.515              0.076
#> 2                 0.811                 0.063     23.515              0.076
#> 3                 0.811                 0.063     23.515              0.076
#> 4                 0.811                 0.063     23.515              0.076
#> 5                 0.811                 0.063     23.515              0.076
#> 6                 0.811                 0.063     23.515              0.076
#>   earlystopfuti_mean_h1 earlystopsupe_mean_h1 ss_mean_h1 rejectnull_mean_h1
#> 1                 0.095                   0.8     22.645              0.834
#> 2                 0.095                   0.8     22.645              0.834
#> 3                 0.095                   0.8     22.645              0.834
#> 4                 0.095                   0.8     22.645              0.834
#> 5                 0.095                   0.8     22.645              0.834
#> 6                 0.095                   0.8     22.645              0.834
#>   lambda gamma  eta
#> 1 0.9625     1 0.90
#> 2 0.9625     1 0.92
#> 3 0.9625     1 0.94
#> 4 0.9625     1 0.96
#> 5 0.9625     1 0.98
#> 6 0.9625     1 1.00
```

However, these functions only use a single core and the process might be
slow. If one has access to multiple cores, we can directly call the
compute_power function using the parallel backend back end as follows

``` r
# Load necessary libraries
library(foreach)
library(doParallel)
library(BOP2FE)

# Define the parameter ranges
lambda_range <- seq(0, 1, by = 0.01)
gamma_range <- seq(0, 1, by = 0.01)
eta_range <- seq(0, 1, by = 0.01)

# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than the total number of cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run the function in parallel
results <- foreach(lambda = lambda_range, .combine = 'rbind') %:%
  foreach(gamma = gamma_range, .combine = 'rbind') %:%
  foreach(eta = eta_range, .combine = 'rbind') %dopar% {
    BOP2FE::get_boundary_oc_nested(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                            lambda = l, gamma = g, eta = e, method = method, seed=seed)
          
  }

# Stop the parallel backend
stopCluster(cl)

# View the results
results%>%
  filter(reject_mean <= 0.1)%>%
  arrange(desc(power_mean))%>%
  head()
```
