#' Boundary values for Nested Endpoint
#' @param H0 Response rate under the null , specified in the following order:
#' - `H0[1]`: Complete Remission (CR) rate,
#' - `H0[2]`: Partial Remission (PR) rate,
#' - `H0[3]`: No Complete Remission or Partial Remission rate, calculated as `1 - (CR + PR)`.
#' @param a alpha values for the beta prior (i.e. usually set to the null response rate)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "power" type function. method=OF is an alternative for "O'Brien-Fleming" 
#' @param seed seed number
#' @importFrom dplyr if_else mutate filter
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{cn11f_max: }{Futility boundry for CR}
#'  \item{cn11s_min: }{superiority boundry for CR}
#'   \item{cn11s_min: }{superiority boundry for CR} 
#'   \item{cn12s_min: }{superiority boundry for CR/PR} 
#'   } 

boundary_nested <- function(H0, a, n, lambda, gamma, eta=NULL,  method = "power",seed = 123) {
  
  
  calculate_postp <- function(H0, beta_a, beta_b) {
    1 - pbeta(H0, beta_a, beta_b)
  }

  # Initializing parameters
  b1 <- c(1, 0, 0)
  b2 <- c(1, 1, 0)
  phi1 <- sum(H0 * b1)
  phi2 <- sum(H0 * b2)

  # number of interim analysis
  nIA<- length(n)
  nsum <- sum(n)


  cf_values <- sapply(seq_along(n), function(i) {
    lambda * (sum(n[1:i]) / nsum)^gamma
  })

  
  if (is.null(eta)& method=="power") {
    stop("eta value should be specified for method-power")
  }
  
  if(method == "OF"){
    cs_values <- sapply(seq_along(n), function(i) {
      2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(sum(n[1:i]) / nsum)) - 1
    })
  }
  else{#power
    cs_values <- sapply(seq_along(n), function(i) {
      1- (1-lambda)  * (sum(n[1:i]) / nsum)^eta
    })
  }
  

  # Function to find boundaries
  find_boundaries <- function(n, H0, a, cf, cs) {
    data <- expand.grid(Y11 = 0:n, Y12 = 0:n) %>%
      dplyr::filter(Y11 + Y12 <= n) %>%
      dplyr::mutate(
        Y13 := n - Y11 - Y12,
        beta_a_11 := b1[1] * (a[1] + Y11) + b1[2] * (a[2] + Y12) + b1[3] * (a[3] + Y13),
        beta_b_11 := (1 - b1[1]) * (a[1] + Y11) + (1 - b1[2]) * (a[2] + Y12) + (1 - b1[3]) * (a[3] + Y13),
        beta_a_12 := b2[1] * (a[1] + Y11) + b2[2] * (a[2] + Y12) + b2[3] * (a[3] + Y13),
        beta_b_12 := (1 - b2[1]) * (a[1] + Y11) + (1 - b2[2]) * (a[2] + Y12) + (1 - b2[3]) * (a[3] + Y13),
        postp11 := calculate_postp(phi1, beta_a_11, beta_b_11),
        postp12 := calculate_postp(phi2, beta_a_12, beta_b_12),
        cn11f := dplyr::if_else(postp11 < cf & postp12 < cf, Y11, NA_integer_),
        cn12f := dplyr::if_else(postp11 < cf & postp12 < cf, Y11 + Y12, NA_integer_),
        cn11s := dplyr::if_else(postp11 >= cs , Y11, NA_integer_),
        cn12s := dplyr::if_else(postp12 >= cs, Y11 + Y12, NA_integer_)
      )

    data.frame(
      cn11f_max = max(data$cn11f, na.rm = TRUE),
      cn12f_max = max(data$cn12f, na.rm = TRUE),
      cn11s_min = min(data$cn11s, na.rm = TRUE),
      cn12s_min = min(data$cn12s, na.rm = TRUE)
    )
  }


  boundaries <- do.call(rbind, lapply(seq_along(n), function(i) {
    cumulative_n <- sum(n[1:i])
    find_boundaries(cumulative_n, H0, a, cf_values[i], cs_values[i])
  }))

  return(boundaries)

}
