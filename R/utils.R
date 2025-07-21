#' Compute Probability Cutoffs for Futility and efficacy Stopping
#'
#' This function computes the probability cutoffs for futility stopping using two methods:
#' power and O'Brien-Fleming type function.
#'
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#'  - 10 is the number of patients enrolled at interim 1,
#'  - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'  - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'  - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param lambda A vector of values for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A vector of values for parameter `gamma` of the cut-off probability for futility 
#' @param eta A vector of values for parameter `eta` of the cut-off probability for efficacy 
#' @param method type of function to be used for the cut off probability for superiority. The default is "power" type function. method=OF 
#' is an alternative for "O'Brien-Fleming" 
#'
#' @return A list containing two elements:
#' \item{cf_values}{A matrix of cutoff values for futility stopping.}
#' \item{cs_values}{A matrix of cutoff values for efficacy.}
#'
#' @importFrom stats qnorm pnorm
#' @keywords internal
#' @export
get_cf_cs_values <- function(n, lambda = NULL, gamma = NULL, eta = NULL, method = "power") {
  nsum <- sum(n)
  # Dimension: lambda * n * gamma
  cf_values <- lambda %o% outer(cumsum(n) / nsum, gamma, '^')
  if(method == "OF") {
    # Dimension: lambda * n
    cs_values <- 2 * pnorm(qnorm((1 + lambda) / 2) %o% sqrt(cumsum(n) / nsum) ^ (-1)) - 1
  } else if(method == 'power') {
    # Dimension: lambda * n * eta
    cs_values <- 1 - (1 - lambda) %o% outer(cumsum(n) / nsum, eta, '^')
  }
  return(list(cf_values = cf_values, cs_values = cs_values))
}

