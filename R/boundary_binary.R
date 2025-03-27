#' Boundary values for Binary Endpoint
#' @param H0 Response rate under the null
#' @param a1 alpha values for the beta prior (i.e. usually set to the null response rate)
#' @param b1 beta values for the beta prior (i.e, usually set to 1-null response rate)
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
#' @param seed for reproducibility
#' @importFrom purrr map_dbl
#' @importFrom dplyr lead lag
#' @importFrom stats pbeta rbinom
#' @importFrom magrittr %>%
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{cnf: }{Futility boundry}
#'  \item{cns: }{superiority boundry}
#'   } 


boundary_binary <- function(H0, a1, b1, n, lambda, gamma, eta=NULL,  method = "power",seed = 123) {
  set.seed(seed)

  #number of interim analysis
  nIA <- length(n)
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

  calculate_postp <- function(x, H0, a1, b1, n) {
    #1 - pbeta(H0, a1 + x, b1 + n - x)
    ifelse(x > n | x < 0, NA, 1 - pbeta(H0, a1 + x, b1 + n - x))
    
  }

  find_boundaries <- function(x_range, H0, a1, b1, n, cf, cs) {
    postp <- purrr::map_dbl(x_range, ~calculate_postp(.x, H0, a1, b1, n))
    next_values <- purrr::map_dbl(x_range + 1, ~calculate_postp(.x, H0, a1, b1, n))
    previous_values <- purrr::map_dbl(x_range - 1, ~calculate_postp(.x, H0, a1, b1, n))

    cnf <- x_range[which(postp < cf & dplyr::lead(postp, default = NA) >= cf)[1]]
    cns <- x_range[which(postp >= cs & dplyr::lag(postp, default = NA) < cs)[1]]

    list(cnf = cnf, cns = cns)
  }

  cnf_values <- numeric(length(n))
  cns_values <- numeric(length(n))

  for (i in seq_along(n)) {
    if (i == 1) {
      x_range <- 0:n[i]
    } else {
      #x_range <- (cnf_values[i-1] + 1):(sum(n[1:i]))
      if (is.na(cnf_values[i-1])) {
        start_value <- -1
      } else {
        start_value <- cnf_values[i-1]
      }
      
      x_range <- (start_value + 1):(sum(n[1:i]))
    
    }

    boundaries <- find_boundaries(x_range, H0, a1, b1, sum(n[1:i]), cf_values[i], cs_values[i])
    cnf_values[i] <- boundaries$cnf
    cns_values[i] <- boundaries$cns
  }

  return(data.frame(cnf = cnf_values, cns = cns_values))
}
