#' Boundary values for Binary Endpoint
#' @param H0 Response rate under the null
#' @param a1 alpha values for the beta prior (i.e. usually set to the null response rate)
#' @param b1 beta values for the beta prior (i.e, usually set to 1-null response rate)
#' @param nIA number of interim analysis
#' @param n sample size at each interim look
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "OBrien-Fleming" function. method=power is an alternative
#' @param seed for reproducibility
#' @importFrom purrr map_dbl
#' @importFrom dplyr lead lag
#' @importFrom stats pbeta rbinom
#' @importFrom magrittr %>%
#' @export

boundary_binary <- function(H0, a1, b1, nIA, n, lambda, gamma, eta=NULL,  method = NULL,seed = 123) {
  set.seed(seed)

  if (length(n) != nIA){
    stop("sample size vector not equal to number of interm analysis")
  }
  nsum <- sum(n)

  # check here if lambda, gamma and eta values are given and are single number
  # if not call the optimal par function
  # if lambda, gamma and eta values are vector call optimal var get the optimal values

  cf_values <- sapply(seq_along(n), function(i) {
    lambda * (sum(n[1:i]) / nsum)^gamma
  })

  if (is.null(eta)|is.null(method)) {
    method <- "OBrien-Fleming"
  }
  

  if(method == "power"){
    cs_values <- sapply(seq_along(n), function(i) {
      1- (1-lambda)  * (sum(n[1:i]) / nsum)^eta
    })
  }
  else{#OBrien-Fleming method
    cs_values <- sapply(seq_along(n), function(i) {
      2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(sum(n[1:i]) / nsum)) - 1
    })
  }

  calculate_postp <- function(x, H0, a1, b1, n) {
    1 - pbeta(H0, a1 + x, b1 + n - x)
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
      x_range <- (cnf_values[i-1] + 1):(sum(n[1:i]))
    }

    boundaries <- find_boundaries(x_range, H0, a1, b1, sum(n[1:i]), cf_values[i], cs_values[i])
    cnf_values[i] <- boundaries$cnf
    cns_values[i] <- boundaries$cns
  }

  return(data.frame(cnf = cnf_values, cns = cns_values))
}
