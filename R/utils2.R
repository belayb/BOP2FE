#' Computes bothe boundry and corresponding operating charactersitcs for binary data 
#' @param H0 A numeric value for the response rate under the null hypothesis
#' @param H1 A numeric value for the response rate under the alternative hypothesis
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result 
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility 
#' @export    
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{fut_boundary_i: }{Futility boundary at the ith analysis}
#' \item{sup_boundary_i: }{Superiority boundary at the ith analysis}
#' \item{earlystopfuti_mean_h0: }{Average number of early stopping due to futility under the null hypothesis}
#'  \item{earlystopsupe_mean_h0: }{Average number of early stopping for futility due to efficacy under the null hypothesis}
#'   \item{ss_mean_h0: }{Average sample size under the null hypothesis} 
#'   \item{rejectnull_mean_h0: }{Average number of hypothesis rejection at the final analysis under the null hypothesis} 
#'   \item{earlystopfuti_mean_h1: }{Average number of early stopping due to futility under the alternative hypothesis} 
#'   \item{earlystopsupe_mean_h1: }{Average number of early early stopping due to efficacy under the alternative hypothesis} 
#'   \item{ss_mean_h1: }{Average sample size under the alternative hypothesis} 
#'   \item{rejectnull_mean_h1: }{Average number of hypothesis rejection at the final analysis under the alternative hypothesis} 
#'   } 
#'   


get_boundary_oc_binary <- function(H0, H1, n, nsim, lambda, gamma, eta=NULL, 
                                   method = "power", seed = NULL) {
  nIA <- length(n)
  a1 <- H0
  b1 <- 1-a1
  # Generate boundary test results
  boundary_tab <- boundary_binary(H0 = H0, a1 = a1, b1 = b1, n = n, 
                                   lambda = lambda, gamma = gamma, 
                                   eta = eta, method = method, 
                                   seed = seed)
  
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- Oc_binary(p = H0, n = n, nsim = nsim, 
                       fb = boundary_tab$cnf, 
                       sb = boundary_tab$cns)
  
  alt_oc <- Oc_binary(p = H1, n = n, nsim = nsim, 
                      fb = boundary_tab$cnf, 
                      sb = boundary_tab$cns)
  
  # Combine results into a single vector
  all_res <- c(boundary_tab[, 1], boundary_tab[, 2], 
               null_oc[1:4], alt_oc[1:4], 
               lambda, gamma, eta)
  
  # names for data frame
  names(all_res) <- c(paste0("fut_boundary", 1:nIA), 
                      paste0("sup_boundary", 1:nIA), 
                      "earlystopfuti_mean_h0", 
                      "earlystopsupe_mean_h0", 
                      "ss_mean_h0", 
                      "rejectnull_mean_h0", 
                      "earlystopfuti_mean_h1", 
                      "earlystopsupe_mean_h1", 
                      "ss_mean_h1", 
                      "rejectnull_mean_h1", 
                      "lambda", "gamma", "eta")
  
  # Convert the results to a data frame
  all_res_df <- data.frame(all_res)
  
  return(all_res_df)
}

################################################################################


#' Computes both the boundry and corresponding operating charactersitcs for nested efficacy endpoint  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Complete Remission (CR) rate,
#' - `H0[2]`: Partial Remission (PR) rate,
#' - `H0[3]`: No Complete Remission or Partial Remission rate, calculated as `1 - (CR + PR)`.
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Complete Remission (CR) rate,
#' - `H1[2]`: Partial Remission (PR) rate,
#' - `H1[3]`: No Complete Remission or Partial Remission rate, calculated as `1 - (CR + PR)`.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility  
#' @export           
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{fut_boundary_CR_i: }{Futility boundary for CR at the ith analysis}
#' \item{fut_boundary_CR.PR_i: }{Futility boundary for CR/PR at the ith analysis}
#' \item{sup_boundary_CR_i: }{Superiority boundary for CR at the ith analysis}
#' \item{sup_boundary_CR.PR_i: }{Superiority boundary for CR/PRat the ith analysis}
#' \item{earlystopfuti_mean_h0: }{Average number of early stopping due to futility under the null hypothesis}
#'  \item{earlystopsupe_mean_h0: }{Average number of early stopping for futility due to efficacy under the null hypothesis}
#'   \item{ss_mean_h0: }{Average sample size under the null hypothesis} 
#'   \item{rejectnull_mean_h0: }{Average number of hypothesis rejection at the final analysis under the null hypothesis} 
#'   \item{earlystopfuti_mean_h1: }{Average number of early stopping due to futility under the alternative hypothesis} 
#'   \item{earlystopsupe_mean_h1: }{Average number of early early stopping due to efficacy under the alternative hypothesis} 
#'   \item{ss_mean_h1: }{Average sample size under the alternative hypothesis} 
#'   \item{rejectnull_mean_h1: }{Average number of hypothesis rejection at the final analysis under the alternative hypothesis} 
#'   } 
#'   

get_boundary_oc_nested <- function(H0, H1, n, nsim, lambda, gamma, eta=NULL, 
                                   method = "power", seed = NULL) {
  a <- H0
  nIA <- length(n)
  # Generate boundary test results
  boundary_tab <- boundary_nested(H0 = H0, a = a, n = n, 
                                   lambda = lambda, gamma = gamma, 
                                   eta = eta, method = method, 
                                   seed = seed)
  
  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))
  
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- Oc_nested(p1 = H0[1],p2=H0[2], p3=H0[3], n = n, nsim = nsim, 
                       fb = fb, 
                       sb = sb)
  
  alt_oc <- Oc_nested(p1 = H1[1],p2=H1[2], p3=H1[3], n = n, nsim = nsim, 
                      fb = fb, 
                      sb = sb)
  
  # Combine results into a single vector
  all_res <- c(boundary_tab[, 1], boundary_tab[, 2], 
               boundary_tab[, 3], boundary_tab[, 4],
               null_oc[1:4], alt_oc[1:4], lambda, gamma, eta)
  
  # names for data frame
  names(all_res) <- c(paste0("fut_boundary_CR_", 1:nIA), 
                      paste0("fut_boundary_CR/PR_", 1:nIA), 
                      paste0("Sup_boundary_CR_", 1:nIA), 
                      paste0("Sup_boundary_CR/PR_", 1:nIA), 
                      "earlystopfuti_mean_h0", 
                      "earlystopsupe_mean_h0", 
                      "ss_mean_h0", 
                      "rejectnull_mean_h0", 
                      "earlystopfuti_mean_h1", 
                      "earlystopsupe_mean_h1", 
                      "ss_mean_h1", 
                      "rejectnull_mean_h1", 
                      "lambda", "gamma", "eta")
  
  # Convert the results to a data frame
  all_res_df <- data.frame(all_res)
  
  return(all_res_df)
}

################################################################################

#' Computes both the boundry and corresponding operating charactersitcs for co-primary efficacy endpoints  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Response - PFS6,
#' - `H0[2]`: Response - no PFS6,
#' - `H0[3]`: No Response - PFS6,
#' - `H0[4]`: No Response - no PFS6.
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Response - PFS6,
#' - `H1[2]`: Response - no PFS6,
#' - `H1[3]`: No Response - PFS6,
#' - `H1[4]`: No Response - no PFS6.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility 
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{fut_boundary_ORR_i: }{Futility boundary for ORR at the ith analysis}
#' \item{fut_boundary_PFS6_i: }{Futility boundary for PFS6 at the ith analysis}
#' \item{sup_boundary_ORR_i: }{Superiority boundary for ORR at the ith analysis}
#' \item{sup_boundary_PFS6_i: }{Superiority boundary for PFS6 at the ith analysis}
#' \item{earlystopfuti_mean_h0: }{Average number of early stopping due to futility under the null hypothesis}
#'  \item{earlystopsupe_mean_h0: }{Average number of early stopping for futility due to efficacy under the null hypothesis}
#'   \item{ss_mean_h0: }{Average sample size under the null hypothesis} 
#'   \item{rejectnull_mean_h0: }{Average number of hypothesis rejection at the final analysis under the null hypothesis} 
#'   \item{earlystopfuti_mean_h1: }{Average number of early stopping due to futility under the alternative hypothesis} 
#'   \item{earlystopsupe_mean_h1: }{Average number of early early stopping due to efficacy under the alternative hypothesis} 
#'   \item{ss_mean_h1: }{Average sample size under the alternative hypothesis} 
#'   \item{rejectnull_mean_h1: }{Average number of hypothesis rejection at the final analysis under the alternative hypothesis} 
#'   } 
#'   

get_boundary_oc_coprimary <- function(H0, H1, n, nsim, lambda, gamma, eta=NULL, 
                                   method = "power", seed = NULL) {
  a <- H0
  nIA <- length(n)
  # Generate boundary test results
  boundary_tab <- boundary_coprimary(H0 = H0, a = a, n = n, 
                                  lambda = lambda, gamma = gamma, 
                                  eta = eta, method = method, 
                                  seed = seed)
  
  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))
  
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- Oc_coprimary(p1 = H0[1],p2=H0[2], p3=H0[3], p4=H0[4],  n = n, nsim = nsim, 
                       fb = fb, 
                       sb = sb,
                       seed = seed)
  
  alt_oc <- Oc_coprimary(p1 = H1[1],p2=H1[2], p3=H1[3], p4=H1[4],  n = n, nsim = nsim, 
                      fb = fb, 
                      sb = sb,
                      seed = seed)
  
  # Combine results into a single vector
  all_res <- c(boundary_tab[, 1], boundary_tab[, 2], 
               boundary_tab[, 3], boundary_tab[, 4],
               null_oc[1:4], alt_oc[1:4], lambda, gamma, eta)
  
  # names for data frame

  names(all_res) <- c(paste0("fut_boundary_ORR_", 1:nIA), 
                      paste0("fut_boundary_PFS6_", 1:nIA), 
                      paste0("Sup_boundary_ORR_", 1:nIA), 
                      paste0("Sup_boundary_PFS6_", 1:nIA), 
                      "earlystopfuti_mean_h0", 
                      "earlystopsupe_mean_h0", 
                      "ss_mean_h0", 
                      "rejectnull_mean_h0", 
                      "earlystopfuti_mean_h1", 
                      "earlystopsupe_mean_h1", 
                      "ss_mean_h1", 
                      "rejectnull_mean_h1", 
                      "lambda", "gamma", "eta")
  
  # Convert the results to a data frame
  all_res_df <- data.frame(all_res)
  
  return(all_res_df)
}

################################################################################

#' Computes both the boundry and corresponding operating charactersitcs for joint efficacy and toxicity endpoints  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Response - toxicity,
#' - `H0[2]`: Response - no toxicity,
#' - `H0[3]`: No Response - toxicity,
#' - `H0[4]`: No Response - no toxicity
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Response - toxicity,
#' - `H1[2]`: Response - no toxicity,
#' - `H1[3]`: No Response - toxicity,
#' - `H1[4]`: No Response - no toxicity.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility 
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{fut_boundary_resp_i: }{Futility boundary for response at the ith analysis}
#' \item{fut_boundary_tox_i: }{Futility boundary for toxicity at the ith analysis}
#' \item{sup_boundary_resp_i: }{Superiority boundary for response at the ith analysis}
#' \item{sup_boundary_tox_i: }{Superiority boundary for toxicity at the ith analysis}
#' \item{earlystopfuti_mean_h0: }{Average number of early stopping due to futility under the null hypothesis}
#'  \item{earlystopsupe_mean_h0: }{Average number of early stopping for futility due to efficacy under the null hypothesis}
#'   \item{ss_mean_h0: }{Average sample size under the null hypothesis} 
#'   \item{rejectnull_mean_h0: }{Average number of hypothesis rejection at the final analysis under the null hypothesis} 
#'   \item{earlystopfuti_mean_h1: }{Average number of early stopping due to futility under the alternative hypothesis} 
#'   \item{earlystopsupe_mean_h1: }{Average number of early early stopping due to efficacy under the alternative hypothesis} 
#'   \item{ss_mean_h1: }{Average sample size under the alternative hypothesis} 
#'   \item{rejectnull_mean_h1: }{Average number of hypothesis rejection at the final analysis under the alternative hypothesis} 
#'   } 
#'   
get_boundary_oc_efftox <- function(H0, H1, n, nsim, lambda, gamma, eta=NULL, 
                                      method = "power", seed = NULL) {
  a <- H0
  nIA <- length(n)
  # Generate boundary test results
  boundary_tab <- boundary_jointefftox(H0 = H0, a = a, n = n, 
                                     lambda = lambda, gamma = gamma, 
                                     eta = eta, method = method, 
                                     seed = seed)
  
  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))
  
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- Oc_jointefftox(p1 = H0[1],p2=H0[2], p3=H0[3], p4=H0[4],  n = n, nsim = nsim, 
                          fb = fb, 
                          sb = sb,
                          seed = seed)
  
  alt_oc <- Oc_jointefftox(p1 = H1[1],p2=H1[2], p3=H1[3], p4=H1[4],  n = n, nsim = nsim, 
                         fb = fb, 
                         sb = sb,
                         seed = seed)
  
  # Combine results into a single vector
  all_res <- c(boundary_tab[, 1], boundary_tab[, 2], 
               boundary_tab[, 3], boundary_tab[, 4],
               null_oc[1:4], alt_oc[1:4], lambda, gamma, eta)
  
  # names for data frame
  names(all_res) <- c(paste0("fut_boundary_resp_", 1:nIA), 
                      paste0("fut_boundary_tox_", 1:nIA), 
                      paste0("Sup_boundary_resp_", 1:nIA), 
                      paste0("Sup_boundary_tox_", 1:nIA), 
                      "earlystopfuti_mean_h0", 
                      "earlystopsupe_mean_h0", 
                      "ss_mean_h0", 
                      "rejectnull_mean_h0", 
                      "earlystopfuti_mean_h1", 
                      "earlystopsupe_mean_h1", 
                      "ss_mean_h1", 
                      "rejectnull_mean_h1", 
                      "lambda", "gamma", "eta")
  
  # Convert the results to a data frame
  all_res_df <- data.frame(all_res)
  
  return(all_res_df)
}



#' Search optimal parameters for binary endpoint
#' 
#' `search_optimal_pars_binary()` is a helper function to run a grid search to find the optimal parameters for futility and efficacy stopping cut-off probabilities. The function runs on 
#' a single core and can take a long time if the search space is too wide. If one has access to multiple cores, one can directly call
#'  `get_boundary_oc_binary()` function and run it in parallel.  
#'  
#' @param H0 A numeric value for the response rate under the null hypothesis
#' @param H1 A numeric value for the response rate under the alternative hypothesis
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param t1e Desired Type - I error rate. If specified it will only return results with type I error rate less the specified value 
#' @param lambda1 starting value for `lambda` values to search
#' @param lambda2 ending value for `lambda` values to search
#' @param grid1 number of `lambda` values to consider between lambda1 and lambda2
#' @param gamma1 starting value for `gamma` values to search
#' @param gamma2 ending value for `gamma` values to search
#' @param grid2 number of `gamma` values to consider between gamma1 and gamma2
#' @param eta1 starting value for `eta` values to search
#' @param eta2 ending value for `eta` values to search
#' @param grid3 number of eta values to consider between eta1 and eta2
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
search_optimal_pars_binary <- function(H0, H1, n, nsim, t1e = NULL, method ="power",
                                       lambda1, lambda2, grid1, gamma1, gamma2, 
                                       grid2, eta1, eta2, grid3, seed=NULL) {
  
  
  # Initialize result a data frame
  nIA <- length(n)
  column_names <- c(paste0("fut_boundary", 1:nIA), 
                    paste0("sup_boundary", 1:nIA), 
                    "earlystopfuti_mean_h0", 
                    "earlystopsupe_mean_h0", 
                    "ss_mean_h0", 
                    "rejectnull_mean_h0", 
                    "earlystopfuti_mean_h1", 
                    "earlystopsupe_mean_h1", 
                    "ss_mean_h1", 
                    "rejectnull_mean_h1", 
                    "lambda", "gamma", "eta")
  
  results <- data.frame(matrix(0, nrow = 0, ncol = length(column_names)))
  colnames(results) <- column_names
  
  
  # Calculate total iterations
  total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
  pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration <- 0
  
  # Iterate over grid values
  for (i in 1:(grid1 + 1)) {
    for (j in 1:(grid2 + 1)) {
      if (method == "power") {
        for (et in 1:(grid3 + 1)) {
          # Calculate lambda, gamma, and eta values for this iteration
          l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
          g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
          e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
          
          # Call the compute_power_nested function
          result <-  get_boundary_oc_binary(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                    lambda = l, gamma = g, eta = e, method = method, seed=seed)
          
          # Append the results
          results <- rbind(results, result)
          
          # Update progress bar
          iteration <- iteration + 1
          utils::setTxtProgressBar(pb, iteration)
        }
      } else {
        # Calculate lambda and gamma values for this iteration
        l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
        g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
        
        # Call the compute_power_nested function without eta
        result <-  get_boundary_oc_binary(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                  lambda = l, gamma = g, eta = NULL, method = method, seed = seed)
        # Append the results
        results <- rbind(results, result)
        
        # Update progress bar
        iteration <- iteration + 1
        utils::setTxtProgressBar(pb, iteration)
      }
    }
  }
  close(pb)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- results %>% dplyr::filter(!is.na(rejectnull_mean_h0) & rejectnull_mean_h0  < t1e)
  
  # Sort results by power
  results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(rejectnull_mean_h1 ))
  
  return(results_sorted)
}

#' Search optimal parameters for nested efficacy endpoint
#' 
#' `search_optimal_pars_nested()` is a helper function to run a grid search to find the optimal parameters for futility and efficacy stopping cut-off probabilities. The function runs on 
#' a single core and can take a long time if the search space is too wide. If one has access to multiple cores, one can directly call
#'  `get_boundary_oc_nested()` function and run it in parallel. 
#'  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Complete Remission (CR) rate,
#' - `H0[2]`: Partial Remission (PR) rate,
#' - `H0[3]`: No Complete Remission or Partial Remission rate, calculated as `1 - (CR + PR)`.
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Complete Remission (CR) rate,
#' - `H1[2]`: Partial Remission (PR) rate,
#' - `H1[3]`: No Complete Remission or Partial Remission rate, calculated as `1 - (CR + PR)`.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param t1e Desired Type - I error rate. If specified it will only return results with type I error rate less the specified value 
#' @param lambda1 starting value for `lambda` values to search
#' @param lambda2 ending value for `lambda` values to search
#' @param grid1 number of `lambda` values to consider between lambda1 and lambda2
#' @param gamma1 starting value for `gamma` values to search
#' @param gamma2 ending value for `gamma` values to search
#' @param grid2 number of `gamma` values to consider between gamma1 and gamma2
#' @param eta1 starting value for `eta` values to search
#' @param eta2 ending value for `eta` values to search
#' @param grid3 number of eta values to consider between eta1 and eta2
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
search_optimal_pars_nested <- function(H0, H1, n, nsim, t1e=NULL, method ="power",
                                       lambda1, lambda2, grid1, gamma1, gamma2, 
                                       grid2, eta1, eta2, grid3, seed=NULL) {
  
  # Initialize result a data frame
  nIA <- length(n)
  column_names <-  c(paste0("fut_boundary_CR_", 1:nIA), 
                     paste0("fut_boundary_CR/PR_", 1:nIA), 
                     paste0("Sup_boundary_CR_", 1:nIA), 
                     paste0("Sup_boundary_CR/PR_", 1:nIA), 
                     "earlystopfuti_mean_h0", 
                     "earlystopsupe_mean_h0", 
                     "ss_mean_h0", 
                     "rejectnull_mean_h0", 
                     "earlystopfuti_mean_h1", 
                     "earlystopsupe_mean_h1", 
                     "ss_mean_h1", 
                     "rejectnull_mean_h1", 
                     "lambda", "gamma", "eta")
  
  results <- data.frame(matrix(0, nrow = 0, ncol = length(column_names)))
  colnames(results) <- column_names
  
  
  # Calculate total iterations
  total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
  pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration <- 0
  
  # Iterate over grid values
  for (i in 1:(grid1 + 1)) {
    for (j in 1:(grid2 + 1)) {
      if (method == "power") {
        for (et in 1:(grid3 + 1)) {
          # Calculate lambda, gamma, and eta values for this iteration
          l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
          g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
          e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
          
          # Call the compute_power_nested function
          result <-  get_boundary_oc_nested(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                            lambda = l, gamma = g, eta = e, method = method, seed=seed)
          
          # Append the results
          results <- rbind(results, result)
          
          # Update progress bar
          iteration <- iteration + 1
          utils::setTxtProgressBar(pb, iteration)
        }
      } else {
        # Calculate lambda and gamma values for this iteration
        l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
        g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
        
        # Call the compute_power_nested function without eta
        result <-  get_boundary_oc_nested(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                  lambda = l, gamma = g, eta = NULL, method = method, seed = seed)
        # Append the results
        results <- rbind(results, result)
        
        # Update progress bar
        iteration <- iteration + 1
        utils::setTxtProgressBar(pb, iteration)
      }
    }
  }
  close(pb)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- results %>% dplyr::filter(!is.na(rejectnull_mean_h0) & rejectnull_mean_h0  < t1e)
  
  # Sort results by power
  results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(rejectnull_mean_h1 ))
  
  return(results_sorted)
}


#' Search optimal parameters for co-primary efficacy endpoint
#' 
#' `search_optimal_pars_coprimary()` is a helper function to run a grid search to find the optimal parameters for futility and efficacy stopping cut-off probabilities. The function runs on 
#' a single core and can take a long time if the search space is too wide. If one has access to multiple cores, one can directly call
#'  `get_boundary_oc_coprimary()` function and run it in parallel.   
#'  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Response - PFS6,
#' - `H0[2]`: Response - no PFS6,
#' - `H0[3]`: No Response - PFS6,
#' - `H0[4]`: No Response - no PFS6.
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Response - PFS6,
#' - `H1[2]`: Response - no PFS6,
#' - `H1[3]`: No Response - PFS6,
#' - `H1[4]`: No Response - no PFS6.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param t1e Desired Type - I error rate. If specified it will only return results with type I error rate less the specified value 
#' @param lambda1 starting value for `lambda` values to search
#' @param lambda2 ending value for `lambda` values to search
#' @param grid1 number of `lambda` values to consider between lambda1 and lambda2
#' @param gamma1 starting value for `gamma` values to search
#' @param gamma2 ending value for `gamma` values to search
#' @param grid2 number of `gamma` values to consider between gamma1 and gamma2
#' @param eta1 starting value for `eta` values to search
#' @param eta2 ending value for `eta` values to search
#' @param grid3 number of eta values to consider between eta1 and eta2
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
search_optimal_pars_coprimary <- function(H0, H1, n, nsim, t1e=NULL, method ="power",
                                       lambda1, lambda2, grid1, gamma1, gamma2, 
                                       grid2, eta1, eta2, grid3, seed=NULL) {
  
  # Initialize result a data frame
  nIA <- length(n)
  column_names <-  c(paste0("fut_boundary_ORR_", 1:nIA), 
                     paste0("fut_boundary_PFS6_", 1:nIA), 
                     paste0("Sup_boundary_ORR_", 1:nIA), 
                     paste0("Sup_boundary_PFS6_", 1:nIA), 
                     "earlystopfuti_mean_h0", 
                     "earlystopsupe_mean_h0", 
                     "ss_mean_h0", 
                     "rejectnull_mean_h0", 
                     "earlystopfuti_mean_h1", 
                     "earlystopsupe_mean_h1", 
                     "ss_mean_h1", 
                     "rejectnull_mean_h1", 
                     "lambda", "gamma", "eta")
  
  results <- data.frame(matrix(0, nrow = 0, ncol = length(column_names)))
  colnames(results) <- column_names
  
  
  # Calculate total iterations
  total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
  pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration <- 0
  
  # Iterate over grid values
  for (i in 1:(grid1 + 1)) {
    for (j in 1:(grid2 + 1)) {
      if (method == "power") {
        for (et in 1:(grid3 + 1)) {
          # Calculate lambda, gamma, and eta values for this iteration
          l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
          g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
          e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
          
          # Call the compute_power_nested function
          result <-  get_boundary_oc_coprimary(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                            lambda = l, gamma = g, eta = e, method = method, seed=seed)
          
          # Append the results
          results <- rbind(results, result)
          
          # Update progress bar
          iteration <- iteration + 1
          utils::setTxtProgressBar(pb, iteration)
        }
      } else {
        # Calculate lambda and gamma values for this iteration
        l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
        g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
        
        # Call the compute_power_nested function without eta
        result <-  get_boundary_oc_coprimary(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                  lambda = l, gamma = g, eta = NULL, method = method, seed=seed)
        # Append the results
        results <- rbind(results, result)
        
        # Update progress bar
        iteration <- iteration + 1
        utils::setTxtProgressBar(pb, iteration)
      }
    }
  }
  close(pb)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- results %>% dplyr::filter(!is.na(rejectnull_mean_h0) & rejectnull_mean_h0  < t1e)
  
  # Sort results by power
  results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(rejectnull_mean_h1 ))
  
  return(results_sorted)
}

#' Search optimal parameters for joint efficacy and toxicity endpoint
#' 
#' `search_optimal_pars_efftox()` is a helper function to run a grid search to find the optimal parameters for futility and efficacy stopping cut-off probabilities. The function runs on 
#' a single core and can take a long time if the search space is too wide. If one has access to multiple cores, one can directly call
#'  `get_boundary_oc_efftox()` function and run it in parallel. 
#'   
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Response - toxicity,
#' - `H0[2]`: Response - no toxicity,
#' - `H0[3]`: No Response - toxicity,
#' - `H0[4]`: No Response - no toxicity
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Response - toxicity,
#' - `H1[2]`: Response - no toxicity,
#' - `H1[3]`: No Response - toxicity,
#' - `H1[4]`: No Response - no toxicity.
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param t1e Desired Type - I error rate. If specified it will only return results with type I error rate less the specified value 
#' @param lambda1 starting value for `lambda` values to search
#' @param lambda2 ending value for `lambda` values to search
#' @param grid1 number of `lambda` values to consider between lambda1 and lambda2
#' @param gamma1 starting value for `gamma` values to search
#' @param gamma2 ending value for `gamma` values to search
#' @param grid2 number of `gamma` values to consider between gamma1 and gamma2
#' @param eta1 starting value for `eta` values to search
#' @param eta2 ending value for `eta` values to search
#' @param grid3 number of eta values to consider between eta1 and eta2
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stoping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
search_optimal_pars_efftox <- function(H0, H1, n, nsim, t1e=NULL, method ="power",
                                          lambda1, lambda2, grid1, gamma1, gamma2, 
                                          grid2, eta1, eta2, grid3, seed=NULL) {
  
  # Initialize result a data frame
  nIA <- length(n)
  column_names <-  c(paste0("fut_boundary_resp_", 1:nIA), 
                     paste0("fut_boundary_tox_", 1:nIA), 
                     paste0("Sup_boundary_resp_", 1:nIA), 
                     paste0("Sup_boundary_tox_", 1:nIA), 
                     "earlystopfuti_mean_h0", 
                     "earlystopsupe_mean_h0", 
                     "ss_mean_h0", 
                     "rejectnull_mean_h0", 
                     "earlystopfuti_mean_h1", 
                     "earlystopsupe_mean_h1", 
                     "ss_mean_h1", 
                     "rejectnull_mean_h1", 
                     "lambda", "gamma", "eta")
  
  results <- data.frame(matrix(0, nrow = 0, ncol = length(column_names)))
  colnames(results) <- column_names
  
  
  # Calculate total iterations
  total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
  pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
  iteration <- 0
  
  # Iterate over grid values
  for (i in 1:(grid1 + 1)) {
    for (j in 1:(grid2 + 1)) {
      if (method == "power") {
        for (et in 1:(grid3 + 1)) {
          # Calculate lambda, gamma, and eta values for this iteration
          l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
          g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
          e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
          
          # Call the compute_power_nested function
          result <-  get_boundary_oc_efftox(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                               lambda = l, gamma = g, eta = e, method = method, seed=seed)
          
          # Append the results
          results <- rbind(results, result)
          
          # Update progress bar
          iteration <- iteration + 1
          utils::setTxtProgressBar(pb, iteration)
        }
      } else {
        # Calculate lambda and gamma values for this iteration
        l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
        g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
        
        # Call the compute_power_nested function without eta
        result <-  get_boundary_oc_efftox(H0 = H0, H1 = H1, nsim=nsim, n = n,
                                  lambda = l, gamma = g, eta = NULL, method = method, seed=seed)
        # Append the results
        results <- rbind(results, result)
        
        # Update progress bar
        iteration <- iteration + 1
        utils::setTxtProgressBar(pb, iteration)
      }
    }
  }
  close(pb)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- results %>% dplyr::filter(!is.na(rejectnull_mean_h0) & rejectnull_mean_h0  < t1e)
  
  # Sort results by power
  results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(rejectnull_mean_h1 ))
  
  return(results_sorted)
}