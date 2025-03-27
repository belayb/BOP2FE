
#' Compute Probability Cutoffs for Futility Stopping
#'
#' This function computes the probability cutoffs for futility stopping using two methods:
#' power and O'Brien-Fleming type function.
#'
#' @param n A vector of sample sizes for each interim analysis of the trial.
#' @param lambda A numeric value representing the baseline probability (default is NULL).
#' @param gamma A numeric value used in the calculation of cutoff values (default is NULL).
#' @param eta A numeric value used in the power method calculation (default is NULL).
#' @param method A character string specifying the method to use for calculating cutoff values.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#'
#' @return A list containing two elements:
#' \item{cf_values}{A numeric vector of cutoff values for futility stopping.}
#' \item{cs_values}{A numeric vector of corresponding probability cutoffs.}
#'
#' @importFrom magrittr %>%
#' @importFrom stats qnorm pnorm
#' @export
get_cf_cs_values<- function(n, lambda=NULL, gamma=NULL, eta= NULL, method = "power"){
  nsum<- sum(n)
  cf_values <- sapply(seq_along(n), function(i) {
    lambda * (sum(n[1:i]) / nsum)^gamma
  })
  
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
  return(list(cf_values=cf_values, cs_values=cs_values))
}


# #' Simulate binary endpoint data
# #'
# #' This function simulates binary data across multiple stages based on specified probabilities and sample sizes.
# #'
# #' @param p Probability of success.
# #' @param n_stage Number of interim analysis in the simulation.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param seed An integer seed for reproducibility of the random number generation.
# #'
# #' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
# #'
# #' @importFrom magrittr %>%
# #' 
# simulate_data_binary2 <- function(p, n_stage, n, nsim, seed = 23456) {
#   set.seed(seed)
#   # Initialize empty data frame
#   resp <- data.frame()
#   # Simulate data for each simulation run
#   for (obs in 1:nsim) {
#     Y_stage <- list()  # List to hold results for each stage
#     for (i in 1:n_stage) {
#       Y_stage[[i]] <- rbinom(1, n[i], p)
#     }
#     # Combine all stages into a single data frame row
#     row_data <- unlist(Y_stage)
#     resp <- rbind(resp, c(obs, row_data))
#   }
#   # Assign dynamic column names for the response data frame
#   #colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("Y", i))))
#   colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("x", i))))
#   
#   # Calculate cumulative sums
#   for (i in 1:n_stage) {
#     #Y_cols <- paste0("Y", 1:i)
#     x_cols <- paste0("x", 1:i)
#     
#     # Calculate cumulative sum values for each stage
#     #resp[[paste0("tsum", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols)))
#     resp[[paste0("sumx", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(x_cols)))
#     
#   }
#   return(resp)
# }
# 
# #' Simulate binary endpoint data2
# #'
# #' This function simulates binary data across multiple stages based on specified probabilities and sample sizes.
# #'
# #' @param p Probability of success.
# #' @param n_stage Number of interim analysis in the simulation.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param seed An integer seed for reproducibility of the random number generation.
# #'
# #' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
# #'
# #' @importFrom magrittr %>%
# #' 
# simulate_data_binary <- function(p, n_stage, n, nsim, seed = 23456) {
#   set.seed(seed)
#   # Initialize empty data frame
#   resp <- data.frame()
#   # Simulate data for each simulation run
#   for (obs in 1:nsim) {
#     Y_stage <- list()  # List to hold results for each stage
#     for (i in 1:n_stage) {
#       Y_stage[[i]] <- rbinom(1, n[i], p)
#     }
#     # Combine all stages into a single data frame row
#     row_data <- unlist(Y_stage)
#     resp <- rbind(resp, c(obs, row_data))
#   }
#   # Assign dynamic column names for the response data frame
#   colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("Y", i))))
# 
#   # Calculate cumulative sums
#   for (i in 1:n_stage) {
#     Y_cols <- paste0("Y", 1:i)
# 
#     # Calculate cumulative sum values for each stage
#     resp[[paste0("tsum", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols)))
# 
#   }
#   return(resp)
# }
# 
# 
# #' Simulate nested endpoint data
# #'
# #' This function simulates nested data across multiple stages based on specified probabilities and sample sizes.
# #'
# #' @param p1 Probability of success for the first outcome.
# #' @param p2 Probability of success for the second outcome.
# #' @param p3 Probability of success for the third outcome.
# #' @param n_stage Number of interim analysis in the simulation.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param seed An integer seed for reproducibility of the random number generation.
# #'
# #' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
# #'
# #' @importFrom magrittr %>%
# #' 
# simulate_data_nested <- function(p1, p2, p3, n_stage, n, nsim, seed = 23456) {
#   set.seed(seed)
#   # Initialize empty data frame
#   resp <- data.frame()
#   # Simulate data for each simulation run
#   for (obs in 1:nsim) {
#     Y_stage <- list()  # List to hold results for each stage
#     
#     # Simulate data for each stage
#     for (i in 1:n_stage) {
#       Y1 <- rbinom(1, n[i], p1)
#       Y2 <- rbinom(1, n[i] - Y1, p2 / (p2 + p3))
#       Y3 <- n[i] - Y1 - Y2
#       
#       # Store the values for this stage
#       Y_stage[[i]] <- c(Y1, Y2, Y3)
#     }
#     
#     # Combine all stages into a single data frame row
#     row_data <- unlist(Y_stage)
#     resp <- rbind(resp, c(obs, row_data))
#   }
#   
#   # Assign dynamic column names for the response data frame
#   colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y"), i, 1:3))))
#   
#   # Calculate cumulative sums
#   for (i in 1:n_stage) {
#     Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
#     Y_cols_i2 <- paste0("Y", 1:i, "2")  # Collect all Y12, Y22,... Yi2 columns
#     Y_cols_i3 <- paste0("Y", 1:i, "3")  # Collect all Y13, Y23,... Yi3 columns
#     
#     # Calculate cumulative sum values for each stage
#     resp[[paste0("tsum", i, "1")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i1)))
#     resp[[paste0("tsum", i, "2")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i2)))
#     resp[[paste0("tsum", i, "3")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i3)))
#   }
#   return(resp)
# }
# 
# 
# 
# 
# #' Simulate co-primary endpoint data
# #'
# #' This function simulates co-primary endpoint data across multiple stages based on specified probabilities and sample sizes.
# #'
# #' @param p1 Probability of response and progression-free survival at 6 months (PFS6).
# #' @param p2 Probability of response and no PFS6 .
# #' @param p3 Probability of no response and PFS6.
# #' @param p4 Probability of no response and no PFS6.
# #' @param n_stage Number of interim analysis in the simulation.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param seed An integer seed for reproducibility of the random number generation.
# #'
# #' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
# #'
# #' @importFrom magrittr %>%
# #' 
# simulate_data_coprimary <- function(p1, p2, p3, p4, n_stage, n, nsim, seed = 23456) {
#   set.seed(seed)
#   # Initialize empty data frame
#   resp <- data.frame()
#   # Simulate data for each simulation run
#   for (obs in 1:nsim) {
#     Y_stage <- list()  # List to hold results for each stage
#     
#     # Simulate data for each stage
#     for (i in 1:n_stage) {
#       Y1 <- rbinom(1, n[i], p1)
#       Y2 <- rbinom(1, n[i] - Y1, p2 / (p2 + p3 + p4))
#       Y3 <- rbinom(1, n[i] - Y1 - Y2, p3 / (p3 + p4))
#       Y4 <- n[i] - Y1 - Y2 - Y3
#       
#       # Store the values for this stage
#       Y_stage[[i]] <- c(Y1, Y2, Y3, Y4)
#     }
#     
#     # Combine all stages into a single data frame row
#     row_data <- unlist(Y_stage)
#     resp <- rbind(resp, c(obs, row_data))
#   }
#   
#   # Assign dynamic column names for the response data frame
#   colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y", "Y"), i, 1:4))))
#   
#   # Calculate cumulative sums
#   for (i in 1:n_stage) {
#    Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
#     Y_cols_i2 <- paste0("Y", 1:i, "2")  # Collect all Y12, Y22,... Yi2 columns
#     Y_cols_i3 <- paste0("Y", 1:i, "3")  # Collect all Y13, Y23,... Yi3 columns
#     Y_cols_i4 <- paste0("Y", 1:i, "4")  # Collect all Y14, Y24,... Yi4 columns
#     
#     # Calculate cumulative sum values for each stage
#     resp[[paste0("tsum", i, "1")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i1)))
#     resp[[paste0("tsum", i, "2")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i2)))
#     resp[[paste0("tsum", i, "3")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i3)))
#     resp[[paste0("tsum", i, "4")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i4)))
#     
#     }
#   return(resp)
# }
# 
# 
# 
# #' Simulate joint efficacy and toxicity endpoint data
# #'
# #' This function simulates joint efficacy and toxicity endpoint data across multiple stages based on specified probabilities and sample sizes.
# #'
# #' @param p1 Probability of response and toxicity.
# #' @param p2 Probability of response and no toxicity .
# #' @param p3 Probability of no response and toxicity
# #' @param p4 Probability of no response and no toxicity
# #' @param n_stage Number of interim analysis in the simulation.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param seed An integer seed for reproducibility of the random number generation.
# #'
# #' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
# #'
# #' @importFrom magrittr %>%
# #' 
# simulate_data_jointefftox <- function(p1, p2, p3, p4, n_stage, n, nsim, seed = 23456) {
#   set.seed(seed)
#   # Initialize empty data frame
#   resp <- data.frame()
#   # Simulate data for each simulation run
#   for (obs in 1:nsim) {
#     Y_stage <- list()  # List to hold results for each stage
#     
#     # Simulate data for each stage
#     for (i in 1:n_stage) {
#       Y1 <- rbinom(1, n[i], p1)
#       Y2 <- rbinom(1, n[i] - Y1, p2 / (p2 + p3 + p4))
#       Y3 <- rbinom(1, n[i] - Y1 - Y2, p3 / (p3 + p4))
#       Y4 <- n[i] - Y1 - Y2 - Y3
#       
#       # Store the values for this stage
#       Y_stage[[i]] <- c(Y1, Y2, Y3, Y4)
#     }
#     
#     # Combine all stages into a single data frame row
#     row_data <- unlist(Y_stage)
#     resp <- rbind(resp, c(obs, row_data))
#   }
#   
#   # Assign dynamic column names for the response data frame
#   colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y", "Y"), i, 1:4))))
#   
#   # Calculate cumulative sums
#   for (i in 1:n_stage) {
#     Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
#     Y_cols_i2 <- paste0("Y", 1:i, "2")  # Collect all Y12, Y22,... Yi2 columns
#     Y_cols_i3 <- paste0("Y", 1:i, "3")  # Collect all Y13, Y23,... Yi3 columns
#     Y_cols_i4 <- paste0("Y", 1:i, "4")  # Collect all Y14, Y24,... Yi4 columns
#     
#     # Calculate cumulative sum values for each stage
#     resp[[paste0("tsum", i, "1")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i1)))
#     resp[[paste0("tsum", i, "2")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i2)))
#     resp[[paste0("tsum", i, "3")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i3)))
#     resp[[paste0("tsum", i, "4")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i4)))
#     
#   }
#   return(resp)
# }
# 
# #' Compute Posterior Probability of Observing a given Response Under the null for a binary outcome
# #' @param xdata a given data set
# #' @param H0 null hypothesis
# #' @param H1 alternative hypothesis
# #' @param n_stage the number of interim analysis
# #' @param n A vector of sample sizes for each interim analysis.
# #'
# 
# calculate_posterior_binary <- function(xdata, H0, H1, n_stage, n,...) {
#   .Defunct()
#   
#   H0 <- as.numeric(H0)
#   H1 <- as.numeric(H1)
#   
#   
#   a1 <- H0
#   b1 <- 1 - a1
#   
#  for (i in 1:n_stage) {
#       tsum_col <- paste0("tsum", i)
#       postp_col <- paste0("postp", i)
#       xdata[[postp_col]] <- 1 - pbeta(H0, a1 + xdata[[tsum_col]], b1 + sum(n[1:i]) - xdata[[tsum_col]])
#     }
#     return(xdata)
#   }
#   
# 
# #' Compute Posterior Probability of Observing a given Response Under the null for a nested outcome
# #'
# #' This function computes the posterior probability of observing a given response under the null hypothesis
# #' across multiple stages of a trial.
# #'
# #' @param xdata A data frame containing the data for which posterior probabilities are to be calculated.
# #' @param H0 A numeric vector representing the null hypothesis probabilities.
# #' @param H1 A numeric vector representing the alternative hypothesis probabilities.
# #' @param n_stage An integer indicating the number of interim analysis in the trial.
#
# #' @return A data frame with additional columns for the posterior probabilities calculated for each stage.
# #'
# #' @importFrom magrittr %>%
# #' @importFrom stats pbeta
# #' @importFrom dplyr select
# #' @importFrom tidyselect all_of
# #' @importFrom purrr pmap_dbl
# #'
# calculate_posterior_nested <- function(xdata, H0, H1, n_stage,...) {
#   .Defunct()
#   
#   H0 <- as.numeric(H0)
#   H1 <- as.numeric(H1)
#   
#   b1 <- c(1, 0, 0)
#   b2 <- c(1, 1, 0)
#   a <- H0
#   
#   phi1 <- sum(H0 * b1)
#   phi2 <- sum(H0 * b2)
#   
#   post_probs <- function(x, phi, b, a) {
#     beta_a <- sum(b * (a + x))
#     beta_b <- sum((1 - b) * (a + x))
#     1 - pbeta(phi, beta_a, beta_b)
#   }
#   
#   phi_list <- list(phi1, phi2)
#   b_list <- list(b1, b2)
#   
#   for (i in 1:n_stage) {
#     tsum_cols <- paste0("tsum", i, 1:3)
#     for (j in 1:2) {
#       phi <- phi_list[[j]]
#       b <- b_list[[j]]
#       postp_col <- paste0("postp", i, j)
#       
#       xdata <- xdata %>%
#         dplyr::mutate(!!postp_col := purrr::pmap_dbl(dplyr::select(., dplyr::all_of(tsum_cols)), ~ post_probs(c(..1, ..2, ..3), phi, b, a)))
#     }
#   }
#   
#   return(xdata)
# }
# 
# 
# #' Compute Posterior Probability of Observing a given Response Under the null for a co-primary outcome
# #'
# #' This function computes the posterior probability of observing a given response under the null hypothesis
# #' across multiple stages of a trial.
# #'
# #' @param xdata A data frame containing the data for which posterior probabilities are to be calculated.
# #' @param H0 A numeric vector representing the null hypothesis probabilities.
# #' @param H1 A numeric vector representing the alternative hypothesis probabilities.
# #' @param n_stage An integer indicating the number of interim analysis in the trial.
# #'
# #' @return A data frame with additional columns for the posterior probabilities calculated for each stage.
# #'
# #' @importFrom magrittr %>%
# #' @importFrom stats pbeta
# #' @importFrom dplyr select
# #' @importFrom tidyselect all_of
# #' @importFrom purrr pmap_dbl
# #'
# #' calculate_posterior_coprimary <- function(xdata, H0, H1, n_stage,...) {
#   .Defunct()
#   
#   H0 <- as.numeric(H0)
#   H1 <- as.numeric(H1)
#   
#   b1 <- c(1, 1, 0, 0)
#   b2 <- c(1, 0, 1, 0)
#   a <- H0
#   
#   phi1 <- sum(H0 * b1)
#   phi2 <- sum(H0 * b2)
#   
#   post_probs <- function(x, phi, b, a) {
#     beta_a <- sum(b * (a + x))
#     beta_b <- sum((1 - b) * (a + x))
#     1 - pbeta(phi, beta_a, beta_b)
#   }
#   
#   phi_list <- list(phi1, phi2)
#   b_list <- list(b1, b2)
#   
#   for (i in 1:n_stage) {
#     tsum_cols <- paste0("tsum", i, 1:4)
#     for (j in 1:2) {
#       phi <- phi_list[[j]]
#       b <- b_list[[j]]
#       postp_col <- paste0("postp", i, j)
#       
#       xdata <- xdata %>%
#         dplyr::mutate(!!postp_col := purrr::pmap_dbl(dplyr::select(., dplyr::all_of(tsum_cols)), ~ post_probs(c(..1, ..2, ..3, ..4), phi, b, a)))
#     }
#   }
#   
#   return(xdata)
# }
# 
# 
# 
# #' Compute Posterior Probability of Observing a given Response Under the null for a joint efficacy and toxicity outcome
# #'
# #' This function computes the posterior probability of observing a given response under the null hypothesis
# #' across multiple stages of a trial.
# #'
# #' @param xdata A data frame containing the data for which posterior probabilities are to be calculated.
# #' @param H0 A numeric vector representing the null hypothesis probabilities.
# #' @param H1 A numeric vector representing the alternative hypothesis probabilities.
# #' @param n_stage An integer indicating the number of interim analysis in the trial.
# #'
# #' @return A data frame with additional columns for the posterior probabilities calculated for each stage.
# #'
# #' @importFrom magrittr %>%
# #' @importFrom stats pbeta
# #' @importFrom dplyr select
# #' @importFrom tidyselect all_of
# #' @importFrom purrr pmap_dbl
# #'
# calculate_posterior_jointefftox <- function(xdata, H0, H1, n_stage, ...) {
#   
#   .Defunct()
#   
#   H0 <- as.numeric(H0)
#   H1 <- as.numeric(H1)
#   
#   b1 <- c(1, 1, 0, 0)
#   b2 <- c(1, 0, 1, 0)
#   a <- H0
#   
#   phi1 <- sum(H0 * b1)
#   phi2 <- sum(H0 * b2)
#   
#   post_probs <- function(x, phi, b, a) {
#     beta_a <- sum(b * (a + x))
#     beta_b <- sum((1 - b) * (a + x))
#     1 - pbeta(phi, beta_a, beta_b)
#   }
#   post_probs2 <- function(x, phi, b, a) {
#     beta_a <- sum(b * (a + x))
#     beta_b <- sum((1 - b) * (a + x))
#     pbeta(phi, beta_a, beta_b)
#   }
#   
#   phi_list <- list(phi1, phi2)
#   b_list <- list(b1, b2)
#   
#   for (i in 1:n_stage) {
#     tsum_cols <- paste0("tsum", i, 1:4)
#     for (j in 1:2) {
#       phi <- phi_list[[j]]
#       b <- b_list[[j]]
#       postp_col <- paste0("postp", i, j)
#       
#       if (j == 1) {
#         xdata <- xdata %>%
#           dplyr::mutate(!!postp_col := purrr::pmap_dbl(dplyr::select(., dplyr::all_of(tsum_cols)), ~ post_probs(c(..1, ..2, ..3, ..4), phi, b, a)))
#       } else {
#         xdata <- xdata %>%
#           dplyr::mutate(!!postp_col := purrr::pmap_dbl(dplyr::select(., dplyr::all_of(tsum_cols)), ~ post_probs2(c(..1, ..2, ..3, ..4), phi, b, a)))
#       }
#     }
#   }
#   
#   return(xdata)
# }
# 
# 
# #' Compute power and type I error rate for binary endpoint 
# #' These functions are gone, no longer available use get_boundary_oc_binary() instead.
# #' @param H0 the null hypothesis
# #' @param H1 the alternative hypothesis
# #' @param n vector of sample sizes
# #' @param nsim number of simulation to run
# #' @param lambda lambda value for
# #' @param gamma gamma value for
# #' @param eta eta value for
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming". 
# #' @param seed for reproducability 
# #' @importFrom magrittr %>%
# #' @importFrom dplyr select mutate case_when
# #' @importFrom purrr reduce
# #'
# #'
# compute_power_binary <- function(H0, H1, n, nsim, lambda = NULL, gamma = NULL, eta = NULL, method = "power", seed = NULL,...){
#   
#   .Defunct("get_boundary_oc_binary", msg="use get_boundary_oc_binary() instead")
#   
#   n_stage <- length(n)
#   
#   y_null <- simulate_data_binary(p=H0, n=n, n_stage=n_stage , nsim = nsim)
#   y_alter <- simulate_data_binary(p=H1, n=n, n_stage=n_stage , nsim = nsim)
#   
#   cfcs<-get_cf_cs_values(n, lambda, gamma, eta, method)
#   cf_values<-cfcs$cf_values
#   cs_values<-cfcs$cs_values
#   
#   
#   y_null <- y_null%>%
#     calculate_posterior_binary(., H0, H1, n_stage, n)%>%
#     mutate(reject = case_when(
#       reduce(1:n_stage, ~ .x | (get(paste0("postp", .y)) < cf_values[.y]), .init = FALSE) ~ 0,
#       reduce(1:n_stage, ~ .x | (get(paste0("postp", .y)) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   y_alter <- y_alter%>%
#     calculate_posterior_binary(., H0, H1, n_stage, n)%>%
#     mutate(power = case_when(
#       reduce(1:n_stage, ~ .x | (get(paste0("postp", .y)) < cf_values[.y] ), .init = FALSE) ~ 0,
#       reduce(1:n_stage, ~ .x | (get(paste0("postp", .y)) >= cs_values[.y] ), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   
#   # Summarizing results
#   reject_mean <- mean(y_null$reject)
#   power_mean <- mean(y_alter$power)
#   
#   if(is.null(eta)) eta=NA
#     
#   result <- data.frame(
#     H0 = H0,
#     H1 = H1,
#     Lambda = lambda,
#     Gamma = gamma,
#     eta = eta,
#     reject_mean = reject_mean,
#     power_mean = power_mean
#   )
#   
#   return(result)
#   
# }
# 
# #' Compute power and type I error rate for nested or ordinal endpoint
# #' These functions are gone, no longer available use get_boundary_oc_nested() instead. 
# #' @param H0 the null hypothesis
# #' @param H1 the alternative hypothesis
# #' @param n vector of sample sizes
# #' @param nsim number of simulation to run
# #' @param lambda lambda value for
# #' @param gamma gamma value for
# #' @param eta eta value for
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming". 
# #' @param seed for reproducability 
# #' @importFrom magrittr %>%
# #' @importFrom dplyr select mutate case_when
# #' @importFrom purrr reduce
# #' 
# #' 
# compute_power_nested <- function(H0, H1, n, nsim, lambda = NULL, gamma = NULL, eta = NULL, method = "power", seed = NULL,...) {
#   
#   .Defunct("get_boundary_oc_nested", msg="use get_boundary_oc_nested() instead")
#   
#   
#   n_stage <- length(n)
#   H0 <- H0
#   H1 <- H1
#   nsim <- nsim
#   
#   y_null <- simulate_data_nested(p1 = H0[1], p2 = H0[2], p3 = H0[3], n = n, n_stage = n_stage, nsim = nsim)
#   y_alter <- simulate_data_nested(p1 = H1[1], p2 = H1[2], p3 = H1[3], n = n, n_stage = n_stage, nsim = nsim)
#   
#   cfcs <- get_cf_cs_values(n, lambda, gamma, eta, method)
#   cf_values <- cfcs$cf_values
#   cs_values <- cfcs$cs_values
#   
#   y_null <- y_null %>%
#     calculate_posterior_nested(., H0, H1, n_stage) %>%
#     dplyr::mutate(reject = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   y_alter <- y_alter %>%
#     calculate_posterior_nested(., H0, H1, n_stage) %>%
#     dplyr::mutate(power = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   # Summarizing results
#   reject_mean <- mean(y_null$reject)
#   power_mean <- mean(y_alter$power)
#   
#   if(is.null(eta)) eta=NA
#   
#   result <- dplyr::tibble(
#     H0 = paste(H0[1], H0[2], sep = ", "),
#     H1 = paste(H1[1], H1[2], sep = ", "),
#     Lambda = lambda,
#     Gamma = gamma,
#     eta = eta,
#     reject_mean = reject_mean,
#     power_mean = power_mean
#   )
#   
#   return(result)
# }
# 
# 
# 
# #' Compute power and type I error rate for co-primary endpoint
# #' These functions are gone, no longer available use get_boundary_oc_coprimary() instead.
# #'  
# #' @param H0 the null hypothesis
# #' @param H1 the alternative hypothesis
# #' @param n vector of sample sizes
# #' @param nsim number of simulation to run
# #' @param lambda lambda value for
# #' @param gamma gamma value for
# #' @param eta eta value for
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming". 
# #' @param seed for reproducability 
# #' @importFrom magrittr %>%
# #' @importFrom dplyr select mutate case_when
# #' @importFrom purrr reduce
# #' 
# #' 
# compute_power_coprimary <- function(H0, H1, n, nsim, lambda = NULL, gamma = NULL, eta = NULL, method = "power", seed = NULL, ...) {
#   
#   .Defunct("get_boundary_oc_coprimary", msg="use get_boundary_oc_coprimary() instead")
#   
#   n_stage <- length(n)
#   H0 <- H0
#   H1 <- H1
#   nsim <- nsim
#   
#   y_null <- simulate_data_coprimary(p1 = H0[1], p2 = H0[2], p3 = H0[3], p4 = H0[4], n = n, n_stage = n_stage, nsim = nsim)
#   y_alter <- simulate_data_coprimary(p1 = H1[1], p2 = H1[2], p3 = H1[3], p4 = H1[4], n = n, n_stage = n_stage, nsim = nsim)
#   
#   cfcs <- get_cf_cs_values(n, lambda, gamma, eta, method)
#   cf_values <- cfcs$cf_values
#   cs_values <- cfcs$cs_values
#   
#   y_null <- y_null %>%
#     calculate_posterior_coprimary(., H0, H1, n_stage) %>%
#     dplyr::mutate(reject = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   y_alter <- y_alter %>%
#     calculate_posterior_coprimary(., H0, H1, n_stage) %>%
#     dplyr::mutate(power = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   # Summarizing results
#   reject_mean <- mean(y_null$reject)
#   power_mean <- mean(y_alter$power)
#   
#   if(is.null(eta)) eta=NA
#   
#   result <- dplyr::tibble(
#     H0 = paste(H0[1], H0[2], sep = ", "),
#     H1 = paste(H1[1], H1[2], sep = ", "),
#     Lambda = lambda,
#     Gamma = gamma,
#     eta = eta,
#     reject_mean = reject_mean,
#     power_mean = power_mean
#   )
#   
#   return(result)
# }
# 
# 
# 
# #' Compute power and type I error rate for joint efficacy and toxicity endpoint 
# #' These functions are gone, no longer available use get_boundary_oc_efftox() instead.
# #' @param H0 the null hypothesis
# #' @param H1 the alternative hypothesis
# #' @param n vector of sample sizes
# #' @param nsim number of simulation to run
# #' @param lambda lambda value for
# #' @param gamma gamma value for
# #' @param eta eta value for
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming". 
# #' @param seed for reproducability 
# #' @importFrom magrittr %>%
# #' @importFrom dplyr select mutate case_when
# #' @importFrom purrr reduce
# #' 
# #' 
# compute_power_jointefftox <- function(H0, H1, n, nsim, lambda = NULL, gamma = NULL, eta = NULL, method = "power", seed = NULL, ...) {
#   
#   .Defunct("get_boundary_oc_efftox", msg = "use get_boundary_oc_efftox() instead")
#   
#   n_stage <- length(n)
#   H0 <- H0
#   H1 <- H1
#   nsim <- nsim
#   
#   y_null <- simulate_data_jointefftox(p1 = H0[1], p2 = H0[2], p3 = H0[3], p4 = H0[4], n = n, n_stage = n_stage, nsim = nsim)
#   y_alter <- simulate_data_jointefftox(p1 = H1[1], p2 = H1[2], p3 = H1[3], p4 = H1[4], n = n, n_stage = n_stage, nsim = nsim)
#   
#   cfcs <- get_cf_cs_values(n, lambda, gamma, eta, method)
#   cf_values <- cfcs$cf_values
#   cs_values <- cfcs$cs_values
#   
#   y_null <- y_null %>%
#     calculate_posterior_jointefftox(., H0, H1, n_stage) %>%
#     dplyr::mutate(reject = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] | get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] & get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   y_alter <- y_alter %>%
#     calculate_posterior_jointefftox(., H0, H1, n_stage) %>%
#     dplyr::mutate(power = dplyr::case_when(
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] | get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
#       purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] & get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
#       TRUE ~ 1
#     ))
#   
#   # Summarizing results
#   reject_mean <- mean(y_null$reject)
#   power_mean <- mean(y_alter$power)
#   
#   if(is.null(eta)) eta=NA
#   
#   result <- dplyr::tibble(
#     H0 = paste(H0[1], H0[2], sep = ", "),
#     H1 = paste(H1[1], H1[2], sep = ", "),
#     Lambda = lambda,
#     Gamma = gamma,
#     eta = eta,
#     reject_mean = reject_mean,
#     power_mean = power_mean
#   )
#   
#   return(result)
# }
# 
# 
# #' Optimal parameters for binary outcome
# #' These functions are gone, no longer available use search_optimal_pars_binary(). 
# #' This function runs a grid search to find the optimal parameters for binary outcome that leads to 
# #' the highest power while controlling the type - I error at the specified level.
# #'
# #' @param H0 Probability of success under the null hypothesis.
# #' @param H1 Probability of success under the alternative hypothesis.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param t1e Desired Type - I error rate
# #' @param t2e Desired Type - II error rate
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
# #' @param lambda1 starting value for lambda values to search
# #' @param lambda2 ending value for lambda values to search
# #' @param grid1 number of lambda values to consider between lambda1 and lambda2
# #' @param gamma1 starting value for gamma values to search
# #' @param gamma2 ending value for gamma values to search
# #' @param grid2 number of gamma values to consider between gamma1 and gamma2
# #' @param eta1 starting value for eta values to search
# #' @param eta2 ending value for eta values to search
# #' @param grid3 number of eta values to consider between eta1 and eta2
# #'
# #' @return A data frame
# #'
# #' @importFrom magrittr %>%
# #' @importFrom utils setTxtProgressBar txtProgressBar
# #' 
# Optimal_pars_binary <- function(H0, H1, n, nsim, t1e, t2e, method ="power",
#                                 lambda1, lambda2, grid1, gamma1, gamma2, grid2, eta1, eta2, grid3,...) {
#  
#   .Defunct("search_optimal_pars_binary", msg = "use search_optimal_pars_binary() instead")
#   
#   # Initialize result a data frame
#   results <- data.frame(
#     reject = numeric(0),
#     power = numeric(0),
#     H0 = character(0),
#     H1 = character(0),
#     lambda = numeric(0),
#     gamma = numeric(0),
#     eta = numeric(0)
#   )
#   
#   
#   # Calculate total iterations
#   total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
#   pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
#   iteration <- 0
#   
#   # Iterate over grid values
#   for (i in 1:(grid1 + 1)) {
#     for (j in 1:(grid2 + 1)) {
#       if (method == "power") {
#         for (et in 1:(grid3 + 1)) {
#           # Calculate lambda, gamma, and eta values for this iteration
#           l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#           g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#           e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
#           
#           # Call the compute_power_nested function
#           result <- compute_power_binary(H0=H0, H1=H1, n=n,nsim=nsim, lambda=l, gamma=g, eta=e, method=method)
#           
#           # Append the results
#           results <- rbind(results, result)
#           
#           # Update progress bar
#           iteration <- iteration + 1
#           utils::setTxtProgressBar(pb, iteration)
#         }
#       } else {
#         # Calculate lambda and gamma values for this iteration
#         l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#         g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#         
#         # Call the compute_power_nested function without eta
#         result <- compute_power_binary(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=NULL, method=method)
#         
#         # Append the results
#         results <- rbind(results, result)
#         
#         # Update progress bar
#         iteration <- iteration + 1
#         utils::setTxtProgressBar(pb, iteration)
#       }
#     }
#   }
#   close(pb)
#   
#   # Filter results by t1e
#   results_filtered <- results %>% dplyr::filter(!is.na(reject_mean) & reject_mean  < t1e)
#   
#   # Sort results by power
#   results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(power_mean))
#   
#   return(results_sorted)
# }
# 
# 
# 
# #' Optimal parameters for nested outcome
# #' These functions are gone, no longer available, use search_optimal_pars_nested() instead.
# #' This function runs a grid search to find the optimal parameters for nested or an ordinal outcome that leads to the highest power while controlling the type - I error at the specified level.
# #'
# #' @param CR0 Probability of success for the first outcome under the null hypothesis.
# #' @param CRPR0 Probability of success for the first or second outcome under the null hypothesis.
# #' @param CR1 Probability of success for the first outcome under the alternative hypothesis.
# #' @param CRPR1 Probability of success for the first or second outcome under the alternative hypothesis.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param t1e Desired Type - I error rate
# #' @param t2e Desired Type - II error rate
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
# #' @param lambda1 starting value for lambda values to search
# #' @param lambda2 ending value for lambda values to search
# #' @param grid1 number of lambda values to consider between lambda1 and lambda2
# #' @param gamma1 starting value for gamma values to search
# #' @param gamma2 ending value for gamma values to search
# #' @param grid2 number of gamma values to consider between gamma1 and gamma2
# #' @param eta1 starting value for eta values to search
# #' @param eta2 ending value for eta values to search
# #' @param grid3 number of eta values to consider between eta1 and eta2
# #'
# #' @return A data frame
# #'
# #' @importFrom magrittr %>%
# #' @importFrom utils setTxtProgressBar txtProgressBar
# #' 
# Optimal_pars_nested <- function(CR0, CRPR0, CR1, CRPR1, n, nsim, t1e, t2e, method ="power",
#                                 lambda1, lambda2, grid1, gamma1, gamma2, grid2, eta1, eta2, grid3, ...) {
#   .Defunct("search_optimal_pars_nested", msg = "use search_optimal_pars_nested() instead")
#   
#    # Define H0 and H1
#   H0 <- c(CR0, CRPR0 - CR0, 1 - CRPR0)
#   H1 <- c(CR1, CRPR1 - CR1, 1 - CRPR1)
#   nsum <- sum(n)
# 
#   # Initialize result a data frame
#   results <- data.frame(
#     reject = numeric(0),
#     power = numeric(0),
#     H0 = character(0),
#     H1 = character(0),
#     lambda = numeric(0),
#     gamma = numeric(0),
#     eta = numeric(0)
#   )
#   
# 
#   # Calculate total iterations
#   total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
#   pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
#   iteration <- 0
# 
#   # Iterate over grid values
#   for (i in 1:(grid1 + 1)) {
#     for (j in 1:(grid2 + 1)) {
#       if (method == "power") {
#         for (et in 1:(grid3 + 1)) {
#           # Calculate lambda, gamma, and eta values for this iteration
#           l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#           g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#           e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
# 
#           # Call the compute_power_nested function
#           result <- compute_power_nested(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=e, method=method)
# 
#           # Append the results
#           results <- rbind(results, result)
# 
#           # Update progress bar
#           iteration <- iteration + 1
#           utils::setTxtProgressBar(pb, iteration)
#         }
#       } else {
#         # Calculate lambda and gamma values for this iteration
#         l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#         g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
# 
#         # Call the compute_power_nested function without eta
#         result <- compute_power_nested(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=NULL, method=method)
# 
#         # Append the results
#         results <- rbind(results, result)
# 
#         # Update progress bar
#         iteration <- iteration + 1
#         utils::setTxtProgressBar(pb, iteration)
#       }
#    }
#   }
#   close(pb)
#   
#   # Filter results by t1e
#   results_filtered <- results %>% dplyr::filter(!is.na(reject_mean) & reject_mean  < t1e)
#   
#   # Sort results by power
#   results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(power_mean))
#   
#   return(results_sorted)
# }
# 
# 
# 
# 
# #' Optimal parameters for co-primary outcome
# #' These functions are gone, no longer available, use search_optimal_pars_coprimary() instead.  
# #' This function runs a grid search to find the optimal parameters for co-primary outcome that leads to the highest power while controlling the type - I error at the specified level.
# #'
# #' @param theta01 Probability of response and progression-free survival at 6 months (PFS6) under the null hypothesis.
# #' @param theta02 Probability of response and no PFS6 under the null hypothesis.
# #' @param theta03 Probability of no response and PFS6 under the null hypothesis.
# #' @param theta04 Probability of no response and no PFS6 under the null hypothesis.
# #' @param theta11 Probability of response and PFS6 under the alternative hypothesis.
# #' @param theta12 Probability of response and no PFS6 under the alternative hypothesis.
# #' @param theta13 Probability of no response and PFS6 under the alternative hypothesis.
# #' @param theta14 Probability of no response and no PFS6 under the alternative hypothesis.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param t1e Desired Type - I error rate
# #' @param t2e Desired Type - II error rate
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
# #' @param lambda1 starting value for lambda values to search
# #' @param lambda2 ending value for lambda values to search
# #' @param grid1 number of lambda values to consider between lambda1 and lambda2
# #' @param gamma1 starting value for gamma values to search
# #' @param gamma2 ending value for gamma values to search
# #' @param grid2 number of gamma values to consider between gamma1 and gamma2
# #' @param eta1 starting value for eta values to search
# #' @param eta2 ending value for eta values to search
# #' @param grid3 number of eta values to consider between eta1 and eta2
# #'
# #' @return A data frame
# #'
# #' @importFrom magrittr %>%
# #' @importFrom utils setTxtProgressBar txtProgressBar
# #' 
# Optimal_pars_coprimary <- function(theta01, theta02, theta03, theta04, 
#                                    theta11, theta12, theta13, theta14, 
#                                    n, nsim, t1e, t2e, method ="power", 
#                                    lambda1, lambda2, grid1, 
#                                    gamma1, gamma2, grid2, 
#                                    eta1, eta2, grid3, ...) {
#   
#   .Defunct("search_optimal_pars_coprimary", msg = "use search_optimal_pars_coprimary() instead")
#   
#   # Define H0 and H1
#   H0 <- c(theta01, theta02, theta03, theta04)
#   H1 <- c(theta11, theta12, theta13, theta14)
#   nsum <- sum(n)
#   
#   # Initialize result a data frame
#  results <- data.frame(
#     reject = numeric(0),
#     power = numeric(0),
#     H0 = character(0),
#     H1 = character(0),
#     lambda = numeric(0),
#     gamma = numeric(0),
#     eta = numeric(0)
#   )
#   
#   
#   # Calculate total iterations
#   total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
#   pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
#   iteration <- 0
#   
#   # Iterate over grid values
#   for (i in 1:(grid1 + 1)) {
#     for (j in 1:(grid2 + 1)) {
#       if (method == "power") {
#         for (et in 1:(grid3 + 1)) {
#           # Calculate lambda, gamma, and eta values for this iteration
#           l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#           g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#           e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
#           
#           # Call the compute_power_coprimary function
#           result <- compute_power_coprimary(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=e, method=method)
#           
#           # Append the results
#           results <- rbind(results, result)
#           
#           # Update progress bar
#           iteration <- iteration + 1
#           utils::setTxtProgressBar(pb, iteration)
#         }
#       } else {
#         # Calculate lambda and gamma values for this iteration
#         l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#         g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#         
#         # Call the compute_power_coprimary function without eta
#         result <- compute_power_coprimary(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=NULL, method=method)
#         
#         # Append the results
#         results <- rbind(results, result)
#         
#         # Update progress bar
#         iteration <- iteration + 1
#         utils::setTxtProgressBar(pb, iteration)
#       }
#     }
#   }
#   close(pb)
#   
#   # Filter results by t1e
#   results_filtered <- results %>% dplyr::filter(!is.na(reject_mean) & reject_mean  < t1e)
#   
#   # Sort results by power
#   results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(power_mean))
#   
#   return(results_sorted)
# }
# 
# 
# 
# 
# #' Optimal parameters for joint efficacy and toxicity outcome
# #'
# #' This function runs a grid search to find the optimal parameters for joint efficacy and toxicity outcome that leads to the highest power while controlling the type - I error at the specified level.
# #'
# #' @param theta01 Probability of response and toxicity under the null hypothesis.
# #' @param theta02 Probability of response and no toxicity under the null hypothesis.
# #' @param theta03 Probability of no response and toxicity under the null hypothesis.
# #' @param theta04 Probability of no response and no toxicity under the null hypothesis.
# #' @param theta11 Probability of response and toxicity under the alternative hypothesis.
# #' @param theta12 Probability of response and no toxicity under the alternative hypothesis.
# #' @param theta13 Probability of no response and toxicity under the alternative hypothesis.
# #' @param theta14 Probability of no response and no toxicity under the alternative hypothesis.
# #' @param n A vector of sample sizes for each interim analysis.
# #' @param nsim Number of simulation runs to perform.
# #' @param t1e Desired Type - I error rate
# #' @param t2e Desired Type - II error rate
# #' @param method A character string specifying the method to use for calculating cutoff values.
# #'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
# #' @param lambda1 starting value for lambda values to search
# #' @param lambda2 ending value for lambda values to search
# #' @param grid1 number of lambda values to consider between lambda1 and lambda2
# #' @param gamma1 starting value for gamma values to search
# #' @param gamma2 ending value for gamma values to search
# #' @param grid2 number of gamma values to consider between gamma1 and gamma2
# #' @param eta1 starting value for eta values to search
# #' @param eta2 ending value for eta values to search
# #' @param grid3 number of eta values to consider between eta1 and eta2
# #'
# #' @return A data frame
# #'
# #' @importFrom magrittr %>%
# #' @importFrom utils setTxtProgressBar txtProgressBar
# #' 
# Optimal_pars_jointefftox <- function(theta01, theta02, theta03, theta04, 
#                                    theta11, theta12, theta13, theta14, 
#                                    n, nsim, t1e, t2e, method ="power", 
#                                    lambda1, lambda2, grid1, 
#                                    gamma1, gamma2, grid2, 
#                                    eta1, eta2, grid3, ...) {
#   .Defunct("search_optimal_pars_efftox", msg = "use search_optimal_pars_efftox() instead")
#   
#   # Define H0 and H1
#   H0 <- c(theta01, theta02, theta03, theta04)
#   H1 <- c(theta11, theta12, theta13, theta14)
#   nsum <- sum(n)
#   
#   # Initialize result a data frame
#   results <- data.frame(
#     reject = numeric(0),
#     power = numeric(0),
#     H0 = character(0),
#     H1 = character(0),
#     lambda = numeric(0),
#     gamma = numeric(0),
#     eta = numeric(0)
#   )
#   
#   
#   # Calculate total iterations
#   total_iterations <- (grid1 + 1) * (grid2 + 1) * (if (method == "power") (grid3 + 1) else 1)
#   pb <- utils::txtProgressBar(min = 0, max = total_iterations, style = 3)
#   iteration <- 0
#   
#   # Iterate over grid values
#   for (i in 1:(grid1 + 1)) {
#     for (j in 1:(grid2 + 1)) {
#       if (method == "power") {
#         for (et in 1:(grid3 + 1)) {
#           # Calculate lambda, gamma, and eta values for this iteration
#           l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#           g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#           e <- if (grid3 == 1) eta1 else ((eta2 - eta1) / grid3) * (et-1) + eta1
#           
#           # Call the compute_power_coprimary function
#           result <- compute_power_jointefftox(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=e, method=method)
#           
#           # Append the results
#           results <- rbind(results, result)
#           
#           # Update progress bar
#           iteration <- iteration + 1
#           utils::setTxtProgressBar(pb, iteration)
#         }
#       } else {
#         # Calculate lambda and gamma values for this iteration
#         l <- if (grid1 == 1) lambda1 else ((lambda2 - lambda1) / grid1) * (i-1) + lambda1
#         g <- if (grid2 == 1) gamma1 else ((gamma2 - gamma1) / grid2) * (j-1) + gamma1
#         
#         # Call the compute_power_coprimary function without eta
#         result <- compute_power_jointefftox(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=NULL, method=method)
#         
#         # Append the results
#         results <- rbind(results, result)
#         
#         # Update progress bar
#         iteration <- iteration + 1
#         utils::setTxtProgressBar(pb, iteration)
#       }
#     }
#   }
#   close(pb)
#   
#   # Filter results by t1e
#   results_filtered <- results %>% dplyr::filter(!is.na(reject_mean) & reject_mean  < t1e)
#   
#   # Sort results by power
#   results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(power_mean))
#   
#   return(results_sorted)
# }
# 
