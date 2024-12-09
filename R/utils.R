#' Simulate nested endpoint data
#'
#' This function simulates nested data across multiple stages based on specified probabilities and sample sizes.
#'
#' @param p1 Probability of success for the first outcome.
#' @param p2 Probability of success for the second outcome.
#' @param p3 Probability of success for the third outcome.
#' @param n_stage Number of interim analysis in the simulation.
#' @param n A vector of sample sizes for each interim analysis.
#' @param nsim Number of simulation runs to perform.
#' @param seed An integer seed for reproducibility of the random number generation.
#'
#' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
#'
#' @examples
#' # Example usage of the simulate_data_nested function
#' result <- simulate_data_nested(p1 = 0.5, p2 = 0.3, p3 = 0.2,
#'  n_stage = 3, n = c(100, 150, 200), nsim = 10)
#' print(result)
#' @importFrom magrittr %>%
#' @export
simulate_data_nested <- function(p1, p2, p3, n_stage, n, nsim, seed = 23456) {
  set.seed(seed)
  # Initialize empty data frame
  resp <- data.frame()
  # Simulate data for each simulation run
  for (obs in 1:nsim) {
    Y_stage <- list()  # List to hold results for each stage
    
    # Simulate data for each stage
    for (i in 1:n_stage) {
      Y1 <- rbinom(1, n[i], p1)
      Y2 <- rbinom(1, n[i] - Y1, p2 / (p2 + p3))
      Y3 <- n[i] - Y1 - Y2
      
      # Store the values for this stage
      Y_stage[[i]] <- c(Y1, Y2, Y3)
    }
    
    # Combine all stages into a single data frame row
    row_data <- unlist(Y_stage)
    resp <- rbind(resp, c(obs, row_data))
  }
  
  # Assign dynamic column names for the response data frame
  colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y"), i, 1:3))))
  
  # Calculate cumulative sums
  for (i in 1:n_stage) {
    Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
    Y_cols_i2 <- paste0("Y", 1:i, "2")  # Collect all Y12, Y22,... Yi2 columns
    Y_cols_i3 <- paste0("Y", 1:i, "3")  # Collect all Y13, Y23,... Yi3 columns
    
    # Calculate cumulative sum values for each stage
    resp[[paste0("tsum", i, "1")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i1)))
    resp[[paste0("tsum", i, "2")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i2)))
    resp[[paste0("tsum", i, "3")]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i3)))
  }
  return(resp)
}


#' Simulate binary endpoint data
#'
#' This function simulates binary data across multiple stages based on specified probabilities and sample sizes.
#'
#' @param p Probability of success.
#' @param n_stage Number of interim analysis in the simulation.
#' @param n A vector of sample sizes for each interim analysis.
#' @param nsim Number of simulation runs to perform.
#' @param seed An integer seed for reproducibility of the random number generation.
#'
#' @return A data frame containing the simulated data for each observation, including cumulative sums for each stage.
#'
#' @examples
#' # Example usage of the simulate_data_nested function
#' result <- simulate_data_binary(0.2, 4, n=c(4,2,2,2), nsim=5)
#' print(result)
#' @importFrom magrittr %>%
#' @export
simulate_data_binary <- function(p, n_stage, n, nsim, seed = 23456) {
  set.seed(seed)
  # Initialize empty data frame
  resp <- data.frame()
  # Simulate data for each simulation run
  for (obs in 1:nsim) {
    Y_stage <- list()  # List to hold results for each stage
    for (i in 1:n_stage) {
      Y_stage[[i]] <- rbinom(1, n[i], p)
    }
    # Combine all stages into a single data frame row
    row_data <- unlist(Y_stage)
    resp <- rbind(resp, c(obs, row_data))
  }
  # Assign dynamic column names for the response data frame
  colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("Y", i))))
  
  # Calculate cumulative sums
  for (i in 1:n_stage) {
    Y_cols <- paste0("Y", 1:i)
    # Calculate cumulative sum values for each stage
    resp[[paste0("tsum", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols)))
  }
  return(resp)
}


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
#'               Options are "power" (default) or "obrien-fleming".
#'
#' @return A list containing two elements:
#' \item{cf_values}{A numeric vector of cutoff values for futility stopping.}
#' \item{cs_values}{A numeric vector of corresponding probability cutoffs.}
#'
#' @examples
#' # Example usage of the get_cf_cs_values function
#' result <- get_cf_cs_values(n = c(100, 150, 200), lambda = 0.05,
#' gamma = 1, eta = 0.5, method = "power")
#' print(result)
#' @importFrom magrittr %>%
#' @importFrom stats qnorm pnorm
#' @export
get_cf_cs_values<- function(n, lambda=NULL, gamma=NULL, eta= NULL, method = NULL){
  nsum<- sum(n)
  cf_values <- sapply(seq_along(n), function(i) {
    lambda * (sum(n[1:i]) / nsum)^gamma
  })

  if(method == "power"){
    cs_values <- sapply(seq_along(n), function(i) {
      1- (1-lambda)  * (sum(n[1:i]) / nsum)^eta
    })
  }
  else{
    cs_values <- sapply(seq_along(n), function(i) {
      2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(sum(n[1:i]) / nsum)) - 1
    })
  }
  return(list(cf_values=cf_values, cs_values=cs_values))
}

#' Compute Posterior Probability of Observing a Given Response Under the Null
#'
#' This function computes the posterior probability of observing a given response under the null hypothesis
#' across multiple stages of a trial.
#'
#' @param xdata A data frame containing the data for which posterior probabilities are to be calculated.
#' @param H0 A numeric vector representing the null hypothesis probabilities.
#' @param H1 A numeric vector representing the alternative hypothesis probabilities.
#' @param n_stage An integer indicating the number of interim analysis in the trial.
#'
#' @return A data frame with additional columns for the posterior probabilities calculated for each stage.
#'
#' @examples
#' # Example usage of the calculate_posterior_nested function
#' #xdata <- data.frame(tsum11 = c(10, 20), tsum12 = c(5, 15), tsum13 = c(2, 8))
#' #h0 <- c(0.5, 0.3, 0.2)
#' #h1 <- c(0.6, 0.4, 0.3)
#' #result <- calculate_posterior_nested(xdata, h0, h1, n_stage = 3)
#' #print(result)
#' @importFrom magrittr %>%
#' @importFrom stats pbeta
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom purrr pmap_dbl
#' @export
#'
calculate_posterior_nested <- function(xdata, H0, H1, n_stage) {
  H0 <- as.numeric(H0)
  H1 <- as.numeric(H1)
  
  b1 <- c(1, 0, 0)
  b2 <- c(1, 1, 0)
  a <- H0
  
  phi1 <- sum(H0 * b1)
  phi2 <- sum(H0 * b2)
  
  post_probs <- function(x, phi, b, a) {
    beta_a <- sum(b * (a + x))
    beta_b <- sum((1 - b) * (a + x))
    1 - pbeta(phi, beta_a, beta_b)
  }
  
  phi_list <- list(phi1, phi2)
  b_list <- list(b1, b2)
  
  for (i in 1:n_stage) {
    tsum_cols <- paste0("tsum", i, 1:3)
    for (j in 1:2) {
      phi <- phi_list[[j]]
      b <- b_list[[j]]
      postp_col <- paste0("postp", i, j)
      
      xdata <- xdata %>%
        dplyr::mutate(!!postp_col := purrr::pmap_dbl(dplyr::select(., dplyr::all_of(tsum_cols)), ~ post_probs(c(..1, ..2, ..3), phi, b, a)))
    }
  }
  
  return(xdata)
}


#' function title - under construction
#' @param H0 the null hypothesis
#' @param H1 the alternative hypothesis
#' @param n vector of sample sizes
#' @param nsim number of simulation to run
#' @param lambda lambda value for
#' @param gamma gamma value for
#' @param eta eta value for
#' @param method method to use for cutoff probability 
#' @param seed for reproducability 
#' 
#add documentation
compute_power_nested <- function(H0, H1, n, nsim, lambda = NULL, gamma = NULL, eta = NULL, method = NULL, seed = NULL) {
  n_stage <- length(n)
  H0 <- H0
  H1 <- H1
  nsim <- nsim
  
  y_null <- simulate_data_nested(p1 = H0[1], p2 = H0[2], p3 = H0[3], n = n, n_stage = n_stage, nsim = nsim)
  y_alter <- simulate_data_nested(p1 = H1[1], p2 = H1[2], p3 = H1[3], n = n, n_stage = n_stage, nsim = nsim)
  
  cfcs <- get_cf_cs_values(n, lambda, gamma, eta, method)
  cf_values <- cfcs$cf_values
  cs_values <- cfcs$cs_values
  
  y_null <- y_null %>%
    calculate_posterior_nested(., H0, H1, n_stage) %>%
    dplyr::mutate(reject = dplyr::case_when(
      purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
      purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
      TRUE ~ 1
    ))
  
  y_alter <- y_alter %>%
    calculate_posterior_nested(., H0, H1, n_stage) %>%
    dplyr::mutate(power = dplyr::case_when(
      purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
      purrr::reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
      TRUE ~ 1
    ))
  
  # Summarizing results
  reject_mean <- mean(y_null$reject)
  power_mean <- mean(y_alter$power)
  
  result <- dplyr::tibble(
    H0 = paste(H0[1], H0[2], sep = ", "),
    H1 = paste(H1[1], H1[2], sep = ", "),
    Lambda = lambda,
    Gamma = gamma,
    eta = eta,
    reject_mean = reject_mean,
    power_mean = power_mean
  )
  
  return(result)
}


#' Optimal parameters for nested outcome
#'
#' This function runs a grid search to find the optimal parameters for nested or an ordinal outcome that leads to the highest power while controlling the type - I error at the specified level.
#'
#' @param CR0 Probability of success for the first outcome under the null hypothesis.
#' @param CRPR0 Probability of success for the first or second outcome under the null hypothesis.
#' @param CR1 Probability of success for the first outcome under the alternative hypothesis.
#' @param CRPR1 Probability of success for the first or second outcome under the alternative hypothesis.
#' @param n A vector of sample sizes for each interim analysis.
#' @param nsim Number of simulation runs to perform.
#' @param t1e Desired Type - I error rate
#' @param t2e Desired Type - II error rate
#' @param method method to use for cutoff
#' @param lambda1 starting value for lambda values to search
#' @param lambda2 ending value for lambda values to search
#' @param grid1 number of lambda values to consider between lambda1 and lambda2
#' @param gamma1 starting value for gamma values to search
#' @param gamma2 ending value for gamma values to search
#' @param grid2 number of gamma values to consider between gamma1 and gamma2
#' @param eta1 starting value for eta values to search
#' @param eta2 ending value for eta values to search
#' @param grid3 number of eta values to consider between eta1 and eta2
#'
#' @return A data frame
#'
#' @importFrom magrittr %>%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
Optimal_pars_nested <- function(CR0, CRPR0, CR1, CRPR1, n, nsim, t1e, t2e, method,
                                lambda1, lambda2, grid1, gamma1, gamma2, grid2, eta1, eta2, grid3) {
  # Define H0 and H1
  H0 <- c(CR0, CRPR0 - CR0, 1 - CRPR0)
  H1 <- c(CR1, CRPR1 - CR1, 1 - CRPR1)
  nsum <- sum(n)

  # Initialize result a data frame
  results <- data.frame(
    reject = numeric(0),
    power = numeric(0),
    H0 = character(0),
    H1 = character(0),
    lambda = numeric(0),
    gamma = numeric(0),
    eta = numeric(0)
  )
  

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
          result <- compute_power_nested(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=e, method=method)

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
        result <- compute_power_nested(H0=H0, H1=H1, n=n, nsim=nsim, lambda=l, gamma=g, eta=NULL, method=method)

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
  results_filtered <- results %>% dplyr::filter(!is.na(reject_mean) & reject_mean  < t1e)
  
  # Sort results by power
  results_sorted <- results_filtered %>% dplyr::arrange(dplyr::desc(power_mean))
  
  return(results_sorted)
}

#' function title- under construction
#' @param xdata a given data set
#' @param H0 null hypothesis
#' @param H1 alternative hypothesis
#' @param n_stage the number of interm analysis
#'

calculate_posterior_binary <- function(xdata, H0, H1, n_stage) {
  H0 <- as.numeric(H0)
  H1 <- as.numeric(H1)

  b1 <- c(1, 0, 0)
  b2 <- c(1, 1, 0)
  a <- H0

  phi1 <- sum(H0 * b1)
  phi2 <- sum(H0 * b2)

  post_probs <- function(x, phi, b, a) {
    beta_a <- sum(b * (a + x))
    beta_b <- sum((1 - b) * (a + x))
    1 - pbeta(phi, beta_a, beta_b)
  }

  for (i in 1:n_stage){
    tsum_cols <- paste0("tsum", i)
    xdata <- xdata%>%
      mutate(!!postp_col := pmap_dbl(select(., all_of(tsum_cols)), ~ post_probs(c(..1, ..2, ..3), phi, b, a)))
  }
  return(xdata)
}

#' function title- under construction
#' @param H0 null hypothesis
#' @param H1 alternative hypothesis
#' @param n sample size vector
#' @param nsim number of simulation
#' @param n_stage the number of interm analysis
#'
#'
compute_power_binary <- function(H0, H1, n, nsim, n_stage){

  y_null <- simulate_data_binary(p=H0, n=n, n_stage=n_stage , nsim = nsim)
  y_alter <- simulate_data_binary(p=H1, n=n, n_stage=n_stage , nsim = nsim)

  cfcs<-get_cf_cs_values(n, lambda, gamma, eta, method)
  cf_values<-cfcs$cf_values
  cs_values<-cfcs$cs_values


  y_null <- y_null%>%
    calculate_posterior_binary(., H0, H1, n_stage)%>%
    mutate(reject = case_when(
      reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
      reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
      TRUE ~ 1
    ))

  y_alter <- y_alter%>%
    calculate_posterior_binary(., H0, H1, n_stage)%>%
    mutate(power = case_when(
      reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) < cf_values[.y] & get(paste0("postp", .y, "2")) < cf_values[.y]), .init = FALSE) ~ 0,
      reduce(1:n_stage, ~ .x | (get(paste0("postp", .y, "1")) >= cs_values[.y] | get(paste0("postp", .y, "2")) >= cs_values[.y]), .init = FALSE) ~ 1,
      TRUE ~ 1
    ))


  # Summarizing results
  reject_mean <- mean(y_null$reject)
  power_mean <- mean(y_alter$power)

  result <- data.frame(
    H0 = H0,
    H1 = H1,
    Lambda = lambda,
    Gamma = gamma,
    eta = eta,
    reject_mean = reject_mean,
    power_mean = power_mean
  )

  return(result)

}
