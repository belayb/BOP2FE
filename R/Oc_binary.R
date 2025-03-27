#' Operating characteristics for for Binary Endpoint
#' @param p Response rate
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation
#' @param fb vector of futility boundary at each interim
#' @param sb vector of superiority boundary at each interim
#' @param seed  for responsibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{earlystopfuti_mean: }{Average number of early stopping due to futility}
#'  \item{earlystopsupe_mean: }{Average number of early stopping for futility due to efficacy}
#'   \item{ss_mean: }{Average sample size"} 
#'   \item{rejectnull_mean: }{"Average number of hypothesis rejection at the final analysis (aka Type-I error 
#'   if the response rate is the null rate or Power if the response rate is the alternative rate.} 
#'   \item{earlystopfuti_sum: }{Total number of early stopping due to futility} 
#'   \item{earlystopsupe_sum: }{Total number of early early stopping due to efficacy} 
#'   \item{ss_sum: }{Sum of sample sizes across simulation"} 
#'   \item{rejectnull_sum: }{Total number of hypothesis rejection at the final analysis} 
#'   } 
#'   

Oc_binary <- function(p, n, nsim, fb, sb, seed = NULL) {
  set.seed(seed)

  num_interims <- length(n)  # Determine the number of interim analyses

  # # Initialize empty data frame with appropriate columns
  # resp <- data.frame(matrix(ncol = num_interims + 1, nrow = 0))
  # 
  # # Name the columns accordingly (e.g., obs, x1, x2, x3, ...)
  # colnames(resp) <- c("obs", paste0("x", 1:num_interims))
  # 
  # # Simulate data
  # for (obs in 1:nsim) {
  #   # Generate binomial random variables for each interim analysis
  #   x <- rbinom(num_interims, n, p)
  #   resp[obs, ] <- c(obs, x)  # Assign values directly to the appropriate row
  # }
  # 
  # # Convert columns to numeric to avoid type issues
  # for (i in 2:ncol(resp)) {
  #   resp[[i]] <- as.numeric(resp[[i]])
  # }
  # 
  # resp<- simulate_data_binary(p, num_interims, n, nsim)
  # resp<- resp[,1:(num_interims+1)]
  # colnames(resp) <- c("obs", paste0("x", 1:num_interims))
  # 
  # 
  # # Calculate cumulative sums for each stage
  # resp[[paste0("sumx", 1)]] <- resp[, paste0("x", 1)]
  # 
  # for (i in 2:num_interims) {
  #   # resp[[paste0("sumx", 1:i)]] <- rowSums(resp[, paste0("x", 1:i)])
  #   resp[[paste0("sumx", i)]] <- rowSums(resp[, paste0("x", 1:i)])
  # 
  # }
  simulate_data_binary2 <- function(p, n_stage, n, nsim, seed = 23456) {
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
    #colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("Y", i))))
    colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0("x", i))))
    
    # Calculate cumulative sums
    for (i in 1:n_stage) {
      #Y_cols <- paste0("Y", 1:i)
      x_cols <- paste0("x", 1:i)
      
      # Calculate cumulative sum values for each stage
      #resp[[paste0("tsum", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols)))
      resp[[paste0("sumx", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(x_cols)))
      
    }
    return(resp)
  }
  # Simulate data 
  resp <- simulate_data_binary2(p=p, n_stage=num_interims, n=n, nsim=nsim, seed = seed)

  # Initialize early stopping criteria columns
  resp$earlystopfuti <- 0
  resp$earlystopsupe <- 0

  # Determine early stopping based on futility and superiority boundaries
  # for (i in 1:(num_interims-1)) {
  #   if (i == 1) {
  #     # Early stopping at the first stage
  #     resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("x", i)]] <= fb[i])
  #     resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("x", i)]] >= sb[i])
  #   } else {
  #     # Early stopping at subsequent stages
  #     resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("sumx", i)]] <= fb[i])
  #     resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("sumx", i)]] >= sb[i])
  #   }
  # }
  # for (i in 1:(num_interims-1)) {
  #   if (i == 1) {
  #     # Early stopping at the first stage
  #     resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("x", i)]] <= fb[i])
  #     resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("x", i)]] >= sb[i])
  #   } else {
  #     # Early stopping at subsequent stages
  #     if (!any(resp$earlystopsupe)) {
  #       resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("sumx", i)]] <= fb[i])
  #     }
  #     if (!any(resp$earlystopfuti)) {
  #       resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("sumx", i)]] >= sb[i])
  #     }
  #   }
  # }

  for (i in 1:(num_interims - 1)) {
    if (i == 1) {
      # Early stopping at the first stage
      resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("x", i)]] <= fb[i])
      resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("x", i)]] >= sb[i])
    } else {
    resp <- resp %>%
      dplyr::mutate(
        earlystopfuti = dplyr::case_when(
          (earlystopsupe == 0)&(!!rlang::sym(paste0("sumx", i)) <= fb[i] ) ~ 1,
          TRUE ~ earlystopfuti
        )
      )
    
    
    resp <- resp %>%
      dplyr::mutate(
        earlystopsupe = dplyr::case_when(
          earlystopsupe == 1 ~ 1,  # if superiority was already triggered, keep it
          (earlystopfuti == 0)&(!!rlang::sym(paste0("sumx", i)) >= sb[i] ) ~ 1,
          TRUE ~ earlystopsupe
        )
      )
  }
  }
  # Calculate sample size based on stopping criteria
  resp$ss <- sum(n)

  for (i in seq_along(n)) {
    # Update ss conditionally, ensuring it updates to the lowest value matching the first condition
    condition_met <- resp[[paste0("sumx", i)]] <= fb[i] | resp[[paste0("sumx", i)]] >= sb[i]
    resp$ss <- ifelse(condition_met & resp$ss == sum(n), sum(n[1:i]), resp$ss)
  }

  # Determine rejection of null hypothesis
  resp$rejectnull <- resp$earlystopsupe | (!resp$earlystopfuti & resp[[paste0("sumx", num_interims)]] >= sb[num_interims])

  # Summary statistics
  summary <- data.frame(
    earlystopfuti_mean = mean(resp$earlystopfuti),
    earlystopsupe_mean = mean(resp$earlystopsupe),
    ss_mean = mean(resp$ss),
    rejectnull_mean = mean(resp$rejectnull),
    earlystopfuti_sum = sum(resp$earlystopfuti),
    earlystopsupe_sum = sum(resp$earlystopsupe),
    ss_sum = sum(resp$ss),
    rejectnull_sum = sum(resp$rejectnull)
  )

  return(summary)
}
