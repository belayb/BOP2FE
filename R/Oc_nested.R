#' Operating characteristics for for nested Endpoint
#' @param p1 Response rate (CR: Complete remission)
#' @param p2 Response rate (PR: Partial remission)
#' @param p3 Response rate (1-(CR+PR))
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
#' @param fb vector of futility boundary at each interim analysis specified 
#' in the following order: c(CR_1, CR/PR_1, CR_2, CR/PR_2, CR_3, CR/PR_3, ...)
#' @param sb vector of superiority boundary at each interim analysis specified 
#' in the following order: c(CR_1, CR/PR_1, CR_2, CR/PR_2, CR_3, CR/PR_3, ...)
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
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

Oc_nested <- function(p1, p2, p3, n, nsim, fb, sb, seed = 12345) {
  set.seed(seed)
  
  # number of interim analysis 
  n_stage <- length(n)
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
  
  # Assign column names for the response data frame
  colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y"), i, 1:3))))
  
  # Calculate CRIA and CRPRIA
  for (i in 1:n_stage) {
    Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
    Y_cols_i12 <- paste0("Y", rep(1:i, each = 2), rep(c("1", "2"), i))  # Collect all Y11, Y12, Y21, Y22,... Yi1, Yi2 columns
    
    # Calculate cumulative CRIA and CRPRIA values for each stage 
    resp[[paste0("CRIA", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i1)))
    resp[[paste0("CRPRIA", i)]] <- rowSums(dplyr::select(resp, dplyr::all_of(Y_cols_i12)))
  }
  
  # Add early stopping criteria for futility and superiority
  resp <- resp %>%
    dplyr::mutate(
      earlystopfuti = 0,  # Initialize the earlystopfuti column
      earlystopsupe = 0   # Initialize the earlystopsupe column
    )
  
  # apply early stopping conditions for all stages
  # for (i in 1:(n_stage - 1)) {
  #   resp <- resp %>%
  #     dplyr::mutate(
  #       earlystopfuti = dplyr::case_when(
  #         earlystopfuti == 1 ~ 1,  # If futility was already triggered, keep it
  #         !!rlang::sym(paste0("CRIA", i)) <= fb[2 * (i - 1) + 1] & !!rlang::sym(paste0("CRPRIA", i)) <= fb[2 * (i - 1) + 2] ~ 1,
  #         TRUE ~ earlystopfuti
  #       ),
  #       earlystopsupe = dplyr::case_when(
  #         earlystopsupe == 1 ~ 1,  # If superiority was already triggered, keep it
  #         !!rlang::sym(paste0("CRIA", i)) >= sb[2 * (i - 1) + 1] | !!rlang::sym(paste0("CRPRIA", i)) >= sb[2 * (i - 1) + 2] ~ 1,
  #         TRUE ~ earlystopsupe
  #       )
  #     )
  # }
  # 
  for (i in 1:(n_stage - 1)) {
      resp <- resp %>%
        dplyr::mutate(
          earlystopfuti = dplyr::case_when(
            earlystopfuti == 1 ~ 1,  # if futility was already triggered, keep it
            (earlystopsupe == 0)&(!!rlang::sym(paste0("CRIA", i)) <= fb[2 * (i - 1) + 1] & 
              !!rlang::sym(paste0("CRPRIA", i)) <= fb[2 * (i - 1) + 2]) ~ 1,
            TRUE ~ earlystopfuti
          )
        )

      resp <- resp %>%
        dplyr::mutate(
          earlystopsupe = dplyr::case_when(
            earlystopsupe == 1 ~ 1,  # if superiority was already triggered, keep it
            (earlystopfuti == 0)&(!!rlang::sym(paste0("CRIA", i)) >= sb[2 * (i - 1) + 1] | 
              !!rlang::sym(paste0("CRPRIA", i)) >= sb[2 * (i - 1) + 2]) ~ 1,
            TRUE ~ earlystopsupe
          )
        )
    }
  
  
  # Add sample size
  resp$ss <- sum(n) # initialize sample size
  
  # Sample size calculation
  for (i in 1:n_stage) {
    idx1 <- 2 * (i - 1) + 1
    idx2 <- 2 * (i - 1) + 2
    # Update sample size
    resp <- resp %>%
      dplyr::mutate(
        cond_1 = !!rlang::sym(paste0("CRIA", i)) <= fb[idx1] & !!rlang::sym(paste0("CRPRIA", i)) <= fb[idx2],
        cond_2 = !!rlang::sym(paste0("CRIA", i)) >= sb[idx1] | !!rlang::sym(paste0("CRPRIA", i)) >= sb[idx2],
        condition_met = cond_1 | cond_2,
        ss = ifelse(condition_met & ss == sum(n), sum(n[1:i]), ss)
      ) %>%
      dplyr::select(-cond_1, -cond_2, -condition_met)  # Remove temporary columns
  }
  
  # Reject null hypothesis
  idx1 <- 2 * (n_stage - 1) + 1
  idx2 <- 2 * (n_stage - 1) + 2
  
  resp <- resp %>%
    dplyr::mutate(rejectnull = dplyr::case_when(
      earlystopsupe == 1 | (earlystopfuti == 0 &
                              (!!rlang::sym(paste0("CRIA", n_stage)) >= sb[idx1] | !!rlang::sym(paste0("CRPRIA", n_stage)) >= sb[idx2])) ~ 1,
      TRUE ~ 0
    ))


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

