#' Operating characteristics for for joint Endpoint
#' @param p1 Response rate (toxicity - OR)
#' @param p2 Response rate (no toxicity - OR)
#' @param p3 Response rate (toxicity - no OR)
#' @param p4 Response rate (no toxicity - No OR)
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param nsim number of simulation
#' @param fb vector of futility boundary at each interim
#' @param sb vector of superiority boundary at each interim
#' @param n_stage number of interim analysis
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @export
#'
Oc_jointefftox <- function(p1, p2, p3, p4, n_stage, n, nsim, fb, sb, seed = 12345) {
  set.seed(seed)

  # Initialize empty data frame
  resp <- data.frame()

  # Simulate data for each simulation run
  for (obs in 1:nsim) {
    Y_stage <- list()  # List to hold results for each stage

    # Simulate data for each stage
    for (i in 1:n_stage) {
      Y1 <- rbinom(1, n[i], p1)
      Y2 <- rbinom(1, n[i] - Y1, p2 / (p2 + p3 + p4))
      Y3 <- rbinom(1, n[i] - Y1 - Y2, p3 / (p3 + p4))
      Y4 <- n[i] - Y1 - Y2 - Y3

      # Store the values for this stage
      Y_stage[[i]] <- c(Y1, Y2, Y3, Y4)
    }

    # Combine all stages into a single data frame row
    row_data <- unlist(Y_stage)
    resp <- rbind(resp, c(obs, row_data))
  }

  # Assign dynamic column names for the response data frame
  colnames(resp) <- c("obs", unlist(lapply(1:n_stage, function(i) paste0(c("Y", "Y", "Y", "Y"), i, 1:4))))

  # Calculate CRIA and CRPRIA dynamically
  for (i in 1:n_stage) {
    #Y_cols_i1 <- paste0("Y", 1:i, "1")  # Collect all Y11, Y21,... Yi1 columns
    Y_cols_i12 <- paste0("Y", rep(1:i, each=2), rep(c("1", "2"), i))  # Collect all Y11, Y12, Y21, Y22,... Yi1, Yi2 columns
    Y_cols_i13 <- paste0("Y", rep(1:i, each=2), rep(c("1", "3"), i)) # Collect all Y11, Y13, Y21, Y23,... Yi1, Yi3 columns
    # Calculate cumulative CRIA and CRPRIA values for each stage dynamically
    resp[[paste0("ORIA", i)]] <- rowSums(resp %>% select(all_of(Y_cols_i12)))
    resp[[paste0("TOXIA", i)]] <- rowSums(resp %>% select(all_of(Y_cols_i13)))
  }

  # Add dynamic early stopping criteria for futility and superiority
  resp <- resp %>%
    mutate(
      earlystopfuti = 0,  # Initialize the earlystopfuti column
      earlystopsupe = 0   # Initialize the earlystopsupe column
    )

  # Dynamically apply early stopping conditions for all stages
  for (i in 1:(n_stage-1)) {
    resp <- resp %>%
      mutate(
        earlystopfuti = case_when(
          earlystopfuti == 1 ~ 1,  # If futility was already triggered, keep it
          !!sym(paste0("ORIA", i)) <= fb[2 * (i - 1) + 1] | !!sym(paste0("TOXIA", i)) >= fb[2 * (i - 1) + 2] ~ 1,
          TRUE ~ earlystopfuti
        ),
        earlystopsupe = case_when(
          earlystopsupe == 1 ~ 1,  # If superiority was already triggered, keep it
          !!sym(paste0("ORIA", i)) >= sb[2 * (i - 1) + 1] & !!sym(paste0("TOXIA", i)) <= sb[2 * (i - 1) + 2] ~ 1,
          TRUE ~ earlystopsupe
        )
      )
  }

  # Add dynamic sample size
  resp$ss <- sum(n) # initialize sample size

  # Sample size calculation
  for (i in 1:n_stage) {
    idx1 <- 2 * (i - 1) + 1
    idx2 <- 2 * (i - 1) + 2
    # Update sample size
    resp <- resp %>%
      mutate(
        cond_1 = !!sym(paste0("ORIA", i)) <= fb[idx1] | !!sym(paste0("TOXIA", i)) >= fb[idx2],
        cond_2 = !!sym(paste0("ORIA", i)) >= sb[idx1] & !!sym(paste0("TOXIA", i)) <= sb[idx2],
        condition_met = cond_1 | cond_2,
        ss = ifelse(condition_met & ss == sum(n), sum(n[1:i]), ss)
      ) %>%
      select(-cond_1, -cond_2, -condition_met)  # Remove temporary columns

  }

  # Reject null hypothesis
  idx1 <- 2 * (n_stage - 1) + 1
  idx2 <- 2 * (n_stage - 1) + 2

  resp <- resp%>%
    mutate(rejectnull := case_when(earlystopsupe == 1 | (earlystopfuti == 0 &
                                                           (!!sym(paste0("ORIA", n_stage)) >= sb[idx1] & !!sym(paste0("TOXIA", n_stage)) <= sb[idx2])) ~ 1,
                                   TRUE ~ 0))


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

