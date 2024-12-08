#' Operating characteristics for for Binary Endpoint
#' @param p Response rate
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param nsim number of simulation
#' @param fb vector of futility boundary at each interim
#' @param sb vector of superiority boundary at each interim
#' @param seed  for responsibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @export

Oc_binary <- function(p, n, nsim, fb, sb, seed = 12345) {
  set.seed(seed)

  num_interims <- length(n)  # Determine the number of interim analyses

  # Initialize empty data frame with appropriate columns
  resp <- data.frame(matrix(ncol = num_interims + 1, nrow = 0))

  # Name the columns accordingly (e.g., obs, x1, x2, x3, ...)
  colnames(resp) <- c("obs", paste0("x", 1:num_interims))

  # Simulate data
  for (obs in 1:nsim) {
    # Generate binomial random variables for each interim analysis
    x <- rbinom(num_interims, n, p)
    resp[obs, ] <- c(obs, x)  # Assign values directly to the appropriate row
  }

  # Convert columns to numeric to avoid type issues
  for (i in 2:ncol(resp)) {
    resp[[i]] <- as.numeric(resp[[i]])
  }

  # Calculate cumulative sums for each stage
  resp[[paste0("sumx", 1)]] <- resp[, paste0("x", 1)]

  for (i in 2:num_interims) {
    # resp[[paste0("sumx", 1:i)]] <- rowSums(resp[, paste0("x", 1:i)])
    resp[[paste0("sumx", i)]] <- rowSums(resp[, paste0("x", 1:i)])

  }

  # Initialize early stopping criteria columns
  resp$earlystopfuti <- 0
  resp$earlystopsupe <- 0

  # Determine early stopping based on futility and superiority boundaries
  for (i in 1:(num_interims-1)) {
    if (i == 1) {
      # Early stopping at the first stage
      resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("x", i)]] <= fb[i])
      resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("x", i)]] >= sb[i])
    } else {
      # Early stopping at subsequent stages
      resp$earlystopfuti <- resp$earlystopfuti | (resp[[paste0("sumx", i)]] <= fb[i])
      resp$earlystopsupe <- resp$earlystopsupe | (resp[[paste0("sumx", i)]] >= sb[i])
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
