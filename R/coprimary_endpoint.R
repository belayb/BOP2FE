#' Boundary values for co-primary Endpoint
#' @param H0 Response rate under the null (Response - PFS6, Response - no PFS6, No response - PFS6, No response - No PFS6)
#' @param a alpha values for the beta prior (i.e. usually set to the null response rate)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis.
#' The value at index \code{i} indicates the number of new patients added at interim analysis \code{i}. 
#' The total sample size at interim \code{i} is the cumulative sum of the values in \code{n} up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as \code{n = c(10, 5, 5, 10)}, where:
#' \itemize{
#'   \item 10 is the number of patients enrolled at interim 1,
#'   \item 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'   \item 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'   \item 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' }
#' @param lambda A vector of values for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A vector of values for parameter `gamma` of the cut-off probability for futility 
#' @param eta A vector of values for parameter `eta` of the cut-off probability for efficacy 
#' @param method type of function to be used for the cut off probability for superiority. The default is "power" type function. method=OF is an alternative for "O'Brien-Fleming" 
#' @param seed seed number
#' @importFrom stats pbeta rbinom
#' @export
#' 
#' @returns A list with the first element corresponding to futility and the second for efficacy boundaries
#' @examples
#' \dontrun{
#' H0=c(0.05,0.05, 0.15, 0.75)
#' a <- H0
#' seed <- 123
#' n <- c(10, 5, 5, 5, 5, 5, 5)
#' method <- "power"
#' test1<- get_boundary_coprimary(H0=H0, a=a, n =n,
#'                                   lambda = seq(0, 1, l = 11),
#'                                   gamma  = seq(0, 1, l = 11),
#'                                   eta    = seq(0, 3, l = 31),
#'                                   method = "power")
#' }
#' 
get_boundary_coprimary <- function(H0, a, n, lambda, gamma, eta = NULL, method = "power", seed = NULL) {
  set.seed(seed)
  
  # beta function 
  calculate_postp <- function(H0, beta_a, beta_b) {
    1 - pbeta(H0, beta_a, beta_b)
  }
  
  # Initializing parameters
  cum_n <- cumsum(n)
  
  b1 <- c(1, 1, 0, 0)
  b2 <- c(1, 0, 1, 0)
  phi1 <- sum(H0 * b1)
  phi2 <- sum(H0 * b2)
  
  
  cf_cs_val <- get_cf_cs_values(n, lambda = lambda, gamma = gamma, eta = eta, method = method)
  cf_values <- cf_cs_val[["cf_values"]]
  cs_values <- cf_cs_val[["cs_values"]]
  
  # Function to generate nested endpoint sequences
  generate_nested_sequences <- function(n, a, b1, b2, phi1, phi2) {
    data_list <- lapply(seq_along(n), function(i) {
      data <- expand.grid(Y1 = 0:n[i], Y2 = 0:n[i], Y3 = 0:n[i])
      data <- subset(data, data$Y1 + data$Y2 + data$Y3 <= n[i])
      data[["Y4"]] <- n[i] - data[["Y1"]] - data[["Y2"]] - data[["Y3"]]
      data[["beta_a_11"]] <- b1[1] * (a[1] + data[["Y1"]]) + b1[2] * (a[2] + data[["Y2"]]) + b1[3] * (a[3] + data[["Y3"]]) + b1[4] * (a[4] + data[['Y4']])
      data[["beta_b_11"]] <- (1 - b1[1]) * (a[1] + data[["Y1"]]) + (1 - b1[2]) * (a[2] + data[["Y2"]]) + (1 - b1[3]) * (a[3] + data[["Y3"]]) + (1 - b1[4]) * (a[4] + data[['Y4']])
      data[["beta_a_12"]] <- b2[1] * (a[1] + data[["Y1"]]) + b2[2] * (a[2] + data[["Y2"]]) + b2[3] * (a[3] + data[["Y3"]]) + b2[4] * (a[4] + data[['Y4']])
      data[["beta_b_12"]] <- (1 - b2[1]) * (a[1] + data[["Y1"]]) + (1 - b2[2]) * (a[2] + data[["Y2"]]) + (1 - b2[3]) * (a[3] + data[["Y3"]]) + (1 - b2[4]) * (a[4] + data[['Y4']])
      
      data[["postp11"]] <- mapply(calculate_postp, phi1, data[["beta_a_11"]], data[["beta_b_11"]])
      data[["postp12"]] <- mapply(calculate_postp, phi2, data[["beta_a_12"]], data[["beta_b_12"]])
      data[["Y12"]] <-  data[["Y1"]] +  data[["Y2"]]
      data[["Y13"]] <-  data[["Y1"]] +  data[["Y3"]]
      
      #data[!duplicated(data[, c('Y12', 'Y13','postp11','postp12')]), ]
      return(data)
    })
    return(data_list)
  }
  
  # Generate sequences for each interim analysis
  nested_sequences <- generate_nested_sequences(cum_n, a, b1, b2, phi1, phi2)
  
  
  combine_data <- function(nested_sequences, var) {
    # Remove duplicated rows based on Y12, Y13, postp11, postp12
    nested_sequences <- lapply(nested_sequences, function(data) {
      data[!duplicated(data[, c('Y12', 'Y13', 'postp11', 'postp12')]), ]
    })
    
    do.call(rbind, lapply(nested_sequences, function(data) {
      max_len <- length(nested_sequences[[(length(nested_sequences))]][[var]])
      vec <- data[[var]]
      length(vec) <- max_len
      return(vec)
    }))
  }
  
  Y12_range <- combine_data(nested_sequences, "Y12")
  Y13_range <- combine_data(nested_sequences, "Y13")
  postp11_range <- combine_data(nested_sequences, "postp11")
  postp12_range <- combine_data(nested_sequences, "postp12")
  
  
  # Maximum values of cnf1
  cnf_max1 <- do.call(
    rbind,
    lapply(seq(lambda), function(i) {
      t(
        sapply(seq(gamma), function(j) {
          apply(
            cbind(0, Y12_range * NA ^ (!((postp11_range < cf_values[i, , j])&(postp12_range < cf_values[i, , j])))), 
            1, 
            FUN = max, 
            na.rm = TRUE
          )
        })
      )
    })
  )
  
  # Maximum values of cnf2
  cnf_max2 <- do.call(
    rbind,
    lapply(seq(lambda), function(i) {
      t(
        sapply(seq(gamma), function(j) {
          apply(
            cbind(0, Y13_range * NA ^ (!((postp11_range < cf_values[i, , j])&(postp12_range < cf_values[i, , j])))), 
            1, 
            FUN = max, 
            na.rm = TRUE
          )
        })
      )
    })
  )
  
  # cnf result
  cnf <- data.frame(
    lambda = rep(lambda, each = length(gamma)),
    gamma  = rep(gamma, times = length(lambda)),
    cnf_max1,cnf_max2
  )
  colnames(cnf) <- c('lambda', 'gamma', paste0('f1', seq(n)),paste0('f2', seq(n)))
  
  # OF boundary
  if(method == "OF") {
    # Minimum values of cns1
    # Dimension: lambda * n
    cns_min1 <- t(
      sapply(seq(lambda), function(i) {
        apply(
          cbind(Y12_range * NA ^ (postp11_range < cs_values[i, ]), cum_n), 
          1, 
          FUN = min, 
          na.rm = TRUE
        )
      })
    )
    # Minimum values of cns2
    cns_min2 <- t(
      sapply(seq(lambda), function(i) {
        apply(
          cbind(Y13_range * NA ^ (postp12_range < cs_values[i, ]), cum_n), 
          1, 
          FUN = min, 
          na.rm = TRUE
        )
      })
    )
    # cns result
    cns <- data.frame(
      lambda = lambda,
      eta    = NA,
      cns_min1, 
      cns_min2
    )
    # Power boundary
  } else if(method == 'power') {
    # Minimum values of cns1
    # Dimension: (lambda * eta) * n
    cns_min1 <- do.call(
      rbind,
      lapply(seq(lambda), function(i) {
        t(
          sapply(seq(eta), function(j) {
            apply(
              cbind(Y12_range * NA ^ (postp11_range < cs_values[i, , j]), cum_n), 
              1, 
              FUN = min, 
              na.rm = TRUE
            )
          })
        )
      })
    )
    # Minimum values of cns2
    cns_min2 <- do.call(
      rbind,
      lapply(seq(lambda), function(i) {
        t(
          sapply(seq(eta), function(j) {
            apply(
              cbind(Y13_range * NA ^ (postp12_range < cs_values[i, , j]), cum_n), 
              1, 
              FUN = min, 
              na.rm = TRUE
            )
          })
        )
      })
    )
    # cns result
    cns <- data.frame(
      lambda = rep(lambda, each = length(eta)),
      eta    = rep(eta, times = length(lambda)),
      cns_min1, 
      cns_min2
    )
  }
  colnames(cns) <- c('lambda', 'eta', paste0('s1', seq(n)), paste0('s2', seq(n)))
  
  # Results of cnf and cns boundaries
  boundaries <- list(
    cnf= cnf,
    cns= cns
  )
  # Return results
  return(boundaries)
}


#' Operating characteristics for for coprimary Endpoint
#' @param p1 Response rate (Response - PFS6)
#' @param p2 Response rate (Response - no PFS6)
#' @param p3 Response rate (No Response - PFS6)
#' @param p4 Response rate (No Response - no PFS6)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis.
#' The value at index \code{i} indicates the number of new patients added at interim analysis \code{i}. 
#' The total sample size at interim \code{i} is the cumulative sum of the values in \code{n} up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as \code{n = c(10, 5, 5, 10)}, where:
#' \itemize{
#'   \item 10 is the number of patients enrolled at interim 1,
#'   \item 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'   \item 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'   \item 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' }
#' @param nsim number of simulation
#' @param fb vector/matrix of futility boundary at each interim analysis specified 
#' in the following order: c(resp_1,..., resp_length(n), PFS6_1, ..., PFS6_length(n))
#' @param sb vector/matrix  of superiority boundary at each interim analysis specified 
#' in the following order: c(resp_1,..., resp_length(n), PFS6_1, ..., PFS6_length(n))
#' @param seed for reproducibility
#' @importFrom stats rmultinom
#' @export
#' 
#' @returns A data frame with the following columns
#' \describe{
#' \item{lambda: }{lambda values for cut-off probability}
#' \item{gamma: }{gamma values for cut-off probability}
#' \item{eta: }{eta values for cut-off probability}
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
#' @examples
#' \dontrun{
#' H0=c(0.05,0.05, 0.15, 0.75)
#' a <- H0
#' seed <- 123
#' n <- c(10, 5, 5, 5, 5, 5, 5)
#' method <- "power"
#' boundary_tab<- get_boundary_coprimary(H0=H0, a=a, n =n,
#'                                lambda = seq(0, 1, l = 11),
#'                                gamma  = seq(0, 1, l = 11),
#'                                eta    = seq(0, 3, l = 31),
#'                                method = "power",
#'                                seed=seed)
#' test_oc<-get_oc_coprimary(
#'   p1 = 0.05,
#'   p2 = 0.05,
#'   p3 = 0.15,
#'   p4 = 0.75,
#'   n = c(10, 5, 5, 5, 5, 5, 5),
#'   nsim = 1000,
#'   fb = boundary_tab$cnf,
#'   sb = boundary_tab$cns,
#'   seed = seed
#' )
#' }
#'   
get_oc_coprimary <- function(p1, p2, p3, p4, n, nsim, fb, sb, seed = NULL) {
  set.seed(seed)
  #Generate data
  num_interims <- length(n)
  simulate_run <- function(i) {
    probs <- c(p1, p2, p3, p4)
    Y <- rmultinom(nsim, n[i], probs)
    t(Y)
  }
  
  resp <- do.call(cbind, lapply(1:num_interims, simulate_run))
  resp <- cbind(1:nsim, resp)
  resp <- as.data.frame(resp)
  # Assign column names for the response data frame
  colnames(resp) <- c("obs", unlist(lapply(1:num_interims, function(i) paste0(c("Y", "Y", "Y", "Y"), i, 1:4))))
  
  sum_oria <- t(apply(resp[,paste0("Y", 1:(num_interims), "1")], 1,cumsum)) + t(apply(resp[,paste0("Y", 1:(num_interims), "2")], 1,cumsum))
  sum_efs6ia <- t(apply(resp[,paste0("Y", 1:(num_interims), "1")], 1,cumsum)) + t(apply(resp[,paste0("Y", 1:(num_interims), "3")], 1,cumsum))
  colnames(sum_oria)<- paste0("ORIA", seq(num_interims))
  colnames(sum_efs6ia) <- paste0("PFS6IA", seq(num_interims))
  
  resp<- as.data.frame(cbind(resp, sum_oria, sum_efs6ia))
  
  # Determine if early stopping for superiority criteria are met
  # Unique combinations of f_{11},...,f_{1length(n)}, f_{21},...,f_{2length(n)} 
  uniq_fb <- fb[!duplicated(t((apply(fb[, c(paste0('f1', seq(n)),paste0('f2', seq(n)))], 1, sort)))),]
  # Unique combinations of s_{1},...,s_{length(n)}
  uniq_sb <- sb[!duplicated(t((apply(sb[, c(paste0('s1', seq(n)),paste0('s2', seq(n)))], 1, sort)))),]
  # Merge two data sets of fb and sb
  uniq_fb_and_sb = merge(uniq_fb, uniq_sb)
  
  # temp
  sb2 = uniq_fb_and_sb[, c('lambda', 'eta', paste0('s1', seq(num_interims - 1)), paste0('s2', seq(num_interims - 1)))]
  earlystopsupe_temp <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, c(paste0('ORIA', seq(num_interims - 1)), paste0('PFS6IA', seq(num_interims - 1)))]), 
        row(as.matrix(resp[, c(paste0('ORIA', seq(num_interims - 1)), paste0('PFS6IA', seq(num_interims - 1)))]))
      ), 
      function(x) t(x >= t(sb2[, -(1:2)]))
    )
  )
  
  earlystopsupe <- '+'(earlystopsupe_temp[, 1:(num_interims - 1)],
                       earlystopsupe_temp[, num_interims:(2*(num_interims-1))])
  
  fb2 = uniq_fb_and_sb[, c('lambda', 'gamma', paste0('f1', seq(num_interims - 1)), paste0('f2', seq(num_interims - 1)))]
  
  earlystopfuti_temp <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, c(paste0('ORIA', seq(num_interims - 1)), paste0('PFS6IA', seq(num_interims - 1)))]), 
        row(as.matrix(resp[, c(paste0('ORIA', seq(num_interims - 1)), paste0('PFS6IA', seq(num_interims - 1)))]))
      ), 
      function(x) t(x <= t(fb2[, -(1:2)]))
    )
  )
  
  earlystopfuti <- '*'(earlystopfuti_temp[, 1:(num_interims - 1)],
                       earlystopfuti_temp[, num_interims:(2*(num_interims-1))])
  
  # Get the time of futility and efficacy stopping 
  sup_time = do.call(
    function(...) pmin(..., na.rm = TRUE),
    as.data.frame(
      t(t(earlystopsupe>0) * seq(num_interims - 1)) * (NA ^ (1 - (earlystopsupe>0)))
    )
  )
  sup_time[is.na(sup_time)] <- num_interims
  
  futi_time = do.call(
    function(...) pmin(..., na.rm = TRUE),
    as.data.frame(
      t(t(earlystopfuti) * seq(num_interims - 1)) * (NA ^ (1 - earlystopfuti))
    )
  )
  futi_time[is.na(futi_time)] <- num_interims
  
  
  earlystopfuti <- matrix(as.double((rowSums(as.matrix(earlystopfuti)) > 0) & (futi_time < sup_time)), ncol = nsim)
  earlystopsupe <- matrix(as.double((rowSums(as.matrix(earlystopsupe)) > 0) & (futi_time > sup_time)), ncol = nsim)
  
  # Check if the stopping conditions are met for sample size calculation
  
  fb21 = uniq_fb_and_sb[, c('lambda', 'gamma', paste0('f1', seq(num_interims)), paste0('f2', seq(num_interims)))]
  conditionfuti_temp <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, c(paste0('ORIA', seq(num_interims)), paste0('PFS6IA', seq(num_interims)))]), 
        row(as.matrix(resp[, c(paste0('ORIA', seq(num_interims)), paste0('PFS6IA', seq(num_interims)))]))
      ), 
      function(x) t(x <= t(fb21[, -(1:2)]))
    )
  )
  conditionfuti <- '*'(conditionfuti_temp[, 1:num_interims],
                       conditionfuti_temp[, (num_interims+1):(2*num_interims)])
  
  sb21 = uniq_fb_and_sb[, c('lambda', 'eta', paste0('s1', seq(num_interims)), paste0('s2', seq(num_interims)))]
  conditionsupe_temp <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, c(paste0('ORIA', seq(num_interims)), paste0('PFS6IA', seq(num_interims)))]), 
        row(as.matrix(resp[, c(paste0('ORIA', seq(num_interims)), paste0('PFS6IA', seq(num_interims)))]))
      ), 
      function(x) t(x >= t(sb21[, -(1:2)]))
    )
  )
  
  conditionsupe <- '+'(conditionsupe_temp[, 1:num_interims],
                       conditionsupe_temp[, (num_interims+1):(2*num_interims)])
  
  condition_met = '+'(
    conditionfuti,
    conditionsupe>0
  )
  # Calculate the sample size at which stopping criteria are met
  ss <- do.call(
    function(...) pmin(..., na.rm = TRUE), 
    as.data.frame(
      t(cumsum(n) * t(condition_met)) * (NA ^ (1 - condition_met))
    )
  )
  ss[is.na(ss)] <- sum(n)
  ss <- matrix(ss, ncol = nsim)
  
  # Determine if the null hypothesis can be rejected
  rejectnull = matrix(
    '+'(
      earlystopsupe,
      '*'(
        (!earlystopfuti),
        t(
          outer(
            resp[[paste0('ORIA', num_interims)]], sb21[[paste0('s1', num_interims)]], '>=') |
            outer(
              resp[[paste0('PFS6IA', num_interims)]], sb21[[paste0('s2', num_interims)]], '>=')
        )
      )
    ) > 0,
    ncol = nsim
  )
  
  # summarize simulation 
  summary = merge(
    data.frame(
      lambda = uniq_fb_and_sb[['lambda']],
      gamma = uniq_fb_and_sb[['gamma']],
      eta = uniq_fb_and_sb[['eta']],
      earlystopfuti_mean = rowSums(earlystopfuti) / nsim, 
      earlystopsupe_mean = rowSums(earlystopsupe) / nsim, 
      ss_mean = rowSums(ss) / nsim, 
      rejectnull_mean = rowSums(rejectnull) / nsim, 
      earlystopfuti_sum = rowSums(earlystopfuti), 
      earlystopsupe_sum = rowSums(earlystopsupe), 
      ss_sum = rowSums(ss), 
      rejectnull_sum = rowSums(rejectnull)
    ),
    uniq_fb_and_sb
  )
  return(summary)
  
}


#' Computes both the boundary and corresponding operating characteristics for co primary endpoints  
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#'  - `H0[1]`: Response - PFS6,
#'  - `H0[2]`: Response - no PFS6,
#'  - `H0[3]`: No Response - PFS6,
#'  - `H0[4]`: No Response - no PFS6
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#'  - `H1[1]`: Response - PFS6,
#'  - `H1[2]`: Response - no PFS6,
#'  - `H1[3]`: No Response - PFS6,
#'  - `H1[4]`: No Response - no PFS6
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis.
#' The value at index \code{i} indicates the number of new patients added at interim analysis \code{i}. 
#' The total sample size at interim \code{i} is the cumulative sum of the values in \code{n} up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as \code{n = c(10, 5, 5, 10)}, where:
#' \itemize{
#'   \item 10 is the number of patients enrolled at interim 1,
#'   \item 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'   \item 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'   \item 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' }
#' @param nsim number of simulation. A value at least 1000 for better result
#' @param lambda A vector of values for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A vector of values for parameter `gamma` of the cut-off probability for futility 
#' @param eta A vector of values for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility 
#' @keywords internal
#' @export
#' 
#' @returns A data frame with the following columns
#' \describe{
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
#'   \item{lambda: }{lambda values for cut-off probability}
#'   \item{gamma: }{gamma values for cut-off probability}
#'   \item{eta: }{eta values for cut-off probability}} 
#' @examples
#' \dontrun{
#' oc_coprimary<-get_boundary_oc_coprimary(
#'   H0=c(0.15,0.30, 0.15, 0.40), 
#'   H1= c(0.18,0.42, 0.02, 0.38),
#'   n = c(10, 5, 5, 5, 5, 5, 5),
#'   nsim = 1000,
#'   lambda = seq(0, 1, l = 11),
#'   gamma  = seq(0, 1, l = 11),
#'   eta    = seq(0, 3, l = 31),
#'   method = "power",
#'   seed = 1
#' )
#' }
#' 
get_boundary_oc_coprimary <- function(
    H0, H1, n, nsim, lambda, gamma, eta = NULL, method = "power", seed = NULL  
) {
  # Generate boundary test results
  boundary_tab <- get_boundary_coprimary(
    H0 = H0, a = H0, n = n, lambda = lambda, gamma = gamma, eta = eta, method = method, seed=seed
  )
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- get_oc_coprimary(
    p1 = H0[1], p2=H0[2], p3=H0[3], p4=H0[4], n = n, nsim = nsim, fb = boundary_tab[['cnf']], sb = boundary_tab[['cns']], seed=seed
  )
  null_oc2 <- null_oc[, colnames(null_oc) %in% c(
    'lambda', 'gamma', 'eta', 'earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean',
    paste0('f1', seq(n)), paste0('f2', seq(n)), paste0('s1', seq(n)),paste0('s2', seq(n))
  )]
  colnames(null_oc2) <- c(
    'lambda', 'gamma', 'eta', paste0(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), '_h0'),
    paste0('f1', seq(n)), paste0('f2', seq(n)), paste0('s1', seq(n)),paste0('s2', seq(n))
  )
  alt_oc <- get_oc_coprimary(
    p1 = H1[1], p2=H1[2], p3=H1[3], p4=H1[4], n = n, nsim = nsim, fb = boundary_tab[['cnf']], sb = boundary_tab[['cns']], seed=seed
  )
  alt_oc2 <- alt_oc[, colnames(alt_oc) %in% c(
    'lambda', 'gamma', 'eta', 'earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean',
    paste0('f1', seq(n)), paste0('f2', seq(n)), paste0('s1', seq(n)),paste0('s2', seq(n))
  )]
  colnames(alt_oc2) <- c(
    'lambda', 'gamma', 'eta', paste0(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), '_h1'),
    paste0('f1', seq(n)), paste0('f2', seq(n)), paste0('s1', seq(n)),paste0('s2', seq(n))
  )
  # Combine results into a single data.frame
  all_res_df <- merge(null_oc2, alt_oc2)
  all_res_df <- all_res_df[, c(
    paste0('f1', seq(n)), paste0('f2', seq(n)), paste0('s1', seq(n)),paste0('s2', seq(n)),
    paste0(rep(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), 2), rep(c('_h0', '_h1'), each = 4)),
    'lambda', 'gamma', 'eta'
  )]
  colnames(all_res_df) <- ifelse(
    colnames(all_res_df) %in% c(paste0('f1', seq(n)),paste0('f2', seq(n))),
    c(paste0('fut_boundary_OR', seq(n)), paste0('fut_boundary_PFS6', seq(n))),
    colnames(all_res_df)
  )
  colnames(all_res_df) <- ifelse(
    colnames(all_res_df) %in% c(paste0('s1', seq(n)),paste0('s2', seq(n))),
    c(paste0('sup_boundary_OR', seq(n)),paste0('sup_boundary_PFS6', seq(n))),
    colnames(all_res_df)
  )
  return(all_res_df)
}

#' Search optimal parameters for co primary endpoint
#' 
#' `search_optimal_pars_coprimary()` is a helper function and calls `get_boundary_oc_coprimary()` to obtain the
#' optimal parameter combinations with type I error less than t1e and sorted according to power. 
#'   
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H0[1]`: Response - PFS6,
#' - `H0[2]`: Response - no PFS6,
#' - `H0[3]`: No Response - PFS6,
#' - `H0[4]`: No Response - no PFS6
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#' - `H1[1]`: Response - PFS6,
#' - `H1[2]`: Response - no PFS6,
#' - `H1[3]`: No Response - PFS6,
#' - `H1[4]`: No Response - no PFS6
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis.
#' The value at index \code{i} indicates the number of new patients added at interim analysis \code{i}. 
#' The total sample size at interim \code{i} is the cumulative sum of the values in \code{n} up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as \code{n = c(10, 5, 5, 10)}, where:
#' \itemize{
#'   \item 10 is the number of patients enrolled at interim 1,
#'   \item 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'   \item 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'   \item 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' }
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
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @returns A data frame with the following columns
#' \describe{
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
#'   \item{lambda: }{lambda values for cut-off probability}
#'   \item{gamma: }{gamma values for cut-off probability}
#'   \item{eta: }{eta values for cut-off probability}} 
#'
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#'test_comprimary <- search_optimal_pars_coprimary(
#'  H0=c(0.05,0.05, 0.15, 0.75),
#'  H1= c(0.15,0.15, 0.20, 0.50),
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000,
#'  t1e = 0.1,
#'  method = "power",
#'  lambda1 = 0,
#'  lambda2 = 1,
#'  grid1 = 11,
#'  gamma1 = 0,
#'  gamma2 = 1,
#'  grid2 = 11,
#'  eta1 = 0,
#'  eta2 = 3,
#'  grid3 = 31,
#'  seed = 123
#')
#'}
#' 
search_optimal_pars_coprimary <- function(
    H0, H1, n, nsim, t1e = NULL, method = "power", 
    lambda1, lambda2, grid1, 
    gamma1, gamma2, grid2, 
    eta1 = NULL, eta2 = NULL, grid3 = NULL, 
    seed = NULL  
) {
  # Set lambda, gamma and eta
  lambda = seq(lambda1, lambda2, l = grid1)
  gamma  = seq(gamma1,  gamma2,  l = grid2)
  if(method == "power") {
    eta    = seq(eta1,    eta2,    l = grid3)
  } else if(method == "OF") {
    eta = NA
  }
  
  # Initialize result a data frame
  all_results = get_boundary_oc_coprimary(H0, H1, n, nsim, lambda, gamma, eta, method, seed)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- all_results[!is.na(all_results$rejectnull_mean_h0) & all_results$rejectnull_mean_h0 < t1e, ]  
  # Sort results by power
  results_sorted <- results_filtered[order(-results_filtered$rejectnull_mean_h1), ]  
  return(results_sorted)
}
