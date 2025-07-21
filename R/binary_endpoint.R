
#' Boundary values for binary Endpoint
#' @param H0 Response rate under the null 
#' @param a1 alpha values for the beta prior (i.e. usually set to the null response rate)
#' @param b1 beta values for the beta prior (i.e. usually set to 1 - the null response rate)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#'  - 10 is the number of patients enrolled at interim 1,
#'  - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'  - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'  - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
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
#' H0 <- 0.2
#' a1 <- H0
#' b1 <- 1- a1
#' seed <- 123
#' n <- c(10, 5, 5, 5, 5, 5, 5)
#' method <- "power"
#' 
#' boundary_binary<- get_boundary_binary(H0=H0, a=a, n =n,
#'                      lambda = seq(0, 1, l = 101),
#'                      gamma  = seq(0, 1, l = 101),
#'                      eta    = seq(0, 3, l = 301),
#'                      method = method,
#'                     seed = seed)
#' }
#' 
get_boundary_binary <- function(H0, a1, b1, n, lambda, gamma, eta = NULL, method = "power", seed=NULL) {
  set.seed(seed)
  # Obtain cf and cs values
  cf_cs_val <- get_cf_cs_values(n, lambda = lambda, gamma = gamma, eta = eta, method = method)
  cf_values <- cf_cs_val[["cf_values"]]
  cs_values <- cf_cs_val[["cs_values"]]
  # Cumulative sample size for each stage (stage = 1,...,length(n))
  cum_n <- cumsum(n)
  # A matrix of patients
  # Dimension: length(n) * sum(n)
  x_range <- do.call(
    rbind, 
    lapply(
      cum_n, 
      function(i) c(seq(i), rep(NA, sum(n) - i))
    )
  )
  # Posterior probability
  postp <- 1 - pbeta(H0, a1 + x_range, b1 + cum_n - x_range)
  # Maximum values of cnf
  # Dimension: (lambda * gamma) * n
  cnf_max <- do.call(
    rbind,
    lapply(seq(lambda), function(i) {
      t(
        sapply(seq(gamma), function(j) {
          apply(
            cbind(0, x_range * NA ^ (postp >= cf_values[i, , j])), 
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
    cnf_max
  )
  colnames(cnf) <- c('lambda', 'gamma', paste0('f', seq(n)))
  # OF boundary
  if(method == "OF") {
    # Minimum values of cns
    # Dimension: lambda * n
    cns_min <- t(
      sapply(seq(lambda), function(i) {
        apply(
          cbind(x_range * NA ^ (postp <  cs_values[i, ]), cum_n), 
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
      cns_min
    )
    # Power boundary
  } else if(method == 'power') {
    # Minimum values of cns
    # Dimension: (lambda * eta) * n
    cns_min <- do.call(
      rbind,
      lapply(seq(lambda), function(i) {
        t(
          sapply(seq(eta), function(j) {
            apply(
              cbind(x_range * NA ^ (postp <  cs_values[i, , j]), cum_n), 
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
      cns_min
    )
  }
  colnames(cns) <- c('lambda', 'eta', paste0('s', seq(n)))
  # Results of cnf and cns boundaries
  boundaries <- list(
    cnf= cnf,
    cns= cns
  )
  # Return results
  return(boundaries)
}


#' Operating characteristics for binary Endpoint
#' @param p Response rate 
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#'  - 10 is the number of patients enrolled at interim 1,
#'  - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'  - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'  - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param nsim number of simulation
#' @param fb vector/matrix of futility boundary at each interim analysis specified 
#' in the following order: c(f_1,..., f_{length(n)})
#' @param sb vector/matrix  of superiority boundary at each interim analysis specified 
#' in the following order: c(s_1,..., s_{length(n)})
#' @param seed for reproducibility
#' @importFrom stats rbinom
#' @export
#' 
#' @returns A data frame with the following columns
#' \itemize{
#' \item{lambda: }{lambda values for cut-off probability}
#' \item{gamma: }{gamma valuesfor cut-off probability}
#' \item{eta: }{eta valuesfor cut-off probability}
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
#'  H0 <- 0.2
#'  a1 <- H0
#'  b1 <- 1-a1
#'  seed <- 123
#'  n <- c(10, 5, 5, 5, 5, 5, 5)
#'  method <- "power"
#'  boundary_tab<- get_boundary_binary(H0=H0, a1=a1, b1=b1, n =n,
#'                                lambda = seq(0, 1, l = 11),
#'                                gamma  = seq(0, 1, l = 11),
#'                                eta    = seq(0, 3, l = 31),
#'                                method = method,
#'                                seed=seed)
#' test_oc<-get_oc_binary(
#'   p = 0.2,
#'   n = c(10, 5, 5, 5, 5, 5, 5),
#'   nsim = 1000,
#'   fb = boundary_tab$cnf,
#'   sb = boundary_tab$cns,
#'   seed = seed
#' )
#' }
#'   
get_oc_binary <- function(p, n, nsim, fb, sb, seed = NULL) {
  set.seed(seed)
  num_interims <- length(n)
  # Simulate binary outcomes for each interim analysis
  # y_stage is a matrix where each row represents a simulation
  # and each column represents the number of successes at each interim
  Y_stage <- matrix(
    rbinom(num_interims * nsim, n, p), 
    ncol = num_interims, 
    byrow = TRUE
  )
  # Calculate cumulative sums of successes for each simulation
  sumx <- t(apply(Y_stage, 1, cumsum))
  resp <- as.data.frame(cbind(Y_stage, sumx))
  
  # Rename columns for clarity
  colnames(resp) <- paste0(
    rep(c('x', 'sumx'), each = num_interims), 
    seq(num_interims)
  )
  # Determine if early stopping for superiority criteria are met
  # Unique combinations of f_{1},...,f_{length(n)}
  uniq_fb <- fb[!duplicated(t((apply(fb[, paste0('f', seq(n))], 1, sort)))),]
  # Unique combinations of s_{1},...,s_{length(n)}
  uniq_sb <- sb[!duplicated(t((apply(sb[, paste0('s', seq(n))], 1, sort)))),]
  # Merge two datasets of fb and sb
  uniq_fb_and_sb = merge(uniq_fb, uniq_sb)
  
  sb2 = uniq_fb_and_sb[, c('lambda', 'eta', paste0('s', seq(num_interims-1)))]
  earlystopsupe <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, paste0('sumx', seq(num_interims - 1))]), 
        row(as.matrix(resp[, paste0('sumx', seq(num_interims - 1))]))
      ), 
      function(x) t(x >= t(sb2[, -(1:2)]))
    )
  )
  
  
  fb2 = uniq_fb_and_sb[, c('lambda', 'gamma', paste0('f', seq(num_interims-1)))]
  earlystopfuti <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, paste0('sumx', seq(num_interims - 1))]), 
        row(as.matrix(resp[, paste0('sumx', seq(num_interims - 1))]))
      ), 
      function(x) t(x <= t(fb2[, -(1:2)]))
    )
  )
  
  # Get the time of futility and efficacy stopping 
  sup_time = do.call(
    function(...) pmin(..., na.rm = TRUE), 
    as.data.frame(
      t(t(earlystopsupe) * seq(num_interims - 1)) * (NA ^ (1 - earlystopsupe))
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
  fb21 = uniq_fb_and_sb[, c('lambda', 'gamma', paste0('f', seq(n)))]
  sb21 = uniq_fb_and_sb[, c('lambda', 'eta', paste0('s', seq(n)))]
  
  condition_met <- do.call(
    rbind,
    lapply(
      split(
        as.matrix(resp[, paste0('sumx', seq(num_interims))]), 
        row(as.matrix(resp[, paste0('sumx', seq(num_interims))]))
      ), 
      function(x) t((x <= t(fb21[, -(1:2)])) | (x >= t(sb21[, -(1:2)])))
    )
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
        t(outer(resp[[paste0('sumx', num_interims)]], sb21[[paste0('s', num_interims)]], '>='))
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


#' Computes both the boundary and corresponding operating characteristics for binary endpoints  
#' @param H0 Null response rates 
#' @param H1 Alternative response rates 
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#'  - 10 is the number of patients enrolled at interim 1,
#'  - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#'  - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#'  - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
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
#'   \item{lambda: }{lambda values for cut-off probability}
#'   \item{gamma: }{gamma valuesfor cut-off probability}
#'   \item{eta: }{eta valuesfor cut-off probability}} 
#' @examples
#' \dontrun{
#' oc_binary<-get_boundary_oc_binary(
#'   H0=0.2, 
#'   H1= 0.4,
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
get_boundary_oc_binary <- function(
    H0, H1, n, nsim, lambda, gamma, eta = NULL, method = "power", seed = NULL  
) {
  # Generate boundary test results
  boundary_tab <- get_boundary_binary(
    H0 = H0, a1 = H0, b1 = 1 - H0, n = n, lambda = lambda, gamma = gamma, eta = eta, method = method, seed=seed
  )
  # Calculate operating characteristics for null and alternative hypotheses
  null_oc <- get_oc_binary(
    p = H0, n = n, nsim = nsim, fb = boundary_tab[['cnf']], sb = boundary_tab[['cns']], seed=seed
  )
  null_oc2 <- null_oc[, colnames(null_oc) %in% c(
    'lambda', 'gamma', 'eta', 'earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean',
    paste0('f', seq(n)), paste0('s', seq(n))
  )]
  colnames(null_oc2) <- c(
    'lambda', 'gamma', 'eta', paste0(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), '_h0'),
    paste0('f', seq(n)), paste0('s', seq(n))
  )
  alt_oc <- get_oc_binary(
    p = H1, n = n, nsim = nsim, fb = boundary_tab[['cnf']], sb = boundary_tab[['cns']], seed=seed
  )
  alt_oc2 <- alt_oc[, colnames(alt_oc) %in% c(
    'lambda', 'gamma', 'eta', 'earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean',
    paste0('f', seq(n)), paste0('s', seq(n))
  )]
  colnames(alt_oc2) <- c(
    'lambda', 'gamma', 'eta', paste0(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), '_h1'),
    paste0('f', seq(n)), paste0('s', seq(n))
  )
  # Combine results into a single data.frame
  all_res_df <- merge(null_oc2, alt_oc2)
  all_res_df <- all_res_df[, c(
    paste0('f', seq(n)), paste0('s', seq(n)),
    paste0(rep(c('earlystopfuti_mean', 'earlystopsupe_mean', 'ss_mean', 'rejectnull_mean'), 2), rep(c('_h0', '_h1'), each = 4)),
    'lambda', 'gamma', 'eta'
  )]
  colnames(all_res_df) <- ifelse(
    colnames(all_res_df) %in% paste0('f', seq(n)),
    paste0('fut_boundary', seq(n)),
    colnames(all_res_df)
  )
  colnames(all_res_df) <- ifelse(
    colnames(all_res_df) %in% paste0('s', seq(n)),
    paste0('sup_boundary', seq(n)),
    colnames(all_res_df)
  )
  return(all_res_df)
}


#' Search optimal parameters for binary endpoint
#' 
#' `search_optimal_pars_binary()` is a helper function and calls `get_boundary_oc_binary()` to obtain the
#'  optimal parameter combinations with type I error less than t1e and sorted according to power. 

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
#'   \item{lambda: }{lambda values for cut-off probability}
#'   \item{gamma: }{gamma valuesfor cut-off probability}
#'   \item{eta: }{eta valuesfor cut-off probability}} 
#'
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#'test_binary <- search_optimal_pars_binary(
#'  H0=0.2,
#'  H1= 0.4,
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
#' 
search_optimal_pars_binary <- function(
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
  all_results = get_boundary_oc_binary(H0, H1, n, nsim, lambda, gamma, eta, method, seed)
  
  # Filter results by t1e
  t1e <- ifelse(is.null(t1e), 1, t1e)
  results_filtered <- all_results[!is.na(all_results$rejectnull_mean_h0) & all_results$rejectnull_mean_h0 < t1e, ]  
  # Sort results by power
  results_sorted <- results_filtered[order(-results_filtered$rejectnull_mean_h1), ]  
  return(results_sorted)
}
