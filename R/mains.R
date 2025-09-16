#' BOP2-FE design for binary endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a binary endpoint
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
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @returns An S3 object of class 'bop2fe' 
#' 
#' @export
#' @examples
#'  \donttest{
#' test_binary <- BOP2FE_binary(
#'  H0=0.2, H1= 0.4,
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000, t1e = 0.1, method = "power",
#'  lambda1 = 0, lambda2 = 1, grid1 = 11,
#'  gamma1 = 0, gamma2 = 1, grid2 = 11,
#'  eta1 = 0, eta2 = 3, grid3 = 31,
#'  seed = 123
#' )
#' summary(test_binary)
#' plot(test_binary)
#' }
#' 
#' 
#' @references 
#' Xu, X., Hashimoto, A., Yimer, B., & Takeda, K. (2025). BOP2-FE: Bayesian optimal phase II design with futility and efficacy stopping boundaries. Journal of Biopharmaceutical Statistics \doi{10.1080/10543406.2025.2558142}.

BOP2FE_binary <- function(H0, H1, n, nsim, t1e = NULL, method = "power", 
                          lambda1, lambda2, grid1, 
                          gamma1, gamma2, grid2, 
                          eta1 = NULL, eta2 = NULL, grid3 = NULL, 
                          seed = NULL) {
  
  design_pars <- list(H0 = H0,
                      H1 = H1, 
                      n = n,
                      cum_n =cumsum(n),
                      nsim = nsim,
                      t1e = t1e,
                      method = method)
  
  all_res <- search_optimal_pars_binary(
    H0, H1, n, nsim, t1e, method, 
    lambda1, lambda2, grid1, 
    gamma1, gamma2, grid2, 
    eta1, eta2, grid3, 
    seed)
  
  search_result <- as.data.frame(all_res)
  

  result <- list(
    end_point = "binary",
    design_pars = design_pars,
    search_result = search_result
  )
  
  class(result) <- "bop2fe"
  
  return(result)
}

#' BOP2-FE design for nested (ordinal) endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a nested (ordinal) endpoint
#' 
#' @param H0 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#'  - `H0[1]`: CR: Complete remission,
#'  - `H0[2]`: PR: Partial remission,
#'  - `H0[3]`: 1-(CR+PR)
#' @param H1 A numeric vector representing the null response rates for different outcomes, specified in the following order:
#'  - `H1[1]`: CR: Complete remission,
#'  - `H1[2]`: PR: Partial remission,
#'  - `H1[3]`: 1-(CR+PR)
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
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#' 
#' @returns An S3 object of class 'bop2fe' 
#' 
#' @export
#' 
#' @examples
#'  \donttest{
#' test_nested <- BOP2FE_nested(
#'  H0=c(0.15,0.15, 0.70), 
#'  H1= c(0.25,0.25, 0.50),
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000, t1e = 0.1, method = "power",
#'  lambda1 = 0, lambda2 = 1, grid1 = 11,
#'  gamma1 = 0, gamma2 = 1, grid2 = 11,
#'  eta1 = 0, eta2 = 3, grid3 = 31,
#'  seed = 123
#')
#'summary(test_nested)
#'plot(test_nested)
#'}
#' 
#' @references 
#' Xu, X., Hashimoto, A., Yimer, B., & Takeda, K. (2025). BOP2-FE: Bayesian optimal phase II design with futility and efficacy stopping boundaries. Journal of Biopharmaceutical Statistics \doi{10.1080/10543406.2025.2558142}.

BOP2FE_nested <- function(H0, H1, n, nsim, t1e = NULL, method = "power", 
                          lambda1, lambda2, grid1, 
                          gamma1, gamma2, grid2, 
                          eta1 = NULL, eta2 = NULL, grid3 = NULL, 
                          seed = NULL) {
  
  
  design_pars <- list(H0 = H0,
                      H1 = H1, 
                      n = n,
                      cum_n =cumsum(n),
                      nsim = nsim,
                      t1e = t1e,
                      method = method)
  
  all_res <- search_optimal_pars_nested(
    H0, H1, n, nsim, t1e, method, 
    lambda1, lambda2, grid1, 
    gamma1, gamma2, grid2, 
    eta1, eta2, grid3, 
    seed)
  
  search_result <- as.data.frame(all_res)
  

  result <- list(
    end_point = "nested",
    design_pars = design_pars,
    search_result = search_result
  )
  
  class(result) <- "bop2fe"
  
  return(result)
}


#' BOP2-FE design for co-primary endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a co-primary endpoint.
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
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#'
#' @returns An S3 object of class 'bop2fe' 
#' 
#' @export
#' @examples
#'  \donttest{
#' test_coprimary <- BOP2FE_coprimary(
#'  H0=c(0.05,0.05, 0.15, 0.75),
#'  H1= c(0.15,0.15, 0.20, 0.50),
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000, t1e = 0.1, method = "power",
#'  lambda1 = 0, lambda2 = 1, grid1 = 11,
#'  gamma1 = 0, gamma2 = 1, grid2 = 11,
#'  eta1 = 0, eta2 = 3, grid3 = 31,
#'  seed = 123
#')
#'summary(test_coprimary)
#'plot(test_coprimary)
#'}
#' 
#' 
#' @references 
#' Xu, X., Hashimoto, A., Yimer, B., & Takeda, K. (2025). BOP2-FE: Bayesian optimal phase II design with futility and efficacy stopping boundaries. Journal of Biopharmaceutical Statistics \doi{10.1080/10543406.2025.2558142}.
#' 
BOP2FE_coprimary <- function(H0, H1, n, nsim, t1e = NULL, method = "power", 
                             lambda1, lambda2, grid1, 
                             gamma1, gamma2, grid2, 
                             eta1 = NULL, eta2 = NULL, grid3 = NULL, 
                             seed = NULL) {
  
  design_pars <- list(H0 = H0,
                      H1 = H1, 
                      n = n,
                      cum_n =cumsum(n),
                      nsim = nsim,
                      t1e = t1e,
                      method = method)
  
  all_res <- search_optimal_pars_coprimary(
    H0, H1, n, nsim, t1e, method, 
    lambda1, lambda2, grid1, 
    gamma1, gamma2, grid2, 
    eta1, eta2, grid3, 
    seed)
  
  search_result <- as.data.frame(all_res)
  

  result <- list(
    end_point = "coprimary",
    design_pars = design_pars,
    search_result = search_result
  )
  
  class(result) <- "bop2fe"
  
  return(result)
}


#' BOP2-FE design for joint efficacy and toxicity endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a joint efficacy and toxicity endpoint.  
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
#' @param method A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param seed for reproducibility             
#' @export
#' 
#' @returns An S3 object of class 'bop2fe' 
#' 
#' @examples
#'  \donttest{
#' test_joint <- BOP2FE_jointefftox(
#'  H0=c(0.15,0.30, 0.15, 0.40),
#'  H1= c(0.18,0.42, 0.02, 0.38),
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000, t1e = 0.1, method = "power",
#'  lambda1 = 0, lambda2 = 1, grid1 = 11,
#'  gamma1 = 0, gamma2 = 1, grid2 = 11,
#'  eta1 = 0, eta2 = 3, grid3 = 31,
#'  seed = 123
#')
#'summary(test_joint)
#'plot(test_joint)
#'}
#' 
#' @references 
#' Xu, X., Hashimoto, A., Yimer, B., & Takeda, K. (2025). BOP2-FE: Bayesian optimal phase II design with futility and efficacy stopping boundaries. Journal of Biopharmaceutical Statistics \doi{10.1080/10543406.2025.2558142} .
#' 
BOP2FE_jointefftox <- function(H0, H1, n, nsim, t1e = NULL, method = "power", 
                               lambda1, lambda2, grid1, 
                               gamma1, gamma2, grid2, 
                               eta1 = NULL, eta2 = NULL, grid3 = NULL, 
                               seed = NULL) {
  

  
  design_pars <- list(H0 = H0,
                      H1 = H1, 
                      n = n,
                      cum_n =cumsum(n),
                      nsim = nsim,
                      t1e = t1e,
                      method = method)
  
  all_res <- search_optimal_pars_efftox(
    H0, H1, n, nsim, t1e, method, 
    lambda1, lambda2, grid1, 
    gamma1, gamma2, grid2, 
    eta1, eta2, grid3, 
    seed)
  
  search_result <- as.data.frame(all_res)
  

  result <- list(
    end_point = "efftox",
    design_pars = design_pars,
    search_result = search_result
  )
  
  class(result) <- "bop2fe"
  
  return(result)
}


#' summarize main results for a given BOP2FE designs
#'
#' @param object the object returned by BOP2FE_xx
#' @param ... additional parameters 
#'
#' @return \code{summary()} returns a list depending on the object entered including design parameters,
#' boundary, operating characteristics. 
#' @rdname summary.bop2fe
#' @method summary bop2fe
#' @export
#' @importFrom utils packageVersion
#' 
summary.bop2fe <- function(object, ...) {
  selected_res <- object$search_result[1,]
  
  grab <- function(data, prefix) {
    cols <- grep(paste0("^", prefix), names(data), value = TRUE)
    list(values = unlist(data[, cols]), nc = length(cols))
  }
  
  if(object$end_point == "binary"){
    fut_boundary <- grab(selected_res, 'fut_boundary')
    sup_boundary <- grab(selected_res, 'sup_boundary')
    nc <- fut_boundary$nc # Number of columns
    summary_tab1 <- matrix(0, 2, nc)
    summary_tab1[1, ] <- fut_boundary$values
    summary_tab1[2, ] <- sup_boundary$values
    row.names(summary_tab1) <- c("Futility boundary", "Efficacy boundary")
    colnames(summary_tab1) <- c(paste0("IA", seq(nc-1)), "FA")
  }else if(object$end_point == "nested"){
    fut_boundary_CR <- grab(selected_res, 'fut_boundary_CR_')
    fut_boundary_CRPR <- grab(selected_res, 'fut_boundary_CR/PR_')
    sup_boundary_CR <- grab(selected_res, 'sup_boundary_CR_')
    sup_boundary_CRPR <- grab(selected_res, 'sup_boundary_CR/PR_')
    nc <- fut_boundary_CR$nc # Number of columns
    summary_tab1 <- matrix(0, 4, nc)
    summary_tab1[1, ] <- fut_boundary_CR$values
    summary_tab1[2, ] <- fut_boundary_CRPR$values
    summary_tab1[3, ] <- sup_boundary_CR$values
    summary_tab1[4, ] <- sup_boundary_CRPR$values
    row.names(summary_tab1) <- c("Futility boundary (CR)", "Futility boundary (CR/PR)",
                                 "Efficacy boundary (CR)", "Efficacy boundary (CR/PR)")
    colnames(summary_tab1) <- c(paste0("IA", seq(nc-1)), "FA")
  }else if(object$end_point == "coprimary"){
    fut_boundary_OR <- grab(selected_res, 'fut_boundary_OR')
    fut_boundary_PFS6 <- grab(selected_res, 'fut_boundary_PFS6')
    sup_boundary_OR <- grab(selected_res, 'sup_boundary_OR')
    sup_boundary_PFS6 <- grab(selected_res, 'sup_boundary_PFS6')
    nc <- fut_boundary_OR$nc # Number of columns
    summary_tab1 <- matrix(0, 4, nc)
    summary_tab1[1, ] <- fut_boundary_OR$values
    summary_tab1[2, ] <- fut_boundary_PFS6$values
    summary_tab1[3, ] <- sup_boundary_OR$values
    summary_tab1[4, ] <- sup_boundary_PFS6$values
    row.names(summary_tab1) <- c("Futility boundary (OR)", "Futility boundary (PFS6)",
                                 "Efficacy boundary (OR)", "Efficacy boundary (PFS6)")
    colnames(summary_tab1) <- c(paste0("IA", seq(nc-1)), "FA")
  } else if(object$end_point == "efftox"){
    fut_boundary_OR <- grab(selected_res, 'fut_boundary_OR')
    fut_boundary_Tox <- grab(selected_res, 'fut_boundary_Tox')
    sup_boundary_OR <- grab(selected_res, 'sup_boundary_OR')
    sup_boundary_Tox <- grab(selected_res, 'sup_boundary_Tox')
    nc <- fut_boundary_OR$nc # Number of columns
    summary_tab1 <- matrix(0, 4, nc)
    summary_tab1[1, ] <- fut_boundary_OR$values
    summary_tab1[2, ] <- fut_boundary_Tox$values
    summary_tab1[3, ] <- sup_boundary_OR$values
    summary_tab1[4, ] <- sup_boundary_Tox$values
    row.names(summary_tab1) <- c("Futility boundary (OR)", "Futility boundary (Tox)",
                                 "Efficacy boundary (OR)", "Efficacy boundary (Tox)")
    colnames(summary_tab1) <- c(paste0("IA", seq(nc-1)), "FA")
  }

  summary_tab2 <- matrix(0, 4, 2)
  summary_tab2[, 1] <- unlist(selected_res[, c('earlystopfuti_mean_h0', 'earlystopsupe_mean_h0', 
                                               'ss_mean_h0', 'rejectnull_mean_h0')])
  summary_tab2[, 2] <- unlist(selected_res[, c('earlystopfuti_mean_h1', 'earlystopsupe_mean_h1',
                                               'ss_mean_h1', 'rejectnull_mean_h1')])
  row.names(summary_tab2) <- c("Early stop for futility", "Early stop for efficacy",
                               "Average sample size", "Reject null")
  colnames(summary_tab2) <- c('Under H0', 'Under H1')
  
  opt_pars <- selected_res[c('lambda', 'gamma', 'eta')]
  row.names(opt_pars) <- NULL
  
  list(
    design_pars = object$design_pars,
    opt_pars = opt_pars,
    boundary = summary_tab1,
    oc = summary_tab2
  )
}


#' Plot the cut-off probability and simulation results for BOP2FE designs
#'
#' Plot the objects returned by other functions, including (1) cut-off probability;
#' (2) boundary values;
#' (3) operating characteristics
#'
#'
#' @param x the object returned by BOP2FE_xx
#' @param ... additional parameters 
#'
#' @return \code{plot()} returns a figure depending on the object entered
#' @rdname plot.bop2fe
#' @method plot bop2fe
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom patchwork plot_annotation wrap_plots 
#' @import ggplot2
#' @export
#' 
plot.bop2fe <- function(x, ...) {
  summary_data <- summary(x)
  summary_tab1 <- summary_data$boundary
  summary_tab2 <- summary_data$oc
  opt_pars <- summary_data$opt_pars
  
  nsum <- sum(x$design_pars$n)
  plot_dat <- data.frame(n = 0:nsum)
  
  plot_dat$postprob_fut <- opt_pars$lambda * (plot_dat$n / nsum)^opt_pars$gamma * 100
  if (x$design_pars$method == "OF") {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + opt_pars$lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  } else { # "power"
    plot_dat$postprob_sup <- (1 - (1 - opt_pars$lambda) * (plot_dat$n / nsum)^opt_pars$eta) * 100
  }
  
  IAs <- cumsum(x$design_pars$n)
  
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = .data$n)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = .data$postprob_sup), color = "red", linewidth = 1, alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$postprob_sup, ymax = 100), fill = "red", alpha = 0.8) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = IAs) +
    ggplot2::scale_y_continuous(name = paste("Cut-off", "\n", "Probability (%)"), breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = IAs, linetype = "dashed") +
    ggplot2::theme_bw()

  
  Oc_tabs2 <- gridExtra::tableGrob(summary_tab2, theme = gridExtra::ttheme_minimal())
  boundary_tab2 <- gridExtra::tableGrob(summary_tab1, theme = gridExtra::ttheme_minimal())
  
  run_date <- Sys.Date()
 
  if(x$end_point=="binary"){
    Info <- paste0("Efficacy cutoff probabilities method - ", x$design_pars$method, ":", " ", "lambda=", opt_pars$lambda, " ", "gamma=", opt_pars$gamma, " ", "eta=", opt_pars$eta)
    Info2 <- paste0("Resp =", x$design_pars$H0)
    Info3 <- paste0("Resp =", x$design_pars$H1)
    title <- "BOP2 FE for binary endpoint"
    
  } else if (x$end_point =="nested"){
    Info <- paste0("Efficacy cutoff probabilities method - ", x$design_pars$method, ":", " ", "lambda=", opt_pars$lambda, " ", "gamma=", opt_pars$gamma, " ", "eta=", opt_pars$eta)
    Info2 <- paste0("CR =", x$design_pars$H0[1], " ", "CR/PR =", x$design_pars$H0[1] + x$design_pars$H0[2])
    Info3 <- paste0("CR =", x$design_pars$H1[1], " ", "CR/PR =", x$design_pars$H1[1] + x$design_pars$H1[2])
    title <- "BOP2 FE for nested endpoint"
  } else if (x$end_point =="coprimary"){
    Info <- paste0("Efficacy cutoff probabilities method - ", x$design_pars$method, ":", " ", "lambda=", opt_pars$lambda, " ", "gamma=", opt_pars$gamma, " ", "eta=", opt_pars$eta)
    Info2 <- paste0("ORR-PFS6 =", x$design_pars$H0[1], " ", "ORR-no PFS6 =", x$design_pars$H0[2], " ", "no ORR-PFS6 =", x$design_pars$H0[3], " ", "no ORR- no PFS6 =", x$design_pars$H0[4])
    Info3 <- paste0("ORR-PFS6 =", x$design_pars$H1[1], " ", "ORR-no PFS6 =", x$design_pars$H1[2], " ", "no ORR-PFS6 =", x$design_pars$H1[3], " ", "no ORR- no PFS6 =", x$design_pars$H1[4])
    title <- "BOP2 FE for coprimary endpoint"
  } else if (x$end_point =="efftox"){
    Info <- paste0("Efficacy cutoff probabilities method - ", x$design_pars$method, ":", " ", "lambda=", opt_pars$lambda, " ", "gamma=", opt_pars$gamma, " ", "eta=", opt_pars$eta)
    Info2 <- paste0("Resp-tox =", x$design_pars$H0[1], " ", "Resp-no tox =", x$design_pars$H0[2], " ", "no Resp-tox =", x$design_pars$H0[3], " ", "no Resp- no tox =", x$design_pars$H0[4])
    Info3 <- paste0("Resp-tox =", x$design_pars$H1[1], " ", "Resp-no tox =", x$design_pars$H1[2], " ", "no Resp-tox =", x$design_pars$H1[3], " ", "no Resp- no tox =", x$design_pars$H1[4])
    title <- "BOP2 FE for joint efficacy and toxicity endpoint"
  } 
  

  layout <- patchwork::wrap_plots(p1, Oc_tabs2, boundary_tab2, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = title,
      subtitle = paste(Info, "\n", "Run Date:", run_date, "\n", "H0:", Info2, "\n", "H1:", Info3),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  print(layout)
}


#' Compute operating characteristics at the optimal boundary 
#'
#' After identifying the optimal boundary that controls the Type I error rate 
#' less than or equal to 0.1 under H0 and maximize the power under H1, it might
#' be of interest to compute the operating characteristics of the optimal boundary 
#' under a different H1 values. This function accepts a single or multiple values of 
#' additional H1 values and compute the operating characteristics for each of them. 
#'
#'
#' @param object the object returned by BOP2FE_xx
#' @param p a single vector or a list of vector for which the operating characteristics is desired.
#' @param endpoint the type of endpoint. Possible options are 'binary', 'nested', 'coprimary', and 'joint'.
#' @param seed for reproducibility 
#'
#' @return \code{simulate_oc()} returns a data frame with the optimal pars and boundary from the given object as well 
#' as the operating characteristics. If a single p vector is supplied the result will be a data frame 
#' with a single row. If multiple p vectors are supplied the data frame will be have multiple rows each 
#' corresponding to the p vectors in the order of their specification 
#' @export
#' 
#' @examples
#'  \donttest{
#' test_nested <- BOP2FE_nested(
#'  H0=c(0.15,0.15, 0.70), 
#'  H1= c(0.25,0.25, 0.50),
#'  n = c(10, 5, 5, 5, 5, 5, 5),
#'  nsim = 1000, t1e = 0.1, method = "power",
#'  lambda1 = 0, lambda2 = 1, grid1 = 11,
#'  gamma1 = 0, gamma2 = 1, grid2 = 11,
#'  eta1 = 0, eta2 = 3, grid3 = 31,
#'  seed = 123
#')
#'
#' # Compute operating characteristics for a single p vector 
#' simulate_oc(test_nested, p=c(0.30,0.30,0.40), 
#'              endpoint = 'nested', seed=123)
#' 
#' # Compute operating characteristics for multiple p vector
#'  simulate_oc(test_nested, p=list(c(0.30,0.30,0.40),c(0.35,0.35,0.30)), 
#'              endpoint = 'nested', seed=123)
#' 
#'}
#' 
simulate_oc <- function(object, p, endpoint, seed = NULL) {
  boundary <- summary(object)$boundary
  opt_pars <- summary(object)$opt_pars
  design_par <- object$design_pars
  
  n <- design_par$n
  nsim <- design_par$nsim
  
  # Helper function to validate and simulate for a single vector
  simulate_single <- function(p_vec) {
    # Validate length of p_vec based on endpoint
    expected_length <- switch(endpoint,
                              binary = 1,
                              nested = 3,
                              coprimary = 4,
                              joint = 4,
                              stop("Invalid endpoint specified. Must be one of: 'binary', 'nested', 'coprimary', 'joint'."))
    
    if (length(p_vec) != expected_length) {
      stop(sprintf("For endpoint '%s', expected a vector of length %d, but got length %d.",
                   endpoint, expected_length, length(p_vec)))
    }
    
    if (!inherits(object,"bop2fe")) {
      stop("Object should be a class of bop2fe")
    }
    # Proceed with simulation
    if (endpoint == "binary") {
      cnf <- cbind(opt_pars[, 1:2], t(boundary[1, ]))
      colnames(cnf) <- c('lambda', 'gamma', paste0('f', seq(n)))
      cns <- cbind(opt_pars[, c(1, 3)], t(boundary[2, ]))
      colnames(cns) <- c('lambda', 'eta', paste0('s', seq(n)))
      return(get_oc_binary(p = p_vec, n = n, nsim = nsim, fb = cnf, sb = cns, seed = seed)[-c(8:11)])
    } else if (endpoint == "nested") {
      p1 <- p_vec[1]; p2 <- p_vec[2]; p3 <- p_vec[3]
      cnf <- cbind(opt_pars[, 1:2], t(boundary[1, ]), t(boundary[2, ]))
      colnames(cnf) <- c('lambda', 'gamma', paste0('f1', seq(n)), paste0('f2', seq(n)))
      cns <- cbind(opt_pars[, c(1, 3)], t(boundary[3, ]), t(boundary[4, ]))
      colnames(cns) <- c('lambda', 'eta', paste0('s1', seq(n)), paste0('s2', seq(n)))
      return(get_oc_nested(p1 = p1, p2 = p2, p3 = p3, n = n, nsim = nsim, fb = cnf, sb = cns, seed = seed)[-c(8:11)])
    } else if (endpoint == "coprimary") {
      p1 <- p_vec[1]; p2 <- p_vec[2]; p3 <- p_vec[3]; p4 <- p_vec[4]
      cnf <- cbind(opt_pars[, 1:2], t(boundary[1, ]), t(boundary[2, ]))
      colnames(cnf) <- c('lambda', 'gamma', paste0('f1', seq(n)), paste0('f2', seq(n)))
      cns <- cbind(opt_pars[, c(1, 3)], t(boundary[3, ]), t(boundary[4, ]))
      colnames(cns) <- c('lambda', 'eta', paste0('s1', seq(n)), paste0('s2', seq(n)))
      return(get_oc_coprimary(p1 = p1, p2 = p2, p3 = p3, p4 = p4, n = n, nsim = nsim, fb = cnf, sb = cns, seed = seed)[-c(8:11)])
    } else if (endpoint == "joint") {
      p1 <- p_vec[1]; p2 <- p_vec[2]; p3 <- p_vec[3]; p4 <- p_vec[4]
      cnf <- cbind(opt_pars[, 1:2], t(boundary[1, ]), t(boundary[2, ]))
      colnames(cnf) <- c('lambda', 'gamma', paste0('f1', seq(n)), paste0('f2', seq(n)))
      cns <- cbind(opt_pars[, c(1, 3)], t(boundary[3, ]), t(boundary[4, ]))
      colnames(cns) <- c('lambda', 'eta', paste0('s1', seq(n)), paste0('s2', seq(n)))
      return(get_oc_jointefftox(p1 = p1, p2 = p2, p3 = p3, p4 = p4, n = n, nsim = nsim, fb = cnf, sb = cns, seed = seed)[-c(8:11)])
    }
  }
  
  # Check if p is a list or a single vector
  if (is.list(p)) {
    results <- lapply(p, simulate_single)
    return(do.call(rbind, results))
  } else {
    return(simulate_single(p))
  }
}

