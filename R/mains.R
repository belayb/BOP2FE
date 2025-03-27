#' BOP2-FE design for binary endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a binary endpoint given optimal
#' design parameter values that are identified to optimize power while controlling
#' type I error rate at a specified value. The optimal parameters (i.e., lambda, gamma, and eta)
#' can be obtained by calling `search_optimal_pars_binary()` function. The function produce a plot
#' that can be saved as a pdf file for documentation.  
#' 
#' 
#' @param H0 the probability of response under the null hypothesis 
#' @param H1 the probability of response under the alternative hypothesis 
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param nsim number of simulation
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom ggplot2  ggplot geom_ribbon geom_line scale_x_continuous scale_y_continuous geom_vline theme_minimal
#' @importFrom tibble tibble
#' @return A list of class \code{"BOP2FE"} containing the following elements:
#' \item{boundary}{A table of stopping boundaries.}
#' \item{Oc}{A table of operating characteristics.}
#' \item{plot}{A plot that shows the boundaries graphically along with tables of 
#' boundary values and operating characteristics}
#' @examples
#' \dontrun{
#' # Example with 7 interim looks
#' BOP2FE_binary(H0 = 0.2, H1=0.4, n = c(10, 5, 5, 5, 5, 5, 5), 
#' lambda = 0.909, gamma = 1, nsim = 10000, seed = 1234)
#' }
#' @export
#'
BOP2FE_binary <- function(H0, H1, n, lambda = NULL, gamma = NULL, eta = NULL, method = "power", nsim = NULL, seed = NULL) {
  
  a1 <- H0
  b1 <- 1 - H0
  nIA <- length(n)
  nsum <- sum(n)
  
  # Check total sample size
  if (nsum == 0) {
    stop("The total (final) sample size can not be zero.")
  }
  
  # Set default method to "OBrien-Fleming" if eta is NULL
  #if (is.null(eta) | is.null(method)) {
  #  method <- "OF"
  #}
  if (is.null(gamma)) {
    gamma <- 0.95
    message("gamma value should be provided. The defult gamma = 0.95 used")
  }
  if (is.null(lambda)) {
    lambda <- 0.95
    message("lambda value should be provided. The defult lambda=0.95 used")
  }
  if (is.null(eta) & method == "power") {
    eta <- 0.95
    message("eta value should be provided. The defult eta=0.95 used")
  }
  
  if (is.null(nsim)) {
    nsim <- 10000
    message("The defult 10000 simulation used")
    
  }
  if (is.null(seed)) {
    seed <- 1234
  }
  
  boundary_tab <- boundary_binary(H0 = H0, a1 = a1, b1 = b1, n = n,
                                  lambda = lambda, gamma = gamma, eta = eta,
                                  method = method, seed = seed)
  
  fb <- boundary_tab$cnf
  sb <- boundary_tab$cns
  
  Oc_tab_null <- Oc_binary(p = H0, n = n, nsim = nsim, fb = fb, sb = sb, seed = seed)
  Oc_tab_alt <- Oc_binary(p = H1, n = n, nsim = nsim, fb = fb, sb = sb, seed = seed)
  
  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)
  
  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100
  
  # Calculate postprob_sup based on method
  if (method == "OF") {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  } else { # "power"
    plot_dat$postprob_sup <- (1 - (1 - lambda) * (plot_dat$n/nsum)^eta)*100
  }
  
  # Plot
  IAs<- cumsum(n)
  
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1, alpha=0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.8) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = paste("Cut-off", "\n", "Probability (%)"), breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = IAs, linetype = "dashed") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Under_H0 = c(Oc_tab_null$earlystopfuti_mean*100, Oc_tab_null$earlystopsupe_mean*100, Oc_tab_null$ss_mean, Oc_tab_null$rejectnull_mean*100),
                            Under_H1 = c(Oc_tab_alt$earlystopfuti_mean*100, Oc_tab_alt$earlystopsupe_mean*100, Oc_tab_alt$ss_mean, Oc_tab_alt$rejectnull_mean*100))
  
  Oc_tabs2 <- dplyr::as_tibble(Oc_tabs,.name_repair = "minimal") %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  names(boundary_tab) <- c("Futility boundary", "Superiority boundary")
  boundary_tab <- tibble::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab2 <- dplyr::as_tibble(boundary_tab,.name_repair = "minimal") %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  run_date <- Sys.Date()
  Info <- paste0("Efficacy cutoff probabilities method - ", method, ":", " ", "lambda=", lambda, " ", "gamma=", gamma, " ",  "eta=", ifelse(is.null(eta)&(method=="OF"), NA, eta))
  layout <- patchwork::wrap_plots(p1, Oc_tabs2, boundary_tab2, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for binary outcome",
      subtitle = paste(Info, "\n", "Run Date:", run_date, "\n", "H0:", H0, ", ", "H1:", H1),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5, face = "bold"))
    )
  
  # Create a list to return
  result <- list(
    boundary = boundary_tab,
    Oc = Oc_tabs,
    plot = layout
  )
  
  # Assign the S3 class
  class(result) <- "BOP2FE"
  
  return(result)
}

#' BOP2-FE design for nested (ordinal) endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a nested (ordinal) endpoint given optimal
#' design parameter values that are identified to optimize power while controlling
#' type I error rate at a specified value. The optimal parameters (i.e., lambda, gamma, and eta)
#' can be obtained by calling `search_optimal_pars_nested()` function. The function produce a plot
#' that can be saved as a pdf file for documentation.  
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
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param nsim number of simulation
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom ggplot2  ggplot geom_ribbon geom_line scale_x_continuous scale_y_continuous geom_vline theme_minimal
#' @importFrom tibble tibble
#' 
#' @return A list of class \code{"BOP2FE"} containing the following elements:
#' \item{boundary}{A table of stopping boundaries.}
#' \item{Oc}{A table of operating characteristics.}
#' \item{plot}{A plot that shows the boundaries graphically along with tables of 
#' boundary values and operating characteristics}
#' @examples
#' \dontrun{
#' # Example with 7 interim looks
#' BOP2FE_nested(H0 = c(0.15,0.15,0.70), H1 = c(0.25,0.25,0.50), n=c(10,5,5,5,5,5,5), 
#' lambda = 0.95, gamma=1, seed = 123)
#' }
#' @export
#'

BOP2FE_nested <- function(H0, H1, n, lambda = NULL, gamma=NULL, eta=NULL,  method = "power", nsim = NULL, seed = NULL){
  
  a <- H0
  p11 = H0[1]; p21 = H1[1];
  p12 = H0[2]; p22 = H1[2];
  p13 = H0[3]; p23 = H1[3];
  
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("The total (final) sample size can not be zero.")
  }

  # Set default method to "OBrien-Fleming" if eta is NULL
  #if (is.null(eta) | is.null(method)) {
  #  method <- "OF"
  #}
  if (is.null(gamma)) {
    gamma <- 0.95
    message("gamma value should be provided. The defult gamma = 0.95 used")
  }
  if (is.null(lambda)) {
    lambda <- 0.95
    message("lambda value should be provided. The defult lambda=0.95 used")
  }
  if (is.null(eta) & method == "power") {
    eta <- 0.95
    message("eta value should be provided. The defult eta=0.95 used")
  }
  
  if (is.null(nsim)) {
    nsim <- 10000
    message("The defult 10000 simulation used")
    
  }
  if (is.null(seed)) {
    seed <- 1234
  }


  boundary_tab <- boundary_nested(H0 = H0, a=a,  n=n,
                                     lambda=lambda, gamma=gamma, eta = eta,
                                     method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab_null <- Oc_nested(p1 = p11, p2 = p12, p3 = p13, n = n, 
                           nsim = nsim, fb = fb, sb =sb, seed = seed)

  Oc_tab_alt <- Oc_nested(p1 = p21, p2 = p22, p3 = p23, n = n,
                           nsim = nsim, fb = fb, sb =sb, seed = seed)
  
  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "OF") {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  } else { # "power"
    plot_dat$postprob_sup <- (1 - (1 - lambda) * (plot_dat$n/nsum)^eta)*100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = paste("Cut-off", "\n", "Probability (%)"), breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Under_H0 = c(Oc_tab_null$earlystopfuti_mean*100, Oc_tab_null$earlystopsupe_mean*100, 
                                      Oc_tab_null$ss_mean, Oc_tab_null$rejectnull_mean*100),
                            Under_H1 = c(Oc_tab_alt$earlystopfuti_mean*100, Oc_tab_alt$earlystopsupe_mean*100, 
                                      Oc_tab_alt$ss_mean, Oc_tab_alt$rejectnull_mean*100))
  
  Oc_tabs2 <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary CR", "Futility boundary CR/PR", "Superiority boundary CR", "Superiority boundary CR/PR")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab2 <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  Info <- paste0("Efficacy cutoff probabilities method - ", method, ":", " ", "lambda=", lambda, " ", "gamma=", gamma, " ",  "eta=", ifelse(is.null(eta)&(method=="OF"), NA, eta))
  Info2 <- paste0("CR=", H0[1], " ", "CR/PR=", H0[1] + H0[2])
  Info3 <- paste0("CR=", H1[1], " ", "CR/PR=", H1[1] + H1[2])
  run_date <- Sys.Date()
  layout <- patchwork::wrap_plots(p1, Oc_tabs2, boundary_tab2, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for nested outcome",
      subtitle = paste(Info, "\n", "Run Date:", run_date, "\n", "H0:", Info2, "\n", "H1:", Info3),
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  # Create a list to return
  result <- list(
    boundary = boundary_tab,
    Oc = Oc_tabs,
    plot = layout
  )
  
  # Assign the S3 class
  class(result) <- "BOP2FE"
  
  return(result)  
}


#' BOP2-FE design for co-primary endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a co-primary endpoint given optimal
#' design parameter values that are identified to optimize power while controlling
#' type I error rate at a specified value. The optimal parameters (i.e., lambda, gamma, and eta)
#' can be obtained by calling `search_optimal_pars_coprimary()` function. The function produce a plot
#' that can be saved as a pdf file for documentation.  
#' 
#' @param H0 Response rate under the null (Response - PFS6, Response - no PFS6, No response - PFS6, No response - No PFS6)
#' @param H1 Response rate under the alternative (Response - PFS6, Response - no PFS6, No response - PFS6, No response - No PFS6)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param nsim number of simulation
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom ggplot2  ggplot geom_ribbon geom_line scale_x_continuous scale_y_continuous geom_vline theme_minimal
#' @importFrom tibble tibble
#' 
#' @return A list of class \code{"BOP2FE"} containing the following elements:
#' \item{boundary}{A table of stopping boundaries.}
#' \item{Oc}{A table of operating characteristics.}
#' \item{plot}{A plot that shows the boundaries graphically along with tables of 
#' boundary values and operating characteristics}
#' @examples
#' \dontrun{
#' # Example with 7 interim looks
#' BOP2FE_coprimary(H0 = c(0.05,0.05,0.15,0.75), H1 =c(0.15,0.15,0.20,0.50), n=c(10,5,5,5,5,5,5), 
#' lambda = 0.95, gamma=1,seed = 123)
#' }
#' @export
#'

BOP2FE_coprimary <- function(H0, H1, n, lambda = NULL, gamma=NULL, eta=NULL,  method = "power", nsim = NULL, seed = NULL){

  a <- H0
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("The total (final) sample size can not be zero.")
  }

  # Set default method to "OBrien-Fleming" if eta is NULL
  #if (is.null(eta) | is.null(method)) {
  #  method <- "OF"
  #}
  if (is.null(gamma)) {
    gamma <- 0.95
    message("gamma value should be provided. The defult gamma = 0.95 used")
  }
  if (is.null(lambda)) {
    lambda <- 0.95
    message("lambda value should be provided. The defult lambda=0.95 used")
  }
  if (is.null(eta) & method == "power") {
    eta <- 0.95
    message("eta value should be provided. The defult eta=0.95 used")
  }
  
  if (is.null(nsim)) {
    nsim <- 10000
    message("The defult 10000 simulation used")
    
  }
  if (is.null(seed)) {
    seed <- 1234
  }

  boundary_tab <- boundary_coprimary(H0 = H0, a=a, n=n,
                                  lambda=lambda, gamma=gamma, eta = eta,
                                  method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab_null <- Oc_coprimary(p1 = H0[1], p2 = H0[2], p3 = H0[3],p4 = H0[4], n = n,
                              nsim = nsim, fb = fb, sb =sb, seed = seed)
  Oc_tab_alt <- Oc_coprimary(p1 = H1[1], p2 = H1[2], p3 = H1[3],p4 = H0[4], n = n,
                              nsim = nsim, fb = fb, sb =sb, seed = seed)
  

  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "OF") {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  } else { # "power"
    plot_dat$postprob_sup <- (1 - (1 - lambda) * (plot_dat$n/nsum)^eta)*100
  }
  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = paste("Cut-off", "\n", "Probability (%)"), breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Under_H0 = c(Oc_tab_null$earlystopfuti_mean*100, Oc_tab_null$earlystopsupe_mean*100, 
                                         Oc_tab_null$ss_mean, Oc_tab_null$rejectnull_mean*100),
                            Under_H1 = c(Oc_tab_alt$earlystopfuti_mean*100, Oc_tab_alt$earlystopsupe_mean*100, 
                                         Oc_tab_alt$ss_mean, Oc_tab_alt$rejectnull_mean*100))
  
  Oc_tabs2 <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary ORR", "Futility boundary PFS6", "Superiority boundary ORR", "Superiority boundary PFS6")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab2 <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  run_date <- Sys.Date()
  
  Info <- paste0("Efficacy cutoff probabilities method - ", method, ":", " ", "lambda=", lambda, " ", "gamma=", gamma, " ",  "eta=", ifelse(is.null(eta)&(method=="OF"), NA, eta))
  Info2 <- paste0("ORR-PFS6=", H0[1], " ", "ORR-no PFS6=", H0[2], " ", "no ORR-PFS6=", H0[3], " ", "no ORR- no PFS6=", H0[4])
  Info3 <- paste0("ORR-PFS6=", H1[1], " ", "ORR-no PFS6=", H1[2], " ", "no ORR-PFS6=", H1[3], " ", "no ORR- no PFS6=", H1[4])
  
  layout <- patchwork::wrap_plots(p1, Oc_tabs2, boundary_tab2, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for co-primary outcome",
      subtitle = paste(Info, "\n", "Run Date:", run_date, "\n", "H0:", Info2, "\n", "H1:", Info3),
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  # Create a list to return
  result <- list(
    boundary = boundary_tab,
    Oc = Oc_tabs,
    plot = layout
  )
  
  # Assign the S3 class
  class(result) <- "BOP2FE"
  
  return(result)  

}



#' BOP2-FE design for joint efficacy and toxicity endpoint 
#' 
#' Computes stopping boundaries and operating characteristics of Bayesian optimal phase II
#' design with efficacy and futility stopping for a joint efficacy and toxicity endpoint given optimal
#' design parameter values that are identified to optimize power while controlling
#' type I error rate at a specified value. The optimal parameters (i.e., lambda, gamma, and eta)
#' can be obtained by calling `search_optimal_pars_efftox()` function. The function produce a plot
#' that can be saved as a pdf file for documentation.  
#' 
#' @param H0 Response rate under the null (toxicity - OR, no toxicity - OR, toxicity - no OR, no toxicity - No OR)
#' @param H1 Response rate under the alternative (toxicity - OR, no toxicity - OR, toxicity - no OR, no toxicity - No OR)
#' @param n A numeric vector representing the additional patients enrolled at each interim analysis. 
#' The value at index `i` indicates the number of new patients added at interim analysis `i`. 
#' The total sample size at interim `i` is the cumulative sum of the values in `n` up to that index. 
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30, 
#' the vector would be represented as `n = c(10, 5, 5, 10)`, where:
#' - 10 is the number of patients enrolled at interim 1,
#' - 5 (15 - 10) is the additional number of patients enrolled at interim 2,
#' - 5 (20 - 15) is the additional number of patients enrolled at interim 3,
#' - 10 (30 - 20) is the additional number of patients enrolled at interim 4.
#' @param lambda A numeric value for parameter `lambda` of the cut-off probability (i.e common for both efficacy and futility cut-off probability)
#' @param gamma A numeric value for parameter `gamma` of the cut-off probability for futility 
#' @param eta A numeric value for parameter `eta` of the cut-off probability for efficacy 
#' @param method A character string specifying the method to use for calculating cutoff values.
#'               Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param nsim number of simulation
#' @param seed for reproducibility
#' @importFrom dplyr mutate case_when
#' @importFrom stats pbeta rbinom
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom ggplot2  ggplot geom_ribbon geom_line scale_x_continuous scale_y_continuous geom_vline theme_minimal
#' @importFrom tibble tibble
#' @return A list of class \code{"BOP2FE"} containing the following elements:
#' \item{boundary}{A table of stopping boundaries.}
#' \item{Oc}{A table of operating characteristics.}
#' \item{plot}{A plot that shows the boundaries graphically along with tables of 
#' boundary values and operating characteristics}
#' @examples
#' \dontrun{
#' # Example with 7 interim looks
#' BOP2FE_jointefftox(H0 = c(0.15, 0.30, 0.15, 0.40), H1=c(0.18, 0.42, 0.02, 0.38), 
#' n=c(10,5,5,5,5,5,5), 
#' lambda = 0.7, gamma=1, seed = 123)
#' }
#' @export
#'
BOP2FE_jointefftox <- function(H0, H1, n, lambda = NULL, gamma=NULL, eta=NULL,  method = "power", nsim = NULL, seed = NULL){

  a <- H0
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("The total (final) sample size can not be zero.")
  }

  # Set default method to "OBrien-Fleming" if eta is NULL
  #if (is.null(eta) | is.null(method)) {
  #  method <- "OF"
  #}
  if (is.null(gamma)) {
    gamma <- 0.95
    message("gamma value should be provided. The defult gamma = 0.95 used")
  }
  if (is.null(lambda)) {
    lambda <- 0.95
    message("lambda value should be provided. The defult lambda=0.95 used")
  }
  if (is.null(eta) & method == "power") {
    eta <- 0.95
    message("eta value should be provided. The defult eta=0.95 used")
  }
  
  if (is.null(nsim)) {
    nsim <- 10000
    message("The defult 10000 simulation used")
    
  }
  if (is.null(seed)) {
    seed <- 1234
  }


  boundary_tab <- boundary_jointefftox(H0 = H0, a=a, n=n,
                                     lambda=lambda, gamma=gamma, eta = eta,
                                     method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab_null <- Oc_jointefftox(p1 = H0[1], p2 = H0[2], p3 = H0[3],p4 = H0[4], n = n,
                         nsim = nsim, fb = fb, sb =sb, seed = seed)
  Oc_tab_alt <- Oc_jointefftox(p1 = H1[1], p2 = H1[2], p3 = H1[3],p4 = H1[4], n = n,
                         nsim = nsim, fb = fb, sb =sb, seed = seed)
  

  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "OF") {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  } else { # "power"
    plot_dat$postprob_sup <- (1 - (1 - lambda) * (plot_dat$n/nsum)^eta)*100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = paste("Cut-off", "\n", "Probability (%)"), breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Under_H0 = c(Oc_tab_null$earlystopfuti_mean*100, Oc_tab_null$earlystopsupe_mean*100,
                                         Oc_tab_null$ss_mean, Oc_tab_null$rejectnull_mean*100),
                            Under_H1 = c(Oc_tab_alt$earlystopfuti_mean*100, Oc_tab_alt$earlystopsupe_mean*100,
                                         Oc_tab_alt$ss_mean, Oc_tab_alt$rejectnull_mean*100))
  
  Oc_tabs2 <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary response", "Futility boundary toxicity", "Superiority boundary response", "Superiority boundary toxicity")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab2 <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  Info <- paste0("Efficacy cutoff probabilities method - ", method, ":", " ", "lambda=", lambda, " ", "gamma=", gamma, " ",  "eta=", ifelse(is.null(eta)&(method=="OF"), NA, eta))
  run_date <- Sys.Date()
  
  Info2 <- paste0("Resp-tox=", H0[1], " ", "Resp-no tox=", H0[2], " ", "no Resp-tox=", H0[3], " ", "no Resp- no tox=", H0[4])
  Info3 <- paste0("Resp-tox=", H1[1], " ", "Resp-no tox=", H1[2], " ", "no Resp-tox=", H1[3], " ", "no Resp- no tox=", H1[4])
  
  layout <- patchwork::wrap_plots(p1, Oc_tabs2, boundary_tab2, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for joint efficacy and toxicity",
      subtitle = paste(Info, "\n", "Run Date:", run_date, "\n", "H0:", Info2, "\n", "H1:", Info3),
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  # Create a list to return
  result <- list(
    boundary = boundary_tab,
    Oc = Oc_tabs,
    plot = layout
  )
  
  # Assign the S3 class
  class(result) <- "BOP2FE"
  
  return(result)  

}

