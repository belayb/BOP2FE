#' BOP2 design for binary endpoint 
#' @param H0 Response rate 
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "OBrien-Fleming" function. method=power is an alternative
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
#' @export
#'
BOP2FE_binary <- function(H0, n, lambda = NULL, gamma=NULL, eta=NULL,  method = NULL, nsim = NULL, seed = NULL){

  a1 <- H0
  b1 <- 1 - H0
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("nsum cannot be zero.")
  }

  # Set default method to "obrain" if eta is NULL
  if (is.null(eta)|is.null(method)) {
    method <- "OBrien-Fleming"
  }
  if(is.null (gamma)) {gamma = 1; message("gamma value should be provided")}
  if(is.null (lambda)) {lambda = 0.95; message("lambda value should be provided")}
  if(is.null (nsim)) {nsim = 10000; }
  if(is.null (seed)) {seed = 1234; }

  boundary_tab <- boundary_binary(H0 = H0, a1 = a1, b1=b1, nIA=nIA, n=n,
                                  lambda=lambda, gamma=gamma, eta = eta,
                                  method = method, seed = seed)

  fb<- boundary_tab$cnf
  sb <- boundary_tab$cns

  Oc_tab <- Oc_binary(p = H0, n = n, nsim = nsim, fb = fb, sb =sb, seed = seed)


  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "power") {
    plot_dat$postprob_sup <- 1 - (1 - lambda) * (plot_dat$n / nsum)^eta
  } else {#"OBrien-Fleming"
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = "Cut-off Probability (%)", breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  


  Oc_tabs <- tibble::tibble(Statistic = c("Early stoping for Futility (%)",
                                  "Early stoping for superiority (%)",
                                  "Average sample size",
                                  "Null rejection (%)"),
                    Value = c(Oc_tab$earlystopfuti_mean, Oc_tab$earlystopsupe_mean, Oc_tab$ss_mean, Oc_tab$rejectnull_mean))

  Oc_tabs<-dplyr::as_tibble(Oc_tabs)%>%gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  names(boundary_tab) <- c("Futility boundary", "Superiority boundary")
  boundary_tab<- tibble::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interm analysis", 1:nIA)
  boundary_tab<-dplyr::as_tibble(boundary_tab)%>%gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)

  Info<- paste0("lambda=", lambda, " ", "gamma=" ,gamma )
  layout <- patchwork::wrap_plots(p1,Oc_tabs, boundary_tab, nrow=3)
  layout<-layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for binary outcome",
      subtitle = Info,
      #caption = paste0("Design Pars: n= ",n, "lambda=", lambda, "gamma=" gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5, face = "bold"))
    )

  return(layout)

}

#' BOP2 design for nested (ordinal) endpoint 
#' @param CR0  complete remission rate under the null (CR)
#' @param CRPR0 complete remission or partial remission rate under the null (CR or PR)
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "OBrien-Fleming" function. method=power is an alternative
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
#' @export
#'

BOP2FE_nested <- function(CR0, CRPR0, n, lambda = NULL, gamma=NULL, eta=NULL,  method = NULL, nsim = NULL, seed = NULL){
  H0 <- c(CR0, CRPR0 - CR0, 1 - CRPR0)
  a <- H0
  p1 = H0[1]
  p2 = H0[2]
  p3 = H0[3]
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("nsum cannot be zero.")
  }

  # Set default method to "obrain" if eta is NULL
  if (is.null(eta)|is.null(method)) {
    method <- "OBrien-Fleming"
  }
  if(is.null (gamma)) {gamma = 1; message("gamma value should be provided")}
  if(is.null (lambda)) {lambda = 0.95; message("lambda value should be provided")}
  if(is.null (nsim)) {nsim = 10000; }
  if(is.null (seed)) {seed = 1234; }


  boundary_tab <- boundary_nested(H0 = H0, a=a, nIA=nIA, n=n,
                                     lambda=lambda, gamma=gamma, eta = eta,
                                     method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab <- Oc_nested(p1 = p1, p2 = p2, p3 = p3, n = n,
                         n_stage = nIA, nsim = nsim, fb = fb, sb =sb, seed = seed)


  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "power") {
    plot_dat$postprob_sup <- 1 - (1 - lambda) * (plot_dat$n / nsum)^eta
  } else {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = "Cut-off Probability (%)", breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Value = c(Oc_tab$earlystopfuti_mean, Oc_tab$earlystopsupe_mean, Oc_tab$ss_mean, Oc_tab$rejectnull_mean))
  
  Oc_tabs <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary CR", "Futility boundary CR/PR", "Superiority boundary CR", "Superiority boundary CR/PR")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  Info <- paste0("lambda=", lambda, " ", "gamma=", gamma)
  layout <- patchwork::wrap_plots(p1, Oc_tabs, boundary_tab, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for nested outcome",
      subtitle = Info,
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  return(layout)
  
}


#' BOP2 design for co-primary endpoint 
#' @param H0 Response rate under the null (Response - PFS6, Response - no PFS6, No response - PFS6, No response - No PFS6)
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "OBrien-Fleming" function. method=power is an alternative
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
#' @export
#'

BOP2FE_coprimary <- function(H0, n, lambda = NULL, gamma=NULL, eta=NULL,  method = NULL, nsim = NULL, seed = NULL){

  a <- H0
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("nsum cannot be zero.")
  }

  # Set default method to "obrain" if eta is NULL
  if (is.null(eta)|is.null(method)) {
    method <- "OBrien-Fleming"
  }
  if(is.null (gamma)) {gamma = 1; message("gamma value should be provided")}
  if(is.null (lambda)) {lambda = 0.95; message("lambda value should be provided")}
  if(is.null (nsim)) {nsim = 10000; }
  if(is.null (seed)) {seed = 1234; }


  boundary_tab <- boundary_coprimary(H0 = H0, a=a, nIA=nIA, n=n,
                                  lambda=lambda, gamma=gamma, eta = eta,
                                  method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab <- Oc_coprimary(p1 = H0[1], p2 = H0[2], p3 = H0[3],p4 = H0[4], n = n,
                         n_stage = nIA, nsim = nsim, fb = fb, sb =sb, seed = seed)


  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "power") {
    plot_dat$postprob_sup <- 1 - (1 - lambda) * (plot_dat$n / nsum)^eta
  } else {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = "Cut-off Probability (%)", breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Value = c(Oc_tab$earlystopfuti_mean, Oc_tab$earlystopsupe_mean, Oc_tab$ss_mean, Oc_tab$rejectnull_mean))
  
  Oc_tabs <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary ORR", "Futility boundary PFS6", "Superiority boundary ORR", "Superiority boundary PFS6")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  Info <- paste0("lambda=", lambda, " ", "gamma=", gamma)
  layout <- patchwork::wrap_plots(p1, Oc_tabs, boundary_tab, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for co-primary outcome",
      subtitle = Info,
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  return(layout)
  

}



#' BOP2 design for joint efficacy and toxicity endpoint 
#' @param H0 Response rate under the null (toxicity - OR, no toxicity - OR, toxicity - no OR, no toxicity - No OR)
#' @param n vector of sample size at each interim look (Note that the sample size at interim i is the difference between sample size at interim i and at interim i-1)
#' @param lambda optimal value for lambda of the cut-off probability
#' @param gamma optimal value for gamma of the cut-off probability
#' @param eta optimal value for eta of the cut-off probability
#' @param method type of function to be used for the cut off probability for superiority. The default is "OBrien-Fleming" function. method=power is an alternative
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
#' @export
#'
BOP2FE_jointefftox <- function(H0, n, lambda = NULL, gamma=NULL, eta=NULL,  method = NULL, nsim = NULL, seed = NULL){

  a <- H0
  nIA <- length(n)
  nsum<- sum(n)
  # Check total sample size
  if (nsum == 0) {
    stop("nsum cannot be zero.")
  }

  # Set default method to "obrain" if eta is NULL
  if (is.null(eta)|is.null(method)) {
    method <- "OBrien-Fleming"
  }
  if(is.null (gamma)) {gamma = 1; message("gamma value should be provided")}
  if(is.null (lambda)) {lambda = 0.95; message("lambda value should be provided")}
  if(is.null (nsim)) {nsim = 10000; }
  if(is.null (seed)) {seed = 1234; }


  boundary_tab <- boundary_jointefftox(H0 = H0, a=a, nIA=nIA, n=n,
                                     lambda=lambda, gamma=gamma, eta = eta,
                                     method = method, seed = seed)

  fb<- c(rbind(boundary_tab$cn11f_max, boundary_tab$cn12f_max))
  sb <- c(rbind(boundary_tab$cn11s_min, boundary_tab$cn12s_min))

  Oc_tab <- Oc_jointefftox(p1 = H0[1], p2 = H0[2], p3 = H0[3],p4 = H0[4], n = n,
                         n_stage = nIA, nsim = nsim, fb = fb, sb =sb, seed = seed)


  # Create data frame
  plot_dat <- data.frame(n = 0:nsum)

  # Calculate postprob_fut
  plot_dat$postprob_fut <- lambda * (plot_dat$n / nsum)^gamma * 100

  # Calculate postprob_sup based on method
  if (method == "power") {
    plot_dat$postprob_sup <- 1 - (1 - lambda) * (plot_dat$n / nsum)^eta
  } else {
    plot_dat$postprob_sup <- (2 * pnorm(qnorm((1 + lambda) / 2) / sqrt(plot_dat$n / nsum)) - 1) * 100
  }

  # Plot
  p1 <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = n)) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_fut), color = "blue", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = postprob_fut), fill = "blue", alpha = 0.7) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_fut, ymax = 100), fill = "gray", alpha = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = postprob_sup), color = "red", linewidth = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = postprob_sup, ymax = 100), fill = "red", alpha = 0.7) +
    ggplot2::scale_x_continuous(name = "Number of Enrolled Participants", breaks = cumsum(n)) +
    ggplot2::scale_y_continuous(name = "Cut-off Probability (%)", breaks = seq(0, 100, by = 20)) +
    ggplot2::geom_vline(xintercept = c(10, 20, 30), linetype = "dashed") +
    #ggplot2::ggtitle("Binary Efficacy Endpoint") +
    ggplot2::theme_minimal()
  
  Oc_tabs <- tibble::tibble(Statistic = c("Early stopping for Futility (%)",
                                          "Early stopping for Superiority (%)",
                                          "Average sample size",
                                          "Null rejection (%)"),
                            Value = c(Oc_tab$earlystopfuti_mean, Oc_tab$earlystopsupe_mean, Oc_tab$ss_mean, Oc_tab$rejectnull_mean))
  
  Oc_tabs <- dplyr::as_tibble(Oc_tabs) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  names(boundary_tab) <- c("Futility boundary response", "Futility boundary toxicity", "Superiority boundary response", "Superiority boundary toxicity")
  boundary_tab <- dplyr::as_tibble(cbind(Pars = names(boundary_tab), t(boundary_tab)))
  colnames(boundary_tab) <- c("Interim analysis", 1:nIA)
  boundary_tab <- dplyr::as_tibble(boundary_tab) %>% gridExtra::tableGrob(theme = gridExtra::ttheme_minimal(), rows = NULL)
  
  Info <- paste0("lambda=", lambda, " ", "gamma=", gamma)
  layout <- patchwork::wrap_plots(p1, Oc_tabs, boundary_tab, nrow = 3)
  layout <- layout +
    patchwork::plot_annotation(
      title = "BOP2 FE for joint efficacy and toxicity",
      subtitle = Info,
      #caption = paste0("Design Pars: n= ", n, "lambda=", lambda, "gamma=", gamma),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = "bold"))
    )
  
  return(layout)
  

}

