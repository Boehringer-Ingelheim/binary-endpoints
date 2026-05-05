compute_mh_sato <- function(pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  alpha <- pars$alpha

  # Perform analysis
  out <- RobinCar::robincar_mh(
    data,
    treatment,
    response,
    covars,
    estimand = "MH",
    ci_type = "Sato"
  )$result

  # Calculate CI
  ucl <- out$estimate + qnorm(1 - alpha / 2) * out$se
  lcl <- out$estimate - qnorm(1 - alpha / 2) * out$se

  # Prepare output
  method <- "MH (Sato)"
  converged <- separation <- sep_z <- NA
  rd_est <- out$estimate
  rd_se <- out$se
  rd_lcl <- lcl
  rd_ucl <- ucl
  rd_pval <- out$`pval (2-sided)`

  # Collect results
  get_results()
}

compute_mh_mgr <- function(pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  alpha <- pars$alpha

  # Perform analysis
  out <- RobinCar::robincar_mh(
    data,
    treatment,
    response,
    covars,
    estimand = "ATE",
    ci_type = "mGR"
  )$result

  # Calculate CI
  ucl <- out$estimate + qnorm(1 - alpha / 2) * out$se
  lcl <- out$estimate - qnorm(1 - alpha / 2) * out$se

  # Prepare output
  method <- "MH (mGR)"
  converged <- separation <- sep_z <- NA
  rd_est <- out$estimate
  rd_se <- out$se
  rd_lcl <- lcl
  rd_ucl <- ucl
  rd_pval <- out$`pval (2-sided)`

  # Collect results
  get_results()
}
