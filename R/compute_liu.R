compute_liu <- function (pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  type_Liu <- pars$type_Liu
  alpha <- pars$alpha

  marg_res <- prepare_margins(data, treatment, response, covars, type_Liu)

  converged <- marg_res$converged
  separation <- marg_res$separation
  sep_z <- marg_res$sep_z
  method = paste0("Liu (", type_Liu, ")")

  if (converged) {

    m <- marg_res$m

    # Prepare output
    rd_est <- summary(m)$AME
    rd_se <- sqrt(summary(m)$SE ^ 2 + var(m$dydx) / length(m$dydx)) # Partial matching
    rd_lcl <- rd_est - rd_se * qnorm(1 - alpha / 2)
    rd_ucl <- rd_est + rd_se * qnorm(1 - alpha / 2)
    rd_pval <- 2 * (1 - pnorm (abs(rd_est / rd_se)))

  } else {

    rd_est <- rd_se <- rd_lcl <- rd_ucl <- rd_pval <- NA_real_

  }

  # Collect results
  get_results()

}
