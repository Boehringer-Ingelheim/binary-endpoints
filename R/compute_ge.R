compute_ge <- function (pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  type_Ge <- pars$type_Ge
  alpha <- pars$alpha

  # Prepare margins
  marg_res <- prepare_margins(data, treatment, response, covars, type_Ge)

  converged <- marg_res$converged
  separation <- marg_res$separation
  sep_z <- marg_res$sep_z
  method = paste0("Ge (", type_Ge, ")")

  if (converged) {

    m <- marg_res$m

    # Prepare output
    rd_est <- summary(m)$AME
    rd_se <- summary(m)$SE
    rd_lcl <- rd_est - qnorm(1 - alpha / 2) * rd_se
    rd_ucl <- rd_est + qnorm(1 - alpha / 2) * rd_se
    rd_pval <- 2 * (1 - pnorm (abs(rd_est / rd_se)))

  } else {

    rd_est <- rd_se <- rd_lcl <- rd_ucl <- rd_pval <- NA_real_

  }

  # Collect results
  get_results()

}
