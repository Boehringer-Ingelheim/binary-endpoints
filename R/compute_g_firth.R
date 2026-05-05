compute_g_firth <- function (pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  flic <- pars$flic
  alpha <- pars$alpha

  # Step 1: Fit LR
  f <- paste(response, "~", treatment, paste("+", covars, collapse = " "))
  fit <- logistf::logistf(
    formula = as.formula(f),
    data = data,
    flic = flic
  )

  # Step 2: Counterfactual Predictions
  m1 <- model.matrix(fit, data = data)
  m1[, paste0(treatment, 1)] <- 1
  pt <- t(m1 %*% coef(fit)) |> plogis()
  m0 <- model.matrix(fit, data = data)
  m0[, paste0(treatment, 1)] <- 0
  pc <- t(m0 %*% coef(fit)) |> plogis()

  # Step 3: Variances
  At <- pt * (1 - pt) / nrow(data)
  Ac <- pc * (1 - pc) / nrow(data)
  V <- vcov(fit)

  # Create model matrices
  dt <- (At %*% m1)
  dc <- (Ac %*% m0)

  # Step 4: Combine (RD estimate + variance estimate by delta method)
  est <- (sum(pt - pc) / nrow(data)) |> as.vector()
  se <- (sqrt((dt %*% V %*% t(dt) + dc %*% V %*% t(dc) - 2 * dc %*% V %*% t(dt)))) |> as.vector()

  # Calculate CI
  ucl <- est + qnorm(1 - alpha / 2) * se
  lcl <- est - qnorm(1 - alpha / 2) * se

  # Prepare output
  method <- ifelse(flic, "Firth (FLIC)", "Firth")
  converged <- separation <- sep_z <- NA
  rd_est <- est
  rd_se <- se
  rd_lcl <- lcl
  rd_ucl <- ucl
  rd_pval <- 2 * (1 - pnorm (abs(est / se)))

  # Collect results
  get_results()

}
