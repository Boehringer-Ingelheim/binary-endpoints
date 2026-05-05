compute_ye <- function (pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  alpha <- pars$alpha
  rand <- pars$rand

  # Perform analysis
  base::suppressWarnings(
    out <- RobinCar::robincar_glm(
      df = data,
      treat_col = treatment,
      response_col = response,
      formula = paste(response, "~", treatment, paste("+", covars, collapse = " ")),
      car_strata_cols = covars,
      car_scheme = ifelse(rand == "simple", "simple", "permuted-block"),
      g_family = binomial,
      contrast_h = "diff"
    )
  )

  # Check for separation
  tmp <- detectseparation::detect_separation(model.matrix(out$main$mod), data[[response]], family = binomial())
  separation <- tmp$outcome
  sep_z <- is.infinite(tmp$coefficients[2]) |> unname()

  method <- "Ye"
  converged <- out$main$mod$converged

  if (converged) {

    out <- out$contrast$result

    # Calculate CI
    ucl <- out$estimate + qnorm(1 - alpha / 2) * out$se
    lcl <- out$estimate - qnorm(1 - alpha / 2) * out$se

    # Prepare output
    rd_est <- out$estimate |> unname()
    rd_se <- out$se
    rd_lcl <- lcl |> unname()
    rd_ucl <- ucl |> unname()
    rd_pval <- out$`pval (2-sided)` |> unname()

  } else {

    rd_est <- rd_se <- rd_lcl <- rd_ucl <-  rd_pval <- NA_real_

  }

  # Collect results
  get_results()

}
