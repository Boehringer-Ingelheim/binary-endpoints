compute_score_test <- function (pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  alpha <- pars$alpha

  method <- "Score Test"

  # Step 0: Adjust data frame
  data <- data |>
    mutate(A = as.numeric(.data[[treatment]]))

  # Step 1: Fit LR
  f <- paste(response, "~ 0 + factor(A)", paste("+", covars, collapse = " "))
  base::suppressWarnings(
    fit <- glm(f, family = binomial, data = data, x = TRUE)
  )

  if (fit$converged) {

    # Step 2: Counterfactual Predictions
    g_predict <- sapply(
      1:2,
      \(a) predict(fit, newdata = data |> mutate(A = a), type = "response")
    )

    # Step 3: Estimate Treatment Effect
    g_prob <- colMeans(g_predict)
    est <- g_prob[2] - g_prob[1]

    # Step 4: Variances
    g_deriv <- sapply(
      1:2,
      \(a) {
        x <- fit$x
        x[, 1:2] <- 0
        x[, a] <- 1
        colMeans(g_predict[, a] * (1 - g_predict[, a]) * x)
      }
    )

    g_var <- var_zhang_etal(g_prob, t(g_predict), fit, t(g_deriv))

    # Step 5: Perform score test and calculate CI
    out <- score_test_diff(g_prob, g_var, fit, 0, alpha)

    # Step 6: Prepare output
    converged <- separation <- sep_z <- NA
    rd_est <- est
    rd_se <- out$se
    rd_lcl <- out$ci[1]
    rd_ucl <- out$ci[2]
    rd_pval <- out$pval

  } else {

    rd_est <- rd_se <- rd_lcl <- rd_ucl <-  rd_pval <- NA_real_

  }

  # Collect results
  get_results()

}
