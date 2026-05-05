compute_boot <- function(pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars
  alpha <- pars$alpha
  R <- pars$R
  stratified <- pars$boot_strat

  method <- "g-comp bootstrap"

  # Step 0.1: Adjust data frame
  data <- data |>
    mutate(A = as.numeric(.data[[treatment]]))

  # Step 0.2: Precompute formula
  f <- paste(response, "~ 0 + factor(A)", paste("+", covars, collapse = " ")) |>
    as.formula()

  # Step 1: Fit LR and get estimate
  est <- get_ate(data, f, treatment, response)

  # Step 2: Bootstrap sampling distribution
  boot_samples <- numeric(R)
  for (i in seq_len(R)) {
    data_boot <- single_bootstrap(data, treatment, stratified)
    boot_samples[i] <- get_ate(data_boot, f, treatment, response)
  }

  # Remove NA values
  boot_samples <- boot_samples[!is.na(boot_samples)]

  # Compute summary stats
  converged <- separation <- sep_z <- NA
  rd_est <- est
  rd_se  <- sd(boot_samples)
  rd_lcl <- quantile(boot_samples, probs = alpha / 2)
  rd_ucl <- quantile(boot_samples, probs = 1 - (alpha / 2))
  rd_pval <- 2 * (1 - pnorm(abs(rd_est / rd_se)))

  # Collect results
  get_results()

}
