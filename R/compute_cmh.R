compute_cmh <- function(pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response
  covars <- pars$covars

  # Convert individual patient data to array
  dat_array <- data |>
    mutate(stratum = interaction(!!!syms(covars))) %>%
    {
      table(
        .[[treatment]],
        .[[response]],
        .$stratum
      )
    }

  # Check for minimum sample size in 2x2 tables
  zero_tables <- any(apply(dat_array, 3L, sum) < 2)

  # Perform analysis
  if (!zero_tables) {
    try(
      cmh_test <- mantelhaen.test(x = dat_array,
                                  alternative = "two.sided",
                                  correct = FALSE),
      silent = TRUE
    )
  }

  # Prepare output
  method <- "cmh test"
  converged <- separation <- sep_z <- NA
  rd_est <- rd_se <- rd_lcl <- rd_ucl <- NA_real_
  rd_pval <- ifelse(zero_tables, NA_real_, ifelse(exists("cmh_test"), cmh_test$p.value, NA_real_))

  # Collect results
  get_results()
}
