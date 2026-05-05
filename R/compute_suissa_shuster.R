compute_suissa_shuster <- function(pars) {

  # Extract arguments
  data <- pars$data
  treatment <- pars$treatment
  response <- pars$response

  # Perform test
  sst_data <- data |>
    mutate(across({{ response }}, \(x) factor(x, levels = c(0,1)))) |>
    dplyr::select({{ treatment }}, {{ response }}) |>
    table()
  sst <- Exact::exact.test(sst_data, method = "z-pooled", to.plot = FALSE)

  # Prepare output
  method <- "suissa-shuster"
  converged <- separation <- sep_z <- NA
  rd_est = ifelse(exists("sst"), sst$estimate, NA_real_)
  rd_se <- rd_lcl <- rd_ucl <- NA_real_
  rd_pval = ifelse(exists("sst"), sst$p.value, NA_real_)

  # Collect results
  get_results()

}
