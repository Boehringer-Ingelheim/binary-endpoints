gen_data <- function (n, pop_coefs, cov_rng, cov_pars, cov_coefs, rand, treat_alloc) {

  # Create covariates
  cov <- create_covs(n, cov_rng, cov_pars)

  # Randomized treatment allocation
  if (rand == "simple") {
    z <- randomize(n, alloc = treat_alloc)
  } else if (rand == "stratified") {
    z <- data.frame(stratum = interaction(cov)) |>
      mutate(
        z = randomize(n(), alloc = treat_alloc),
        .by = stratum
      ) |>
      purrr::pluck("z")
  }

  # Make predictions
  y <- lr_predict(n, z, c(pop_coefs, cov_coefs), cov)

  # Bind together and coerce
  out <- cbind(z, cov, y) |>
    mutate(z = factor(z))

  return (out)

}
