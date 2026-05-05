single_scenario <- function(n_sim,
                            params,
                            methods,
                            alpha,
                            type_Ge, type_Liu, flic, R, boot_strat) {

  n <- params$n
  prev <- params$prev
  diff <- params$diff
  cov_rng <- params$cov_rng[[1]]
  cov_pars <- params$cov_pars[[1]]
  cov_coefs <- params$cov_coefs[[1]]
  rand <- params$rand
  treat_alloc <- params$treat_alloc[[1]]
  pop_coefs = params$pop_coefs[[1]]

  data_params <- list(
    n = n,
    pop_coefs = pop_coefs,
    cov_rng = cov_rng,
    cov_pars = cov_pars,
    cov_coefs = cov_coefs,
    rand = rand,
    treat_alloc = treat_alloc
  )
  analysis_params <- list(
    methods = methods,
    alpha = alpha,
    type_Ge = type_Ge,
    type_Liu = type_Liu,
    flic = flic,
    R = R,
    boot_strat = boot_strat
  )

  # Sequential execution (parallel handled by Slurm)
  out <- purrr::map(seq_len(n_sim), ~ single_trial(data_params, analysis_params))

  out |>
    dplyr::bind_rows(.id = "run") |>
    dplyr::bind_cols(
      tibble::tibble(
        n, prev, diff,
        list(cov_coefs),
        rand,
        list(treat_alloc)
      )
    )
}
