single_trial <- function (data_params, analysis_params) {

  data <- do.call(gen_data, data_params)

  full_analysis_params <- list(
    data = data,
    treatment = "z",
    response = "y",
    covars = paste0("x", 1:length(data_params$cov_rng)),
    methods= analysis_params$methods,
    alpha = analysis_params$alpha,
    type_Ge = analysis_params$type_Ge,
    type_Liu = analysis_params$type_Liu,
    flic = analysis_params$flic,
    rand = data_params$rand,
    R = analysis_params$R,
    boot_strat = analysis_params$boot_strat
  )

  res <- run_analyses(full_analysis_params)
  return (res)

}
