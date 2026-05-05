find_coefs <- function (prev, diff, cov, coefs, tol = 1e-3, max_iter = 100, int_tol = 1e-3) {

  lprobs <- as.matrix(cov) %*% coefs
  b0 <- find_b0_cpp(lprobs, prev, tol, max_iter, int_tol)
  lprobs <- as.matrix(cbind(1, cov)) %*% c(b0, coefs)
  bz <- find_bz_cpp(lprobs, diff, tol, max_iter, int_tol)
  return (c(b0 = b0, bz = bz))

}

randomize <- function (n, alloc) {

  ni <- floor(n / sum(alloc))
  z <- rep(1:0, times = (ni * alloc))
  nrest <- n - length(z)
  z <- c(z, sample(rep(1:0, times = alloc), nrest))
  sample(z)

}

create_covs <- function (n, cov_rng, cov_pars) {

  cov_list <- list()
  for (i in 1:length(cov_rng)) {
    cov_list[[paste0("x", i)]] <- do.call(cov_rng[i], c(n, as.list(cov_pars[[i]])))
  }

  return (as.data.frame(cov_list))

}

lr_predict <- function (n, z, coefs, cov) {

  lprobs <-  as.matrix(cbind(1, z, cov)) %*% coefs
  probs <- plogis(lprobs)
  rbinom(n, 1, probs)

}

def_population_coefs <- function (prev, diff, cov_rng, cov_pars, cov_coefs, n_pop = 1e5) {

  # Create covariates
  cov <- create_covs(n_pop, cov_rng, cov_pars)

  # Determine true regression coefficients for specified PATE
  find_coefs(prev, diff, cov, cov_coefs)

}

run_analyses <- function (pars) {

  lapply(
    pars$methods,
    \(x) do.call(x, list(pars))
  ) |>
    bind_rows()

}

get_results <- function(...) {
  env <- parent.frame()
  tibble(
    "method" = env$method,
    "converged" = env$converged,
    "separation" = env$separation,
    "separation_treatment" = env$sep_z,
    "rd_est" = env$rd_est,
    "rd_se" = env$rd_se,
    "rd_lcl" = env$rd_lcl,
    "rd_ucl" = env$rd_ucl,
    "rd_pval" = env$rd_pval
  )
}

prepare_margins <- function (data,
                             treatment,
                             response,
                             covars,
                             type) {

  # Fit logistic regression model
  base::suppressWarnings(
    out <- glm(
      paste(response, "~", treatment, paste("+", covars, collapse = " ")),
      family = binomial,
      data = data
    )
  )
  converged <- out$converged

  # Check for separation
  tmp <- detectseparation::detect_separation(model.matrix(out), data[[response]], family = binomial())
  separation <- tmp$outcome
  sep_z <- is.infinite(tmp$coefficients[2])

  if (converged) {

    # Calculate marginal effects
    m <- margins::margins(
      model = out,
      variables = treatment,
      vcov = if (type == "orig") vcov(out) else sandwich::vcovHC(out, type = type)
    )

  } else {

    m <- NULL

  }

  return (
    list(
      converged = converged,
      separation = separation,
      sep_z = sep_z,
      m = m,
      type = type
    )
  )

}

var_zhang_etal <- function(mu, g, wm, muDbeta) {

  ## R code is based on code provided as supporting information for the
  ## article 'A Robust Score Test in G‐Computation for Covariate Adjustment
  ## in Randomized Clinical Trials Leveraging Different Variance Estimators
  ## via Influence Functions' by Zhang et al. (2025, Stat. Med)

  # Variance estimator based on M-estimation
  # mu:      g-comp estimate per arm (k arms)
  # g:       predicted values of each subject for each arm (dim = c(k, n))
  # wm:      fitted working model
  # muDbeta: derivative of mu per arm w.r.t. beta (dim = c(k, p))

  # influence function for beta
  betaPhi <- (vcov(wm) * nobs(wm)) %*% t(wm$x * (wm$y - wm$fitted.values))

  # influence function for mu
  muPhi <- muDbeta %*% betaPhi + g - mu  # dim = c(k, n)

  # variance for mu
  Sigma <- var(t(muPhi)) / nobs(wm)

  return (Sigma)

}

score_test_diff <- function(mu, Sigma, wm, h0, alpha) {

  ## R code is based on code provided as supporting information for the
  ## article 'A Robust Score Test in G‐Computation for Covariate Adjustment
  ## in Randomized Clinical Trials Leveraging Different Variance Estimators
  ## via Influence Functions' by Zhang et al. (2025, Stat. Med)

  # mu:    g-comp estimates per arm
  # Sigma: Variance of per-arm g-comp estimators
  # wm:    working model
  # h0:    difference under H0
  # alpha: Type I error rate

  level <- (1 - alpha)
  v <- Sigma[1, 1] - 2 * Sigma[1, 2] + Sigma[2, 2]
  z <- (mu[2] - mu[1] - h0) / sqrt(v + (mu[2] - mu[1] - h0)^2 / nobs(wm))
  pval <- 2 * (1 - pnorm(abs(z))) # 2-sided
  ci <- mu[2] - mu[1] + sqrt(v * qchisq(level, df = 1) /
                               (1 - qchisq(level, df = 1) / nobs(wm))) * c(-1, 1)

  return(
    list(
      se = sqrt(v + (mu[2] - mu[1] - h0)^2 / nobs(wm)),
      pval = pval,
      ci = ci
    )
  )

}

# Fast calculation of g-computation estimator for bootstrap analysis
get_ate <- function(data, formula, treatment, response) {

  # Build design matrix and fit LR
  X <- model.matrix(formula, data)
  y <- data[[response]]

  base::suppressWarnings(
    out <- tryCatch(fastglm::fastglm(X, y, family = binomial()), error = \(x) return(NULL))
  )
  if (is.null(out)) return(NA_real_)

  # Change design matrix to produce counterfactual predictions
  X[,1] <- 1
  X[,2] <- 0
  pc <- tryCatch(fastglm:::predict.fastglm(out, X, type = "response"), error = function(e) return(rep(NA_real_, nrow(data))))

  X[,1] <- 0
  X[,2] <- 1
  pt <- tryCatch(fastglm:::predict.fastglm(out, X, type = "response"), error = function(e) return(rep(NA_real_, nrow(data))))

  # Estimate and return ATE
  sum(pt - pc) / nrow(data)

}

single_bootstrap <- function(data, treatment, stratified) {

  idx <- 0

  while (length(unique(data[idx, treatment])) < 2) {
    if (stratified) {
      idx <- lapply(
        X = split(seq_len(nrow(data)), data[[treatment]]),
        FUN = \(x) sample(x, length(x), replace = TRUE)
      ) |>
        unlist()
    } else {
      idx <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    }
  }

  return (data[idx, , drop = FALSE])

}
