#-------------------------------------------------------------------------------
# required packages
require(MASS)
require(glmnet)

#-------------------------------------------------------------------------------
# functions

#' Inverse logit
#' 
#' Computes the inverse of the logit function (expit)
#' 
#' @param x argument (numeric vector).
#' 
#' @return 
#' Vector of the same length as x.
expit <- function(x) 1 / (1 + exp(-x))

#' Get e-value from numerator and denominator likelihood and observation
#' 
#' @param numerator numerator probability.
#' @param denominator denominator probability.
#' @param y observations (binary).
#' @param unbiased use the method guaranteed to have expected value 1 as 
#'     described in the paper (only for E-CRT).
#' @param M number of simulations for unbiased method.
#' 
#' @return 
#' Vector of e-values, same length as input.
get_e <- function(numerator, denominator, y, unbiased = FALSE, M = NULL) {
  if (unbiased) {
    denominator <- numerator / (M + 1) + M / (M + 1) * denominator
  }
  out <- (1 - y - numerator) / (1 - y - denominator)
  out[numerator == denominator] <- 1
  out
}

#' Find first time where an e-process exceeds the level 1 / alpha
#' 
#' @param e e-process (vector).
#' @param alphas vector/matrix of alphas of interest. If matrix: number of 
#'     columns should match the length of e, and each column give the 
#'     respective level for rejecting at the given stage.     
#' @param tmax maximal time to return if the process does not exceed the level.
#' 
#' @return 
#' A data.frame containing the sample sizes and whether the level was crossed.
get_first_rejection <- function(e, alphas, tmax) {
  m <- NROW(alphas)
  if (!is.matrix(alphas)) alphas <- matrix(alphas, ncol = 1)
  size <- integer(m)
  rejected <- logical(m)
  for (j in seq_len(m)) {
    if (max(e) >= 1 / alphas[j, ]) {
      rejected[j] <- TRUE
      size[j] <- min(which(e >= 1 / alphas[j, ]))
    } else {
      size[j] <- tmax
    }
  }
  data.frame(size = size, rejected = rejected, alphas = alphas)
}

#' Find worst case and average sample size for rejecting with given power
#' 
#' @param times vector of times where a sequential test has achieved a maximum
#'     test statistics.
#' @param rejected whether the hypothesis was rejected at the times.
#' @param beta desired power (1 minus type II error).
#' @param tmax maximum time to use if hypothesis was not rejected.
#' 
#' @return 
#' A data.frame containing the minimum sample size required to reject with the
#' given power, the average sample size, and whether the power is achieved.
find_sample_size <- function(times, rejected, beta, tmax) {
  times[!rejected] <- tmax
  worst <- quantile(times, probs = beta, type = 1)
  power_achieved <- (mean(rejected[times <= worst]) >= beta)
  average <- mean(pmin(times, worst))
  data.frame(
    worst = worst,
    average = average,
    power_achieved = power_achieved,
    beta = beta
  )
}

#' Try to fit a logistic regression model
#' 
#' Fits a logistic regression model with glm to the data and returns FALSE if
#' this creates errors.
#' 
#' @param formula model formula.
#' @param data data for the formula (containing all variables).
#' 
#' @return
#' The fitted model, or FALSE if an error occurred.
try_logistic_glm <- function(formula, data) {
  tryCatch(
    expr = glm(formula = formula, family = binomial, data = data),
    error = function(out) FALSE
  )
}

#' Try to fit a penalized logistic regression model
#' 
#' Fits a penalized logistic regression model with glmnet. Returns FALSE if
#' this results in an error.
#' 
#' @param x matrix of covariates.
#' @param y vector of observations (binary).
#' @param lambda the regularization parameter.
#' @param not_penalize column positions of variables in x which should NOT be
#'     penalized. None if this is NULL.
#' @param fifb1 return FALSE if the first coefficient (beta[1,1]) is zero
#'     (only for running MLE)
#' 
#' @return
#' The fitted model, or FALSE if an error occurred (or if the coefficient of
#' interest equals zero and fifb1 is TRUE).
try_logistic_glm_penalized <- function(
    x,
    y,
    lambda,
    not_penalize = NULL,
    fifb1 = FALSE
  ) {
  penalty.factor <- rep(1, ncol(x))
  if (!is.null(not_penalize)) penalty.factor[not_penalize] <- 0
  out <- tryCatch(
    expr =  glmnet::glmnet(
      x = x,
      y = y,
      family = "binomial",
      lambda = lambda,
      penalty.factor = penalty.factor
    ),
    error = function(out) FALSE
  )
  if (!isFALSE(out)) {
    if (fifb1 && out$beta[1, 1] == 0) return(FALSE)
  }
  out
}

#' Try to apply cv.glmnet for logistic regression
#'
#' Applies cv.glmnet to given data. Returns FALSE if this resulted in an error,
#' and otherwise the lambda.1se parameter.
#'
#' @param x matrix of covariates.
#' @param y vector of observations (binary).
#' @param not_penalize column positions of variables in x which should NOT be
#'   penalized. None if this is NULL.
#'
#' @return The fitted model, or FALSE if an error occurred.
try_cv_glmnet <- function(x, y, not_penalize = NULL) {
  penalty.factor <- rep(1, ncol(x))
  if (!is.null(not_penalize)) penalty.factor[not_penalize] <- 0
  tryCatch(
    expr = glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      penalty.factor = penalty.factor
    )$lambda.min,
    error = function(out) FALSE
  )
}

#' Predict from a logistic glm
#' 
#' Predict from the model with the new data and formats the output as an
#' unnamed vector.
#' 
#' @param logistic regression model (fitted with glm).
#' @param newdata data.frame containing all relevant variables.
#' 
#' @return 
#' A vector of length nrow(newdata) containing the predictions (probability of
#' outcome 1, not the linear predictor).
predict_logistic_glm <- function(model, newdata) {
  c(predict(model, newdata, type = "response"), use.names = FALSE)
}

#' Simulate from conditional of multivariate normal distribution
#' 
#' @param z the values on which the distribution is conditioned (a vector or a
#'     matrix, with the columns of z corresponding to the values to condition
#'     on).
#' @param Mu the joint mean vector. The last entries correspond to the mean of
#'     the vector (or columns of the matrix) z.
#' @param Sigma the joint covariance matrix, with the lower right square being
#'     the variance matrix corresponding to z.
#'     
#' @return 
#' Either a vector of length nsim, or a matrix with nsim columns and ncol(z)
#' rows.
sim_normal_covariate <- function(z, Mu, Sigma, nsim) {
  d1 <- seq_len(NROW(z))
  d2 <- (NROW(z) + 1):nrow(Sigma)
  mu1 <- Mu[d1]
  mu2 <- Mu[d2]
  sigma11 <- Sigma[d1, d1]
  sigma12 <- Sigma[d1, d2]
  sigma21 <- Sigma[d2, d1]
  sigma22 <- Sigma[d2, d2]
  newmu <- c(mu2 + sigma12 %*% solve(sigma11, z - mu1))
  newsigma <- c(sigma22 - sigma12 %*% solve(sigma11, sigma12))
  out <- rnorm(nsim * NCOL(z), mean = newmu, sd = sqrt(newsigma))
  if (is.matrix(z)) out <- matrix(out, nrow = NCOL(z))
  out
}

#' Sequential conditional randomization test for logistic regression
#' 
#' Sequentially computes conditional randomization e-variables
#' 
#' @param x covariate of interest (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param sim_fun a function taking as input a one (vector of) nuisance 
#'     covariate(s) and simulates the corresponding x covariates from it.
#' @param n0 minimal sample size for which the e-variables are computed.
#' @param eps small positive number. Probabilities are truncated to lie in the
#'     interval [eps, 1-eps]. Default is 0 (no truncation).
#' 
#' @return
#' A data.frame containing the numerator and the denominator of the e-variables,
#' and the observations.
e_crt_logistic <- function(x, y, z, sim_fun, n0 = NULL, eps = 0) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  n <- nrow(data)
  numerator <- rep(1, n)
  denominator <- rep(1, n)
  formula <- formula(paste("y ~", paste(colnames(data)[-1], collapse = "+")))
  for (k in (n0 + 1):n) {
    model <- try_logistic_glm(formula, data[seq_len(k - 1), ])
    if (isFALSE(model)) next
    numerator[k] <- 
      pmin(pmax(eps, predict_logistic_glm(model, data[k, ])), 1 - eps)
    sims <- sim_fun(z[k, ])
    newdata <- cbind(x = sims, as.data.frame(z[k, , drop = FALSE]))
    sim_pred <- pmin(pmax(eps, predict_logistic_glm(model, newdata)), 1 - eps)
    denominator[k] <- mean(sim_pred)
  }
  data.frame(numerator = numerator, denominator = denominator, y = y)
}

#' Sequential conditional randomization test for logistic regression with
#' Gaussian distribution for X|Z
#' 
#' Sequentially computes conditional randomization e-variables. The 
#' distribution of X|Z is either given, or estimated on available data.
#' It is assumed that X|Z is normal, with m
#' 
#' @param x covariate of interest (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param sim_fun a function taking as input a one (vector of) nuisance 
#'     covariate(s) and simulates the corresponding x covariates from it.
#' @param n0 minimal sample size for which the e-variables are computed.
#' @param eps small positive number. Probabilities are truncated to lie in the
#'     interval [eps, 1-eps]. Default is 0 (no truncation).
#' @param coef vector containing coefficients and variance for X|Z (can be
#'     omitted if xz_sample is given).
#' @param xz_sample sample of X and Z to estimate the conditional distribution
#'     of X|Z. Can be omitted if coef is given.
#' @param update_sample if TRUE, at each step the available data is added to
#'     xz_sample to re-estimate the conditional distribution.
#' @param bias_fun a bias function to modify the conditional mean of X|Z.
#' @param nsim number of simulations for conditional randomization.
#' 
#' @return
#' A data.frame containing the numerator and the denominator of the e-variables,
#' and the observations.
e_crt_logistic_gaussian <- function(
    x,
    y,
    z,
    n0 = NULL,
    eps = 0,
    coef = NULL,
    xz_sample = NULL,
    update_sample = FALSE,
    bias_fun = identity,
    nsim = 500
  ) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  n <- nrow(data)
  numerator <- rep(1, n)
  denominator <- rep(1, n)
  formula <- formula(paste("y ~", paste(colnames(data)[-1], collapse = "+")))
  if (!is.null(coef) || !update_sample) {
    if (!is.null(coef)) {
      sdev <- sqrt(coef[length(coef)])
      beta <- coef[-length(coef)]
    } else {
      colnames(xz_sample) <- c("x", colnames(z))
      xz_sample <- as.data.frame(xz_sample)
      model <- lm(x ~ . - 1, data = xz_sample)
      beta <- coef(model)
      sdev <- sqrt(mean(residuals(model)^2))
    }
    sim_fun <- function(zz) rnorm(nsim, mean = bias_fun(sum(beta * zz)), sd = sdev)
    for (k in (n0 + 1):n) {
      model <- try_logistic_glm(formula, data[seq_len(k - 1), ])
      if (isFALSE(model)) next
      numerator[k] <- 
        pmin(pmax(eps, predict_logistic_glm(model, data[k, ])), 1 - eps)
      sims <- sim_fun(z[k, ])
      newdata <- cbind(x = sims, as.data.frame(z[k, , drop = FALSE]))
      sim_pred <- pmin(pmax(eps, predict_logistic_glm(model, newdata)), 1 - eps)
      denominator[k] <- mean(sim_pred)
    }
  } else {
    n_xz <- nrow(xz_sample) - 1
    xz_sample <- rbind(xz_sample, cbind(x, z))
    colnames(xz_sample) <- c("x", colnames(z))
    xz_sample <- as.data.frame(xz_sample)
    xz_formula <- formula("x ~ . - 1")
    for (k in (n0 + 1):n) {
      model <- try_logistic_glm(formula, data[seq_len(k - 1), ])
      if (isFALSE(model)) next
      model_xz <- lm(xz_formula, data = xz_sample[seq_len(n_xz + k), ], )
      sdev <- sqrt(mean(residuals(model_xz)^2))
      numerator[k] <- 
        pmin(pmax(eps, predict_logistic_glm(model, data[k, ])), 1 - eps)
      sims <- rnorm(nsim, mean = bias_fun(predict(model_xz, data[k, ])), sd = sdev)
      newdata <- cbind(x = sims, as.data.frame(z[k, , drop = FALSE]))
      sim_pred <- pmin(pmax(eps, predict_logistic_glm(model, newdata)), 1 - eps)
      denominator[k] <- mean(sim_pred)
    }
  }
  data.frame(numerator = numerator, denominator = denominator, y = y)
}


#' Oracle CRT e-values for loistic regression
#' 
#' Computes CRT e-values with the true event probabilities.
#' 
#' @param x covariate (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param sim_fun a function taking as input a one (vector of) nuisance 
#'     covariate(s) and simulates the corresponding x covariates from it.
#' @param par parameter vector (including intercept). The parameters should be
#'     in the order (1, z, x), where "1" refers to the intercept, "z" to all
#'     parameters for the z covariates, and "x" to the parameter for x.
#' 
#' @return
#' A data.frame containing the numerator and the denominator of the e-variables,
#' and the observations.
e_crt_logistic_oracle <- function(x, y, z, sim_fun, par) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
  }
  covariate_mat <- cbind(1, z, x)
  numerator <- c(expit(covariate_mat %*% par))
  sims <- sim_fun(t(z))
  linear_predictor <- sims * par[length(par)] + 
     c(covariate_mat[, -ncol(covariate_mat)] %*% par[-length(par)])
  denominator <- rowMeans(expit(linear_predictor))
  data.frame(numerator = numerator, denominator = denominator, y = y)
}

#' Running MLE logistic
#' 
#' Testing in a logistic regression model with the running maximum likelihood
#' method.
#' 
#' @param x covariate (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param n0 minimal sample size for which the e-variables are computed.
#' 
#' @return
#' A data.frame containing the likelihood ratio process M.
running_mle_logistic <- function(x, y, z, n0) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  n <- nrow(data)
  M <- rep(1, n)
  lnumerator <- rep(0, n)
  formula_full <- 
    formula(paste("y ~", paste(colnames(data)[-1], collapse = "+")))
  formula_null <- formula(paste("y ~", paste(colnames(z), collapse = "+")))
  q <- ncol(data)
  for (k in (n0+1):(n0 + q)) {
    # when the number of variables is smaller than the sample size, the MLE in
    # the denominator is a perfect fit, so the denominator equals 1
    model_full <- try_logistic_glm(formula_full, data[seq_len(k - 1), ])
    if (isFALSE(model_full)) next
    lnumerator[k] <- lnumerator[k - 1] + 
      log(abs(predict_logistic_glm(model_full, data[k, ]) - 1 + y[k]))
    M[k] <- exp(lnumerator[k])
  }
  for (k in (n0 + q):n) {
    model_full <- try_logistic_glm(formula_full, data[seq_len(k - 1), ])
    model_null <- try_logistic_glm(formula_null, data[(n0 + 1):k, ])
    if (isFALSE(model_null) || isFALSE(model_full)) {
      lnumerator[k] <- lnumerator[k - 1]
      M[k] <- M[k - 1]
    }
    lnumerator[k] <- lnumerator[k - 1] +
      log(abs(predict_logistic_glm(model_full, data[k, ]) - 1 + y[k]))
    M[k] <- exp(lnumerator[k] - 
                  sum(log(abs(fitted(model_null) + y[(n0 + 1):k] - 1))))
  }
  data.frame(M = M) 
}

#' Running MLE logistic with penalized alternative
#' 
#' Testing in a logistic regression model with the running maximum likelihood
#' method. The estimator under the alternative uses L1 penalization.
#' 
#' @param x covariate (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param n0 minimal sample size for which the e-variables are computed.
#' @param not_penalized indices of covariates not to penalize in the model. The
#'     covariates will be combined as cbind(x, z) (first x, then z).
#' 
#' @details 
#' The penalization parameter is estimated by cross-validation (cv.glmnet) and
#' only updated every 10 observations to reduce the computation time.
#' 
#' @return
#' A data.frame containing the likelihood ratio process M.
running_mle_logistic_penalized <- function(
    x,
    y,
    z,
    n0,
    not_penalize = NULL
  ) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  covariates <- cbind(x, z)
  n <- nrow(data)
  M <- rep(1, n)
  lnumerator <- rep(0, n)
  formula_null <- formula(paste("y ~", paste(colnames(z), collapse = "+")))
  cv_lambda <- try_cv_glmnet(
    x = covariates[seq_len(n0), ],
    y = y[seq_len(n0)]
  )
  q <- ncol(data)
  for (k in (n0 + 1):(n0 + q)) {
    if (isFALSE(cv_lambda) || ((k  - 1L) %% 10L == 0L)) {
      # update lambda every 10th observation, or try until it works if something
      # went wrong in the previous update
      cv_lambda <- try_cv_glmnet(
        x = covariates[seq_len(k - 1), ],
        y = y[seq_len(k - 1)],
        not_penalize = not_penalize
      )
    }
    if (isFALSE(cv_lambda)) {
      last_false <- k
      lnumerator[k] <- lnumerator[k - 1]
      M[k] <- M[k - 1]
      next
    }
    model_full <- try_logistic_glm_penalized(
      x = covariates[seq_len(k - 1), ],
      y = y[seq_len(k - 1)],
      lambda = cv_lambda,
      not_penalize = not_penalize
    )
    if (isFALSE(model_full)) {
      lnumerator[k] <- lnumerator[k - 1]
      M[k] <- M[k - 1]
      next
    }
    lnumerator[k] <- lnumerator[k - 1] +
      log(abs(predict_logistic_glm(model_full, covariates[k, , drop = FALSE]) - 1 + y[k]))
    M[k] <- exp(lnumerator[k])
  }
  for (k in (n0 + q):n) {
    if (isFALSE(cv_lambda) || ((k  - 1L) %% 10L == 0L)) {
      cv_lambda <- try_cv_glmnet(
        x = covariates[seq_len(k - 1), ],
        y = y[seq_len(k - 1)],
        not_penalize = not_penalize
      )
    }
    if (isFALSE(cv_lambda)) {
      lnumerator[k] <- lnumerator[k - 1]
      M[k] <- M[k - 1]
      next
    }
    model_full <- try_logistic_glm_penalized(
      x = covariates[seq_len(k - 1), ],
      y = y[seq_len(k - 1)],
      lambda = cv_lambda,
      not_penalize = not_penalize
    )
    model_null <- try_logistic_glm(formula_null, data[(n0 + 1):k, ])
    if (isFALSE(model_full) || isFALSE(model_null)) {
      lnumerator[k] <- lnumerator[k - 1]
      M[k] <- M[k - 1]
      next
    }
    lnumerator[k] <- lnumerator[k - 1] +
      log(abs(predict_logistic_glm(model_full, covariates[k, , drop = FALSE]) - 1 + y[k]))
    M[k] <- exp(lnumerator[k] -
                  sum(log(abs(fitted(model_null) + y[(n0 + 1):k] - 1))))
  }
  data.frame(M = M) 
}

#' Conditional randomization test with Gaussian distribution for X|Z
#' 
#' Computes the conditional randomization test for logistic regression. Also
#' computes the p-value when the test statistic is absolute correlation between
#' Y and X, instead of the likelihood of the logistic model.
#' 
#' @param x covariate (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param stops sample sizes at which the conditional randomization test is
#'     performed.
#' @param coef vector containing coefficients and variance for X|Z (can be
#'     omitted if xz_sample is given).
#' @param xz_sample sample of X and Z to estimate the conditional distribution
#'     of X|Z. Can be omitted if coef is given.
#' @param reuse_data if TRUE, the conditional distribution of X|Z is estimated
#'     on the same data which is used for testing.
#' @param bias_fun a bias function to modify the conditional mean of X|Z.
#' @param nsim number of simulations for conditional randomization.
#'     
#' @return
#' A data.frame containing the p-values of the tests, for all desired sample 
#' sizes.
crt_logistic_gaussian <- function(
    x,
    y,
    z,
    stops,
    sim_fun,
    coef = NULL,
    xz_sample = NULL,
    reuse_data = NULL,
    bias_fun = identity,
    nsim = 500
  ) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  p_crt <- rep(1, length(stops))
  p_crt_cor <- rep(1, length(stops))
  form <- formula(paste("y ~", paste(colnames(data)[-1], collapse = "+")))
  n <- length(x)
  if (!is.null(coef) | !is.null(xz_sample)) {
    if (!is.null(coef)) {
      sdev <- sqrt(coef[length(coef)])
      beta <- coef[-length(coef)]
    } else {
      model <- lm(xz_sample[, 1] ~ xz_sample[, -1] - 1)
      beta <- coef(model)
      sdev <- sqrt(mean(residuals(model)^2))
    }
    means <- bias_fun(z %*% beta)
    sims <- matrix(nrow = n, ncol = nsim)
    for (j in seq_len(nsim)) sims[, j] <- rnorm(n, means, sdev)
    for (k in seq_along(stops)) {
      nn <- seq_len(stops[k])
      model <- try_logistic_glm(form, data[nn, ])
      if (isFALSE(model)) next
      llikelihood <- sum(log(abs(1 - y[nn] - fitted(model))))
      llikelihood_sim <- rep(NA, nsim)
      cor_true <- abs(cor(y[nn], x[nn]))
      cor_sim <- rep(NA, nsim)
      for (s in seq_len(ncol(sims))) {
        cor_sim[s] <- abs(cor(y[nn], sims[nn, s]))
        df_new <- data[seq_len(stops[k]), ]
        df_new$x <- sims[nn, s]
        sim_model <- try_logistic_glm(form, df_new)
        if (isFALSE(sim_model)) next
        llikelihood_sim[s] <- sum(log(abs(1 - y[nn] - fitted(sim_model))))
      }
      p_crt[k] <- 
        (1 + sum(llikelihood_sim >= llikelihood, na.rm = TRUE)) / (nsim + 1)
      p_crt_cor[k] <- 
        (1 + sum(cor_sim >= cor_true, na.rm = TRUE)) / (nsim + 1)
    }
  } else {
    form_xz <- formula(paste0("x ~", paste0(colnames(z), collapse = "+"), -1))
    for (k in seq_along(stops)) {
      nn <- seq_len(stops[k])
      model <- try_logistic_glm(form, data[nn, ])
      if (isFALSE(model)) next
      model_xz <- lm(form_xz, data = data[nn, ])
      means <- bias_fun(fitted(model_xz))
      sdev <- sqrt(mean(residuals(model_xz)^2))
      llikelihood <- sum(log(abs(1 - y[nn] - fitted(model))))
      llikelihood_sim <- rep(NA, nsim)
      cor_true <- abs(cor(y[nn], x[nn]))
      cor_sim <- rep(NA, nsim)
      for (s in seq_len(ncol(sims))) {
        df_new <- data[nn, ]
        sim <- sim
        df_new$x <- sim
        cor_sim[s] <- abs(cor(y[nn], sim))
        sim_model <- try_logistic_glm(form, df_new)
        if (isFALSE(sim_model)) next
        llikelihood_sim[s] <- sum(log(abs(1 - y[nn] - fitted(sim_model))))
      }
      p_crt[k] <- 
        (1 + sum(llikelihood_sim >= llikelihood, na.rm = TRUE)) / (nsim + 1)
      p_crt_cor[k] <- 
        (1 + sum(cor_sim >= cor_true, na.rm = TRUE)) / (nsim + 1)
    }
  }
  data.frame(stops = stops, p_crt = p_crt, p_crt_cor = p_crt_cor) 
}


#' Conditional randomization test and asymptotic test for logistic regression
#' 
#' Computes the conditional randomization test and the likelihood ratio test for
#' logistic regression.
#' 
#' @param x covariate (vector)
#' @param y response variable (binary vector).
#' @param z nuisance covariate(s) (either a single vector or a matrix).
#' @param stops sample sizes at which the conditional randomization test is
#'     performed.
#' @param sim_fun a function taking as input a one (vector of) nuisance 
#'     covariate(s) and simulates the corresponding x covariates from it.
#'     
#' @return
#' A data.frame containing the p-values of the tests, for all desired sample 
#' sizes.
crt_lrt_logistic <- function(x, y, z, stops, sim_fun) {
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
    colnames(z) <- "z1"
  } else {
    colnames(z) <- paste0("z", seq_len(ncol(z)))
  }
  data <- as.data.frame(cbind(y = y, x = x, z))
  p_crt <- rep(1, length(stops))
  p_asymptotic <- rep(1, length(stops))
  form <- formula(paste("y ~", paste(colnames(data)[-1], collapse = "+")))
  sims <- sim_fun(t(z))
  nsim <- ncol(sims)
  for (k in seq_along(stops)) {
    nn <- seq_len(stops[k])
    model <- try_logistic_glm(form, data[nn, ])
    if (isFALSE(model)) next
    p_asymptotic[k] <- summary(model)$coefficients["x", 4]
    llikelihood <- sum(log(abs(1 - y[nn] - fitted(model))))
    llikelihood_sim <- rep(NA, nsim)
    for (s in seq_len(ncol(sims))) {
      df_new <- data[seq_len(stops[k]), ]
      df_new$x <- sims[seq_len(stops[k]), s]
      sim_model <- try_logistic_glm(form, df_new)
      if (isFALSE(sim_model)) next
      llikelihood_sim[s] <- sum(log(abs(1 - y[nn] - fitted(sim_model))))
    }
    p_crt[k] <- 
      (1 + sum(llikelihood_sim >= llikelihood, na.rm = TRUE)) / (nsim + 1)
  }
  data.frame(stops = stops, p_crt = p_crt, p_asymptotic = p_asymptotic) 
}